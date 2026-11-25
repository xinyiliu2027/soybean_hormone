# 完整独立续跑v2：拆分内皮层胁迫/激素图，新增全部细胞类型图。
# 已有内皮层结果直接读取；全细胞图使用相同regulon、RDS和已保存P1/P2/P3。
# 不重跑GENIE3，不重新聚类，不重新拟时序，不改变已有结果目录。
# source("/data01/liuxy/Project/TF_split_stress_hormones_allcelltypes_v2.R",encoding="UTF-8")
# 本脚本在Source / Rscript时自动运行一次，不必再输入main()。
# 如只想加载函数：options(tf_layout_v2.skip_main=TRUE) 后 source(...)。

PARENT_DIR <- "/data01/liuxy/Project/CK2_added_v3_Results/Ensemble_PI_Heat_Cold_PEG_with_CK2"
# 已知上次结果目录时可直接填写，例如 file.path(PARENT_DIR,"TF_activity_rerun_20261028_123456")。
# NULL：优先采用当前会话已有且含结果的OUT_DIR，否则查找最近完成的一次结果；不新建目录。
SOURCE_RUN_DIR <- NULL
SESSION_RUN_DIR <- if(exists("OUT_DIR",inherits=FALSE)&&is.character(OUT_DIR)&&length(OUT_DIR)==1L)OUT_DIR else NULL
SOURCE_PREFIX <- "TF_exactCK_P123_all8_"
PREFIX <- "TF_layout_split_alltypes_v2_"
INPUT_RDS <- "/data01/liuxy/Project/muscat_results_ck_ref_paper/CK_reference/CKref_merged_seurat_annotated.rds"
OUT_DIR <- NULL # main()中固定为已有SOURCE_RUN_DIR
CELLTYPE_COL <- "cell_type"
TREATMENT_COL <- "treatment"
SAMPLE_COL <- "orig.ident"
STATES <- c("P1","P2","P3")
STRESS <- c("Heat","Cold","PEG")
HORMONES <- c("ABA","IAA","MeJA","eBL","SA")
TREATMENTS <- c(STRESS,HORMONES)
GROUPS <- c("CK",TREATMENTS)
ENDO_TYPES <- paste0("Endodermis_",STATES)
TYPE_ORDER <- c("Root Cap","Epidermis","Cortex","Dividing Cell","Meristem","QC",
                ENDO_TYPES,"Pericycle","Procambium","Phloem","Xylem")
MIN_CELLS_PER_SAMPLE_STAGE <- 10L # 如上次configuration有保存，以保存值为准
MIN_REPLICATES <- 2L
TOP_HEATMAP_TFS <- 60L # 没找到上一张图的列清单时才使用
FOCUS_TFS <- character()
UCELL_MAX_RANK <- 1500L # 以保存结果的maxRank为准
SCORE_CHUNK <- 1000L
REUSE_ALLTYPE_SUMMARY <- TRUE # 对输入文件内容、regulon、标签、软件版本等校验后才复用
RUN_ALL_CELLTYPES <- TRUE
# 轻量可复现绘图，不要求ggplot2/pheatmap；全细胞评分需要SeuratObject/Matrix/UCell/BiocParallel。
CELL_WIDTH_MM <- 2.6
CELL_HEIGHT_MM <- 4.7
TF_FONT_PT <- 10.5
CELLTYPE_FONT_PT <- 12
LEGEND_FONT_PT <- 10.5
TITLE_FONT_PT <- 14
DENDRO_HEIGHT_MM <- 18
PNG_DPI <- 250
TF_LABEL_ANGLE <- 90
FIGURES_WRITTEN <- character()
out <- function(x)file.path(OUT_DIR,paste0(PREFIX,x))
srcpath <- function(x)file.path(OUT_DIR,paste0(SOURCE_PREFIX,x))
tabout <- function(x,name)utils::write.csv(x,out(name),row.names=FALSE,na="",fileEncoding="UTF-8")
msg <- function(...)message(format(Sys.time(),"%H:%M:%S")," ",paste0(...,collapse=""))
need <- function(p) {
  z<-p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)]
  if(length(z))stop("缺少R包：",paste(z,collapse=", "),"。全细胞评分需要这些包；内皮层分图已保存时不会丢失。")
}

canon <- function(x) tolower(gsub("[[:space:]_.-]+","",trimws(as.character(x))))
normalize_groups <- function(x) {
  aliases <- list(CK=c("CK","REF","CKref","Control","Mock","CTRL"),
                  Heat=c("Heat","High_Temp","highTemp","High_temperature","HT","HS"),
                  Cold=c("Cold","Low_Temp","lowTemp","Low_temperature","LT","CS"),PEG=c("PEG"),
                  ABA=c("ABA","Abscisic acid"),IAA=c("IAA","Indole-3-acetic acid"),
                  MeJA=c("MeJA","MEJA","Methyl jasmonate"),eBL=c("eBL","24-epibrassinolide","Epibrassinolide"),
                  SA=c("SA","Salicylic acid"))
  lk <- do.call(rbind,lapply(names(aliases),function(g)data.frame(key=unique(canon(aliases[[g]])),group=g)))
  if(anyDuplicated(lk$key))stop("条件别名冲突。")
  unname(lk$group[match(canon(x),lk$key)])
}
canon_gene <- function(x) {
  y <- toupper(trimws(as.character(x))); hit <- grepl("^GLYMA[._-]?[0-9]{2}G[0-9]{6}(\\.[0-9]+)*$",y)
  y[hit] <- sub("^GLYMA[._-]?([0-9]{2}G[0-9]{6}).*$","Glyma.\\1",y[hit])
  y[!hit] <- trimws(as.character(x[!hit])); y
}
read_counts <- function(obj,cells) {
  a <- methods::slot(obj,"assays")[["RNA"]]; if(is.null(a))stop("缺少RNA assay。")
  if(inherits(a,"Assay5")) {
    ls <- grep("^counts($|[.])",SeuratObject::Layers(a,search=NA),value=TRUE)
    if(!length(ls))stop("没有RNA counts。")
    if(length(ls)>1L)a <- SeuratObject::JoinLayers(a,layers="counts",new="counts")
    ls <- grep("^counts($|[.])",SeuratObject::Layers(a,search=NA),value=TRUE)
    if(length(ls)!=1L)stop("counts layers未唯一合并。")
    x <- SeuratObject::LayerData(a,layer=ls[1])
  } else x <- SeuratObject::GetAssayData(a,slot="counts")
  if(!all(cells %in% colnames(x)))stop("RNA counts缺少细胞。")
  x <- methods::as(x[,cells,drop=FALSE],"dgCMatrix")
  if(anyDuplicated(rownames(x))||is.null(rownames(x)))stop("基因名缺失/重复。")
  if(any(Matrix::colSums(x)<=0))stop("存在零文库细胞。")
  if(length(x@x))for(k in seq.int(1L,length(x@x),by=1000000L)) {
    v<-x@x[k:min(length(x@x),k+999999L)]
    if(any(!is.finite(v)|v<0|abs(v-round(v))>1e-6))stop("需要原始非负整数RNA counts。")
  }
  x
}

score_cells <- function(x,reg) {
  original<-names(reg);safe<-sprintf("R%05d",seq_along(reg));names(reg)<-safe
  scores<-matrix(NA_real_,ncol(x),length(reg),dimnames=list(colnames(x),original))
  for(start in seq.int(1L,ncol(x),by=SCORE_CHUNK)) {
    ix<-start:min(ncol(x),start+SCORE_CHUNK-1L)
    z<-as.matrix(UCell::ScoreSignatures_UCell(x[,ix,drop=FALSE],features=reg,
                                              maxRank=min(UCELL_MAX_RANK,nrow(x)),BPPARAM=BiocParallel::SerialParam()))
    cols<-paste0(safe,"_UCell")
    if(!all(cols %in% colnames(z))||!all(colnames(x)[ix] %in% rownames(z)))stop("UCell返回名称与输入不匹配。")
    scores[ix,]<-z[colnames(x)[ix],cols,drop=FALSE]
    msg("UCell已评分 ",max(ix),"/",ncol(x)," 个细胞。")
  }
  if(any(!is.finite(scores))||any(scores<0|scores>1))stop("UCell评分超出[0,1]或有缺失。")
  scores
}

sample_average <- function(scores,md,design) {
  if(!identical(rownames(scores),md$cell))stop("评分与细胞元数据未对齐。")
  meta<-do.call(rbind,lapply(seq_len(nrow(design)),function(i)
    data.frame(sample=design$sample[i],group=design$group[i],state=STATES,n_cells=0L)))
  mat<-matrix(NA_real_,nrow(meta),ncol(scores),dimnames=list(paste(meta$sample,meta$state,sep="|"),colnames(scores)))
  for(i in seq_len(nrow(meta))) {
    ix<-which(md$sample==meta$sample[i]&md$state==meta$state[i]);meta$n_cells[i]<-length(ix)
    if(length(ix)>=MIN_CELLS_PER_SAMPLE_STAGE)mat[i,]<-colMeans(scores[ix,,drop=FALSE])
  }
  meta$Status<-ifelse(meta$n_cells==0,"No_cells",ifelse(meta$n_cells<MIN_CELLS_PER_SAMPLE_STAGE,"Below_minimum_cells","Available"))
  list(mat=mat,meta=meta)
}
condition_average <- function(av) {
  meta<-do.call(rbind,lapply(GROUPS,function(g)data.frame(group=g,state=STATES)))
  mat<-matrix(NA_real_,nrow(meta),ncol(av$mat),dimnames=list(paste(meta$group,meta$state,sep="|"),colnames(av$mat)))
  meta$N_samples_available<-0L;meta$N_samples_total<-0L
  for(i in seq_len(nrow(meta))) {
    ix<-which(av$meta$group==meta$group[i]&av$meta$state==meta$state[i])
    use<-ix[av$meta$Status[ix]=="Available"]
    meta$N_samples_total[i]<-length(ix);meta$N_samples_available[i]<-length(use)
    if(length(use))mat[i,]<-colMeans(av$mat[use,,drop=FALSE])
  }
  meta$Status<-ifelse(meta$N_samples_available==0,"Missing",
                      ifelse(meta$N_samples_available<meta$N_samples_total,"Incomplete_samples",
                             ifelse(meta$N_samples_available<MIN_REPLICATES,"Single_sample_descriptive","Complete")))
  list(mat=mat,meta=meta)
}
make_delta <- function(cm) {
  z<-do.call(rbind,lapply(TREATMENTS,function(g) {
    a<-cm[paste(g,STATES,sep="|"),,drop=FALSE]-cm[paste("CK",STATES,sep="|"),,drop=FALSE]
    rownames(a)<-paste(g,STATES,sep="|");a
  }));z
}

text_width_mm <- function(labels,size,bold=FALSE) {
  max(vapply(labels,function(s)grid::convertWidth(grid::grobWidth(grid::textGrob(s,
                                                                                 gp=grid::gpar(fontsize=size,fontface=if(bold)'bold' else 'plain'))),'mm',valueOnly=TRUE),numeric(1)),0)
}
short_celltype <- function(x) {
  x[x=='Endodermis_Unassigned_or_other_lineage']<-'Endodermis (unassigned/other)'
  x
}
family_palette <- function(fams) {
  known<-c(AP2='#FF8A80',ARF='#56B4E9','ARR-B'='#C6AC00',bHLH='#FF80DF',bZIP='#00BFEF',
           C2H2='#57CA28',C3H='#C5A1FF','CO-like'='#67BCF0',Dof='#D88EF0',ERF='#36BDEB',
           'G2-like'='#00CBCD',GATA='#F18CE8',GRAS='#A0C800','HD-ZIP'='#90C500',HSF='#00C5DB',
           MYB='#00CB79',MYB_related='#C6C000','NF-YA'='#FFA14B',Trihelix='#00CEA0',Unknown='#999999')
  fams<-sort(unique(fams));extra<-setdiff(fams,names(known))
  if(length(extra))known<-c(known,setNames(grDevices::hcl.colors(length(extra),'Dark 3'),extra))
  known[fams]
}

layout_metrics <- function(z,tf_labels,cell_labels,fams) {
  # 在临时无文件图形设备上按字体实际宽度测量，避免TF标签/右侧细胞名被截断。
  grDevices::pdf(file=NULL,width=14,height=8)
  on.exit(grDevices::dev.off(),add=TRUE)
  # 垂直标签的横向占用取决于字体高度；字体放大时自动留足间距，避免相邻字母重叠。
  actual_cw<-max(CELL_WIDTH_MM,TF_FONT_PT*25.4/72*1.15)
  body_w<-ncol(z)*actual_cw
  body_h<-nrow(z)*CELL_HEIGHT_MM
  if(any(grepl('\n',cell_labels,fixed=TRUE)))body_h<-max(body_h,nrow(z)*5.1)
  right_names<-text_width_mm(cell_labels,CELLTYPE_FONT_PT)+5
  tf_bottom<-text_width_mm(tf_labels,TF_FONT_PT)+6
  family_w<-text_width_mm(c('Family',fams),LEGEND_FONT_PT,TRUE)+9
  top<-17+DENDRO_HEIGHT_MM+5
  legend_h<-7+length(fams)*4.8
  list(left=4,top=top,body_w=body_w,body_h=body_h,right_names=right_names,
       legend_w=family_w,tf_bottom=tf_bottom,width=4+body_w+right_names+14+family_w+5,
       height=top+max(body_h+tf_bottom,legend_h)+5)
}

render_dendrogram <- function(hc,w,h) {
  if(is.null(hc))return(invisible(NULL))
  n<-length(hc$order);positions<-numeric(n);positions[hc$order]<-seq_len(n)-.5
  heights<-hc$height;maxh<-max(heights);if(!is.finite(maxh)||maxh<=0)maxh<-1
  center<-numeric(n-1L)
  for(i in seq_len(n-1L)) {
    nodes<-hc$merge[i,]
    xx<-vapply(nodes,function(k)if(k<0)positions[-k] else center[k],numeric(1))
    yy<-vapply(nodes,function(k)if(k<0)0 else heights[k],numeric(1))
    center[i]<-mean(xx)
    grid::grid.segments(x0=grid::unit(xx/n*w,'mm'),x1=grid::unit(xx/n*w,'mm'),
                        y0=grid::unit(yy/maxh*h,'mm'),y1=grid::unit(rep(heights[i]/maxh*h,2),'mm'),gp=grid::gpar(lwd=.65))
    grid::grid.segments(x0=grid::unit(xx[1]/n*w,'mm'),x1=grid::unit(xx[2]/n*w,'mm'),
                        y0=grid::unit(heights[i]/maxh*h,'mm'),y1=grid::unit(heights[i]/maxh*h,'mm'),gp=grid::gpar(lwd=.65))
  }
}

render_panel <- function(z,labels,cell_labels,families,colors,tree,co,geom,bounds,subtitle,gaps=integer()) {
  grid::grid.newpage()
  # 使用mm原点定位，保证导出的PDF/PNG中字体与色块尺寸一致。
  w<-geom$width;h<-geom$height;left<-geom$left;top<-geom$top
  bw<-geom$body_w;bh<-geom$body_h
  ytop<-h-top;nr<-nrow(z);nc<-ncol(z);cw<-bw/nc;ch<-bh/nr
  grid::grid.rect(gp=grid::gpar(fill='white',col=NA))
  grid::grid.text(co,x=grid::unit(w/2,'mm'),y=grid::unit(h-5,'mm'),
                  gp=grid::gpar(fontsize=TITLE_FONT_PT,fontface='bold'))
  grid::grid.text(subtitle,x=grid::unit(w/2,'mm'),y=grid::unit(h-11,'mm'),
                  gp=grid::gpar(fontsize=10.5))
  grid::pushViewport(grid::viewport(x=grid::unit(left,'mm'),y=grid::unit(ytop+5,'mm'),
                                    width=grid::unit(bw,'mm'),height=grid::unit(DENDRO_HEIGHT_MM,'mm'),just=c('left','bottom')))
  render_dendrogram(tree,bw,DENDRO_HEIGHT_MM-1)
  grid::popViewport()
  xx<-left+(seq_len(nc)-.5)*cw
  grid::grid.rect(x=grid::unit(xx,'mm'),y=grid::unit(ytop+2.5,'mm'),width=grid::unit(cw,'mm'),height=grid::unit(3,'mm'),
                  gp=grid::gpar(fill=unname(colors[families]),col=NA))
  grid::grid.text('Family',x=grid::unit(left+bw+1.5,'mm'),y=grid::unit(ytop+2.5,'mm'),just='left',gp=grid::gpar(fontsize=10,fontface='bold'))
  palette<-grDevices::colorRampPalette(c('#2166AC','#67A9CF','#FFFFBF','#FDAE61','#B2182B'))(100)
  v<-as.vector(z);lo<-bounds[1];hi<-bounds[2]
  bins<-pmin(100,pmax(1,floor((pmax(lo,pmin(hi,v))-lo)/(hi-lo)*100)+1))
  fill<-rep('#D9D9D9',length(v));ok<-is.finite(v);fill[ok]<-palette[bins[ok]]
  grid::grid.rect(x=grid::unit(rep(xx,each=nr),'mm'),y=grid::unit(rep(ytop-(seq_len(nr)-.5)*ch,nc),'mm'),
                  width=grid::unit(cw,'mm'),height=grid::unit(ch,'mm'),gp=grid::gpar(fill=fill,col=NA))
  if(length(gaps))grid::grid.segments(x0=grid::unit(left,'mm'),x1=grid::unit(left+bw,'mm'),
                                      y0=grid::unit(ytop-gaps*ch,'mm'),y1=grid::unit(ytop-gaps*ch,'mm'),gp=grid::gpar(col='white',lwd=2))
  grid::grid.text(cell_labels,x=grid::unit(left+bw+1.5,'mm'),y=grid::unit(ytop-(seq_len(nr)-.5)*ch,'mm'),
                  just='left',gp=grid::gpar(fontsize=CELLTYPE_FONT_PT,lineheight=.95))
  grid::grid.text(labels,x=grid::unit(xx,'mm'),y=grid::unit(ytop-bh-1.5,'mm'),rot=TF_LABEL_ANGLE,
                  just=c('right','centre'),gp=grid::gpar(fontsize=TF_FONT_PT))
  # 色标含义由当前图的标题、subtitle与bounds共同指定。
  bx<-left+bw+geom$right_names+2;barh<-min(bh,48)
  grid::grid.rect(x=grid::unit(bx,'mm'),y=grid::unit(ytop-(seq_len(100)-.5)/100*barh,'mm'),
                  width=grid::unit(3.2,'mm'),height=grid::unit(barh/100+.02,'mm'),gp=grid::gpar(fill=rev(palette),col=NA))
  ticks<-seq(lo,hi,length.out=5)
  grid::grid.text(formatC(ticks,format='f',digits=2),x=grid::unit(bx+2.5,'mm'),y=grid::unit(ytop-(hi-ticks)/(hi-lo)*barh,'mm'),
                  just='left',gp=grid::gpar(fontsize=LEGEND_FONT_PT))
  lx<-bx+14
  grid::grid.text('Family',x=grid::unit(lx-1.5,'mm'),y=grid::unit(ytop+2.5,'mm'),just='left',gp=grid::gpar(fontsize=LEGEND_FONT_PT,fontface='bold'))
  fn<-names(colors);yy<-ytop-(seq_along(fn)-.5)*4.8
  grid::grid.rect(x=grid::unit(lx,'mm'),y=grid::unit(yy,'mm'),width=grid::unit(3.5,'mm'),height=grid::unit(4.2,'mm'),gp=grid::gpar(fill=unname(colors),col=NA))
  grid::grid.text(fn,x=grid::unit(lx+3,'mm'),y=grid::unit(yy,'mm'),just='left',gp=grid::gpar(fontsize=LEGEND_FONT_PT))
}

# 每列一个TF，对当前全部处理×三阶段的可用行统一标准化；不填补缺失值。
column_zscore <- function(m) {
  z<-matrix(NA_real_,nrow(m),ncol(m),dimnames=dimnames(m));qc<-list()
  for(j in seq_len(ncol(m))) {
    ok<-is.finite(m[,j]);n<-sum(ok);mu<-if(n)mean(m[ok,j]) else NA_real_
    s<-if(n>=2)stats::sd(m[ok,j]) else NA_real_
    state<-if(n<2)'Insufficient_observations_Z_missing' else if(!is.finite(s)||s<=1e-12)'Constant_delta_Z_zero' else 'Scaled'
    if(n>=2)z[ok,j]<-if(state=='Scaled')(m[ok,j]-mu)/s else 0
    qc[[j]]<-data.frame(TF=colnames(m)[j],N_observed=n,Mean_delta=mu,SD_delta=s,Status=state)
  }
  list(z=z,qc=do.call(rbind,qc))
}

write_matrix <- function(m,file)tabout(data.frame(Group=rownames(m),m,check.names=FALSE),file)
select_tfs <- function(delta,dict) {
  magnitude<-apply(delta,2,function(x)if(any(is.finite(x)))max(abs(x),na.rm=TRUE) else -Inf)
  ranked<-names(magnitude)[order(-magnitude,names(magnitude))]
  ranked<-ranked[is.finite(magnitude[ranked])]
  if(length(setdiff(FOCUS_TFS,dict$TF)))stop("FOCUS_TFS包含未评分TF：",paste(setdiff(FOCUS_TFS,dict$TF),collapse=","))
  sel<-unique(c(FOCUS_TFS,head(ranked,TOP_HEATMAP_TFS)))
  if(!length(sel))stop("没有可展示的阶段匹配比较。")
  tabout(data.frame(TF=names(magnitude),Max_absolute_Delta=as.numeric(magnitude),Displayed=names(magnitude) %in% sel),"TF_display_selection.csv")
  sel
}
plot_heatmap <- function(z,dict,stem,title,subtitle,bounds,gaps=integer(),tree=NULL) {
  if(!nrow(z)||!ncol(z))stop("热图矩阵为空。")
  d<-dict[match(colnames(z),dict$TF),,drop=FALSE]
  if(anyNA(d$Label))stop("热图TF与字典不匹配。")
  ord<-if(is.null(tree))seq_len(ncol(z)) else tree$order
  z<-z[,ord,drop=FALSE];d<-d[ord,,drop=FALSE]
  labels<-gsub("|"," | ",rownames(z),fixed=TRUE)
  colors<-family_palette(d$Family);geom<-layout_metrics(z,d$Label,labels,names(colors))
  render<-function()render_panel(z,d$Label,labels,d$Family,colors,tree,title,geom,bounds,subtitle,gaps)
  grDevices::pdf(out(paste0(stem,".pdf")),width=geom$width/25.4,height=geom$height/25.4,useDingbats=FALSE)
  tryCatch(render(),finally=grDevices::dev.off())
  grDevices::png(out(paste0(stem,".png")),width=geom$width/25.4,height=geom$height/25.4,units="in",res=PNG_DPI)
  tryCatch(render(),finally=grDevices::dev.off())
  write_matrix(z,paste0(stem,"_values.csv"))
  FIGURES_WRITTEN <<- unique(c(FIGURES_WRITTEN,out(paste0(stem,c(".pdf",".png")))))
}
make_tree <- function(z) {
  if(ncol(z)<2L)return(NULL)
  # 不填补缺失进行聚类；只用所有显示TF都有数据的行。无充分共同数据时保持原列顺序。
  a<-z[rowSums(is.finite(z))==ncol(z),,drop=FALSE]
  if(nrow(a)<2L)return(NULL)
  stats::hclust(stats::dist(t(a)),method="average")
}

resolve_output_directory <- function() {
  valid <- function(d)is.character(d)&&length(d)==1L&&!is.na(d)&&file.exists(file.path(d,paste0(SOURCE_PREFIX,"scores_and_stages.rds")))
  candidates<-unique(c(PARENT_DIR,Sys.glob(file.path(PARENT_DIR,"TF_activity_rerun_*"))))
  candidates<-candidates[vapply(candidates,valid,logical(1))]
  if(!is.null(SOURCE_RUN_DIR)) {
    if(!valid(SOURCE_RUN_DIR))stop("指定的SOURCE_RUN_DIR中没有已完成评分文件：",SOURCE_RUN_DIR)
    chosen<-SOURCE_RUN_DIR;reason<-"Explicit SOURCE_RUN_DIR"
  } else if(valid(SESSION_RUN_DIR)) {
    chosen<-SESSION_RUN_DIR;reason<-"Existing session OUT_DIR containing completed scores"
  } else {
    if(!length(candidates))stop("未找到 ",SOURCE_PREFIX,"scores_and_stages.rds。请将SOURCE_RUN_DIR设为截图对应的已有结果目录；不需要重跑GENIE3。")
    complete<-file.exists(file.path(candidates,paste0(SOURCE_PREFIX,"SUCCESS.txt")))
    pool<-if(any(complete))candidates[complete] else candidates
    # 同级时间戳目录中优先最近成功结果；无SUCCESS才用最近已保存评分。
    stamp<-file.info(file.path(pool,paste0(SOURCE_PREFIX,if(any(complete))"SUCCESS.txt" else "scores_and_stages.rds")))$mtime
    chosen<-pool[order(stamp,decreasing=TRUE,na.last=TRUE)[1]]
    reason<-if(any(complete))"Most recently successful existing run" else "Most recent saved scores; no old SUCCESS found"
  }
  chosen<-normalizePath(chosen,winslash="/",mustWork=TRUE)
  candidates<-unique(c(candidates,chosen));candidates<-unique(normalizePath(candidates,winslash="/",mustWork=TRUE))
  audit<-data.frame(Directory=candidates,Previous_success=file.exists(file.path(candidates,paste0(SOURCE_PREFIX,"SUCCESS.txt"))),
                    Selected=candidates==chosen,Selection_rule=reason)
  list(directory=chosen,audit=audit)
}
validate_dictionary <- function(d,ids) {
  if(!is.data.frame(d)||!all(c("TF","Label","Family") %in% names(d))||anyNA(d$TF)||anyDuplicated(d$TF)||!all(ids %in% d$TF))
    stop("原TF名称字典缺失或与评分列不匹配。")
  d<-d[match(ids,d$TF),,drop=FALSE]
  if(anyNA(d$Label)||any(!nzchar(d$Label))||anyDuplicated(d$Label))stop("TF显示标签缺失/重复。")
  if(any(grepl("Glyma[._-]?[0-9]{2}G[0-9]+",d$Label,ignore.case=TRUE)))
    stop("字典Label仍是原始Glyma ID；请使用你P1图中的Family [TF编号]字典。")
  d$Family[is.na(d$Family)|!nzchar(d$Family)]<-"Unknown";d
}
load_completed_endodermis <- function() {
  f<-srcpath("scores_and_stages.rds");b<-readRDS(f)
  if(!is.list(b)||!all(c("scores","metadata","regulons","TF_names","maxRank","UCell_version") %in% names(b)))
    stop("评分RDS不符合你上传的脚本结构：",f)
  m<-b$metadata;s<-b$scores
  if(!is.data.frame(m)||!all(c("cell","sample","group","state") %in% names(m))||anyNA(m$cell)||anyDuplicated(m$cell))stop("保存的阶段元数据不完整。")
  if(!is.matrix(s)||!is.numeric(s)||anyDuplicated(colnames(s))||anyDuplicated(rownames(s))||!setequal(m$cell,rownames(s)))stop("保存的细胞评分矩阵异常。")
  s<-s[m$cell,,drop=FALSE]
  if(any(!is.finite(s))||any(s<0|s>1))stop("保存的UCell分数不是[0,1]有限值。")
  if(anyNA(m$state)||!all(m$state %in% STATES)||anyNA(m$group)||!all(m$group %in% GROUPS))
    stop("需要这次9条件、P1/P2/P3的已完成结果；不混入旧分期。")
  if(anyNA(m$sample)||any(!nzchar(m$sample)))stop("保存样本ID缺失。")
  design<-unique(m[c("sample","group")]);if(anyDuplicated(design$sample))stop("保存的样本ID跨条件重复。")
  # 保存评分的m可能不含完全未分配的样本；用原sample_design补回并保留空阶段。
  if(file.exists(srcpath("sample_design.csv"))) {
    ds<-read.csv(srcpath("sample_design.csv"),check.names=FALSE,stringsAsFactors=FALSE)
    if(!all(c("sample","group") %in% names(ds))||anyDuplicated(ds$sample)||anyNA(ds$sample)||anyNA(ds$group)||!all(ds$group %in% GROUPS))stop("原sample_design.csv结构错误。")
    ix<-match(design$sample,ds$sample)
    if(anyNA(ix)||any(design$group!=ds$group[ix]))stop("保存样本设计与评分不一致。")
    design<-ds[c("sample","group")]
  }
  if(!all(GROUPS %in% design$group))stop("保存样本设计缺少本次请求条件。")
  cfg<-list()
  if(file.exists(srcpath("configuration.rds"))) {
    cfg<-readRDS(srcpath("configuration.rds"))
    if(!is.null(cfg$INPUT_RDS))INPUT_RDS<<-cfg$INPUT_RDS
    if(!is.null(cfg$SAMPLE_COL))SAMPLE_COL<<-cfg$SAMPLE_COL
    if(!is.null(cfg$MIN_CELLS_PER_SAMPLE_STAGE))MIN_CELLS_PER_SAMPLE_STAGE<<-as.integer(cfg$MIN_CELLS_PER_SAMPLE_STAGE)
    if(!is.null(cfg$GROUPS)&&!setequal(cfg$GROUPS,GROUPS))stop("原分析条件与本次8处理+CK不一致。")
  }
  if(length(MIN_CELLS_PER_SAMPLE_STAGE)!=1||!is.finite(MIN_CELLS_PER_SAMPLE_STAGE)||MIN_CELLS_PER_SAMPLE_STAGE<1)stop("保存的最少细胞数参数无效。")
  if(length(b$maxRank)!=1||!is.finite(b$maxRank)||b$maxRank<1)stop("保存的UCell maxRank无效。")
  UCELL_MAX_RANK<<-as.integer(b$maxRank)
  av<-sample_average(s,m,design);cm<-condition_average(av);delta<-make_delta(cm$mat)
  if(!is.null(b$Delta)) {
    if(!identical(dimnames(delta),dimnames(b$Delta))||!isTRUE(all.equal(delta,b$Delta,tolerance=1e-10)))
      stop("从保存评分重算的差值无法复现原图；请核对MIN_CELLS或原配置，不静默改动。")
  }
  dict<-validate_dictionary(b$TF_names,colnames(s))
  # 精确采用你P1图已经展示的TF，不另挑一批产生选择差异。
  plotfile<-srcpath("Figure03_vs_CK_raw_Delta_values.csv")
  if(file.exists(plotfile)) {
    a<-read.csv(plotfile,check.names=FALSE,stringsAsFactors=FALSE)
    if(!"Group" %in% names(a)||anyDuplicated(a$Group))stop("原图导出表Group异常。")
    selected<-setdiff(names(a),"Group")
    if(!length(selected)||!all(selected %in% colnames(delta))||!setequal(a$Group,rownames(delta)))stop("原图的TF/比较行与保存评分不一致。")
    aa<-as.matrix(a[,selected,drop=FALSE])
    if(!is.numeric(aa)||!isTRUE(all.equal(unname(aa),unname(delta[a$Group,selected,drop=FALSE]),tolerance=1e-8)))stop("原图导出值与原始差值不一致。")
  } else selected<-select_tfs(delta,dict)
  if(!is.list(b$regulons)||!all(selected %in% names(b$regulons)))stop("保存regulon缺少所选TF。")
  reg<-b$regulons[selected]
  if(any(vapply(reg,function(v)!is.character(v)||anyNA(v)||!length(v)||any(!nzchar(v)),logical(1))))stop("regulon靶基因集异常。")
  # 后续只评分本次展示的TF，节省全部细胞计算量；不是重训网络。
  reg<-lapply(reg,unique)
  z<-column_zscore(delta[,selected,drop=FALSE]);tree<-make_tree(z$z)
  tabout(dict,"TF_display_dictionary.csv")
  write_matrix(delta[,selected,drop=FALSE],"Endodermis_all8_raw_Delta.csv")
  tabout(cm$meta,"Endodermis_availability.csv")
  tabout(data.frame(TF=selected,Label=dict$Label[match(selected,dict$TF)]),"Displayed_TFs.csv")
  list(bundle=b,metadata=m,design=design,scores=s[,selected,drop=FALSE],sample=av,condition=cm,
       delta=delta[,selected,drop=FALSE],selected=selected,dict=dict,reg=reg,tree=tree,zdelta=z$z)
}
plot_endodermis_split <- function(e) {
  lim<-max(abs(e$delta),na.rm=TRUE);if(!is.finite(lim)||lim<1e-8)lim<-.01
  sets<-list(Stress=STRESS,Hormones=HORMONES)
  for(nm in names(sets)) {
    g<-sets[[nm]];rows<-unlist(lapply(g,function(k)paste(k,STATES,sep="|")))
    a<-e$delta[rows,,drop=FALSE];z<-e$zdelta[rows,,drop=FALSE]
    rownames(a)<-rownames(z)<-sub("[|](P[123])$","|Endodermis_\\1",rows)
    plot_heatmap(a,e$dict,paste0("Figure01_Endodermis_",nm,"_vs_CK_raw_Delta"),
                 paste(paste(g,collapse=" / "),"vs CK | Endodermis P1/P2/P3"),
                 "Raw Delta UCell; same TFs, order and color scale across the two panels",c(-lim,lim),seq(3,nrow(a)-1,3),e$tree)
    plot_heatmap(z,e$dict,paste0("Figure02_Endodermis_",nm,"_vs_CK_Delta_Z"),
                 paste(paste(g,collapse=" / "),"vs CK | Endodermis"),
                 "Column Z of Delta; jointly scaled across all 8 treatments x 3 stages",c(-3,3),seq(3,nrow(z)-1,3),e$tree)
  }
  msg("内皮层原图已拆分：三种胁迫一张，五种激素一张；原始差值不变。")
}

canonical_celltype <- function(x) {
  x<-trimws(as.character(x));key<-canon(x)
  standards<-c("Root Cap","Epidermis","Cortex","Dividing Cell","Meristem","QC","Pericycle","Procambium","Phloem","Xylem")
  j<-match(key,canon(standards));hit<-!is.na(j);x[hit]<-standards[j[hit]]
  x[key %in% c("lateralroot","procambium")]<-"Procambium"
  x
}
make_alltype_metadata <- function(raw,e) {
  if(!all(c(CELLTYPE_COL,TREATMENT_COL,SAMPLE_COL) %in% names(raw)))stop("原始RDS缺少细胞类型/处理/样本列。")
  if(is.null(rownames(raw))||anyDuplicated(rownames(raw)))stop("原始RDS细胞条码缺失/重复。")
  co<-normalize_groups(raw[[TREATMENT_COL]]);ct<-canonical_celltype(raw[[CELLTYPE_COL]])
  selected<-!is.na(co)&co %in% GROUPS
  m<-data.frame(cell=rownames(raw)[selected],sample=as.character(raw[[SAMPLE_COL]][selected]),
                group=co[selected],Celltype=ct[selected],stringsAsFactors=FALSE)
  if(anyNA(m$sample)||any(!nzchar(trimws(m$sample))))stop("原RDS生物样本ID缺失。")
  di<-unique(m[c("sample","group")]);if(anyDuplicated(di$sample))stop("同一生物样本ID跨处理。")
  if(!all(GROUPS %in% di$group))stop("原RDS缺少某些请求处理；不能只补出一部分却宣称全部完成。")
  ix<-match(e$metadata$cell,m$cell)
  if(anyNA(ix)||any(e$metadata$sample!=m$sample[ix])||any(e$metadata$group!=m$group[ix]))
    stop("原RDS细胞/样本/处理与本次已完成内皮层评分不一致；不能混合不同输入。")
  key<-canon(m$Celltype);isendo<-!is.na(key)&grepl("^endodermis|内皮层",key)
  if(any(!isendo[ix]))stop("已保存内皮层细胞在当前RDS被注释为其他类型；检查CELLTYPE_COL。")
  mapping<-match(m$cell,e$metadata$cell)
  m$Stage_source<-ifelse(isendo,"No_accepted_stage_in_saved_scores","Original_celltype")
  ok<-isendo&!is.na(mapping)
  m$Celltype[ok]<-paste0("Endodermis_",e$metadata$state[mapping[ok]])
  m$Stage_source[ok]<-"Exact_saved_P1_P2_P3"
  m$Included<-!is.na(m$Celltype)&nzchar(m$Celltype)&(!isendo|ok)
  m$Exclusion_reason<-ifelse(m$Included,"",ifelse(isendo,"No_accepted_saved_P1_P2_P3","Missing_celltype"))
  # 所有未分配内皮层保持审计，不画一个伪造的第四阶段。
  audit<-m[!m$Included,,drop=FALSE]
  m<-m[m$Included,,drop=FALSE]
  present<-unique(m$Celltype)
  types<-c(intersect(TYPE_ORDER,present),sort(setdiff(present,TYPE_ORDER)))
  if(!all(ENDO_TYPES %in% types))stop("缺少已保存P1/P2/P3之一。")
  di<-di[order(match(di$group,GROUPS),di$sample),,drop=FALSE]
  list(md=m,design=di,types=types,excluded=audit)
}
aggregate_one_sample <- function(scores,md,sample,condition,types) {
  if(!identical(as.character(rownames(scores)),as.character(md$cell))||any(md$sample!=sample)||any(md$group!=condition))stop("聚合时样本/细胞顺序不一致。")
  meta<-data.frame(sample=sample,group=condition,Celltype=types,n_cells=0L)
  mat<-matrix(NA_real_,length(types),ncol(scores),dimnames=list(paste(sample,types,sep="|"),colnames(scores)))
  for(i in seq_along(types)) {
    ix<-which(md$Celltype==types[i]);meta$n_cells[i]<-length(ix)
    if(length(ix)>=MIN_CELLS_PER_SAMPLE_STAGE)mat[i,]<-colMeans(scores[ix,,drop=FALSE])
  }
  meta$Status<-ifelse(meta$n_cells==0,"No_cells",ifelse(meta$n_cells<MIN_CELLS_PER_SAMPLE_STAGE,"Below_minimum_cells","Available"))
  list(meta=meta,mat=mat)
}
alltype_condition_means <- function(av,types) {
  meta<-do.call(rbind,lapply(GROUPS,function(g)data.frame(group=g,Celltype=types)))
  mat<-matrix(NA_real_,nrow(meta),ncol(av$mat),dimnames=list(paste(meta$group,meta$Celltype,sep="|"),colnames(av$mat)))
  meta$N_samples_total<-0L;meta$N_samples_available<-0L
  for(i in seq_len(nrow(meta))) {
    ix<-which(av$meta$group==meta$group[i]&av$meta$Celltype==meta$Celltype[i])
    use<-ix[av$meta$Status[ix]=="Available"]
    meta$N_samples_total[i]<-length(ix);meta$N_samples_available[i]<-length(use)
    if(length(use))mat[i,]<-colMeans(av$mat[use,,drop=FALSE])
  }
  meta$Status<-ifelse(meta$N_samples_available==0,"Missing",
                      ifelse(meta$N_samples_available<meta$N_samples_total,"Incomplete_samples",
                             ifelse(meta$N_samples_available<MIN_REPLICATES,"Single_sample_descriptive","Complete")))
  list(mat=mat,meta=meta)
}
alltype_deltas <- function(cm,types) {
  do.call(rbind,lapply(TREATMENTS,function(g) {
    a<-cm[paste(g,types,sep="|"),,drop=FALSE]-cm[paste("CK",types,sep="|"),,drop=FALSE]
    rownames(a)<-paste(g,types,sep="|");a
  }))
}
fingerprint <- function(x) {
  f<-tempfile(fileext=".rds");on.exit(unlink(f),add=TRUE);saveRDS(x,f,version=2,compress=FALSE)
  unname(tools::md5sum(f))
}
get_alltype_summary <- function(e) {
  need(c("SeuratObject","Matrix","UCell","BiocParallel"))
  if(!file.exists(INPUT_RDS))stop("找不到原始RDS：",INPUT_RDS)
  uv<-as.character(utils::packageVersion("UCell"))
  if(!identical(uv,as.character(e$bundle$UCell_version)))
    stop("UCell版本与已完成内皮层评分不同（当前",uv,"，保存",e$bundle$UCell_version,
         "）。内皮层分图已保存；全细胞图请在原评分环境运行，避免混合不同版本的分数。")
  msg("核对原始RDS与评分文件指纹；随后只补做这批TF的全细胞评分。")
  key<-fingerprint(list(version="split_alltypes_v2_1",input=normalizePath(INPUT_RDS,winslash="/"),
                        input_md5=unname(tools::md5sum(INPUT_RDS)),bundle_md5=unname(tools::md5sum(srcpath("scores_and_stages.rds"))),
                        regulons=e$reg,selected=e$selected,min_cells=MIN_CELLS_PER_SAMPLE_STAGE,maxRank=UCELL_MAX_RANK,UCell=uv,
                        columns=c(CELLTYPE_COL,TREATMENT_COL,SAMPLE_COL),groups=GROUPS,type_order=TYPE_ORDER))
  cache<-out("AllCelltypes_summary_cache.rds")
  if(REUSE_ALLTYPE_SUMMARY&&file.exists(cache)) {
    ans<-readRDS(cache)
    if(is.list(ans)&&identical(ans$key,key)&&all(c("av","cm","types","qc") %in% names(ans))) {
      msg("复用完全匹配的全细胞评分汇总缓存。");return(ans)
    }
    msg("缓存与当前输入不一致，将重新补充全细胞评分。")
  }
  obj<-readRDS(INPUT_RDS);if(!inherits(obj,"Seurat"))stop("输入RDS不是Seurat对象。")
  mm<-make_alltype_metadata(methods::slot(obj,"meta.data"),e)
  tabout(mm$md,"AllCelltypes_cell_metadata.csv");tabout(mm$excluded,"AllCelltypes_excluded_Endodermis_and_unlabelled.csv")
  tabout(mm$design,"AllCelltypes_sample_design.csv")
  x<-read_counts(obj,mm$md$cell);rm(obj);invisible(gc())
  ids<-canon_gene(rownames(x))
  if(anyNA(ids)||anyDuplicated(ids)||any(!nzchar(ids)))stop("原RDS基因名规范化后冲突；不能改变表达背景后继续。")
  rownames(x)<-ids
  if(UCELL_MAX_RANK>nrow(x))stop("原RNA基因背景小于保存的maxRank。")
  missing<-lapply(e$reg,setdiff,y=rownames(x))
  if(any(lengths(missing)>0))stop("原RDS缺少之前实际评分使用的靶基因；本次不悄悄改变regulon。")
  avs<-list();qcs<-list()
  for(i in seq_len(nrow(mm$design))) {
    sample<-mm$design$sample[i];condition<-mm$design$group[i]
    md<-mm$md[mm$md$sample==sample,,drop=FALSE]
    if(nrow(md)) {
      msg("全细胞评分：",sample,"，",nrow(md),"个细胞；",length(e$selected),"个候选regulon。")
      sc<-score_cells(x[,md$cell,drop=FALSE],e$reg)
    } else sc<-matrix(numeric(),nrow=0,ncol=length(e$selected),dimnames=list(character(),e$selected))
    # 检查新评分能复现同一批内皮层细胞，防止模型/基因背景/软件差异被忽略。
    overlap<-intersect(md$cell,rownames(e$scores))
    error<-if(length(overlap))max(abs(sc[overlap,e$selected,drop=FALSE]-e$scores[overlap,e$selected,drop=FALSE])) else NA_real_
    qcs[[i]]<-data.frame(Sample=sample,Condition=condition,Checked_Endodermis_cells=length(overlap),Max_absolute_UCell_difference=error)
    if(is.finite(error)&&error>1e-6) {
      tabout(do.call(rbind,qcs),"Endodermis_rescoring_validation.csv")
      stop("补评分无法复现已保存内皮层UCell：",sample,"，误差=",signif(error,5),"。不把两套分数混为同一图。")
    }
    avs[[i]]<-aggregate_one_sample(sc,md,sample,condition,mm$types)
    rm(sc);invisible(gc())
  }
  av<-list(mat=do.call(rbind,lapply(avs,`[[`,"mat")),meta=do.call(rbind,lapply(avs,`[[`,"meta")))
  cm<-alltype_condition_means(av,mm$types)
  ans<-list(key=key,av=av,cm=cm,types=mm$types,qc=do.call(rbind,qcs))
  saveRDS(ans,cache);ans
}

plot_all_celltypes <- function(a,e) {
  cm<-a$cm$mat[,e$selected,drop=FALSE];types<-a$types
  delta<-alltype_deltas(cm,types);zz<-column_zscore(delta);tree<-make_tree(zz$z)
  lim<-max(abs(delta),na.rm=TRUE);if(!is.finite(lim)||lim<1e-8)lim<-.01
  tabout(a$qc,"Endodermis_rescoring_validation.csv")
  tabout(a$cm$meta,"AllCelltypes_condition_availability.csv")
  tabout(cbind(a$av$meta,as.data.frame(a$av$mat,check.names=FALSE)),"AllCelltypes_sample_mean_UCell.csv")
  write_matrix(cm,"AllCelltypes_condition_raw_UCell.csv")
  write_matrix(delta,"AllCelltypes_all8_vs_CK_raw_Delta.csv")
  write_matrix(zz$z,"AllCelltypes_all8_vs_CK_Delta_Z.csv");tabout(zz$qc,"AllCelltypes_Delta_Z_parameters.csv")
  ck<-cm[paste("CK",types,sep="|"),,drop=FALSE];rownames(ck)<-types
  ckz<-column_zscore(ck);tabout(ckz$qc,"CK_activity_Z_parameters.csv")
  write_matrix(ck,"CK_all_celltypes_raw_UCell.csv")
  plot_heatmap(ckz$z,e$dict,"Figure03_CK_all_celltypes_activity_Z",
               "Endodermis-trained candidate regulons | CK all cell types",
               "Column Z of CK raw UCell across cell types; equal-sample means; grey = unavailable",
               c(-3,3),tree=make_tree(ckz$z))
  sets<-list(Stress=STRESS,Hormones=HORMONES)
  for(nm in names(sets)) {
    g<-sets[[nm]];rows<-unlist(lapply(g,function(k)paste(k,types,sep="|")))
    gaps<-seq(length(types),length(rows)-1,by=length(types))
    plot_heatmap(delta[rows,,drop=FALSE],e$dict,paste0("Figure04_AllCelltypes_",nm,"_vs_CK_raw_Delta"),
                 paste(paste(g,collapse=" / "),"vs CK | all cell types"),
                 "Raw Delta UCell; matched cell identities; equal-sample means; same scale across both panels",c(-lim,lim),gaps,tree)
    plot_heatmap(zz$z[rows,,drop=FALSE],e$dict,paste0("Figure05_AllCelltypes_",nm,"_vs_CK_Delta_Z"),
                 paste(paste(g,collapse=" / "),"vs CK | all cell types"),
                 "Column Z of Delta across all 8 treatments x cell types; red need not mean positive Delta",c(-3,3),gaps,tree)
  }
  # 每种处理也单独输出一张：右边仅为celltype，便于像参考图一样排版。
  for(g in TREATMENTS) {
    z<-delta[paste(g,types,sep="|"),,drop=FALSE];rownames(z)<-types
    plot_heatmap(z,e$dict,paste0("Figure06_",g,"_vs_CK_celltypes_right_raw_Delta"),
                 paste(g,"vs CK | cell-type-matched candidate regulon response"),
                 "Raw Delta UCell; all cell types; Endodermis subdivided by exact saved P1/P2/P3",
                 c(-lim,lim),tree=tree)
  }
}
main <- function() {
  FIGURES_WRITTEN <<- character()
  chosen<-resolve_output_directory();OUT_DIR<<-chosen$directory
  msg("使用已有输出目录：",OUT_DIR,"；不新建时间戳目录。")
  tabout(chosen$audit,"source_directory_candidates.csv")
  for(nm in c("SUCCESS.txt","PARTIAL.txt"))if(file.exists(out(nm)))unlink(out(nm))
  e<-load_completed_endodermis()
  plot_endodermis_split(e)
  writeLines(c("Endodermis split panels finished using saved values; all-cell-type panels not yet completed.",
               paste("Output:",OUT_DIR)),out("PARTIAL.txt"))
  alltypes<-NULL
  if(RUN_ALL_CELLTYPES) {
    alltypes<-get_alltype_summary(e)
    plot_all_celltypes(alltypes,e)
  }
  figures<-FIGURES_WRITTEN
  tabout(data.frame(File=basename(figures),Path=figures),"figure_manifest.csv")
  utils::capture.output(sessionInfo(),file=out("sessionInfo.txt"))
  info<-c("Completed: saved Endodermis results split into Stress and Hormones panels.",
          paste("All-cell-type figures requested/completed:",RUN_ALL_CELLTYPES),
          paste("Output directory (unchanged):",OUT_DIR),paste("New prefix:",PREFIX),
          "Original TF gene sets, display labels and exact saved P1/P2/P3 are retained.",
          "Endodermis raw Delta values are unchanged; only plot layout is split.",
          "The two Endodermis panels share their TF order and raw Delta color scale.",
          "Endodermis Delta Z uses the combined 8-treatment x 3-stage reference, not separate per-panel scaling.",
          "All-cell-type scoring uses the same UCell version and gene sets; overlapping Endodermis cells are numerically validated.",
          "Only the displayed TFs are additionally scored in other cell types; this is not new all-cell-type network training.",
          "Non-assigned Endodermis is audited/excluded, not fabricated as a fourth stage or forced into P1/P2/P3.",
          "CK activity Z and treatment-minus-CK Delta Z are different quantities and have separate filenames.",
          "All-cell-type Delta Z uses all eight treatments x retained cell identities; its reference differs from Endodermis-only Z.",
          "Raw values, missingness and available biological replicate counts are exported. No new significance tests were performed.",
          "TF labels retain verified symbols or Family [TFxxxx]; internal keys are not official gene symbols.",
          "Original input RDS and old results are not overwritten. Re-running this follow-up updates only this new prefix.")
  if(file.exists(out("PARTIAL.txt")))unlink(out("PARTIAL.txt"))
  writeLines(info,out("SUCCESS.txt"))
  msg("全部完成。输出目录：",OUT_DIR,"；前缀：",PREFIX)
  print(data.frame(Figure=basename(figures)),row.names=FALSE)
  invisible(list(out_dir=OUT_DIR,figures=figures,alltypes=alltypes))
}
run_followup <- function() {
  tryCatch(main(),error=function(err) {
    if(!is.null(OUT_DIR)&&dir.exists(OUT_DIR))writeLines(c("PARTIAL / FAILED",conditionMessage(err),
                                                           "Existing completed input results were not modified. Any saved split panels remain available.",paste("Output:",OUT_DIR)),out("PARTIAL.txt"))
    stop(conditionMessage(err),call.=FALSE)
  })
}
if(!isTRUE(getOption("tf_layout_v2.skip_main")))run_followup()