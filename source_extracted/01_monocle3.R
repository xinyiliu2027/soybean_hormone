
# 在RStudio中可 source(..., encoding='UTF-8') 后运行 main_v8()

# -*- coding: UTF-8 -*-
# -*- coding: UTF-8 -*-
# 完整独立脚本 v7c：六基因原UMAP表达、六处理对比、拟时配色及重复表达动态。
# 用法：source('soybean_endodermis_sixgenes_originalUMAP_v7c_ptPlasma.txt', encoding='UTF-8')
#
# 相对 v7b 只改了拟时序配色这一件事，RUN_DIR、输出文件名、数值口径、绘图布局完全不变：
#  (A) PT_COLORS 换成 plasma 连续色带：0 深蓝紫 -> 紫 -> 品红 -> 橙红 -> 亮黄(上限)，
#      与参考图一致（参考图即 matplotlib plasma）。
#  (B) 拟时色阶默认改为 linear 连续映射（PT_COLOR_MODE='linear'），颜色随真实拟时值
#      连续变化，不再把颜色停靠点按分位拉成等间距；原分位版本仍额外输出一份
#      _pooled_quantile 对照文件（PT_ALSO_EMIT_QUANTILE <- FALSE 可完全不输出）。
#  (C) 色条刻度改为整齐圆整值（0/5/10/15/20），linear 模式按真实数值定位；分位模式
#      仍标注真实分位值（PT_LEGEND_NICE_TICKS <- FALSE 可关掉圆整）。
#
# 沿用 v7b 的三处修复：不加载 org.G.max.eg.db / AnnotationDbi；sessionInfo 安全写出；
# smooth.spline 自由度自适应且静默，不再出现 "smoothing parameter too large" 警告。
#
# 只加载已存CDS和CSV，不修改UMAP/轨迹/根/拟时序，不运行降维、graph_test或GO。
# PDF/PNG/CSV直接写回原RUN_DIR，文件名前缀由 OUT_PREFIX 控制（默认仍为 v7）。
RUN_DIR <- '/data01/liuxy/Project/CK2_added_v3_Results/Endodermis_monocle3_originalUMAP_grouped_v3/20261119_163425'
OUT_PREFIX <- 'v7'
TARGET_GENES <- c('Glyma.17G152800','Glyma.06G050000','Glyma.09G139000',
                  'Glyma.06G003800','Glyma.13G094400','Glyma.12G173300')
USER_LABELS <- c('ARR-B candidate','MYB84/RAX3-like','MYB candidate','MYB48-like','MYB36-like','RING-CH candidate')
# 上述短名仅来自用户提供注释，不作为本次验证过的正交同源或功能结论。
GROUPS <- c('CK','ABA','eBL','IAA','MeJA','SA')
GROUP_COLORS <- c(CK='#4D4D4D',ABA='#D55E00',eBL='#009E73',IAA='#0072B2',MeJA='#CC79A7',SA='#E69F00')
PREFERRED_DB_VERSION <- 'v4'
FEATURE_OVERRIDES <- NULL  # 可选命名向量：c('Glyma.17G152800'='实际矩阵行名')
PT_COLOR_MODE <- 'linear' # 主图：颜色随真实拟时值连续变化（参考图口径）；也可改回 'pooled_quantile'
PT_ALSO_EMIT_QUANTILE <- TRUE # 额外输出一份分位间距版本作对照；FALSE 则只出主图
PT_LEGEND_NICE_TICKS <- TRUE # 色条刻度用整齐圆整值（0/5/10/15/20）
PT_PROBS <- c(0,.10,.35,.65,.90,1)
# plasma 连续色带：0=深蓝紫 -> 紫 -> 品红 -> 橙红 -> 亮黄=上限，与参考图色条一致。
PT_COLORS <- c('#0D0887','#6A00A8','#B12A90','#E16462','#FCA636','#F0F921')
EXPRESSION_COLORS <- c('#352A87','#8157B7','#C44C8A','#EC8A55','#F3D55B')
EXPR_UPPER_QUANTILE <- .99 # 同一基因全部细胞正表达值的99%分位；高于阈值饱和但原数值保留
ZERO_COLOR <- '#D9D9D9'
POINT_CEX <- .45
DRAW_GRAPH_ON_GENE_MAP <- TRUE
SELECT_PATHS <- NULL
N_TREND_BINS <- 20L
MIN_CELLS_PER_REPLICATE_BIN <- 5L
SMOOTH_DF <- 3L
PNG_DPI <- 200L

parse_soybean_id <- function(ids) {
  raw<-as.character(ids);s<-trimws(raw)
  s<-sub('^(gene|transcript):','',s,ignore.case=TRUE)
  pat<-'^(Glyma)[._-]?([0-9]{2})[gG]([0-9]{6})(.*)$'
  valid<-!is.na(s)&grepl(pat,s,ignore.case=TRUE)
  base<-s;suffix<-rep('',length(s));ver<-rep(NA_character_,length(s))
  base[valid]<-sub(pat,'Glyma.\\2G\\3',s[valid],ignore.case=TRUE)
  suffix[valid]<-sub(pat,'\\4',s[valid],ignore.case=TRUE)
  vv<-regexpr('(^|[._-])v[0-9]+(?=$|[._-])',suffix,perl=TRUE,ignore.case=TRUE)
  has<-which(vv>0)
  if(length(has))ver[has]<-tolower(sub('^[._-]','',regmatches(suffix,vv)))
  tail<-gsub('(^|[._-])v[0-9]+(?=$|[._-])','',suffix,perl=TRUE,ignore.case=TRUE)
  safe<-grepl('^([._-]([0-9]+|[tTpP][0-9]+))*$',tail)
  data.frame(raw_id=raw,clean_id=ifelse(valid&safe,base,s),version=ver,recognized=valid&safe)
}
match_soybean_gid <- function(features,dbkeys,preferred='v4') {
  f<-parse_soybean_id(features);b<-parse_soybean_id(dbkeys);preferred<-tolower(preferred)
  b$priority<-ifelse(!is.na(b$version)&b$version==preferred,1L,ifelse(is.na(b$version),2L,3L))
  groups<-split(seq_len(nrow(b)),b$clean_id)
  out<-data.frame(feature_id=features,clean_id=f$clean_id,feature_version=f$version,
                  GID=NA_character_,db_version=NA_character_,match_status='unmatched',n_db_candidates=0L)
  for(i in seq_len(nrow(f))) {
    if(!is.na(f$version[i])&&f$version[i]!=preferred) {
      out$match_status[i]<-'source_nonpreferred_version_requires_mapping';next
    }
    ii<-groups[[f$clean_id[i]]]
    if(!length(ii))next
    out$n_db_candidates[i]<-length(ii)
    pr<-min(b$priority[ii]);best<-ii[b$priority[ii]==pr]
    if(pr==3L){out$match_status[i]<-'db_only_nonpreferred_version';next}
    exact<-best[b$raw_id[best]==f$raw_id[i]]
    if(length(exact)==1L)best<-exact
    if(length(best)!=1L){out$match_status[i]<-'ambiguous_same_priority';next}
    out$GID[i]<-b$raw_id[best];out$db_version[i]<-b$version[best]
    out$match_status[i]<-if(pr==1L)'preferred_v4' else 'unique_unversioned_GID'
  }
  out$matched<-!is.na(out$GID);out
}

read_table_v7 <- function(root,name,required=TRUE) {
  f<-file.path(root,name)
  if(!file.exists(f)) {
    if(required)stop('Missing saved file: ',f)
    return(NULL)
  }
  read.csv(f,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding='UTF-8',na.strings=c('','NA'))
}
save_v7 <- function(stem,fun,width=12,height=8) {
  pdf(paste0(stem,'.pdf'),width=width,height=height,useDingbats=FALSE)
  tryCatch(fun(),finally=dev.off())
  png(paste0(stem,'.png'),width=round(width*PNG_DPI),height=round(height*PNG_DPI),res=PNG_DPI,
      type=if(capabilities('cairo'))'cairo' else getOption('bitmapType','Xlib'))
  tryCatch(fun(),finally=dev.off())
}

# 安全包版本查询：DESCRIPTION 缺失/损坏的包返回 NA，绝不抛出
# "Error in x[['Version']] : subscript out of bounds"。
pkg_version_chr <- function(pkg) {
  tryCatch({
    v<-suppressWarnings(utils::packageVersion(pkg))
    if(is.null(v)||!length(v))NA_character_ else as.character(v)
  },error=function(e)NA_character_,warning=function(w)NA_character_)
}
# 安全加载检查：坏包只返回 FALSE，不会让脚本中止在 requireNamespace 内部。
pkg_available <- function(pkg) {
  isTRUE(suppressWarnings(suppressMessages(tryCatch(
    requireNamespace(pkg,quietly=TRUE),error=function(e)FALSE))))
}
write_session_info_v7 <- function(file,root) {
  pk<-sort(unique(c('base','stats','utils','graphics','grDevices','methods','monocle3',
                    'SummarizedExperiment','SingleCellExperiment',loadedNamespaces())))
  L<-c(paste0('R version: ',R.version.string),
       paste0('Platform: ',R.version$platform),
       paste0('Written at: ',format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z')),
       paste0('No OrgDb / AnnotationDbi package is loaded by this script.'),
       '',
       'Loaded namespaces and versions; unknown = DESCRIPTION unreadable (reported, never fatal):')
  for(p in pk) {
    v<-pkg_version_chr(p)
    L<-c(L,sprintf('  %-30s %s',p,if(is.na(v))'unknown' else v))
  }
  tryCatch(writeLines(L,file,useBytes=FALSE),error=function(e)NULL)
  # 完整 sessionInfo 只作附加信息，失败也不影响主流程。
  si<-suppressWarnings(tryCatch(utils::capture.output(utils::sessionInfo()),error=function(e)NULL))
  if(length(si))
    tryCatch(writeLines(si,file.path(root,paste0(OUT_PREFIX,'_sessionInfo_full.txt')),useBytes=FALSE),
             error=function(e)NULL)
  invisible(file)
}

make_pt_scale <- function(pt,mode=PT_COLOR_MODE) {
  v<-pt[is.finite(pt)]
  if(!length(v))stop('No finite pseudotime.')
  if(!mode %in% c('pooled_quantile','linear'))stop('Invalid PT_COLOR_MODE.')
  if(diff(range(v))==0)return(list(knots=v[1],colors=PT_COLORS[1],mode=mode))
  knots<-if(mode=='linear')seq(min(v),max(v),length.out=length(PT_COLORS)) else
    as.numeric(quantile(v,PT_PROBS,names=FALSE,type=7))
  knots<-unique(knots)
  # 若大量拟时并列，则合并重复断点；始终保持颜色映射单调。
  list(knots=knots,colors=colorRampPalette(PT_COLORS)(length(knots)),mode=mode)
}
color_by_scale <- function(values,scale,na_color='#C5C5C5') {
  out<-rep(na_color,length(values));ok<-is.finite(values)
  if(length(scale$knots)==1L) {out[ok]<-scale$colors[1];return(out)}
  u<-approx(scale$knots,seq(0,1,length.out=length(scale$knots)),xout=values[ok],rule=2)$y
  pal<-colorRampPalette(scale$colors)(256)
  out[ok]<-pal[pmin(256L,pmax(1L,1L+floor(255*u)))]
  out
}
make_expr_scale <- function(v) {
  p<-v[is.finite(v)&v>0]
  upper<-if(length(p))as.numeric(quantile(p,EXPR_UPPER_QUANTILE,names=FALSE)) else 1
  if(!is.finite(upper)||upper<=0)upper<-1
  list(knots=seq(0,upper,length.out=length(EXPRESSION_COLORS)),colors=EXPRESSION_COLORS,
       upper=upper,n_clipped=sum(v>upper,na.rm=TRUE),mode='linear')
}
color_legend_v7 <- function(scale,title) {
  par(mar=c(2,.4,3,3.2));plot.new();plot.window(xlim=c(0,1),ylim=c(0,1))
  title(main=title,cex.main=.82)
  if(length(scale$knots)<2) {
    rect(.2,.2,.45,.8,col=scale$colors[1],border=NA)
    text(.5,.5,format(scale$knots[1],digits=3),adj=0,cex=.8);return(invisible(NULL))
  }
  y<-seq(.15,.87,length.out=257)
  rect(.15,head(y,-1),.4,tail(y,-1),col=colorRampPalette(scale$colors)(256),border=NA)
  # 数值 -> 色条纵向位置：linear 模式即连续渐变，分位模式等价于原等间距断点。
  ypos<-function(v)approx(scale$knots,seq(.15,.87,length.out=length(scale$knots)),
                          xout=v,rule=2,ties='ordered')$y
  tv<-scale$knots
  if(isTRUE(PT_LEGEND_NICE_TICKS)&&identical(scale$mode,'linear')) {
    tt<-pretty(range(scale$knots),4)
    tt<-tt[tt>=min(scale$knots)-1e-9&tt<=max(scale$knots)+1e-9]
    if(length(tt)>=2L)tv<-tt
  }
  pp<-as.numeric(ypos(tv))
  segments(.4,pp,.47,pp)
  labels<-vapply(tv,function(z)format(signif(z,3),trim=TRUE),character(1))
  text(.51,pp,labels,adj=0,cex=.75)
  if(title=='Expression')text(.1,.04,'Grey = zero',adj=0,cex=.7)
}

draw_umap_v7 <- function(d,ix,values,scale,title,edges=NULL,nodes=NULL,root_nodes=NULL,
                         expression=FALSE,missing=FALSE) {
  par(mar=c(3.1,3.3,3,1),xaxs='r',yaxs='r')
  plot(d$x,d$y,col='#EEEEEE',pch=16,cex=POINT_CEX*.75,asp=1,
       xlab='Original UMAP 1',ylab='Original UMAP 2',main=title,cex.main=.95,cex.axis=.8,cex.lab=.8)
  if(!is.null(edges)&&nrow(edges))segments(edges$x,edges$y,edges$xend,edges$yend,col='#707070',lwd=.8)
  if(missing) {
    mtext('Gene missing / ambiguous; see mapping audit',side=3,line=.2,cex=.68,col='#A32929')
    return(invisible(NULL))
  }
  vv<-values[ix];cols<-color_by_scale(vv,scale)
  if(expression)cols[is.finite(vv)&vv==0]<-ZERO_COLOR
  # 低值先画，高值后画，避免高表达细胞被零表达覆盖；坐标不变。
  ord<-order(vv,na.last=FALSE)
  points(d$x[ix[ord]],d$y[ix[ord]],pch=16,col=cols[ord],cex=POINT_CEX)
  if(length(root_nodes)&&!is.null(nodes)) {
    rr<-nodes[nodes$node %in% root_nodes,,drop=FALSE]
    points(rr$x,rr$y,pch=23,bg='white',col='black',cex=1.1,lwd=1)
  }
}

# 仅在同一路径内求曲线；各处理/重复使用同一套分箱边界。
# 缺失箱为NA并断线；样条只在至少4个连续有效箱内求值，不跨缺失区域或向外推。
# df 至少比有效点数少2，且拟合全程 suppressWarnings，避免
#"smoothing parameter value too large / setting df = 1" 警告（不改变曲线口径）。
smooth_contiguous <- function(pt,value,valid,df=SMOOTH_DF) {
  out<-rep(NA_real_,length(pt))
  if(!any(valid))return(out)
  runs<-split(which(valid),cumsum(c(TRUE,diff(which(valid))!=1)))
  for(ii in runs) {
    out[ii]<-value[ii]
    if(length(ii)>=4L&&length(unique(pt[ii]))>=4L) {
      df_i<-max(1L,min(as.integer(df),length(ii)-2L))
      fit<-tryCatch(suppressWarnings(smooth.spline(pt[ii],value[ii],df=df_i)),error=function(e)NULL)
      if(!is.null(fit)) {
        pr<-tryCatch(suppressWarnings(predict(fit,pt[ii])$y),error=function(e)NULL)
        if(!is.null(pr)&&length(pr)==length(ii))out[ii]<-pmax(0,pr)
      }
    }
  }
  out
}
path_trends_v7 <- function(d,ex,ix,path) {
  pt<-d$local_pseudotime[ix]
  br<-unique(as.numeric(quantile(pt,seq(0,1,length.out=N_TREND_BINS+1),names=FALSE)))
  if(length(br)<5L)return(list(status='too_few_distinct_time_bins',bins=data.frame(),cor=data.frame()))
  bn<-cut(pt,br,include.lowest=TRUE,labels=FALSE)
  centers<-(head(br,-1)+tail(br,-1))/2;nb<-length(centers)
  rows<-list();corr<-list()
  for(tg in GROUPS) {
    samples<-unique(d$sample[d$treatment==tg])
    for(s in samples) {
      jj<-which(d$treatment[ix]==tg & d$sample[ix]==s)
      for(gene in rownames(ex)) {
        v<-ex[gene,ix[jj]];p<-pt[jj]
        ok<-is.finite(v)&is.finite(p)
        rho<-if(sum(ok)>=10&&length(unique(v[ok]))>1&&length(unique(p[ok]))>1)
          suppressWarnings(cor(p[ok],v[ok],method='spearman')) else NA_real_
        corr[[paste(tg,s,gene)]]<-data.frame(path=path,treatment=tg,sample=s,gene=gene,
                                             n_cells=sum(ok),detected_fraction=if(sum(ok))mean(v[ok]>0) else NA_real_,
                                             rho_spearman_descriptive=rho)
        means<-rep(NA_real_,nb);ns<-integer(nb)
        for(b in seq_len(nb)) {
          kk<-jj[bn[jj]==b];vv<-ex[gene,ix[kk]];vv<-vv[is.finite(vv)];ns[b]<-length(vv)
          if(ns[b]>=MIN_CELLS_PER_REPLICATE_BIN)means[b]<-mean(vv)
        }
        smooth<-smooth_contiguous(centers,means,is.finite(means))
        rows[[paste(tg,s,gene)]]<-data.frame(path=path,treatment=tg,sample=s,gene=gene,
                                             bin=seq_len(nb),left=head(br,-1),right=tail(br,-1),pseudotime=centers,
                                             n_cells=ns,mean_log1p=means,smooth_log1p=smooth)
      }
    }
  }
  list(status='completed',bins=do.call(rbind,rows),cor=do.call(rbind,corr))
}

render_analysis_v7 <- function(d,ex,edges,nodes,roots,paths,path_nodes,root) {
  csv<-function(x,nm)write.csv(x,file.path(root,paste0(OUT_PREFIX,'_',nm,'.csv')),row.names=FALSE,na='')
  settings<-list();expr_scales<-lapply(seq_len(nrow(ex)),function(i)make_expr_scale(ex[i,]))
  names(expr_scales)<-rownames(ex)
  # 1. 第二张图重新着色；每个分量各自定义色阶，六处理共用，绝不按每组拉伸。
  for(cc in unique(d$component)) {
    ids<-which(d$component==cc & is.finite(d$local_pseudotime))
    if(!length(ids))next
    ee<-edges[edges$component==cc,,drop=FALSE]
    rr<-if(!is.null(roots))roots$root_node[roots$component==cc] else character()
    for(mode in unique(c(PT_COLOR_MODE, if(isTRUE(PT_ALSO_EMIT_QUANTILE)) 'pooled_quantile'))) {
      sc<-make_pt_scale(d$local_pseudotime[ids],mode)
      settings[[paste(cc,mode)]]<-data.frame(component=cc,mode=mode,pseudotime=sc$knots,color=sc$colors)
      save_v7(file.path(root,paste0(OUT_PREFIX,'_P2_component_',cc,'_',mode)),function() {
        layout(matrix(c(1,2,3,7,4,5,6,7),2,4,byrow=TRUE),widths=c(1,1,1,.38))
        par(oma=c(2,0,3,0))
        for(tg in GROUPS) {
          ii<-ids[d$treatment[ids]==tg]
          draw_umap_v7(d,ii,d$local_pseudotime,sc,paste(tg,'| local PT'),ee,nodes,rr)
        }
        color_legend_v7(sc,'Local PT')
        mtext(paste('Endodermis | original UMAP | component',cc),3,outer=TRUE,line=1,font=2)
        mtext(if(mode=='pooled_quantile')'Quantile-spaced color stops pooled across treatments; legend = original PT values.' else
          'Linear color mapping; legend = original PT values.',1,outer=TRUE,line=.4,cex=.8)
      },14,9)
    }
  }
  csv(do.call(rbind,settings),'pseudotime_color_knots')
  # 2. 六基因总体图：各基因单独色阶，数值上限明确写出；同一基因所有组范围一致。
  save_v7(file.path(root,paste0(OUT_PREFIX,'_sixgenes_originalUMAP_overview')),function() {
    par(mfrow=c(2,3),oma=c(2,0,2,0))
    for(g in TARGET_GENES) {
      sc<-expr_scales[[g]];miss<-all(!is.finite(ex[g,]))
      draw_umap_v7(d,seq_len(nrow(d)),ex[g,],sc,g,
                   if(DRAW_GRAPH_ON_GENE_MAP)edges else NULL,expression=TRUE,missing=miss)
      if(!miss)mtext(paste0('Expression 0 to ',signif(sc$upper,3),' (positive q',EXPR_UPPER_QUANTILE,')'),3,line=.2,cex=.62)
    }
    mtext('Six candidate genes | original Endodermis UMAP',3,outer=TRUE,font=2)
    mtext('Grey = zero. Colored = log1p(count / saved size factor). Per-gene scale; see individual panels for legends.',1,outer=TRUE,cex=.75)
  },13,9)
  for(g in TARGET_GENES) {
    sc<-expr_scales[[g]];miss<-all(!is.finite(ex[g,]))
    save_v7(file.path(root,paste0(OUT_PREFIX,'_',g,'_six_treatments')),function() {
      layout(matrix(c(1,2,3,7,4,5,6,7),2,4,byrow=TRUE),widths=c(1,1,1,.38))
      par(oma=c(2,0,3,0))
      for(tg in GROUPS)draw_umap_v7(d,which(d$treatment==tg),ex[g,],sc,tg,
                                    if(DRAW_GRAPH_ON_GENE_MAP)edges else NULL,expression=TRUE,missing=miss)
      if(miss) {
        par(mar=c(1,1,2,1));plot.new();text(.5,.5,'No mapped\nexpression',cex=.8)
      } else color_legend_v7(sc,'Expression')
      mtext(paste(g,'| same expression scale across all six treatments'),3,outer=TRUE,line=1,font=2)
      mtext(paste0('log1p(count / size factor); upper=',signif(sc$upper,3),
                   '; values above cap=',sc$n_clipped,'. Grey: zero; pale background: reference cells.'),1,outer=TRUE,cex=.75)
    },14,9)
  }
  # 3. 保存每个生物学重复的表达/检测率，避免用单细胞数代替重复数。
  summary<-list()
  for(g in TARGET_GENES)for(tg in GROUPS)for(s in unique(d$sample[d$treatment==tg])) {
    v<-ex[g,d$treatment==tg & d$sample==s];v<-v[is.finite(v)]
    summary[[paste(g,tg,s)]]<-data.frame(gene=g,treatment=tg,sample=s,n_cells=length(v),
                                         mean_log1p=if(length(v))mean(v) else NA_real_,median_log1p=if(length(v))median(v) else NA_real_,
                                         detected_fraction=if(length(v))mean(v>0) else NA_real_)
  }
  csv(do.call(rbind,summary),'gene_expression_per_replicate')
  csv(do.call(rbind,lapply(TARGET_GENES,function(g)data.frame(gene=g,
                                                              upper=expr_scales[[g]]$upper,positive_quantile=EXPR_UPPER_QUANTILE,n_clipped=expr_scales[[g]]$n_clipped))),
      'expression_color_limits')
  # 4. 每路径、每基因的六处理拟时表达曲线，两个重复各画一条；不跨分支混合。
  status<-list()
  if(!is.null(paths)&&!is.null(path_nodes)) {
    if(!is.null(SELECT_PATHS)) {
      if(!all(SELECT_PATHS %in% paths$path))stop('Unknown SELECT_PATHS.')
      paths<-paths[paths$path %in% SELECT_PATHS,,drop=FALSE]
    }
    for(i in seq_len(nrow(paths))) {
      pa<-paths$path[i];cc<-paths$component[i]
      ix<-which(d$component==cc & d$nearest_node %in% path_nodes$node[path_nodes$path==pa]&
                  is.finite(d$local_pseudotime))
      if(length(ix)<20L) {status[[pa]]<-data.frame(path=pa,status='too_few_cells');next}
      tr<-path_trends_v7(d,ex,ix,pa);status[[pa]]<-data.frame(path=pa,status=tr$status)
      if(tr$status!='completed')next
      csv(tr$bins,paste0(pa,'_replicate_expression_bins'));csv(tr$cor,paste0(pa,'_descriptive_correlations'))
      for(g in TARGET_GENES) {
        dd<-tr$bins[tr$bins$gene==g,,drop=FALSE]
        save_v7(file.path(root,paste0(OUT_PREFIX,'_',g,'_',pa,'_pseudotime_trends')),function() {
          par(mfrow=c(2,3),oma=c(3,0,3,0))
          yy<-c(dd$mean_log1p,dd$smooth_log1p);yy<-yy[is.finite(yy)]
          ylim<-if(length(yy)&&max(yy)>0)c(0,max(yy)*1.08) else c(0,1)
          for(tg in GROUPS) {
            z<-dd[dd$treatment==tg,,drop=FALSE]
            par(mar=c(3.3,3.5,2.5,1))
            plot(NA,xlim=range(d$local_pseudotime[ix]),ylim=ylim,xlab='Saved local pseudotime',
                 ylab='Mean log1p expression',main=tg,cex.lab=.8)
            samples<-unique(z$sample)
            for(j in seq_along(samples)) {
              zz<-z[z$sample==samples[j],,drop=FALSE]
              lines(zz$pseudotime,zz$smooth_log1p,col=GROUP_COLORS[tg],lty=j,lwd=1.5)
              points(zz$pseudotime,zz$mean_log1p,col=GROUP_COLORS[tg],pch=if(j==1)16 else 1,cex=.6)
            }
            tryCatch(legend('topright',legend=samples,col=GROUP_COLORS[tg],lty=seq_along(samples),cex=.6,bty='n'),
                     error=function(e)NULL)
            if(!any(is.finite(z$mean_log1p)))text(mean(range(d$local_pseudotime[ix])),mean(ylim),
                                                  'No supported expression bins',cex=.75,col='grey40')
          }
          mtext(paste(g,'|',pa,'| biological replicate trends'),3,outer=TRUE,line=1,font=2)
          mtext('Points = observed bin means; lines = within-supported-run spline. Gaps stay missing. No treatment-DE test.',
                1,outer=TRUE,line=1,cex=.75)
        },13,8)
      }
    }
  }
  if(length(status))csv(do.call(rbind,status),'path_trend_status')
  invisible(TRUE)
}

main_v7 <- function() {
  if(!dir.exists(RUN_DIR))stop('RUN_DIR not found: ',RUN_DIR)
  root<-normalizePath(RUN_DIR,winslash='/',mustWork=TRUE)
  for(pkg in c('monocle3','SummarizedExperiment','SingleCellExperiment'))
    if(!pkg_available(pkg))stop('Required package missing or unloadable: ',pkg)
  d<-read_table_v7(root,'13_cell_metadata_pseudotime_states.csv')
  req<-c('cell','x','y','treatment','sample','component','nearest_node','local_pseudotime')
  if(!all(req %in% names(d)))stop('Missing metadata fields: ',paste(setdiff(req,names(d)),collapse=', '))
  d$treatment[d$treatment=='Control']<-'CK';d$treatment[d$treatment=='MEJA']<-'MeJA'
  if(anyNA(d$sample)||anyNA(d$treatment)||anyDuplicated(d$cell)||!all(d$treatment %in% GROUPS))stop('Invalid metadata.')
  cds<-readRDS(file.path(root,'endodermis_monocle3_ordered.rds'))
  if(!setequal(colnames(cds),d$cell))stop('Saved CDS and metadata contain different cells.')
  xy<-SingleCellExperiment::reducedDim(cds,'UMAP')[d$cell,1:2,drop=FALSE]
  delta<-max(abs(xy-as.matrix(d[,c('x','y')])))
  write.csv(data.frame(max_abs_coordinate_delta=delta),
            file.path(root,paste0(OUT_PREFIX,'_UMAP_coordinate_audit.csv')),row.names=FALSE)
  if(!is.finite(delta)||delta>1e-7)stop('Saved coordinates mismatch; refusing to replace original UMAP.')
  counts<-SummarizedExperiment::assay(cds,'counts')
  sf<-monocle3::size_factors(cds);names(sf)<-colnames(cds);sf<-sf[d$cell]
  if(any(!is.finite(sf)|sf<=0))stop('Invalid saved size factors.')
  mapping<-match_soybean_gid(TARGET_GENES,rownames(counts),PREFERRED_DB_VERSION)
  names(mapping)[names(mapping)=='GID']<-'feature_id_in_counts'
  mapping$user_supplied_label<-USER_LABELS
  if(!is.null(FEATURE_OVERRIDES)) {
    if(is.null(names(FEATURE_OVERRIDES))||anyDuplicated(names(FEATURE_OVERRIDES))||
       !all(names(FEATURE_OVERRIDES) %in% TARGET_GENES)||!all(FEATURE_OVERRIDES %in% rownames(counts)))
      stop('FEATURE_OVERRIDES must be a named target-to-existing-feature vector.')
    ii<-match(names(FEATURE_OVERRIDES),mapping$feature_id)
    mapping$feature_id_in_counts[ii]<-unname(FEATURE_OVERRIDES);mapping$match_status[ii]<-'explicit_user_override'
  }
  mapping$matched<-!is.na(mapping$feature_id_in_counts)
  write.csv(mapping,file.path(root,paste0(OUT_PREFIX,'_target_gene_mapping_audit.csv')),row.names=FALSE,na='')
  ex<-matrix(NA_real_,length(TARGET_GENES),nrow(d),dimnames=list(TARGET_GENES,d$cell))
  for(i in which(mapping$matched))ex[i,]<-log1p(as.numeric(counts[mapping$feature_id_in_counts[i],d$cell])/sf)
  write.csv(data.frame(cell=d$cell,treatment=d$treatment,sample=d$sample,
                       component=d$component,local_pseudotime=d$local_pseudotime,t(ex),check.names=FALSE),
            file.path(root,paste0(OUT_PREFIX,'_sixgenes_cell_expression.csv')),row.names=FALSE,na='')
  # 保留已有graph_test统计的原始字段，不把其q值解释成处理差异p值。
  gt<-read_table_v7(root,'25_graph_test_pooled_candidate_screen.csv',FALSE)
  if(!is.null(gt)&&'feature_id' %in% names(gt)) {
    stats<-merge(mapping,gt,by.x='feature_id_in_counts',by.y='feature_id',all.x=TRUE,sort=FALSE)
    write.csv(stats,file.path(root,paste0(OUT_PREFIX,'_sixgenes_cached_pooled_graph_test.csv')),row.names=FALSE,na='')
  }
  # 注释数据库分支已整体移除：不加载 org.G.max.eg.db / AnnotationDbi，不做 AnnotationDbi::select。
  # 基因行号只靠 match_soybean_gid() 对 rownames(counts) 的字符串解析确定。
  writeLines(c('OrgDb annotation intentionally skipped in this version (v7c).',
               'org.G.max.eg.db and AnnotationDbi are never loaded, so a broken DESCRIPTION cannot abort the run.',
               'Gene identity in the expression matrices comes only from Glyma-ID string matching against rownames(counts);',
               'see v7_target_gene_mapping_audit.csv for per-gene match status.'),
             file.path(root,paste0(OUT_PREFIX,'_annotation_status.txt')))
  edges<-read_table_v7(root,'09_graph_edges.csv');nodes<-read_table_v7(root,'08_principal_nodes.csv')
  edges$component<-nodes$component[match(edges$from,nodes$node)]
  if(anyNA(edges$component))stop('Saved graph node/edge mismatch.')
  roots<-read_table_v7(root,'11b_roots_for_separate_component_heatmaps.csv',FALSE)
  paths<-read_table_v7(root,'26b_root_to_leaf_paths.csv',FALSE)
  pn<-read_table_v7(root,'26c_path_node_membership.csv',FALSE)
  render_analysis_v7(d,ex,edges,nodes,roots,paths,pn,root)
  writeLines(c('Completed six-gene expression maps and saved-pseudotime recoloring.',
               paste('Matched genes:',sum(mapping$matched),'of',nrow(mapping)),
               'Unmatched/ambiguous genes have explicit placeholders; they are not treated as zero expression.',
               paste('Pseudotime color mode:',PT_COLOR_MODE),
               'Pseudotime palette = plasma-like ramp: 0 dark indigo -> purple -> magenta -> orange -> yellow (upper end).',
               'Linear mode maps color continuously to real pseudotime; a pooled-quantile variant is also written with the _pooled_quantile suffix.',
               'Color stops are pooled within component and shared by six treatments. Values and coordinates are unchanged.',
               'Expression is log1p(raw counts / saved size factor). A per-gene color range is shared by all treatments.',
               'User gene names/homology/function hints were not independently validated; use Glyma IDs for interpretation.',
               'Spearman coefficients are descriptive per replicate; no cell-level p-values are used as biological-replicate evidence.',
               'CK denotes the control treatment in these figures, not a cytokinin treatment.',
               'Different graph components have separate local roots; their absolute pseudotimes must not be pooled.',
               'No OrgDb / AnnotationDbi package is loaded; gene matching is string-based only.',
               'No new UMAP, graph, root, graph_test or GO analysis was performed.'),
             file.path(root,paste0(OUT_PREFIX,'_analysis_README.txt')))
  write_session_info_v7(file.path(root,paste0(OUT_PREFIX,'_sessionInfo.txt')),root)
  if(any(!mapping$matched))warning('Some target genes are missing or ambiguous; see v7_target_gene_mapping_audit.csv.')
  cat('Finished in original directory: ',root,'\n',sep='')
  invisible(root)
}

main_v7()
# 完整独立脚本 v7b：六基因原UMAP表达、六处理对比、拟时配色及重复表达动态。
# 用法：source('soybean_endodermis_sixgenes_originalUMAP_v7b_noOrgDb.R', encoding='UTF-8')
#
# 相对 v7 只有三处改动，其余逻辑、输出目录、数值口径完全不变：
#  (1) 彻底移除 org.G.max.eg.db / AnnotationDbi：不再定义 ORGDB_PACKAGE、不再
#      requireNamespace、不再 AnnotationDbi::select。目标基因 -> counts 行名的对应
#      只依赖字符串解析 match_soybean_gid()（v4 版本优先），不需要任何注释数据库。
#      因此 "DESCRIPTION file is missing or broken" 这条警告不会再出现。
#  (2) 结尾的 sessionInfo() 换成不会抛错的安全版本：逐包 tryCatch 取版本，任何一个包
#      的 DESCRIPTION 损坏都只记为 unknown，不会让脚本以
#      "Error in x[['Version']] : subscript out of bounds" 中止。
#  (3) smooth_contiguous() 的自由度改为 min(SMOOTH_DF, 有效点数 - 2) 并 suppressWarnings，
#      消除 "smoothing parameter value too large / setting df = 1" 警告；曲线口径不变。
#
# 只加载已存CDS和CSV，不修改UMAP/轨迹/根/拟时序，不运行降维、graph_test或GO。
# PDF/PNG/CSV直接写回原RUN_DIR，文件名前缀由 OUT_PREFIX 控制（默认仍为 v7）。
RUN_DIR <- '/data01/liuxy/Project/CK2_added_v3_Results/Endodermis_monocle3_originalUMAP_grouped_v3/20261119_163425'
OUT_PREFIX <- 'v7'
TARGET_GENES <- c('Glyma.17G152800','Glyma.06G050000','Glyma.09G139000',
                  'Glyma.06G003800','Glyma.13G094400','Glyma.12G173300')
USER_LABELS <- c('ARR-B candidate','MYB84/RAX3-like','MYB candidate','MYB48-like','MYB36-like','RING-CH candidate')
# 上述短名仅来自用户提供注释，不作为本次验证过的正交同源或功能结论。
GROUPS <- c('CK','ABA','eBL','IAA','MeJA','SA')
GROUP_COLORS <- c(CK='#4D4D4D',ABA='#D55E00',eBL='#009E73',IAA='#0072B2',MeJA='#CC79A7',SA='#E69F00')
PREFERRED_DB_VERSION <- 'v4'
FEATURE_OVERRIDES <- NULL  # 可选命名向量：c('Glyma.17G152800'='实际矩阵行名')
PT_COLOR_MODE <- 'pooled_quantile' # 或linear；默认同时输出线性对照
PT_PROBS <- c(0,.10,.35,.65,.90,1)
PT_COLORS <- c('#243B8A','#007FA8','#60C2B0','#F4D166','#E47D31','#8C1D40')
EXPRESSION_COLORS <- c('#352A87','#8157B7','#C44C8A','#EC8A55','#F3D55B')
EXPR_UPPER_QUANTILE <- .99 # 同一基因全部细胞正表达值的99%分位；高于阈值饱和但原数值保留
ZERO_COLOR <- '#D9D9D9'
POINT_CEX <- .45
DRAW_GRAPH_ON_GENE_MAP <- TRUE
SELECT_PATHS <- NULL
N_TREND_BINS <- 20L
MIN_CELLS_PER_REPLICATE_BIN <- 5L
SMOOTH_DF <- 3L
PNG_DPI <- 200L

parse_soybean_id <- function(ids) {
  raw<-as.character(ids);s<-trimws(raw)
  s<-sub('^(gene|transcript):','',s,ignore.case=TRUE)
  pat<-'^(Glyma)[._-]?([0-9]{2})[gG]([0-9]{6})(.*)$'
  valid<-!is.na(s)&grepl(pat,s,ignore.case=TRUE)
  base<-s;suffix<-rep('',length(s));ver<-rep(NA_character_,length(s))
  base[valid]<-sub(pat,'Glyma.\\2G\\3',s[valid],ignore.case=TRUE)
  suffix[valid]<-sub(pat,'\\4',s[valid],ignore.case=TRUE)
  vv<-regexpr('(^|[._-])v[0-9]+(?=$|[._-])',suffix,perl=TRUE,ignore.case=TRUE)
  has<-which(vv>0)
  if(length(has))ver[has]<-tolower(sub('^[._-]','',regmatches(suffix,vv)))
  tail<-gsub('(^|[._-])v[0-9]+(?=$|[._-])','',suffix,perl=TRUE,ignore.case=TRUE)
  safe<-grepl('^([._-]([0-9]+|[tTpP][0-9]+))*$',tail)
  data.frame(raw_id=raw,clean_id=ifelse(valid&safe,base,s),version=ver,recognized=valid&safe)
}
match_soybean_gid <- function(features,dbkeys,preferred='v4') {
  f<-parse_soybean_id(features);b<-parse_soybean_id(dbkeys);preferred<-tolower(preferred)
  b$priority<-ifelse(!is.na(b$version)&b$version==preferred,1L,ifelse(is.na(b$version),2L,3L))
  groups<-split(seq_len(nrow(b)),b$clean_id)
  out<-data.frame(feature_id=features,clean_id=f$clean_id,feature_version=f$version,
                  GID=NA_character_,db_version=NA_character_,match_status='unmatched',n_db_candidates=0L)
  for(i in seq_len(nrow(f))) {
    if(!is.na(f$version[i])&&f$version[i]!=preferred) {
      out$match_status[i]<-'source_nonpreferred_version_requires_mapping';next
    }
    ii<-groups[[f$clean_id[i]]]
    if(!length(ii))next
    out$n_db_candidates[i]<-length(ii)
    pr<-min(b$priority[ii]);best<-ii[b$priority[ii]==pr]
    if(pr==3L){out$match_status[i]<-'db_only_nonpreferred_version';next}
    exact<-best[b$raw_id[best]==f$raw_id[i]]
    if(length(exact)==1L)best<-exact
    if(length(best)!=1L){out$match_status[i]<-'ambiguous_same_priority';next}
    out$GID[i]<-b$raw_id[best];out$db_version[i]<-b$version[best]
    out$match_status[i]<-if(pr==1L)'preferred_v4' else 'unique_unversioned_GID'
  }
  out$matched<-!is.na(out$GID);out
}

read_table_v7 <- function(root,name,required=TRUE) {
  f<-file.path(root,name)
  if(!file.exists(f)) {
    if(required)stop('Missing saved file: ',f)
    return(NULL)
  }
  read.csv(f,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding='UTF-8',na.strings=c('','NA'))
}
save_v7 <- function(stem,fun,width=12,height=8) {
  pdf(paste0(stem,'.pdf'),width=width,height=height,useDingbats=FALSE)
  tryCatch(fun(),finally=dev.off())
  png(paste0(stem,'.png'),width=round(width*PNG_DPI),height=round(height*PNG_DPI),res=PNG_DPI,
      type=if(capabilities('cairo'))'cairo' else getOption('bitmapType','Xlib'))
  tryCatch(fun(),finally=dev.off())
}

# 安全包版本查询：DESCRIPTION 缺失/损坏的包返回 NA，绝不抛出
# "Error in x[['Version']] : subscript out of bounds"。
pkg_version_chr <- function(pkg) {
  tryCatch({
    v<-suppressWarnings(utils::packageVersion(pkg))
    if(is.null(v)||!length(v))NA_character_ else as.character(v)
  },error=function(e)NA_character_,warning=function(w)NA_character_)
}
# 安全加载检查：坏包只返回 FALSE，不会让脚本中止在 requireNamespace 内部。
pkg_available <- function(pkg) {
  isTRUE(suppressWarnings(suppressMessages(tryCatch(
    requireNamespace(pkg,quietly=TRUE),error=function(e)FALSE))))
}
write_session_info_v7 <- function(file,root) {
  pk<-sort(unique(c('base','stats','utils','graphics','grDevices','methods','monocle3',
                    'SummarizedExperiment','SingleCellExperiment',loadedNamespaces())))
  L<-c(paste0('R version: ',R.version.string),
       paste0('Platform: ',R.version$platform),
       paste0('Written at: ',format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z')),
       paste0('No OrgDb / AnnotationDbi package is loaded by this script.'),
       '',
       'Loaded namespaces and versions; unknown = DESCRIPTION unreadable (reported, never fatal):')
  for(p in pk) {
    v<-pkg_version_chr(p)
    L<-c(L,sprintf('  %-30s %s',p,if(is.na(v))'unknown' else v))
  }
  tryCatch(writeLines(L,file,useBytes=FALSE),error=function(e)NULL)
  # 完整 sessionInfo 只作附加信息，失败也不影响主流程。
  si<-suppressWarnings(tryCatch(utils::capture.output(utils::sessionInfo()),error=function(e)NULL))
  if(length(si))
    tryCatch(writeLines(si,file.path(root,paste0(OUT_PREFIX,'_sessionInfo_full.txt')),useBytes=FALSE),
             error=function(e)NULL)
  invisible(file)
}

make_pt_scale <- function(pt,mode=PT_COLOR_MODE) {
  v<-pt[is.finite(pt)]
  if(!length(v))stop('No finite pseudotime.')
  if(!mode %in% c('pooled_quantile','linear'))stop('Invalid PT_COLOR_MODE.')
  if(diff(range(v))==0)return(list(knots=v[1],colors=PT_COLORS[1],mode=mode))
  knots<-if(mode=='linear')seq(min(v),max(v),length.out=length(PT_COLORS)) else
    as.numeric(quantile(v,PT_PROBS,names=FALSE,type=7))
  knots<-unique(knots)
  # 若大量拟时并列，则合并重复断点；始终保持颜色映射单调。
  list(knots=knots,colors=colorRampPalette(PT_COLORS)(length(knots)),mode=mode)
}
color_by_scale <- function(values,scale,na_color='#C5C5C5') {
  out<-rep(na_color,length(values));ok<-is.finite(values)
  if(length(scale$knots)==1L) {out[ok]<-scale$colors[1];return(out)}
  u<-approx(scale$knots,seq(0,1,length.out=length(scale$knots)),xout=values[ok],rule=2)$y
  pal<-colorRampPalette(scale$colors)(256)
  out[ok]<-pal[pmin(256L,pmax(1L,1L+floor(255*u)))]
  out
}
make_expr_scale <- function(v) {
  p<-v[is.finite(v)&v>0]
  upper<-if(length(p))as.numeric(quantile(p,EXPR_UPPER_QUANTILE,names=FALSE)) else 1
  if(!is.finite(upper)||upper<=0)upper<-1
  list(knots=seq(0,upper,length.out=length(EXPRESSION_COLORS)),colors=EXPRESSION_COLORS,
       upper=upper,n_clipped=sum(v>upper,na.rm=TRUE),mode='linear')
}
color_legend_v7 <- function(scale,title) {
  par(mar=c(2,.4,3,3.2));plot.new();plot.window(xlim=c(0,1),ylim=c(0,1))
  title(main=title,cex.main=.82)
  if(length(scale$knots)<2) {
    rect(.2,.2,.45,.8,col=scale$colors[1],border=NA)
    text(.5,.5,format(scale$knots[1],digits=3),adj=0,cex=.8);return(invisible(NULL))
  }
  y<-seq(.15,.87,length.out=257)
  rect(.15,head(y,-1),.4,tail(y,-1),col=colorRampPalette(scale$colors)(256),border=NA)
  ticks<-seq(.15,.87,length.out=length(scale$knots))
  segments(.4,ticks,.47,ticks)
  labels<-vapply(scale$knots,function(z)format(signif(z,3),trim=TRUE),character(1))
  text(.51,ticks,labels,adj=0,cex=.75)
  if(title=='Expression')text(.1,.04,'Grey = zero',adj=0,cex=.7)
}

draw_umap_v7 <- function(d,ix,values,scale,title,edges=NULL,nodes=NULL,root_nodes=NULL,
                         expression=FALSE,missing=FALSE) {
  par(mar=c(3.1,3.3,3,1),xaxs='r',yaxs='r')
  plot(d$x,d$y,col='#EEEEEE',pch=16,cex=POINT_CEX*.75,asp=1,
       xlab='Original UMAP 1',ylab='Original UMAP 2',main=title,cex.main=.95,cex.axis=.8,cex.lab=.8)
  if(!is.null(edges)&&nrow(edges))segments(edges$x,edges$y,edges$xend,edges$yend,col='#707070',lwd=.8)
  if(missing) {
    mtext('Gene missing / ambiguous; see mapping audit',side=3,line=.2,cex=.68,col='#A32929')
    return(invisible(NULL))
  }
  vv<-values[ix];cols<-color_by_scale(vv,scale)
  if(expression)cols[is.finite(vv)&vv==0]<-ZERO_COLOR
  # 低值先画，高值后画，避免高表达细胞被零表达覆盖；坐标不变。
  ord<-order(vv,na.last=FALSE)
  points(d$x[ix[ord]],d$y[ix[ord]],pch=16,col=cols[ord],cex=POINT_CEX)
  if(length(root_nodes)&&!is.null(nodes)) {
    rr<-nodes[nodes$node %in% root_nodes,,drop=FALSE]
    points(rr$x,rr$y,pch=23,bg='white',col='black',cex=1.1,lwd=1)
  }
}

# 仅在同一路径内求曲线；各处理/重复使用同一套分箱边界。
# 缺失箱为NA并断线；样条只在至少4个连续有效箱内求值，不跨缺失区域或向外推。
# df 至少比有效点数少2，且拟合全程 suppressWarnings，避免
#"smoothing parameter value too large / setting df = 1" 警告（不改变曲线口径）。
smooth_contiguous <- function(pt,value,valid,df=SMOOTH_DF) {
  out<-rep(NA_real_,length(pt))
  if(!any(valid))return(out)
  runs<-split(which(valid),cumsum(c(TRUE,diff(which(valid))!=1)))
  for(ii in runs) {
    out[ii]<-value[ii]
    if(length(ii)>=4L&&length(unique(pt[ii]))>=4L) {
      df_i<-max(1L,min(as.integer(df),length(ii)-2L))
      fit<-tryCatch(suppressWarnings(smooth.spline(pt[ii],value[ii],df=df_i)),error=function(e)NULL)
      if(!is.null(fit)) {
        pr<-tryCatch(suppressWarnings(predict(fit,pt[ii])$y),error=function(e)NULL)
        if(!is.null(pr)&&length(pr)==length(ii))out[ii]<-pmax(0,pr)
      }
    }
  }
  out
}
path_trends_v7 <- function(d,ex,ix,path) {
  pt<-d$local_pseudotime[ix]
  br<-unique(as.numeric(quantile(pt,seq(0,1,length.out=N_TREND_BINS+1),names=FALSE)))
  if(length(br)<5L)return(list(status='too_few_distinct_time_bins',bins=data.frame(),cor=data.frame()))
  bn<-cut(pt,br,include.lowest=TRUE,labels=FALSE)
  centers<-(head(br,-1)+tail(br,-1))/2;nb<-length(centers)
  rows<-list();corr<-list()
  for(tg in GROUPS) {
    samples<-unique(d$sample[d$treatment==tg])
    for(s in samples) {
      jj<-which(d$treatment[ix]==tg & d$sample[ix]==s)
      for(gene in rownames(ex)) {
        v<-ex[gene,ix[jj]];p<-pt[jj]
        ok<-is.finite(v)&is.finite(p)
        rho<-if(sum(ok)>=10&&length(unique(v[ok]))>1&&length(unique(p[ok]))>1)
          suppressWarnings(cor(p[ok],v[ok],method='spearman')) else NA_real_
        corr[[paste(tg,s,gene)]]<-data.frame(path=path,treatment=tg,sample=s,gene=gene,
                                             n_cells=sum(ok),detected_fraction=if(sum(ok))mean(v[ok]>0) else NA_real_,
                                             rho_spearman_descriptive=rho)
        means<-rep(NA_real_,nb);ns<-integer(nb)
        for(b in seq_len(nb)) {
          kk<-jj[bn[jj]==b];vv<-ex[gene,ix[kk]];vv<-vv[is.finite(vv)];ns[b]<-length(vv)
          if(ns[b]>=MIN_CELLS_PER_REPLICATE_BIN)means[b]<-mean(vv)
        }
        smooth<-smooth_contiguous(centers,means,is.finite(means))
        rows[[paste(tg,s,gene)]]<-data.frame(path=path,treatment=tg,sample=s,gene=gene,
                                             bin=seq_len(nb),left=head(br,-1),right=tail(br,-1),pseudotime=centers,
                                             n_cells=ns,mean_log1p=means,smooth_log1p=smooth)
      }
    }
  }
  list(status='completed',bins=do.call(rbind,rows),cor=do.call(rbind,corr))
}

render_analysis_v7 <- function(d,ex,edges,nodes,roots,paths,path_nodes,root) {
  csv<-function(x,nm)write.csv(x,file.path(root,paste0(OUT_PREFIX,'_',nm,'.csv')),row.names=FALSE,na='')
  settings<-list();expr_scales<-lapply(seq_len(nrow(ex)),function(i)make_expr_scale(ex[i,]))
  names(expr_scales)<-rownames(ex)
  # 1. 第二张图重新着色；每个分量各自定义色阶，六处理共用，绝不按每组拉伸。
  for(cc in unique(d$component)) {
    ids<-which(d$component==cc & is.finite(d$local_pseudotime))
    if(!length(ids))next
    ee<-edges[edges$component==cc,,drop=FALSE]
    rr<-if(!is.null(roots))roots$root_node[roots$component==cc] else character()
    for(mode in unique(c(PT_COLOR_MODE,'linear'))) {
      sc<-make_pt_scale(d$local_pseudotime[ids],mode)
      settings[[paste(cc,mode)]]<-data.frame(component=cc,mode=mode,pseudotime=sc$knots,color=sc$colors)
      save_v7(file.path(root,paste0(OUT_PREFIX,'_P2_component_',cc,'_',mode)),function() {
        layout(matrix(c(1,2,3,7,4,5,6,7),2,4,byrow=TRUE),widths=c(1,1,1,.38))
        par(oma=c(2,0,3,0))
        for(tg in GROUPS) {
          ii<-ids[d$treatment[ids]==tg]
          draw_umap_v7(d,ii,d$local_pseudotime,sc,paste(tg,'| local PT'),ee,nodes,rr)
        }
        color_legend_v7(sc,'Local PT')
        mtext(paste('Endodermis | original UMAP | component',cc),3,outer=TRUE,line=1,font=2)
        mtext(if(mode=='pooled_quantile')'Quantile-spaced color stops pooled across treatments; legend = original PT values.' else
          'Linear color mapping; legend = original PT values.',1,outer=TRUE,line=.4,cex=.8)
      },14,9)
    }
  }
  csv(do.call(rbind,settings),'pseudotime_color_knots')
  # 2. 六基因总体图：各基因单独色阶，数值上限明确写出；同一基因所有组范围一致。
  save_v7(file.path(root,paste0(OUT_PREFIX,'_sixgenes_originalUMAP_overview')),function() {
    par(mfrow=c(2,3),oma=c(2,0,2,0))
    for(g in TARGET_GENES) {
      sc<-expr_scales[[g]];miss<-all(!is.finite(ex[g,]))
      draw_umap_v7(d,seq_len(nrow(d)),ex[g,],sc,g,
                   if(DRAW_GRAPH_ON_GENE_MAP)edges else NULL,expression=TRUE,missing=miss)
      if(!miss)mtext(paste0('Expression 0 to ',signif(sc$upper,3),' (positive q',EXPR_UPPER_QUANTILE,')'),3,line=.2,cex=.62)
    }
    mtext('Six candidate genes | original Endodermis UMAP',3,outer=TRUE,font=2)
    mtext('Grey = zero. Colored = log1p(count / saved size factor). Per-gene scale; see individual panels for legends.',1,outer=TRUE,cex=.75)
  },13,9)
  for(g in TARGET_GENES) {
    sc<-expr_scales[[g]];miss<-all(!is.finite(ex[g,]))
    save_v7(file.path(root,paste0(OUT_PREFIX,'_',g,'_six_treatments')),function() {
      layout(matrix(c(1,2,3,7,4,5,6,7),2,4,byrow=TRUE),widths=c(1,1,1,.38))
      par(oma=c(2,0,3,0))
      for(tg in GROUPS)draw_umap_v7(d,which(d$treatment==tg),ex[g,],sc,tg,
                                    if(DRAW_GRAPH_ON_GENE_MAP)edges else NULL,expression=TRUE,missing=miss)
      if(miss) {
        par(mar=c(1,1,2,1));plot.new();text(.5,.5,'No mapped\nexpression',cex=.8)
      } else color_legend_v7(sc,'Expression')
      mtext(paste(g,'| same expression scale across all six treatments'),3,outer=TRUE,line=1,font=2)
      mtext(paste0('log1p(count / size factor); upper=',signif(sc$upper,3),
                   '; values above cap=',sc$n_clipped,'. Grey: zero; pale background: reference cells.'),1,outer=TRUE,cex=.75)
    },14,9)
  }
  # 3. 保存每个生物学重复的表达/检测率，避免用单细胞数代替重复数。
  summary<-list()
  for(g in TARGET_GENES)for(tg in GROUPS)for(s in unique(d$sample[d$treatment==tg])) {
    v<-ex[g,d$treatment==tg & d$sample==s];v<-v[is.finite(v)]
    summary[[paste(g,tg,s)]]<-data.frame(gene=g,treatment=tg,sample=s,n_cells=length(v),
                                         mean_log1p=if(length(v))mean(v) else NA_real_,median_log1p=if(length(v))median(v) else NA_real_,
                                         detected_fraction=if(length(v))mean(v>0) else NA_real_)
  }
  csv(do.call(rbind,summary),'gene_expression_per_replicate')
  csv(do.call(rbind,lapply(TARGET_GENES,function(g)data.frame(gene=g,
                                                              upper=expr_scales[[g]]$upper,positive_quantile=EXPR_UPPER_QUANTILE,n_clipped=expr_scales[[g]]$n_clipped))),
      'expression_color_limits')
  # 4. 每路径、每基因的六处理拟时表达曲线，两个重复各画一条；不跨分支混合。
  status<-list()
  if(!is.null(paths)&&!is.null(path_nodes)) {
    if(!is.null(SELECT_PATHS)) {
      if(!all(SELECT_PATHS %in% paths$path))stop('Unknown SELECT_PATHS.')
      paths<-paths[paths$path %in% SELECT_PATHS,,drop=FALSE]
    }
    for(i in seq_len(nrow(paths))) {
      pa<-paths$path[i];cc<-paths$component[i]
      ix<-which(d$component==cc & d$nearest_node %in% path_nodes$node[path_nodes$path==pa]&
                  is.finite(d$local_pseudotime))
      if(length(ix)<20L) {status[[pa]]<-data.frame(path=pa,status='too_few_cells');next}
      tr<-path_trends_v7(d,ex,ix,pa);status[[pa]]<-data.frame(path=pa,status=tr$status)
      if(tr$status!='completed')next
      csv(tr$bins,paste0(pa,'_replicate_expression_bins'));csv(tr$cor,paste0(pa,'_descriptive_correlations'))
      for(g in TARGET_GENES) {
        dd<-tr$bins[tr$bins$gene==g,,drop=FALSE]
        save_v7(file.path(root,paste0(OUT_PREFIX,'_',g,'_',pa,'_pseudotime_trends')),function() {
          par(mfrow=c(2,3),oma=c(3,0,3,0))
          yy<-c(dd$mean_log1p,dd$smooth_log1p);yy<-yy[is.finite(yy)]
          ylim<-if(length(yy)&&max(yy)>0)c(0,max(yy)*1.08) else c(0,1)
          for(tg in GROUPS) {
            z<-dd[dd$treatment==tg,,drop=FALSE]
            par(mar=c(3.3,3.5,2.5,1))
            plot(NA,xlim=range(d$local_pseudotime[ix]),ylim=ylim,xlab='Saved local pseudotime',
                 ylab='Mean log1p expression',main=tg,cex.lab=.8)
            samples<-unique(z$sample)
            for(j in seq_along(samples)) {
              zz<-z[z$sample==samples[j],,drop=FALSE]
              lines(zz$pseudotime,zz$smooth_log1p,col=GROUP_COLORS[tg],lty=j,lwd=1.5)
              points(zz$pseudotime,zz$mean_log1p,col=GROUP_COLORS[tg],pch=if(j==1)16 else 1,cex=.6)
            }
            tryCatch(legend('topright',legend=samples,col=GROUP_COLORS[tg],lty=seq_along(samples),cex=.6,bty='n'),
                     error=function(e)NULL)
            if(!any(is.finite(z$mean_log1p)))text(mean(range(d$local_pseudotime[ix])),mean(ylim),
                                                  'No supported expression bins',cex=.75,col='grey40')
          }
          mtext(paste(g,'|',pa,'| biological replicate trends'),3,outer=TRUE,line=1,font=2)
          mtext('Points = observed bin means; lines = within-supported-run spline. Gaps stay missing. No treatment-DE test.',
                1,outer=TRUE,line=1,cex=.75)
        },13,8)
      }
    }
  }
  if(length(status))csv(do.call(rbind,status),'path_trend_status')
  invisible(TRUE)
}

main_v7 <- function() {
  if(!dir.exists(RUN_DIR))stop('RUN_DIR not found: ',RUN_DIR)
  root<-normalizePath(RUN_DIR,winslash='/',mustWork=TRUE)
  for(pkg in c('monocle3','SummarizedExperiment','SingleCellExperiment'))
    if(!pkg_available(pkg))stop('Required package missing or unloadable: ',pkg)
  d<-read_table_v7(root,'13_cell_metadata_pseudotime_states.csv')
  req<-c('cell','x','y','treatment','sample','component','nearest_node','local_pseudotime')
  if(!all(req %in% names(d)))stop('Missing metadata fields: ',paste(setdiff(req,names(d)),collapse=', '))
  d$treatment[d$treatment=='Control']<-'CK';d$treatment[d$treatment=='MEJA']<-'MeJA'
  if(anyNA(d$sample)||anyNA(d$treatment)||anyDuplicated(d$cell)||!all(d$treatment %in% GROUPS))stop('Invalid metadata.')
  cds<-readRDS(file.path(root,'endodermis_monocle3_ordered.rds'))
  if(!setequal(colnames(cds),d$cell))stop('Saved CDS and metadata contain different cells.')
  xy<-SingleCellExperiment::reducedDim(cds,'UMAP')[d$cell,1:2,drop=FALSE]
  delta<-max(abs(xy-as.matrix(d[,c('x','y')])))
  write.csv(data.frame(max_abs_coordinate_delta=delta),
            file.path(root,paste0(OUT_PREFIX,'_UMAP_coordinate_audit.csv')),row.names=FALSE)
  if(!is.finite(delta)||delta>1e-7)stop('Saved coordinates mismatch; refusing to replace original UMAP.')
  counts<-SummarizedExperiment::assay(cds,'counts')
  sf<-monocle3::size_factors(cds);names(sf)<-colnames(cds);sf<-sf[d$cell]
  if(any(!is.finite(sf)|sf<=0))stop('Invalid saved size factors.')
  mapping<-match_soybean_gid(TARGET_GENES,rownames(counts),PREFERRED_DB_VERSION)
  names(mapping)[names(mapping)=='GID']<-'feature_id_in_counts'
  mapping$user_supplied_label<-USER_LABELS
  if(!is.null(FEATURE_OVERRIDES)) {
    if(is.null(names(FEATURE_OVERRIDES))||anyDuplicated(names(FEATURE_OVERRIDES))||
       !all(names(FEATURE_OVERRIDES) %in% TARGET_GENES)||!all(FEATURE_OVERRIDES %in% rownames(counts)))
      stop('FEATURE_OVERRIDES must be a named target-to-existing-feature vector.')
    ii<-match(names(FEATURE_OVERRIDES),mapping$feature_id)
    mapping$feature_id_in_counts[ii]<-unname(FEATURE_OVERRIDES);mapping$match_status[ii]<-'explicit_user_override'
  }
  mapping$matched<-!is.na(mapping$feature_id_in_counts)
  write.csv(mapping,file.path(root,paste0(OUT_PREFIX,'_target_gene_mapping_audit.csv')),row.names=FALSE,na='')
  ex<-matrix(NA_real_,length(TARGET_GENES),nrow(d),dimnames=list(TARGET_GENES,d$cell))
  for(i in which(mapping$matched))ex[i,]<-log1p(as.numeric(counts[mapping$feature_id_in_counts[i],d$cell])/sf)
  write.csv(data.frame(cell=d$cell,treatment=d$treatment,sample=d$sample,
                       component=d$component,local_pseudotime=d$local_pseudotime,t(ex),check.names=FALSE),
            file.path(root,paste0(OUT_PREFIX,'_sixgenes_cell_expression.csv')),row.names=FALSE,na='')
  # 保留已有graph_test统计的原始字段，不把其q值解释成处理差异p值。
  gt<-read_table_v7(root,'25_graph_test_pooled_candidate_screen.csv',FALSE)
  if(!is.null(gt)&&'feature_id' %in% names(gt)) {
    stats<-merge(mapping,gt,by.x='feature_id_in_counts',by.y='feature_id',all.x=TRUE,sort=FALSE)
    write.csv(stats,file.path(root,paste0(OUT_PREFIX,'_sixgenes_cached_pooled_graph_test.csv')),row.names=FALSE,na='')
  }
  # 注释数据库分支已整体移除：不加载 org.G.max.eg.db / AnnotationDbi，不做 AnnotationDbi::select。
  # 基因行号只靠 match_soybean_gid() 对 rownames(counts) 的字符串解析确定。
  writeLines(c('OrgDb annotation intentionally skipped in this version (v7b).',
               'org.G.max.eg.db and AnnotationDbi are never loaded, so a broken DESCRIPTION cannot abort the run.',
               'Gene identity in the expression matrices comes only from Glyma-ID string matching against rownames(counts);',
               'see v7_target_gene_mapping_audit.csv for per-gene match status.'),
             file.path(root,paste0(OUT_PREFIX,'_annotation_status.txt')))
  edges<-read_table_v7(root,'09_graph_edges.csv');nodes<-read_table_v7(root,'08_principal_nodes.csv')
  edges$component<-nodes$component[match(edges$from,nodes$node)]
  if(anyNA(edges$component))stop('Saved graph node/edge mismatch.')
  roots<-read_table_v7(root,'11b_roots_for_separate_component_heatmaps.csv',FALSE)
  paths<-read_table_v7(root,'26b_root_to_leaf_paths.csv',FALSE)
  pn<-read_table_v7(root,'26c_path_node_membership.csv',FALSE)
  render_analysis_v7(d,ex,edges,nodes,roots,paths,pn,root)
  writeLines(c('Completed six-gene expression maps and saved-pseudotime recoloring.',
               paste('Matched genes:',sum(mapping$matched),'of',nrow(mapping)),
               'Unmatched/ambiguous genes have explicit placeholders; they are not treated as zero expression.',
               paste('Pseudotime color mode:',PT_COLOR_MODE),
               'Quantile color stops are pooled within component and shared by six treatments. Values and coordinates are unchanged.',
               'Expression is log1p(raw counts / saved size factor). A per-gene color range is shared by all treatments.',
               'User gene names/homology/function hints were not independently validated; use Glyma IDs for interpretation.',
               'Spearman coefficients are descriptive per replicate; no cell-level p-values are used as biological-replicate evidence.',
               'CK denotes the control treatment in these figures, not a cytokinin treatment.',
               'Different graph components have separate local roots; their absolute pseudotimes must not be pooled.',
               'No OrgDb / AnnotationDbi package is loaded; gene matching is string-based only.',
               'No new UMAP, graph, root, graph_test or GO analysis was performed.'),
             file.path(root,paste0(OUT_PREFIX,'_analysis_README.txt')))
  write_session_info_v7(file.path(root,paste0(OUT_PREFIX,'_sessionInfo.txt')),root)
  if(any(!mapping$matched))warning('Some target genes are missing or ambiguous; see v7_target_gene_mapping_audit.csv.')
  cat('Finished in original directory: ',root,'\n',sep='')
  invisible(root)
}

main_v7()




## ============================================================
## M9 (HormoneProgram09) Top hub 基因 × 细胞类型/胁迫 表达热图
## ============================================================
## ============================================================
## M9 (HormoneProgram09) Top hub 基因 × 细胞类型 × 五种激素 表达热图
## 激素：ABA, IAA, eBL, MeJA, SA
## ============================================================
## M9 (HormoneProgram09) Top hub 基因 × 细胞类型 × 五种激素 表达热图
## 激素：ABA, IAA, eBL, MeJA, SA
## 数据：长表，激素标签在 contrast 列（形如 ABA-REF）
## ============================================================

library(data.table)
library(ggplot2)
library(pheatmap)
library(dplyr)

## ---------- 1. DEG 文件自动查找 ----------
PAPER_DIR <- "/data01/liuxy/Project/muscat_results_ck_ref_paper"
DEG_FILE_CANDIDATES <- c(
  file.path(PAPER_DIR, "All_DEG_filtered_paper_thresholds.csv"),
  file.path(PAPER_DIR, "All_DEG_filtered_combined.csv"),
  file.path(PAPER_DIR, "All_DEG_filtered_paper_thresholds.csv.gz"),
  file.path(PAPER_DIR, "All_DEG_filtered_combined.csv.gz")
)
DEG_FILE <- DEG_FILE_CANDIDATES[file.exists(DEG_FILE_CANDIDATES)][1]
if (is.na(DEG_FILE)) stop("找不到 DEG 文件")
cat("读取 DEG 文件：", DEG_FILE, "\n")

## ---------- 2. 读取 Top50 kME 表，筛 M9 (HormoneProgram09) ----------
top50 <- read.csv(
  "/data01/liuxy/Project/CK2_added_v3_Results/hdWGCNA_Hormones5_CK2_Endo3_v1/Hormones5_Endo3_v1_module_Top50_kME.csv"
)
m9_top50 <- top50[top50$Program == "HormoneProgram09", ]
m9_top50 <- m9_top50[order(-m9_top50$kME), ]
target_genes <- m9_top50$Gene
cat("纳入热图的 M9 基因数：", length(target_genes), "\n")

## ---------- 3. 读取 DEG 表 ----------
deg <- fread(DEG_FILE)

col_gene     <- "gene"
col_cell     <- "cluster_id"
col_contrast <- "contrast"
col_log2fc   <- "logFC"

## ---------- 4. 解析 contrast -> 激素名，只保留五种激素 ----------
## contrast 形如：ABA-REF / eBL-REF / IAA-REF / MEJA-REF / SA-REF ...
parse_hormone <- function(x) {
  x <- sub("-REF$", "", x)          # 去掉尾部 -REF
  x <- sub("_vs_REF$", "", x)       # 兼容 comparison 写法
  x
}
deg[, Hormone := parse_hormone(get(col_contrast))]

## 统一命名（MEJA -> MeJA）
hormone_rename <- c(
  "ABA"  = "ABA",
  "IAA"  = "IAA",
  "eBL"  = "eBL",
  "MEJA" = "MeJA",
  "MeJA" = "MeJA",
  "SA"   = "SA"
)
deg[, Hormone := hormone_rename[Hormone]]

cat("\n=== 全部激素取值计数 ===\n")
print(table(deg$Hormone, useNA = "ifany"))

KEEP_HORMONE <- c("ABA", "IAA", "eBL", "MeJA", "SA")
deg <- deg[Hormone %in% KEEP_HORMONE]
cat("\n保留激素：", paste(unique(deg$Hormone), collapse = ", "),
    "；行数：", nrow(deg), "\n")

## ---------- 5. 筛选目标基因 ----------
sub <- deg[get(col_gene) %in% target_genes]
if (nrow(sub) == 0) stop("目标基因在 DEG 表中未匹配到，请检查基因 ID")
cat("筛选后行数：", nrow(sub), "\n")

## 核对：有多少 M9 基因在 DEG 表里出现过
hit <- target_genes[target_genes %in% unique(sub[[col_gene]])]
cat("在 DEG 中命中的基因数：", length(hit), "/", length(target_genes), "\n")
miss <- setdiff(target_genes, unique(sub[[col_gene]]))
if (length(miss) > 0) cat("未命中：", paste(miss, collapse = ", "), "\n")

## ---------- 6. 长表 -> 宽表：行为基因，列为 细胞类型_激素 ----------
sub[, combo := paste(get(col_cell), Hormone, sep = "|")]

wide <- dcast(
  sub,
  as.formula(paste(col_gene, "~ combo")),
  value.var = col_log2fc,
  fun.aggregate = mean,
  na.rm = TRUE
)

genes_vec <- wide[[col_gene]]
mat <- as.matrix(wide[, -1, with = FALSE])
rownames(mat) <- genes_vec
mat <- mat[target_genes[target_genes %in% rownames(mat)], , drop = FALSE]

## ---------- 7. 顶部注释 ----------
annot_cols <- data.table(combo = colnames(mat))
annot_cols[, c("Cell_Type", "Hormone") := tstrsplit(combo, "\\|")]
annot_col_df <- as.data.frame(annot_cols[, .(Cell_Type, Hormone)])
rownames(annot_col_df) <- annot_cols$combo

## ---------- 8. 配色 ----------
cell_type_colors <- c(
  "Cortex"        = "#7fd4c1",
  "Dividing Cell" = "#fcf9a8",
  "Endodermis"    = "#b39ddb",
  "Epidermis"     = "#ef6c6c",
  "Lateral Root"  = "#8d6e63",
  "Meristem"      = "#6fb3d2",
  "Pericycle"     = "#f5a65b",
  "Phloem"        = "#f8c8c8",
  "QC"            = "#fce4ec",
  "Root Cap"      = "#b088c8",
  "Xylem"         = "#c8e6c9"
)
hormone_colors <- c(
  "ABA"  = "#90c9a8",
  "IAA"  = "#6fb3d2",
  "eBL"  = "#f5a65b",
  "MeJA" = "#b088c8",
  "SA"   = "#ef6c6c"
)
ann_colors <- list(
  Cell_Type = cell_type_colors[names(cell_type_colors) %in% annot_col_df$Cell_Type],
  Hormone   = hormone_colors[names(hormone_colors) %in% annot_col_df$Hormone]
)

## ---------- 9. 画热图 ----------
hm_colors <- colorRampPalette(c("#3949ab", "#9fa8da", "#eeeeee", "#ef9a9a", "#b71c1c"))(100)

out_pdf <- file.path(PAPER_DIR, "M9_hub_genes_hormone_heatmap.pdf")
pdf(out_pdf, width = 14, height = 10)

pheatmap(
  mat,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale        = "none",
  color        = hm_colors,
  annotation_col = annot_col_df,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  show_colnames = FALSE,
  main = "M9 (HormoneProgram09) hub genes under five hormone treatments",
  treeheight_row = 0,
  treeheight_col = 0,
  na_col = "#dddddd",
  fontsize_row = 8,
  annotation_legend = TRUE
)
dev.off()

cat("热图已保存：", out_pdf, "\n")
fwrite(sub, file.path(PAPER_DIR, "M9_hub_genes_hormone_expression.csv"))
cat("筛选明细表已导出：", file.path(PAPER_DIR, "M9_hub_genes_hormone_expression.csv"), "\n")


library(data.table)
library(ggplot2)
library(pheatmap)
library(dplyr)

## ---------- 1. DEG 文件自动查找 ----------
PAPER_DIR <- "/data01/liuxy/Project/muscat_results_ck_ref_paper"
DEG_FILE_CANDIDATES <- c(
  file.path(PAPER_DIR, "All_DEG_filtered_paper_thresholds.csv"),
  file.path(PAPER_DIR, "All_DEG_filtered_combined.csv"),
  file.path(PAPER_DIR, "All_DEG_filtered_paper_thresholds.csv.gz"),
  file.path(PAPER_DIR, "All_DEG_filtered_combined.csv.gz")
)
DEG_FILE <- DEG_FILE_CANDIDATES[file.exists(DEG_FILE_CANDIDATES)][1]
if (is.na(DEG_FILE)) stop("找不到 DEG 文件")
cat("读取 DEG 文件：", DEG_FILE, "\n")

## ---------- 2. 读取 Top50 kME 表，筛 M9 (HormoneProgram09) ----------
top50 <- read.csv(
  "/data01/liuxy/Project/CK2_added_v3_Results/hdWGCNA_Hormones5_CK2_Endo3_v1/Hormones5_Endo3_v1_module_Top50_kME.csv"
)

m9_top50 <- top50[top50$Program == "HormoneProgram09", ]
m9_top50 <- m9_top50[order(-m9_top50$kME), ]
print(m9_top50, row.names = FALSE)

## 基因列名在你的 CSV 里是 Gene（大写 G）
target_genes <- m9_top50$Gene
cat("纳入热图的 M9 基因数：", length(target_genes), "\n")

## ---------- 3. 读取 DEG 表 ----------
deg <- fread(DEG_FILE)
cat("DEG 表列名：\n"); print(colnames(deg))

## 按你实际列名固定映射
col_gene   <- "gene"
col_cell   <- "cluster_id"
col_contrast <- "contrast"
col_log2fc <- "logFC"

## ---------- 4. 解析 contrast -> 激素名，只保留五种激素 ----------
## contrast 形如 ABA_vs_REF / IAA1_vs_REF1 / eBL2_vs_REF2 ...
parse_hormone <- function(x) {
  x <- sub("_vs_.*$", "", x)          # 去掉 _vs_REF...
  x <- sub("[0-9]+$", "", x)          # 去掉尾部重复编号 1/2
  x
}
deg[, Hormone := parse_hormone(get(col_contrast))]

## 五组激素统一命名（MEJA -> MeJA，与图例一致）
hormone_rename <- c(
  "ABA"  = "ABA",
  "IAA"  = "IAA",
  "eBL"  = "eBL",
  "MEJA" = "MeJA",
  "MeJA" = "MeJA",
  "SA"   = "SA"
)
deg[, Hormone := hormone_rename[Hormone]]

KEEP_HORMONE <- c("ABA", "IAA", "eBL", "MeJA", "SA")
deg <- deg[Hormone %in% KEEP_HORMONE]
cat("保留激素：", paste(unique(deg$Hormone), collapse = ", "), "\n")

## ---------- 5. 筛选目标基因 ----------
sub <- deg[get(col_gene) %in% target_genes]
if (nrow(sub) == 0) stop("目标基因在 DEG 表中未匹配到，请检查基因 ID")
cat("筛选后行数：", nrow(sub), "\n")

## ---------- 6. 长表 -> 宽表：行为基因，列为 细胞类型_激素 ----------
sub[, combo := paste(get(col_cell), Hormone, sep = "|")]

wide <- dcast(
  sub,
  as.formula(paste(col_gene, "~ combo")),
  value.var = col_log2fc,
  fun.aggregate = mean,
  na.rm = TRUE
)

genes_vec <- wide[[col_gene]]
mat <- as.matrix(wide[, -1, with = FALSE])
rownames(mat) <- genes_vec
## 按 M9 hub 排名排列行
mat <- mat[target_genes[target_genes %in% rownames(mat)], , drop = FALSE]

## ---------- 7. 顶部注释（拆回 Cell_Type / Hormone） ----------
annot_cols <- data.table(combo = colnames(mat))
annot_cols[, c("Cell_Type", "Hormone") := tstrsplit(combo, "\\|")]

annot_col_df <- as.data.frame(annot_cols[, .(Cell_Type, Hormone)])
rownames(annot_col_df) <- annot_cols$combo

## ---------- 8. 配色 ----------
cell_type_colors <- c(
  "Cortex"        = "#7fd4c1",
  "Dividing Cell" = "#fcf9a8",
  "Endodermis"    = "#b39ddb",
  "Epidermis"     = "#ef6c6c",
  "Meristem"      = "#6fb3d2",
  "Pericycle"     = "#f5a65b",
  "Phloem"        = "#f8c8c8",
  "QC"            = "#fce4ec",
  "Root Cap"      = "#b088c8",
  "Xylem"         = "#c8e6c9"
)
hormone_colors <- c(
  "ABA"  = "#90c9a8",
  "IAA"  = "#6fb3d2",
  "eBL"  = "#f5a65b",
  "MeJA" = "#b088c8",
  "SA"   = "#ef6c6c"
)

ann_colors <- list(
  Cell_Type = cell_type_colors[names(cell_type_colors) %in% annot_col_df$Cell_Type],
  Hormone   = hormone_colors[names(hormone_colors) %in% annot_col_df$Hormone]
)

## ---------- 9. 画热图 ----------
hm_colors <- colorRampPalette(c("#3949ab", "#9fa8da", "#eeeeee", "#ef9a9a", "#b71c1c"))(100)

out_pdf <- file.path(PAPER_DIR, "M9_hub_genes_hormone_heatmap.pdf")
pdf(out_pdf, width = 14, height = 10)

pheatmap(
  mat,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale        = "none",
  color        = hm_colors,
  annotation_col = annot_col_df,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  show_colnames = FALSE,
  main = "M9 (HormoneProgram09) hub genes under five hormone treatments",
  treeheight_row = 0,
  treeheight_col = 0,
  na_col = "#dddddd",
  fontsize_row = 8,
  annotation_legend = TRUE
)
dev.off()

cat("热图已保存：", out_pdf, "\n")

## ---------- 10. 导出筛选明细 ----------
fwrite(sub, file.path(PAPER_DIR, "M9_hub_genes_hormone_expression.csv"))
cat("筛选明细表已导出：", file.path(PAPER_DIR, "M9_hub_genes_hormone_expression.csv"), "\n")

# 读取 Top50 Hub 基因表
# 注意路径中只有一个 hdWGCNA_Hormones5_CK2_Endo3_v1
top50 <- read.csv("/data01/liuxy/Project/CK2_added_v3_Results/hdWGCNA_Hormones5_CK2_Endo3_v1/Hormones5_Endo3_v1_module_Top50_kME.csv")

# 筛选 M9 (HormoneProgram09) 
m9_top50 <- top50[top50$Program == "HormoneProgram09", ]
print(m9_top50, row.names = FALSE)
##
mods <- read.csv("/data01/liuxy/Project/CK2_added_v3_Results/hdWGCNA_Hormones5_CK2_Endo3_v1/Hormones5_Endo3_v1_network_modules.csv")

target_genes <- c("Glyma.11G179600", "Glyma.16G103500")
result <- mods[mods$gene_name %in% target_genes, c("gene_name", "module", "original_module", "color")]

print(result, row.names = FALSE)
# -*- coding: UTF-8 -*-
# 独立完整绘图脚本：六处理拟时山峦图 + 同轴分组热图 + 对应缓存GO。
# source('soybean_endodermis_sixgroup_ridges_v6.txt', encoding='UTF-8')
# 只需R自带包；读取v3/v5已保存CSV。不读取Seurat，不重新降维/拟时/聚类/GO。
# 上方六条密度：CK、ABA、eBL、IAA、MeJA、SA（每组内有效生物学重复等权）。
# 下方按处理分成六块热图；每块沿用原模块、行顺序、row Z-score。
# 同组/同路径模块的GO才可列在旁边；缺少或无显著GO会明确注明。
# 不将不同处理的独立模块当作一套共同基因，不将不同分量/路径混为一条时间轴。

RUN_DIR <- '/data01/liuxy/Project/CK2_added_v3_Results/Endodermis_monocle3_originalUMAP_grouped_v3/20261119_163425'
CACHE_VERSION <- 'auto' # auto:存在非空v5模块表则用v5，否则v3；也可显式'v5'或'v3'
SELECT_PATHS <- NULL    # NULL:绘制全部已保存路径；或c('C01_P01')
HEATMAP_GROUPS <- c('CK','ABA','eBL','IAA','MeJA','SA')
GROUP_ORDER <- c('CK','ABA','eBL','IAA','MeJA','SA')
GROUP_COLORS <- c(CK='#4D4D4D',ABA='#D55E00',eBL='#009E73',IAA='#0072B2',MeJA='#CC79A7',SA='#E69F00')
MAX_ROWS_PER_GROUP <- 150L # 仅压缩绘图行数，原矩阵、模块和GO统计完全不改
MIN_CELLS_PER_REPLICATE_KDE <- 3L
KDE_GRID <- 512L
BANDWIDTH <- NULL       # NULL:每条路径统一带宽，六组共用；不要各组单独调带宽
DRAW_REPLICATE_LINES <- TRUE
GO_FDR <- 0.05
GO_TERMS_PER_MODULE <- 2L
Z_LIMIT <- 2.5
PNG_DPI <- 180L
FIGURE_WIDTH <- 15

read_csv_safe <- function(f,required=TRUE) {
  if(!file.exists(f)) {
    if(required)stop('Missing saved file: ',f)
    return(NULL)
  }
  read.csv(f,stringsAsFactors=FALSE,check.names=FALSE,fileEncoding='UTF-8',na.strings=c('','NA'))
}
integral_trap <- function(x,y)sum(diff(x)*(head(y,-1)+tail(y,-1))/2)

# 对有限支持区间做反射KDE，减少边界处的密度流失，再归一化到面积1。
# 不改变拟时序；所有组/重复共用同一带宽和网格，峰高可在这张图内比较。
reflected_kde <- function(pt,limits,bw,n=KDE_GRID) {
  if(length(pt)<MIN_CELLS_PER_REPLICATE_KDE||length(unique(pt))<2L)return(NULL)
  xx<-seq(limits[1],limits[2],length.out=n)
  aug<-c(pt,2*limits[1]-pt,2*limits[2]-pt)
  fit<-density(aug,bw=bw,from=limits[1],to=limits[2],n=n,cut=0)
  yy<-pmax(0,approx(fit$x,fit$y,xout=xx,rule=1)$y)
  area<-integral_trap(xx,yy)
  if(!is.finite(area)||area<=0)return(NULL)
  data.frame(pseudotime=xx,density=yy/area)
}

prepare_densities <- function(z,all_cells,limits) {
  bw<-BANDWIDTH
  if(is.null(bw))bw<-max(bw.nrd0(z$local_pseudotime),diff(limits)/80)
  if(length(bw)!=1L||!is.finite(bw)||bw<=0)stop('Invalid KDE bandwidth.')
  lines<-list();audit<-list();curves<-list()
  for(tg in GROUP_ORDER) {
    zz<-z[z$treatment==tg,,drop=FALSE]
    ss<-unique(all_cells$sample[all_cells$treatment==tg])
    fits<-list()
    for(s in ss) {
      v<-zz$local_pseudotime[zz$sample==s]
      fit<-reflected_kde(v,limits,bw)
      status<-if(length(v)<MIN_CELLS_PER_REPLICATE_KDE)'too_few_cells' else
        if(length(unique(v))<2L)'all_pseudotimes_tied' else 'used'
      audit[[paste(tg,s,sep='__')]]<-data.frame(treatment=tg,sample=s,n_path=length(v),
                                                n_total=sum(all_cells$treatment==tg & all_cells$sample==s),kde_status=status,bandwidth=bw)
      if(!is.null(fit)) {
        fits[[s]]<-fit$density
        lines[[paste(tg,s,sep='__')]]<-data.frame(treatment=tg,sample=s,fit)
      }
    }
    if(length(fits)) {
      yy<-rowMeans(do.call(cbind,fits))
      curves[[tg]]<-data.frame(treatment=tg,pseudotime=seq(limits[1],limits[2],length.out=KDE_GRID),
                               density=yy,n_cells=nrow(zz),n_replicates_used=length(fits),n_replicates_expected=length(ss))
    }
  }
  list(mean=if(length(curves))do.call(rbind,curves) else data.frame(),
       replicates=if(length(lines))do.call(rbind,lines) else data.frame(),
       audit=do.call(rbind,audit),bandwidth=bw)
}

sample_heatmap_rows <- function(modules,nmax) {
  if(length(modules)<=nmax)return(seq_along(modules))
  selected<-unlist(lapply(split(seq_along(modules),modules),function(ii) {
    n<-max(1L,floor(nmax*length(ii)/length(modules)))
    ii[unique(round(seq(1,length(ii),length.out=min(length(ii),n))))]
  }),use.names=FALSE)
  sort(selected)
}

load_heatmap <- function(root,folder,tg,path) {
  stem<-file.path(root,folder,tg,paste0('B_heatmap_',tg,'_',path))
  f<-paste0(stem,'_rowZ_matrix.csv')
  if(!file.exists(f))return(list(ok=FALSE,reason='No saved supported heatmap for this group/path.'))
  tab<-read_csv_safe(f)
  if(!nrow(tab)||!all(c('treatment','path','feature_id','module') %in% names(tab)))
    stop('Invalid heatmap cache: ',f)
  if(anyNA(tab$treatment)||anyNA(tab$path)||any(tab$treatment!=tg)|any(tab$path!=path))
    stop('Heatmap group/path mismatch: ',f)
  nc<-setdiff(names(tab),c('treatment','path','feature_id','module'))
  mat<-as.matrix(tab[,nc,drop=FALSE]);storage.mode(mat)<-'double'
  ax<-read_csv_safe(paste0(stem,'_column_pseudotime.csv'))
  pt<-ax$local_pseudotime
  if(ncol(mat)!=length(pt)||length(pt)<2L||any(!is.finite(pt))||any(diff(pt)<=0)||any(!is.finite(mat)))
    stop('Invalid cached matrix or pseudotime axis: ',f)
  if(anyNA(tab$module)||anyDuplicated(tab$feature_id))stop('Invalid gene/module labels: ',f)
  ix<-sample_heatmap_rows(as.character(tab$module),MAX_ROWS_PER_GROUP)
  list(ok=TRUE,z=mat[ix,,drop=FALSE],pt=pt,modules=as.character(tab$module[ix]),
       genes=tab$feature_id[ix],n_saved=nrow(mat),file=f)
}

draw_ridges <- function(den,limits,title_text,show_axis=TRUE) {
  par(mar=c(if(show_axis)2.8 else .8,4.8,2.6,1),xaxs='i',yaxs='i')
  plot(NA,xlim=limits,ylim=c(.5,7.3),axes=FALSE,xlab='',ylab='',main=title_text,cex.main=1)
  base<-setNames(rev(seq_along(GROUP_ORDER)),GROUP_ORDER)
  axis(2,at=base,labels=names(base),las=1,tick=FALSE,cex.axis=.9)
  if(show_axis) {axis(1,cex.axis=.8);mtext('Saved local pseudotime',side=1,line=1.8,cex=.85)}
  # 所有山峦共享一个密度高度系数，绝不分别归一到同样峰高。
  maximum<-if(nrow(den$mean))max(den$mean$density) else 1
  if(DRAW_REPLICATE_LINES&&nrow(den$replicates))maximum<-max(maximum,den$replicates$density)
  for(tg in GROUP_ORDER) {
    y0<-base[tg];abline(h=y0,col='grey90',lwd=.5)
    z<-if(nrow(den$mean))den$mean[den$mean$treatment==tg,,drop=FALSE] else data.frame()
    a<-den$audit[den$audit$treatment==tg,,drop=FALSE]
    if(nrow(z)) {
      yy<-y0+.78*z$density/maximum
      polygon(c(z$pseudotime,rev(z$pseudotime)),c(yy,rep(y0,nrow(z))),
              col=adjustcolor(GROUP_COLORS[tg],alpha.f=.55),border=NA)
      lines(z$pseudotime,yy,col=GROUP_COLORS[tg],lwd=1.2)
      if(DRAW_REPLICATE_LINES) {
        rr<-den$replicates[den$replicates$treatment==tg,,drop=FALSE]
        for(s in unique(rr$sample)) {
          one<-rr[rr$sample==s,,drop=FALSE]
          lines(one$pseudotime,y0+.78*one$density/maximum,col=adjustcolor(GROUP_COLORS[tg],alpha.f=.8),lty=3,lwd=.65)
        }
      }
    }
    label<-paste0('n=',sum(a$n_path),'; reps=',sum(a$kde_status=='used'),'/',nrow(a))
    if(!nrow(z))label<-paste(label,'(no KDE)')
    text(limits[2],y0+.83,labels=label,adj=c(1,0),cex=.65,col='grey30')
  }
}

draw_heatmap <- function(h,tg,path,limits,last=FALSE) {
  par(mar=c(if(last)3 else .4,4.8,2.2,1),xaxs='i',yaxs='i')
  if(!h$ok) {
    plot(NA,xlim=limits,ylim=c(0,1),axes=FALSE,xlab='',ylab='',main=paste(tg,'|',path),cex.main=.9)
    text(mean(limits),.5,h$reason,cex=.8,col='grey40')
  } else {
    nr<-nrow(h$z);nc<-ncol(h$z)
    plot(NA,xlim=limits,ylim=c(0,nr),axes=FALSE,xlab='',ylab='',
         main=paste0(tg,' | ',path,' | ',nr,' displayed / ',h$n_saved,' cached genes'),cex.main=.85)
    # 按实际拟时坐标放置各列，不把稀疏/不等间距箱误画成等时间宽度。
    edges<-c(h$pt[1],(head(h$pt,-1)+tail(h$pt,-1))/2,tail(h$pt,1))
    if(min(edges)<limits[1]-1e-6||max(edges)>limits[2]+1e-6)stop('Heatmap exceeds path pseudotime support.')
    pal<-colorRampPalette(c('#5B8E3E','#F7F7F7','#B12A90'))(201)
    zz<-pmax(-Z_LIMIT,pmin(Z_LIMIT,h$z))
    ci<-pmin(201L,pmax(1L,1L+floor((zz+Z_LIMIT)/(2*Z_LIMIT)*200)))
    rect(rep(head(edges,-1),each=nr),rep(nr-seq_len(nr),nc),
         rep(tail(edges,-1),each=nr),rep(nr-seq_len(nr)+1,nc),col=pal[as.vector(ci)],border=NA)
    r<-rle(h$modules);end<-cumsum(r$lengths);beg<-c(0,head(end,-1))
    mc<-setNames(hcl.colors(length(unique(h$modules)),'Dark 3'),unique(h$modules))
    dx<-diff(limits)*.015
    for(i in seq_along(r$values)) {
      if(i>1)abline(h=nr-beg[i],col='white',lwd=2)
      rect(limits[1]-dx,nr-end[i],limits[1]-.003*diff(limits),nr-beg[i],
           col=mc[r$values[i]],border=NA,xpd=NA)
      text(limits[1]-dx*1.5,nr-(beg[i]+end[i])/2,r$values[i],adj=1,cex=.65,xpd=NA)
    }
    box(col='grey75')
  }
  if(last){axis(1,cex.axis=.8);mtext('Saved local pseudotime (same axis as density)',1,line=1.8,cex=.8)}
}

draw_go_text <- function(go,tg,path,available,mods) {
  par(mar=c(.5,.5,2.2,.5))
  plot.new();plot.window(xlim=c(0,1),ylim=c(0,1))
  title(main=paste(tg,'|',path,'| GO BP'),cex.main=.85,adj=0)
  if(!available) {
    text(0,.85,'No cached GO table.\nRun the GO repair script first.',adj=c(0,1),cex=.8);return(invisible(NULL))
  }
  hit<-go[go$treatment==tg & go$path==path & is.finite(go$FDR_across_all_group_modules)&
            go$FDR_across_all_group_modules<GO_FDR,,drop=FALSE]
  if(!nrow(hit)) {
    text(0,.85,'No cached significant GO terms.\nConsult the matching GO audit for\nuntested/failed/non-significant modules.',
         adj=c(0,1),cex=.75);return(invisible(NULL))
  }
  txt<-character()
  for(m in unique(mods)) {
    x<-hit[as.character(hit$module)==m,,drop=FALSE]
    x<-head(x[order(x$FDR_across_all_group_modules,-x$Count),,drop=FALSE],GO_TERMS_PER_MODULE)
    txt<-c(txt,paste0(m,':'))
    if(!nrow(x))txt<-c(txt,'  No significant term') else for(j in seq_len(nrow(x))) {
      s<-paste0(x$Description[j],' (FDR ',formatC(x$FDR_across_all_group_modules[j],digits=2,format='g'),')')
      txt<-c(txt,paste0('  ',strwrap(s,width=43)))
    }
  }
  if(!length(txt))txt<-'No supported heatmap module to annotate.'
  text(0,seq(.94,.04,length.out=max(length(txt),2))[seq_along(txt)],txt,adj=c(0,.5),
       cex=min(.77,13/max(length(txt),1)))
}

save_devices <- function(stem,width,height,fun) {
  pdf(paste0(stem,'.pdf'),width=width,height=height,useDingbats=FALSE)
  tryCatch(fun(),finally=dev.off())
  png(paste0(stem,'.png'),width=round(width*PNG_DPI),height=round(height*PNG_DPI),res=PNG_DPI,
      type=if(capabilities('cairo'))'cairo' else getOption('bitmapType','Xlib'))
  tryCatch(fun(),finally=dev.off())
}

main_ridges <- function() {
  if(!dir.exists(RUN_DIR))stop('Saved output directory does not exist: ',RUN_DIR)
  root<-normalizePath(RUN_DIR,winslash='/',mustWork=TRUE)
  d<-read_csv_safe(file.path(root,'13_cell_metadata_pseudotime_states.csv'))
  req<-c('cell','treatment','sample','component','nearest_node','local_pseudotime')
  if(!all(req %in% names(d)))stop('Cell metadata missing: ',paste(setdiff(req,names(d)),collapse=', '))
  if(anyDuplicated(d$cell)||anyNA(d$sample)||anyNA(d$treatment))stop('Invalid cell/sample/treatment labels.')
  d$treatment[d$treatment=='Control']<-'CK';d$treatment[d$treatment=='MEJA']<-'MeJA'
  if(!all(d$treatment %in% GROUP_ORDER))stop('Unexpected treatments in cached metadata.')
  if(!all(HEATMAP_GROUPS %in% GROUP_ORDER))stop('HEATMAP_GROUPS must use the six treatment names.')
  v5mod<-read_csv_safe(file.path(root,'v5_29_group_path_gene_modules.csv'),FALSE)
  ver<-CACHE_VERSION
  if(ver=='auto')ver<-if(!is.null(v5mod)&&nrow(v5mod)>0)'v5' else 'v3'
  if(!ver %in% c('v3','v5'))stop('CACHE_VERSION must be auto, v3 or v5.')
  folder<-if(ver=='v5')'B_heatmaps_paper_v5' else 'B_heatmaps'
  gofile<-if(ver=='v5')'v5_38_GO_BP_all_groups_all_tested.csv' else '38_GO_BP_all_groups_all_tested.csv'
  go<-read_csv_safe(file.path(root,gofile),FALSE)
  go_ok<-!is.null(go)&&all(c('treatment','path','module','Description','Count','FDR_across_all_group_modules') %in% names(go))
  paths<-read_csv_safe(file.path(root,'26b_root_to_leaf_paths.csv'))
  pn<-read_csv_safe(file.path(root,'26c_path_node_membership.csv'))
  if(anyDuplicated(paths$path))stop('Duplicate path definition.')
  if(!is.null(SELECT_PATHS)) {
    if(!all(SELECT_PATHS %in% paths$path))stop('SELECT_PATHS not in saved path table.')
    paths<-paths[paths$path %in% SELECT_PATHS,,drop=FALSE]
  }
  if(!nrow(paths))stop('No saved paths selected.')
  report<-list()
  for(i in seq_len(nrow(paths))) {
    path<-paths$path[i];component<-paths$component[i]
    cat('Drawing ',path,' from ',ver,' cache...\n',sep='')
    result<-tryCatch({
      z<-d[d$component==component & d$nearest_node %in% pn$node[pn$path==path]&is.finite(d$local_pseudotime),,drop=FALSE]
      if(nrow(z)<3||length(unique(z$local_pseudotime))<2)stop('Insufficient resolved pseudotime for this path.')
      limits<-range(z$local_pseudotime)
      den<-prepare_densities(z,d,limits)
      stem<-file.path(root,paste0('v6_Figure1_six_hormone_groups_',path))
      write.csv(den$audit,paste0(stem,'_density_audit.csv'),row.names=FALSE)
      if(nrow(den$mean))write.csv(den$mean,paste0(stem,'_density_curves.csv'),row.names=FALSE)
      if(nrow(den$replicates))write.csv(den$replicates,paste0(stem,'_replicate_curves.csv'),row.names=FALSE)
      hs<-setNames(lapply(HEATMAP_GROUPS,function(tg)load_heatmap(root,folder,tg,path)),HEATMAP_GROUPS)
      # v5有模块表但该组没有通过筛选时，不自动偷换成v3旧模块。
      if(ver=='v5')for(tg in names(hs)) {
        if(!any(v5mod$treatment==tg & v5mod$path==path))hs[[tg]]<-list(ok=FALSE,reason='No supported v5 module for this group/path.')
      }
      hm_audit<-do.call(rbind,lapply(names(hs),function(tg) {
        h<-hs[[tg]];data.frame(treatment=tg,path=path,cache=ver,
                               status=if(h$ok)'plotted' else h$reason,
                               n_displayed=if(h$ok)nrow(h$z) else 0,n_saved=if(h$ok)h$n_saved else 0)
      }))
      write.csv(hm_audit,paste0(stem,'_heatmap_audit.csv'),row.names=FALSE)
      title<-paste('Endodermis |',path,'| six treatment distributions')
      save_devices(paste0(stem,'_ridges_only'),10,6,function() {
        par(oma=c(1,0,0,0));draw_ridges(den,limits,title)
        mtext('Density area = 1 per treatment; equal weight of usable replicates; dotted lines = replicates.',1,outer=TRUE,cex=.65)
      })
      htotal<-5+2.25*length(hs)
      save_devices(stem,FIGURE_WIDTH,htotal,function() {
        ng<-length(hs)
        layout(rbind(c(1,2),cbind(seq(3,by=2,length.out=ng),seq(4,by=2,length.out=ng))),
               widths=c(2.05,1),heights=c(3.5,rep(2.25,ng)))
        par(oma=c(2,0,1,0))
        draw_ridges(den,limits,title)
        par(mar=c(1,.5,2,.5));plot.new();plot.window(xlim=c(0,1),ylim=c(0,1))
        text(0,.97,paste('Treatment density +',ver,'cached heatmaps'),adj=c(0,1),font=2,cex=.9)
        legend('topleft',inset=c(0,.10),legend=GROUP_ORDER,fill=GROUP_COLORS,bty='n',cex=.85,ncol=2)
        text(0,.54,paste(strwrap(paste('One common pseudotime axis within this saved path.',
                                       'Density area is one per treatment; peak height is not cell count or gene expression.',
                                       'Heatmaps use local group modules. Matching M numbers across groups do not mean the same genes.'),48),collapse='\n'),
             adj=c(0,1),cex=.78)
        pal<-colorRampPalette(c('#5B8E3E','#F7F7F7','#B12A90'))(100)
        rect(seq(0,.79,length.out=100),.12,seq(.008,.798,length.out=100),.17,col=pal,border=NA)
        text(c(0,.4,.8),.08,c(-Z_LIMIT,'row Z-score',Z_LIMIT),cex=.8)
        for(j in seq_along(hs)) {
          tg<-names(hs)[j];h<-hs[[j]]
          draw_heatmap(h,tg,path,limits,last=j==length(hs))
          draw_go_text(go,tg,path,go_ok,if(h$ok)unique(h$modules) else character())
        }
        mtext('Exploratory root; pseudotime is not elapsed treatment time. Blank heatmap regions are outside cached support.',
              side=1,outer=TRUE,line=.4,cex=.7)
      })
      'completed'
    },error=function(e)paste('failed:',conditionMessage(e)))
    report[[path]]<-data.frame(path=path,component=component,cache=ver,status=result)
    write.csv(do.call(rbind,report),file.path(root,'v6_sixgroup_ridge_status.csv'),row.names=FALSE)
    cat(path,': ',result,'\n',sep='')
  }
  writeLines(c(paste('Original output directory:',root),paste('Cache version:',ver),
               'Only CSV caches were read. No embeddings, trajectory, graph_test, expression modules or GO were recomputed.',
               'Every path is plotted separately; shared trunk cells can occur in multiple paths.',
               'Ridges: common bandwidth, bounded reflected KDE, equal averaging of usable replicate densities.',
               'Per-group density integrates to one. Cell counts and replicate availability are explicitly labelled.',
               'Heatmaps: original saved row order and Z-scores. Groups are separate blocks, not a shared module gene set.',
               'GO: only the matching cache version, group and path. Blank/unavailable GO is never replaced by another group.',
               'Plots and numerical curves are saved here; examine v6_sixgroup_ridge_status.csv for completion.'),
             file.path(root,'v6_sixgroup_ridge_README.txt'))
  if(any(vapply(report,function(x)x$status!='completed',logical(1))))
    stop('Some paths failed; see v6_sixgroup_ridge_status.csv. Completed figures remain saved.')
  invisible(root)
}

main_ridges()



# -*- coding: UTF-8 -*-
# v4_dotplot：GO改为气泡点图版本，其余逻辑和原版v4完全一致
# 运行：source('soybean_endodermis_resume_GO_replot_v4_dotplot.txt', encoding='UTF-8')
# 或：Rscript soybean_endodermis_resume_GO_replot_v4_dotplot.txt
# 改动点：draw_go() 从geom_col条形图 → 气泡dot plot（size=Count，color=-log10(FDR)）
# 不读原始Seurat对象，不重新降维/聚类/建图/定根/graph_test/平滑/聚类模块。
# 原输出目录不变，重绘覆盖同名图，已保存的CDS/CSV分析输入不改。

RUN_DIR <- '/data01/liuxy/Project/CK2_added_v3_Results/Endodermis_monocle3_originalUMAP_grouped_v3/20261119_163425'
REDRAW_TRAJECTORY <- TRUE
REDRAW_HEATMAPS <- TRUE
REDRAW_SUMMARIES <- TRUE
RESUME_GO <- TRUE
ORGDB_PACKAGE <- 'org.G.max.eg.db'
KEYTYPE_OVERRIDE <- 'GID'
PREFERRED_DB_VERSION <- 'v4'
ORGDB_SEARCH_DIR <- '/home/liuxy'
ORGDB_TARBALL <- NULL
GENE_ID_MAP_TSV <- NULL                 # 可选：feature_id、GID两列正式一对一映射
MIN_DB_MATCH_RATE <- 0.50
MIN_DETECTED_FRACTION <- 0.10           # 与v3一致，用于恢复各路径的GO检测背景
GO_MIN_GENES <- 5L
GO_MIN_SET_SIZE <- 10L
GO_MAX_SET_SIZE <- 500L
GO_FDR_CUTOFF <- 0.05                   # 不为出图而放宽显著性阈值
GO_TOP_TERMS <- 6L
PT_PALETTE <- 'plasma'                 # 紫→红橙→黄，替换不明显的单一蓝色
PT_COLOR_TRANSFORM <- 'sqrt'           # 只改变颜色映射，拟时序数值不改；图例标原始值
# 所有处理使用同一全局色标，绝不按每组各自归一化。另存linear线性色标版本作对照。
POINT_SIZE <- 0.55
BACKGROUND_ALPHA <- 0.30
SEED <- 20261119L

# __ID_HELPERS__

# 可测试的包加载分支：已加载 > 已安装 > 从本地tar.gz安装。
orgdb_load_route <- function(is_loaded,is_available) {
  if(is_loaded) 'reuse_loaded_namespace' else if(is_available) 'load_installed_package' else 'install_local_package'
}

# 从已保存的counts与graph_test恢复相同的检测背景；不重新做graph_test。
restore_universe <- function(count_matrix,cell_ids,tested_features,min_fraction=0.10) {
  if(!length(cell_ids)||anyDuplicated(cell_ids)) stop('路径细胞列表为空或含重复barcode。')
  if(!all(cell_ids %in% colnames(count_matrix))) stop('保存的路径细胞不全在CDS中，拒绝混用运行结果。')
  cc <- count_matrix[,cell_ids,drop=FALSE]
  n_detect <- if(inherits(cc,'Matrix')) Matrix::rowSums(cc>0) else rowSums(cc>0)
  detectable <- rownames(count_matrix)[n_detect>=max(10L,ceiling(min_fraction*length(cell_ids)))]
  intersect(tested_features,detectable)
}

main_resume <- function() {
  if(!dir.exists(RUN_DIR)) stop('找不到已完成的输出目录：',RUN_DIR)
  run_dir <- normalizePath(RUN_DIR,winslash='/',mustWork=TRUE)
  logcon <- file(file.path(run_dir,'resume_replot_v4_dotplot.log'),'wt',encoding='UTF-8')
  sink(logcon,split=TRUE);on.exit({sink();close(logcon)},add=TRUE)
  on.exit(writeLines(capture.output(sessionInfo()),file.path(run_dir,'sessionInfo_v4_dotplot.txt')),add=TRUE)
  cat('Using saved run, same output directory: ',run_dir,'\n',sep='')
  need <- function(pkgs) {
    bad <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
    if(length(bad)) stop('缺少依赖：',paste(bad,collapse=', '))
  }
  need('ggplot2');suppressPackageStartupMessages(library(ggplot2))
  set.seed(SEED)
  groups <- c('CK','ABA','eBL','IAA','MeJA','SA')
  group_colors <- c(CK='#4D4D4D',ABA='#D55E00',eBL='#009E73',IAA='#0072B2',MeJA='#CC79A7',SA='#E69F00')
  read_table <- function(name,required=TRUE) {
    f <- file.path(run_dir,name)
    if(!file.exists(f)) {
      if(required) stop('缺少已保存文件：',f)
      return(NULL)
    }
    read.csv(f,check.names=FALSE,stringsAsFactors=FALSE,na.strings=c('','NA'),fileEncoding='UTF-8')
  }
  csv <- function(x,name) {
    path <- file.path(run_dir,paste0(name,'.csv'))
    dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
    write.csv(x,path,row.names=FALSE,na='',fileEncoding='UTF-8')
  }
  theme_pub <- theme_classic(base_size=12)+theme(
    plot.title=element_text(size=15,face='bold'),plot.subtitle=element_text(size=10),
    plot.caption=element_text(size=8),strip.background=element_blank(),
    strip.text=element_text(size=13,face='bold'),plot.margin=margin(12,16,12,12))
  save_plot <- function(p,stem,w=11,h=7) {
    path <- file.path(run_dir,stem);dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
    ggsave(paste0(path,'.pdf'),plot=p,width=w,height=h,limitsize=FALSE,bg='white')
    ggsave(paste0(path,'.png'),plot=p,width=w,height=h,dpi=300,limitsize=FALSE,bg='white')
  }
  placeholder <- function(title,reason,stem) {
    p <- ggplot()+annotate('text',x=0,y=0,label=paste(strwrap(reason,78),collapse='\n'),size=4)+
      xlim(-1,1)+ylim(-1,1)+theme_void()+labs(title=title)+
      theme(plot.title=element_text(face='bold',hjust=0.5,size=15))
    save_plot(p,stem,11,5)
  }
  stages <- list()
  run_stage <- function(name,fun) {
    cat('\n--- ',name,' ---\n',sep='')
    out <- tryCatch({fun();list(status='completed',detail='')},error=function(e)
      list(status='failed',detail=conditionMessage(e)))
    stages[[name]] <<- data.frame(stage=name,status=out$status,detail=out$detail)
    csv(do.call(rbind,stages),'resume_replot_v4_dotplot_status')
    cat(name,': ',out$status,' ',out$detail,'\n',sep='')
    invisible(out)
  }
  d <- read_table('13_cell_metadata_pseudotime_states.csv')
  for(nm in c('cell','x','y','treatment','component','pseudotime','local_pseudotime'))
    if(!nm %in% names(d)) stop('保存的细胞表缺少列：',nm)
  d$treatment <- factor(d$treatment,levels=groups)
  if(anyNA(d$treatment)||anyDuplicated(d$cell)) stop('细胞表有未知处理或重复barcode。')
  edges <- read_table('09_graph_edges.csv')
  nodes <- read_table('08_principal_nodes.csv')
  roots <- read_table('11_root_and_interval_definition.csv')
  local_roots <- read_table('11b_roots_for_separate_component_heatmaps.csv',required=FALSE)
  target <- unique(roots$component)
  if(length(target)!=1L) stop('主分量定义不唯一。')
  d$pt_plot <- ifelse(d$component==target & is.finite(d$pseudotime),d$pseudotime,NA_real_)
  edges$component <- nodes$component[match(edges$from,nodes$node)]
  if(anyNA(edges$component)) stop('保存的边表与节点表不匹配。')
  color_scale <- function(limits,transform=PT_COLOR_TRANSFORM,label='Pseudotime') {
    if(any(!is.finite(limits))||diff(limits)<=0) stop('拟时序颜色范围无效。')
    scale_color_viridis_c(option=PT_PALETTE,begin=0,end=0.97,limits=limits,trans=transform,
                          name=label,na.value='#BDBDBD',guide=guide_colorbar(barheight=grid::unit(42,'mm'),
                                                                             barwidth=grid::unit(4,'mm')))
  }
  
  # ==================== 1. 仅用已保存坐标与拟时序重绘 ====================
  if(REDRAW_TRAJECTORY) run_stage('trajectory_replot',function() {
    finite <- d$pt_plot[is.finite(d$pt_plot)]
    if(!length(finite)||min(finite)<0) stop('主拟时序为空或存在负值。')
    limits <- range(finite)
    csv(data.frame(palette=PT_PALETTE,transform=PT_COLOR_TRANSFORM,
                   min_pseudotime=limits[1],max_pseudotime=limits[2],shared_across_groups=TRUE),'v4_pseudotime_color_settings')
    facet_edges <- do.call(rbind,lapply(groups,function(tg) {
      ee<-edges[edges$component %in% unique(d$component[d$treatment==tg]),,drop=FALSE]
      ee$treatment<-factor(rep(tg,nrow(ee)),levels=groups);ee
    }))
    make_a <- function(transform) ggplot(d,aes(x,y))+
      geom_point(data=d[,c('x','y')],aes(x,y),inherit.aes=FALSE,
                 color='grey75',alpha=BACKGROUND_ALPHA,size=0.35)+
      geom_segment(data=facet_edges,aes(x=x,y=y,xend=xend,yend=yend),inherit.aes=FALSE,
                   color='grey30',linewidth=0.35)+
      geom_point(aes(color=pt_plot),size=POINT_SIZE,alpha=0.95)+color_scale(limits,transform)+
      geom_point(data=nodes[nodes$node %in% roots$root_node,,drop=FALSE],aes(x,y),
                 inherit.aes=FALSE,shape=23,fill='white',color='black',size=3,stroke=0.9)+
      coord_equal()+theme_pub+labs(title='Endodermis | original UMAP',
                                   subtitle=paste0('Saved trajectory | primary component ',target,' | ',transform,' color scale'),
                                   x='Original UMAP 1',y='Original UMAP 2',
                                   caption='Same color scale for all groups; legend shows original pseudotime. Root direction remains exploratory unless validated.')
    p<-make_a(PT_COLOR_TRANSFORM)
    save_plot(p+facet_wrap(~treatment,ncol=3,drop=FALSE),'A_pseudotime_6groups',14,9)
    save_plot(p+facet_wrap(~treatment,ncol=1,drop=FALSE),'A_pseudotime_vertical',8,22)
    save_plot(make_a('identity')+facet_wrap(~treatment,ncol=3,drop=FALSE),
              'A_pseudotime_6groups_plasma_linear',14,9)
    # 每个处理单独保存，所有图仍使用同一坐标范围与同一色标。
    for(tg in groups) {
      z<-d[d$treatment==tg,,drop=FALSE]
      ee<-edges[edges$component %in% unique(z$component),,drop=FALSE]
      one<-ggplot(z,aes(x,y,color=pt_plot))+
        geom_point(data=d[,c('x','y')],aes(x,y),inherit.aes=FALSE,
                   color='grey75',alpha=BACKGROUND_ALPHA,size=0.35)+
        geom_segment(data=ee,aes(x=x,y=y,xend=xend,yend=yend),inherit.aes=FALSE,
                     color='grey30',linewidth=0.35)+geom_point(size=POINT_SIZE,alpha=0.95)+
        color_scale(limits)+coord_equal(xlim=range(d$x),ylim=range(d$y))+theme_pub+
        labs(title=paste(tg,'| Endodermis'),subtitle=paste(PT_COLOR_TRANSFORM,'color scale; saved pseudotime'),
             x='Original UMAP 1',y='Original UMAP 2')
      save_plot(one,paste0('A_pseudotime_',tg),8,7)
    }
    if(!is.null(local_roots)) for(cc in local_roots$component) {
      z<-d[d$component==cc & is.finite(d$local_pseudotime),,drop=FALSE]
      if(nrow(z)<10L||diff(range(z$local_pseudotime))<=0) next
      ee<-edges[edges$component==cc,,drop=FALSE]
      p<-ggplot(z,aes(x,y,color=local_pseudotime))+
        geom_segment(data=ee,aes(x=x,y=y,xend=xend,yend=yend),inherit.aes=FALSE,
                     color='grey30',linewidth=0.35)+geom_point(size=POINT_SIZE)+
        color_scale(range(z$local_pseudotime),label='Local PT')+
        facet_wrap(~treatment,ncol=3,drop=FALSE)+coord_equal()+theme_pub+
        labs(title=paste('Endodermis | component',cc),
             subtitle='Saved local pseudotime; compare treatments within this component only',
             x='Original UMAP 1',y='Original UMAP 2')
      save_plot(p,paste0('A_local_component_',cc,'_6groups'),13,8)
    }
  })
  
  # ==================== 2. 从Z-score矩阵重绘，不再平滑或聚类 ====================
  modules <- read_table('29_group_path_gene_modules.csv',required=FALSE)
  if(REDRAW_HEATMAPS) run_stage('heatmap_replot',function() {
    need('pheatmap')
    files<-list.files(file.path(run_dir,'B_heatmaps'),pattern='_rowZ_matrix[.]csv$',recursive=TRUE,full.names=TRUE)
    if(!length(files)) stop('没有保存的_rowZ_matrix.csv，请查看旧版heatmap_manifest.csv。')
    checks<-list()
    for(f in files) {
      status<-tryCatch({
        ztab<-read.csv(f,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding='UTF-8')
        tg<-unique(ztab$treatment);path<-unique(ztab$path)
        if(length(tg)!=1L||!tg %in% groups||length(path)!=1L) stop('矩阵组名/路径不唯一。')
        numeric_cols<-setdiff(names(ztab),c('treatment','path','feature_id','module'))
        z<-as.matrix(ztab[,numeric_cols,drop=FALSE]);storage.mode(z)<-'double'
        rownames(z)<-ztab$feature_id
        if(any(!is.finite(z))) stop('Z-score矩阵有非有限值。')
        stem<-sub('_rowZ_matrix[.]csv$','',f)
        axis<-read.csv(paste0(stem,'_column_pseudotime.csv'),stringsAsFactors=FALSE)
        if(nrow(axis)!=ncol(z)) stop('保存的热图列和拟时序轴长度不一致。')
        if(any(diff(axis$local_pseudotime)<=0)) stop('热图列轴不是递增拟时序。')
        mods<-as.character(ztab$module);lev<-unique(mods)
        gaps<-head(cumsum(rle(mods)$lengths),-1L)
        # 行与列直接沿用保存的聚类顺序；只更新注释与色板。
        z[]<-pmax(-2.5,pmin(2.5,z));colnames(z)<-sprintf('PT%03d',seq_len(ncol(z)))
        annrow<-data.frame(Module=factor(mods,levels=lev),row.names=rownames(z))
        anncol<-data.frame(Treatment=factor(rep(tg,ncol(z)),levels=groups),
                           Pseudotime=axis$local_pseudotime,row.names=colnames(z))
        tick<-unique(round(seq(1,ncol(z),length.out=5)));labels<-rep('',ncol(z))
        labels[tick]<-format(round(axis$local_pseudotime[tick],2),trim=TRUE)
        ph<-pheatmap::pheatmap(z,cluster_rows=FALSE,cluster_cols=FALSE,gaps_row=gaps,
                               annotation_row=annrow,annotation_col=anncol,
                               annotation_colors=list(Treatment=group_colors,
                                                      Module=setNames(colorRampPalette(c('#E69F00','#009E73','#CC79A7'))(length(lev)),lev),
                                                      Pseudotime=c('#0D0887','#CC4778','#F0F921')),
                               color=colorRampPalette(c('#2166AC','#F7F7F7','#B2182B'))(100),
                               breaks=seq(-2.5,2.5,length.out=101),legend_breaks=c(-2,-1,0,1,2),
                               labels_col=labels,angle_col='0',show_rownames=FALSE,show_colnames=TRUE,
                               annotation_names_row=FALSE,border_color=NA,fontsize=10,
                               main=paste0(tg,' | ',path,'\nSaved clustered trends | row Z-score'),silent=TRUE)
        grDevices::pdf(paste0(stem,'.pdf'),width=9,height=9)
        tryCatch({grid::grid.newpage();grid::grid.draw(ph$gtable)},finally=grDevices::dev.off())
        grDevices::png(paste0(stem,'.png'),width=2700,height=2700,res=300,
                       type=if(capabilities('cairo'))'cairo' else getOption('bitmapType','Xlib'))
        tryCatch({grid::grid.newpage();grid::grid.draw(ph$gtable)},finally=grDevices::dev.off())
        'completed'
      },error=function(e)paste('failed:',conditionMessage(e)))
      checks[[f]]<-data.frame(matrix_file=basename(f),status=status)
    }
    csv(do.call(rbind,checks),'v4_heatmap_redraw_audit')
    if(any(vapply(checks,function(x)x$status!='completed',logical(1))))
      stop('部分缓存热图重绘失败，详见v4_heatmap_redraw_audit.csv；GO仍将继续。')
  })
  
  if(REDRAW_SUMMARIES) run_stage('replicate_summary_replot',function() {
    st<-read_table('21_per_replicate_interval_proportions.csv',FALSE)
    if(!is.null(st)) {
      st$treatment<-factor(st$treatment,levels=groups)
      p<-ggplot(st,aes(treatment,proportion_all,color=treatment))+
        geom_point(position=position_jitter(width=.08,height=0,seed=SEED),size=3)+
        facet_wrap(~interval,nrow=1)+scale_color_manual(values=group_colors)+theme_pub+
        theme(legend.position='none')+labs(x=NULL,y='Proportion per biological replicate')
      save_plot(p,'23_replicate_interval_proportions',14,5)
    }
    ss<-read_table('22_per_replicate_pseudotime_summary.csv',FALSE)
    if(!is.null(ss)) {
      ss$treatment<-factor(ss$treatment,levels=groups)
      p<-ggplot(ss,aes(treatment,median_pseudotime,color=treatment))+
        geom_point(position=position_jitter(width=.08,height=0,seed=SEED),size=3)+
        scale_color_manual(values=group_colors)+theme_pub+theme(legend.position='none')+
        labs(x=NULL,y='Saved median pseudotime per biological replicate')
      save_plot(p,'24_replicate_median_pseudotime')
    }
  })
  
  # ==================== 3. 修复OrgDb加载，补做此前未执行的GO ====================
  if(RESUME_GO) run_stage('GO_resume',function() {
    need(c('AnnotationDbi','clusterProfiler','GO.db','monocle3','SummarizedExperiment','Matrix'))
    if(is.null(modules)||!nrow(modules)) stop('没有保存的29_group_path_gene_modules.csv，不能凭空生成GO模块。')
    for(nm in c('job','treatment','path','component','module','feature_id'))
      if(!nm %in% names(modules)) stop('模块表缺少：',nm)
    loaded<-ORGDB_PACKAGE %in% loadedNamespaces()
    lib<-file.path(run_dir,'R_library')
    if(!loaded && dir.exists(lib)) .libPaths(unique(c(lib,.libPaths())))
    available<-if(loaded)TRUE else requireNamespace(ORGDB_PACKAGE,quietly=TRUE)
    route<-orgdb_load_route(loaded,available)
    cat('OrgDb route: ',route,'\n',sep='')
    if(route=='install_local_package') {
      tb<-ORGDB_TARBALL
      if(is.null(tb)) {
        candidates<-list.files(ORGDB_SEARCH_DIR,pattern='^org[.]G[.]max[.]eg[.]db.*[.]tar[.]gz$',full.names=TRUE)
        if(length(candidates)!=1L) stop('OrgDb未安装且安装包不唯一/不存在；请填写ORGDB_TARBALL。')
        tb<-candidates[1]
      }
      if(!file.exists(tb)) stop('安装包不存在：',tb)
      dir.create(lib,showWarnings=FALSE,recursive=TRUE)
      install.packages(tb,repos=NULL,type='source',lib=lib)
      .libPaths(unique(c(lib,.libPaths())))
      if(!requireNamespace(ORGDB_PACKAGE,quietly=TRUE)) stop('本地OrgDb安装后仍无法加载。')
    }
    # 已加载包直接复用，其命名空间路径才是实际生效的包路径。
    ns<-asNamespace(ORGDB_PACKAGE);db<-get(ORGDB_PACKAGE,envir=ns)
    if(!methods::is(db,'OrgDb')) stop('同名对象不是OrgDb。')
    if(!KEYTYPE_OVERRIDE %in% AnnotationDbi::keytypes(db)) stop('该OrgDb没有GID键类型。')
    dbkeys<-AnnotationDbi::keys(db,keytype=KEYTYPE_OVERRIDE)
    csv(AnnotationDbi::metadata(db),'32_orgdb_metadata')
    writeLines(c(paste('Load route:',route),paste('Actual namespace path:',getNamespaceInfo(ns,'path')),
                 paste('Loaded package version:',as.character(getNamespaceVersion(ns))),
                 paste('Keytype:',KEYTYPE_OVERRIDE),paste('Preferred ID marker:',PREFERRED_DB_VERSION)),
               file.path(run_dir,'34_annotation_provenance_v4_dotplot.txt'))
    cds_path<-file.path(run_dir,'endodermis_monocle3_ordered.rds')
    if(!file.exists(cds_path)) stop('缺少已保存的ordered CDS，无法恢复GO背景。')
    cat('Loading saved CDS only: ',cds_path,'\n',sep='')
    cds<-readRDS(cds_path)
    counts<-SummarizedExperiment::assay(cds,'counts')
    if(!setequal(colnames(counts),d$cell)) stop('保存的CDS与细胞结果表来自不同运行。')
    gt<-read_table('25_graph_test_pooled_candidate_screen.csv')
    if(!all(c('feature_id','q_value','morans_I') %in% names(gt))) stop('保存的graph_test表缺少关键列。')
    tested<-gt$feature_id[is.finite(gt$q_value)&is.finite(gt$morans_I)]
    idmap<-match_soybean_gid(rownames(counts),dbkeys,PREFERRED_DB_VERSION)
    if(!is.null(GENE_ID_MAP_TSV)) {
      mp<-read.delim(GENE_ID_MAP_TSV,stringsAsFactors=FALSE,check.names=FALSE)
      if(!all(c('feature_id','GID') %in% names(mp))||anyDuplicated(mp$feature_id)||anyDuplicated(mp$GID))
        stop('映射表必须是feature_id/GID两列的一对一正式映射。')
      idmap$GID<-mp$GID[match(idmap$feature_id,mp$feature_id)];idmap$match_status<-'user_mapping'
    }
    idmap$matched<-!is.na(idmap$GID)&idmap$GID %in% dbkeys
    csv(idmap,'35_feature_to_GID_audit');csv(as.data.frame(table(status=idmap$match_status)),'35b_ID_match_status_counts')
    match_rate<-mean(idmap$matched);cat('GID match rate: ',round(match_rate*100,2),'%\n',sep='')
    if(match_rate<MIN_DB_MATCH_RATE) stop('GID匹配率过低；已保存匹配审计，不能可靠进行GO。')
    map_ids<-function(ids)unique(idmap$GID[idmap$matched & idmap$feature_id %in% ids])
    results<-list();audit_rows<-list();jobs<-unique(modules[,c('job','treatment','path','component')])
    if(anyDuplicated(jobs$job)) stop('同一个job对应多个处理/路径，拒绝混用模块。')
    path_nodes<-read_table('26c_path_node_membership.csv',FALSE)
    for(i in seq_len(nrow(jobs))) {
      info<-jobs[i,,drop=FALSE];job<-info$job;tg<-info$treatment;pa<-info$path
      if(!tg %in% groups) stop('模块表包含未知处理组：',tg)
      cat('GO: ',tg,' | ',pa,'\n',sep='')
      stem<-file.path('B_heatmaps',tg,paste0('B_heatmap_',tg,'_',pa))
      cells<-read_table(paste0(stem,'_cells.csv'),FALSE)
      if(!is.null(cells)) {
        if(any(cells$treatment!=tg)||any(cells$path!=pa)) stop('路径细胞表标签不匹配：',job)
        cell_ids<-cells$cell
      } else {
        if(is.null(path_nodes)) stop('缺少路径细胞表及路径节点表：',job)
        pn<-path_nodes$node[path_nodes$path==pa]
        cell_ids<-d$cell[d$treatment==tg & d$nearest_node %in% pn & is.finite(d$local_pseudotime)]
      }
      if(!all(d$treatment[match(cell_ids,d$cell)]==tg)) stop('路径细胞混入其他处理：',job)
      universe_features<-restore_universe(counts,cell_ids,tested,MIN_DETECTED_FRACTION)
      background<-map_ids(universe_features)
      # 用长度一致的向量写空表，避免零匹配时data.frame报“行数不一致”。
      csv(data.frame(GID=background),file.path('C_GO',tg,paste0('GO_',tg,'_',pa,'_input_universe')))
      job_modules<-modules[modules$job==job,,drop=FALSE]
      for(mm in unique(job_modules$module)) {
        key<-paste(job,mm,sep='__');ids<-job_modules$feature_id[job_modules$module==mm]
        genes<-intersect(map_ids(ids),background)
        ar<-data.frame(job=job,treatment=tg,component=info$component,path=pa,module=mm,
                       n_features=length(ids),n_GID=length(genes),n_input_background=length(background),
                       n_tested_terms=0L,n_significant_global=0L,status='',detail='',stringsAsFactors=FALSE)
        csv(data.frame(GID=genes),file.path('C_GO',tg,paste0('GO_',tg,'_',pa,'_',mm,'_input_genes')))
        if(length(genes)<GO_MIN_GENES||length(background)<10L) {
          ar$status<-'insufficient_mapped_genes';audit_rows[[key]]<-ar;next
        }
        eg<-tryCatch(clusterProfiler::enrichGO(gene=genes,universe=background,OrgDb=db,
                                               keyType=KEYTYPE_OVERRIDE,ont='BP',pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1,
                                               minGSSize=GO_MIN_SET_SIZE,maxGSSize=GO_MAX_SET_SIZE,readable=FALSE),error=function(e)e)
        if(inherits(eg,'error')) {
          ar$status<-'enrichment_error';ar$detail<-conditionMessage(eg);audit_rows[[key]]<-ar;next
        }
        # 保存完整检验表，避免as.data.frame按对象内阈值再次隐藏行。
        rr<-if(is.null(eg))data.frame() else eg@result
        ar$status<-if(nrow(rr))'tested' else 'no_testable_GO_terms';ar$n_tested_terms<-nrow(rr)
        audit_rows[[key]]<-ar
        if(nrow(rr)) {
          rr$job<-job;rr$treatment<-tg;rr$component<-info$component;rr$path<-pa;rr$module<-mm
          results[[key]]<-rr
          saveRDS(eg,file.path(run_dir,'C_GO',tg,paste0('GO_',tg,'_',pa,'_',mm,'_enrichResult.rds')))
          csv(data.frame(GID=eg@universe),file.path('C_GO',tg,paste0('GO_',tg,'_',pa,'_',mm,'_effective_universe')))
        }
      }
      csv(do.call(rbind,audit_rows),'37_GO_group_path_module_audit')
    }
    if(length(results)) {
      allgo<-do.call(rbind,results)
      allgo$FDR_across_all_group_modules<-p.adjust(allgo$pvalue,method='BH')
      allgo$significant_global<-is.finite(allgo$FDR_across_all_group_modules)&
        allgo$FDR_across_all_group_modules<GO_FDR_CUTOFF
    } else allgo<-data.frame(ID=character(),Description=character(),pvalue=numeric(),p.adjust=numeric(),
                             Count=integer(),job=character(),treatment=character(),component=integer(),path=character(),
                             module=character(),FDR_across_all_group_modules=numeric(),significant_global=logical())
    csv(allgo,'38_GO_BP_all_groups_all_tested')
    csv(allgo[allgo$significant_global,,drop=FALSE],'39_GO_BP_all_groups_significant')
    final_audit<-do.call(rbind,audit_rows)
    for(i in seq_len(nrow(final_audit))) {
      ar<-final_audit[i,]
      n<-sum(allgo$job==ar$job & allgo$module==ar$module & allgo$significant_global)
      final_audit$n_significant_global[i]<-n
      if(final_audit$status[i]=='tested' && n==0L) final_audit$status[i]<-'tested_no_global_significance'
    }
    csv(final_audit,'37_GO_group_path_module_audit')
    
    # ========== 这里是【修改后的气泡绘图函数】，替换原来geom_col条形图 ==========
    draw_go<-function(z,tg,pa,title_suffix='',exploratory=FALSE) {
      top<-do.call(rbind,lapply(split(z,z$module),function(x)
        head(x[order(x$FDR_across_all_group_modules,x$pvalue,-x$Count),,drop=FALSE],GO_TOP_TERMS)))
      top$key<-paste(top$module,top$ID,sep='::');top$key<-factor(top$key,levels=rev(unique(top$key)))
      labels<-setNames(vapply(top$Description,function(s)paste(strwrap(s,48),collapse='\n'),character(1)),as.character(top$key))
      top$FDR_strength <- -log10(pmax(top$FDR_across_all_group_modules,1e-300))
      
      ggplot(top,aes(x=Count, y=key))+
        geom_point(aes(size=Count, color=FDR_strength))+
        scale_y_discrete(labels=labels)+
        scale_color_viridis_c(option='magma',direction=-1,name='-log10(global FDR)')+
        scale_size_continuous(range=c(2,9),name='Gene count')+
        facet_wrap(~module,ncol=1,scales='free_y')+
        theme_pub+
        labs(title=paste(tg,'|',pa,'| GO BP',title_suffix),
             subtitle=if(exploratory)'Top tested terms only; NOT significant at global FDR < 0.05' else
               paste('Global FDR <',GO_FDR_CUTOFF,'| soybean GID'),
             x='Mapped gene count',y=NULL)
    }
    # ==========================================================================
    
    for(tg in groups) {
      group_go<-allgo[allgo$treatment==tg,,drop=FALSE]
      csv(group_go,file.path('C_GO',tg,paste0('GO_BP_',tg,'_all_results')))
      jtg<-jobs[jobs$treatment==tg,,drop=FALSE]
      if(!nrow(jtg)) placeholder(paste(tg,'| GO BP'),'No saved heatmap modules for this group.',
                                 file.path('C_GO',tg,paste0('GO_BP_',tg,'_NO_MODULE')))
      for(i in seq_len(nrow(jtg))) {
        job<-jtg$job[i];pa<-jtg$path[i];stem<-file.path('C_GO',tg,paste0('GO_BP_',tg,'_',pa))
        jj<-group_go[group_go$job==job,,drop=FALSE];hit<-jj[jj$significant_global,,drop=FALSE]
        csv(jj,paste0(stem,'_all_results'));csv(hit,paste0(stem,'_significant'))
        if(nrow(hit)) {
          nn<-sum(pmin(table(hit$module),GO_TOP_TERMS))
          save_plot(draw_go(hit,tg,pa),stem,12,max(6,.48*nn+1.3*length(unique(hit$module))))
        } else {
          aa<-final_audit[final_audit$job==job,,drop=FALSE]
          reason<-paste('No significant GO BP result. Module status:',paste(unique(aa$status),collapse='; '))
          placeholder(paste(tg,'|',pa,'| GO BP'),reason,stem)
          # 若真实检验完成但不显著，可另外看明确标注“未达显著”的描述性排序图。
          if(nrow(jj)) save_plot(draw_go(jj,tg,pa,exploratory=TRUE),paste0(stem,'_top_tested_NOT_SIGNIFICANT'),
                                 12,max(6,.48*sum(pmin(table(jj$module),GO_TOP_TERMS))+1.3*length(unique(jj$module))))
        }
        for(mm in unique(modules$module[modules$job==job])) {
          hh<-hit[hit$module==mm,,drop=FALSE]
          if(nrow(hh)) save_plot(draw_go(hh,tg,pa,paste('|',mm)),paste0(stem,'_',mm),
                                 11,max(4,.5*min(nrow(hh),GO_TOP_TERMS)+2)) else {
                                   aa<-final_audit[final_audit$job==job & final_audit$module==mm,,drop=FALSE]
                                   placeholder(paste(tg,'|',pa,'|',mm,'| GO BP'),
                                               paste('No significant GO result. Status:',paste(unique(aa$status),collapse='; ')),paste0(stem,'_',mm))
                                 }
        }
      }
    }
    writeLines(c('The previous GO step stopped before enrichment because the OrgDb was already loaded.',
                 paste('v4_dotplot load route:',route),paste('Tested GO rows:',nrow(allgo)),
                 paste('Global FDR significant rows:',sum(allgo$significant_global)),
                 'Check 37_GO_group_path_module_audit.csv to distinguish errors, mapping shortages and true lack of significance.',
                 'Files ending _NO_SIGNIFICANT from v3 are historical; use the freshly written canonical GO_BP_<group>_<path> plots and v4_dotplot audit.',
                 'No significance threshold was relaxed. Descriptive top-tested plots are explicitly labelled NOT SIGNIFICANT.'),
               file.path(run_dir,'GO_resume_v4_dotplot_summary.txt'))
    if(any(final_audit$status=='enrichment_error')) stop('部分模块GO仍有计算错误；其他结果已保存，详见37审计表。')
    cat('GO enrichment completed: ',nrow(allgo),' tested rows; ',sum(allgo$significant_global),' significant rows.\n',sep='')
  })
  writeLines(c(paste('Existing run directory:',run_dir),
               'Saved coordinates, pseudotimes, counts, graph-test results and gene modules were reused.',
               'No embedding, graph, root, graph_test, smoothing or clustering was recomputed.',
               paste('Primary pseudotime palette:',PT_PALETTE,'; color transform:',PT_COLOR_TRANSFORM),
               'All groups share the same scale. Sqrt transforms colors only; legend ticks remain original pseudotime values.',
               'A_pseudotime_6groups_plasma_linear is the linear-scale comparison plot.',
               'GO was newly computed because it did not run successfully in v3.',
               'Check resume_replot_v4_dotplot_status.csv for actual completion of each stage.',
               'GO figure changed: dot bubble plot (size=gene count, color=-log10(FDR)) instead of bar chart.'),
             file.path(run_dir,'REPlot_v4_dotplot_README.txt'))
  if(any(vapply(stages,function(x)x$status=='failed',logical(1))))
    stop('续跑有部分步骤失败；已完成图表仍保存于原目录，查看resume_replot_v4_dotplot_status.csv。',call.=FALSE)
  cat('\nFinished in the SAME directory: ',run_dir,'\n',sep='')
  invisible(run_dir)
}

main_resume()
下面为从头开始

##激素拟时
# -*- coding: UTF-8 -*-
# 大豆：5天苗，处理24小时；CK + ABA/eBL/IAA/MeJA/SA；每组2个生物学重复。
# 运行：Rscript soybean_endodermis_monocle3_grouped_v3.txt
# v3：直接复制原Seurat UMAP坐标，不重新PCA/UMAP；每组单独计算平滑聚类热图和模块GO。
# 每幅热图/GO及表格注明Treatment、Component、Path、Module；没有结果也输出组名与原因。
# 修复：按提供的日志明确使用完整注释cell_type，不自动优先使用稀疏的new_cell_type。
# 日志中cell_type有7473个Endodermis（全部处理）；六组筛选后的数量以本次输出为准。
# 可在R控制台运行 source('soybean_endodermis_monocle3_grouped_v3.txt', encoding='UTF-8')。
# 如果控制台处于Browse[1]>，先输入Q退出错误调试，再source整个新文件。
# 这里只读取输入RDS，不覆盖它。输出保存到独立目录。
# 图A：六组共同建图，固定坐标与色标，按处理分面；不强行生成Y形。
# 初次没有早期细胞证据时使用探索性起点；P1/P2/P3不是已验证的分化状态。
# 处理24小时不是拟时序24小时；CK也不是自动的早期细胞。
# 同一时间采样的处理效应可能影响轨迹，需结合大豆标记与重复验证。
# 依赖安装（请在服务器匹配的R/Bioconductor环境中预先完成，不固定旧版Bioc）：
# install.packages(c('Seurat','ggplot2','igraph','pheatmap','remotes','BiocManager'))
# BiocManager::install(c('SingleCellExperiment','SummarizedExperiment',
#                       'AnnotationDbi','clusterProfiler','GO.db'))
# monocle3安装依赖与安装方法：https://cole-trapnell-lab.github.io/monocle3/docs/installation/
# remotes::install_github('cole-trapnell-lab/monocle3', upgrade='never')

# ========================== 1. 参数区 ==========================
INPUT_RDS <- '/data01/liuxy/Project/CK2_added_v3_Results/seurat_obj_with_CK2_v3.rds'
OUT_DIR <- '/data01/liuxy/Project/CK2_added_v3_Results/Endodermis_monocle3_originalUMAP_grouped_v3'
ASSAY <- 'RNA'                         # 必须为原始counts，不能用integrated表达值
TREATMENT_COL <- 'treatment'            # 根据本次运行日志确定
CELLTYPE_COL <- 'cell_type'             # 根据日志采用覆盖完整对象的注释列
ENDODERMIS_VALUES <- c('Endodermis')     # 根据日志中的精确标签；不混入其他类型
# 不按new/merged等列名判断最终注释，也不自动混合不同注释列。
# new_cell_type/merged_cell_type在日志中各有192718个NA，仅625个Endodermis。
# 当前分析使用cell_type的注释定义；若以后核实其他列才是最终完整注释，可在这里更改。
# 设CELLTYPE_COL <- NULL恢复严格自动识别；候选不一致时停止而不猜测。
# 如果自动匹配不合适，填写精确标签，例如 c('Endodermis','Endodermis_1')。
SAMPLE_COL <- NULL                      # 自动识别明确的两重复样本列；也可手动填写 'sample_id'
# 仅在每组恰好两个样本且候选划分无歧义时自动识别；仍须核对样本审计表。
# 若SAMPLE_COL只有1/2且各处理重复使用，代码会以 treatment::sample 组合为样本ID。
# 必须确认这两项确实是生物学重复，而不是测序lane、技术重复或处理标签。
EXPECTED_REPLICATES <- 2L
SEED <- 20261119L
CORES <- 4L
MIN_GENE_CELLS <- 10L
SEURAT_REDUCTION <- 'umap'              # 原RDS的UMAP；若实际叫umap.integrated，请修改此处
# 本版没有重新降维的备用分支。缺少指定坐标时停止，列出可用reduction供选择。
# 默认不回归treatment，不按处理分别建图，也不强制连接不相连的partition。
# 已确认有技术批次时，可另行制定批次校正方案；不能把处理差异直接当作批次去掉。

# 起点优先级：ROOT_NODES > ROOT_CELL_IDS > ROOT_STATE_COL/VALUES > 探索性端点。
ROOT_NODES <- character(0)              # 如 c('Y_35')，从本次诊断图读取真实节点ID
ROOT_CELL_IDS <- character(0)           # 已知早期细胞barcode，按最多细胞投影的节点定根
ROOT_STATE_COL <- NULL                  # 已有早期/成熟状态注释列
ROOT_STATE_VALUES <- character(0)       # 如 c('early_endodermis')，必须确有依据
ROOT_BIOLOGICALLY_VALIDATED <- FALSE     # 只有经标记/注释验证后才能改TRUE
TARGET_COMPONENT <- NULL                # NULL=细胞最多的连通分量；编号见诊断表
# 多个不相连分量没有可直接比较的拟时尺度。默认只给一个分量定根，其余标灰并统计。
EXTRA_COMPONENT_ROOTS <- character(0)    # 可为其他分量各指定一个真实节点ID；否则用探索性端点
# 主图A与状态比较仍用共同主分量；扩展对其他分量另定局部起点，局部拟时序不可跨分量比较。
STAGE_REFERENCE <- 'CK'                 # 在同一有根分量中，以CK的1/3、2/3分位数定界
POINT_SIZE <- 0.35
PLOT_MAX_PER_GROUP <- Inf               # 默认画所有细胞；只影响绘图，不影响计算与比例
MARKER_GENES <- character(0)            # 填大豆表达矩阵中的确切基因ID用于趋势图

# 扩展：共同图候选基因；各处理内分别计算起点到末端路径热图和对应模块GO BP。
RUN_GENE_ANALYSIS <- TRUE               # 较慢；先看图A可设FALSE
GRAPH_Q_CUTOFF <- 0.01
GRAPH_MORAN_CUTOFF <- 0.10
MAX_HEATMAP_GENES <- 180L               # 每张热图最多180个有清晰动态模式的基因
HEATMAP_MODULES <- 3L
N_PSEUDOTIME_BINS <- 30L
MIN_PATH_CELLS <- 80L                   # 每个处理在一条路径上的最少细胞
MIN_CELLS_PER_BIN <- 8L                 # 分位数分箱，避免稀疏末端单细胞形成红条
MIN_BINS <- 6L
MIN_CELLS_PER_SAMPLE_BIN <- 3L
SMOOTH_DF <- 4                         # 低自由度样条，用于展示趋势，不当作显著性检验
HEATMAP_GRID_POINTS <- 80L
MIN_DETECTED_FRACTION <- 0.10
MIN_SMOOTH_AMPLITUDE <- 0.15             # log1p归一化表达的平滑趋势振幅
MIN_TREND_R2 <- 0.15                    # 相对于分箱均值的描述性拟合优度
MAX_BIN_GAP_FRACTION <- 0.40             # 大段没有细胞支持时跳过，不用平滑填满
MAX_PATHS_PER_GROUP <- 3L               # 每组按细胞数最多的3条合格路径；Inf=全部
GO_TOP_TERMS_PER_MODULE <- 6L
GO_FDR_CUTOFF <- 0.05
RUN_GO <- TRUE
ORGDB_PACKAGE <- 'org.G.max.eg.db'
KEYTYPE_OVERRIDE <- 'GID'
PREFERRED_DB_VERSION <- 'v4'
ORGDB_SEARCH_DIR <- '/home/liuxy'
ORGDB_TARBALL <- NULL                   # 可指定 '/home/liuxy/org.G.max.eg.db_1.44.0.tar.gz'
# PREFERRED_DB_VERSION用于基因ID版本标记清理和数据库记录优选，不是R包版本。
# 如Glyma.01G000100.v4、Glyma.01G000100_v4先提取基础ID并保留v4标签。
# 同一基础ID优先数据库v4，其次无版本的唯一GID；其他版本或多义匹配明确记录。
# 默认表达矩阵已有非v4版本标签的基因不自动跨版本转换；正式映射表可覆盖此限制。
# counts行名保持原样；清理仅用于注释/GO匹配，传给enrichGO的是数据库原始GID。
GENE_ID_MAP_TSV <- NULL                 # 可选两列：feature_id、GID；严格一对一映射
MIN_DB_MATCH_RATE <- 0.50               # 低于此阈值停止GO扩展，不影响已保存的图A

# 只识别明确的大豆Glyma基础ID及常见后缀，避免gsub删除生物学ID本身。
parse_soybean_id <- function(ids) {
  raw <- as.character(ids)
  s <- trimws(raw)
  s <- sub('^(gene|transcript):', '', s, ignore.case=TRUE)
  pattern <- '^(Glyma)[._-]?([0-9]{2})[gG]([0-9]{6})(.*)$'
  valid <- !is.na(s) & grepl(pattern,s,ignore.case=TRUE)
  base <- s; suffix <- rep('',length(s)); ver <- rep(NA_character_,length(s))
  base[valid] <- sub(pattern,'Glyma.\\2G\\3',s[valid],ignore.case=TRUE)
  suffix[valid] <- sub(pattern,'\\4',s[valid],ignore.case=TRUE)
  # 后缀中独立的v数字是版本；.1等纯数字后缀仅作为转录本后缀清理。
  vv <- regexpr('(^|[._-])v[0-9]+(?=$|[._-])',suffix,perl=TRUE,ignore.case=TRUE)
  has <- which(vv>0)
  if(length(has)) {
    hit <- regmatches(suffix,vv)
    ver[has] <- tolower(sub('^[._-]','',hit))
  }
  tail_without_version <- gsub('(^|[._-])v[0-9]+(?=$|[._-])','',suffix,perl=TRUE,ignore.case=TRUE)
  # 允许空后缀、.1/.1.2、.t1/.p1；其他后缀保留原文以避免错误合并。
  safe_tail <- grepl('^([._-]([0-9]+|[tTpP][0-9]+))*$',tail_without_version)
  canonical <- ifelse(valid & safe_tail,base,s)
  data.frame(raw_id=raw,clean_id=canonical,version=ver,
             recognized=valid & safe_tail,stringsAsFactors=FALSE)
}

match_soybean_gid <- function(features,dbkeys,preferred='v4') {
  f <- parse_soybean_id(features); b <- parse_soybean_id(dbkeys)
  preferred <- tolower(preferred)
  b$priority <- ifelse(!is.na(b$version)&b$version==preferred,1L,
                       ifelse(is.na(b$version),2L,3L))
  by_id <- split(seq_len(nrow(b)),b$clean_id)
  out <- data.frame(feature_id=features,clean_id=f$clean_id,feature_version=f$version,
                    GID=NA_character_,db_version=NA_character_,match_status='unmatched',
                    n_db_candidates=0L,stringsAsFactors=FALSE)
  for(i in seq_len(nrow(f))) {
    # 有明确非首选版本时不伪装成v4。用户可提供真实跨版本映射覆盖。
    if(!is.na(f$version[i]) && f$version[i]!=preferred) {
      out$match_status[i] <- 'source_nonpreferred_version_requires_mapping';next
    }
    ii <- by_id[[f$clean_id[i]]]
    if(!length(ii)) next
    out$n_db_candidates[i] <- length(ii)
    best <- ii[b$priority[ii]==min(b$priority[ii])]
    if(min(b$priority[ii])==3L) {
      out$match_status[i] <- 'db_only_nonpreferred_version';next
    }
    # 同优先级有多个转录本时不随便选第一条；若原始精确匹配存在则采用它。
    exact <- best[b$raw_id[best]==f$raw_id[i]]
    if(length(exact)==1) best <- exact
    if(length(best)!=1) {out$match_status[i]<-'ambiguous_same_priority';next}
    j <- best[1];out$GID[i]<-b$raw_id[j];out$db_version[i]<-b$version[j]
    out$match_status[i] <- if(b$priority[j]==1L)'preferred_v4' else 'unique_unversioned_GID'
  }
  out$matched <- !is.na(out$GID)
  out
}

# 在任何细胞筛选之前，对照全部注释候选的覆盖率和六组Endodermis计数。
audit_celltype_columns <- function(md,treatment,groups,match_endodermis) {
  text_columns <- names(md)[vapply(md,function(x)is.character(x)||is.factor(x),logical(1))]
  candidates <- text_columns[vapply(md[text_columns],function(x)any(match_endodermis(x)),logical(1)) |
                               grepl('cell.*type|celltype|annotation',text_columns,ignore.case=TRUE)]
  if(!length(candidates)) stop('metadata没有可审计的细胞类型注释列。')
  selected <- !is.na(treatment) & treatment %in% groups
  result <- lapply(candidates,function(nm) {
    labels <- as.character(md[[nm]])
    missing <- is.na(labels) | tolower(trimws(labels)) %in% c('','na','nan','none','null','<na>')
    endo <- !missing & match_endodermis(labels)
    nn <- table(factor(treatment[selected & endo],levels=groups))
    rr <- data.frame(candidate_column=nm,n_cells_all=nrow(md),
                     n_missing_all=sum(missing),endodermis_all=sum(endo),
                     n_cells_six_groups=sum(selected),n_missing_six_groups=sum(missing & selected),
                     endodermis_six_groups=sum(endo & selected),groups_with_endodermis=sum(nn>0),
                     missing_endodermis_groups=paste(groups[nn==0],collapse=';'),stringsAsFactors=FALSE)
    for(j in seq_along(groups)) rr[[paste0('Endodermis_',groups[j])]] <- as.integer(nn[j])
    rr
  })
  do.call(rbind,result)
}

# 平滑趋势与聚类是展示步骤，不生成伪装成生物学重复检验的p值。
make_smooth_heatmap <- function(expression,pseudotime,samples=NULL) {
  stopifnot(ncol(expression)==length(pseudotime),!is.null(rownames(expression)))
  if(any(!is.finite(pseudotime))) stop('Non-finite local pseudotime.')
  if(length(unique(pseudotime))<MIN_BINS+1L) stop('Too few distinct pseudotime values.')
  # 等数量分箱避免旧版稀疏末端由单个高表达细胞主导。
  nb <- min(N_PSEUDOTIME_BINS,floor(length(pseudotime)/MIN_CELLS_PER_BIN))
  if(nb<MIN_BINS) stop('Too few cells for adequately populated bins.')
  br <- unique(as.numeric(quantile(pseudotime,seq(0,1,length.out=nb+1L),type=7)))
  if(length(br)<MIN_BINS+1L) stop('Tied pseudotime values leave too few bins.')
  bins <- cut(pseudotime,breaks=br,include.lowest=TRUE,labels=FALSE)
  has_reps <- !is.null(samples) && length(samples)==length(pseudotime) && !anyNA(samples)
  means <- list();bin_rows <- list()
  for(b in sort(unique(bins))) {
    ii <- which(bins==b)
    if(length(ii)<MIN_CELLS_PER_BIN) next
    nrep <- NA_integer_;used <- ii
    if(has_reps) {
      chunks <- split(ii,as.character(samples[ii]))
      chunks <- chunks[lengths(chunks)>=MIN_CELLS_PER_SAMPLE_BIN]
      if(!length(chunks)) next
      # 先算每个重复的均值，再对可用重复等权平均，避免细胞数多的重复占优。
      by_rep <- vapply(chunks,function(jj)rowMeans(expression[,jj,drop=FALSE]),numeric(nrow(expression)))
      by_rep <- matrix(by_rep,nrow=nrow(expression))
      value <- rowMeans(by_rep);nrep <- length(chunks);used <- unlist(chunks,use.names=FALSE)
    } else value <- rowMeans(expression[,ii,drop=FALSE])
    means[[as.character(b)]] <- value
    bin_rows[[as.character(b)]] <- data.frame(bin=b,mean_pt=mean(pseudotime[used]),
                                              n_cells=length(used),n_replicates=nrep,min_pt=min(pseudotime[used]),max_pt=max(pseudotime[used]))
  }
  if(length(means)<MIN_BINS) stop('Too few supported bins after cell/replicate checks.')
  bin_info <- do.call(rbind,bin_rows)
  mat <- do.call(cbind,means);rownames(mat)<-rownames(expression)
  xx <- bin_info$mean_pt
  if(any(diff(xx)<=0)) stop('Non-increasing bin centers.')
  if(max(diff(xx))/diff(range(xx))>MAX_BIN_GAP_FRACTION)
    stop('Large unsupported pseudotime gap; not smoothed across the gap.')
  grid_pt <- seq(min(xx),max(xx),length.out=HEATMAP_GRID_POINTS)
  detected <- rowMeans(expression>0)
  smooth_mat <- matrix(NA_real_,nrow(expression),length(grid_pt),dimnames=list(rownames(expression),NULL))
  quality <- data.frame(feature_id=rownames(expression),detected_fraction=detected,
                        amplitude=NA_real_,trend_R2=NA_real_,score=NA_real_,eligible=FALSE)
  for(i in seq_len(nrow(mat))) {
    yy <- as.numeric(mat[i,]); total <- sum((yy-mean(yy))^2)
    if(!is.finite(total)||total<=1e-10||detected[i]<MIN_DETECTED_FRACTION) next
    fit <- tryCatch(smooth.spline(x=xx,y=yy,df=min(SMOOTH_DF,length(xx)-1L)),error=function(e)NULL)
    if(is.null(fit)) next
    smooth_mat[i,] <- predict(fit,grid_pt)$y
    r2 <- 1-sum((yy-predict(fit,xx)$y)^2)/total
    amp <- diff(range(smooth_mat[i,]))
    quality$amplitude[i] <- amp;quality$trend_R2[i] <- r2;quality$score[i] <- amp*max(r2,0)
    quality$eligible[i] <- is.finite(amp)&&amp>=MIN_SMOOTH_AMPLITUDE&&is.finite(r2)&&r2>=MIN_TREND_R2
  }
  selected <- which(quality$eligible)
  if(length(selected)<max(6,HEATMAP_MODULES*2L)) stop('Too few genes with supported smooth dynamics.')
  selected <- head(selected[order(-quality$score[selected],quality$feature_id[selected])],MAX_HEATMAP_GENES)
  sm <- smooth_mat[selected,,drop=FALSE]
  z <- t(scale(t(sm)))
  if(any(!is.finite(z))) stop('Non-finite row-scaled trend matrix.')
  # Ward聚类等价于按标准化曲线的整体形状分组；聚类在截断色阶之前进行。
  hc <- hclust(dist(z),method='ward.D2');k<-min(HEATMAP_MODULES,nrow(z)-1L)
  raw_module <- cutree(hc,k=k)
  centers <- do.call(rbind,lapply(sort(unique(raw_module)),function(m)colMeans(z[raw_module==m,,drop=FALSE])))
  peak <- apply(centers,1,which.max)
  cluster_order <- order(peak,seq_along(peak))
  remap <- setNames(paste0('M',seq_len(k)),sort(unique(raw_module))[cluster_order])
  modules <- unname(remap[as.character(raw_module)]);names(modules)<-rownames(z)
  # 模块按平均峰值从早到晚排列；模块内部仍用层次聚类，不随机排序。
  ord <- unlist(lapply(paste0('M',seq_len(k)),function(m) {
    ids <- names(modules)[modules==m]
    if(length(ids)<2L) return(ids)
    h <- hclust(dist(z[ids,,drop=FALSE]),method='ward.D2')
    dend <- reorder(as.dendrogram(h),wts=apply(z[ids,,drop=FALSE],1,which.max),agglo.FUN=mean)
    ids[order.dendrogram(dend)]
  }),use.names=FALSE)
  z<-z[ord,,drop=FALSE];sm<-sm[ord,,drop=FALSE];modules<-modules[ord]
  gaps<-head(cumsum(as.integer(table(factor(modules,levels=paste0('M',seq_len(k)))))), -1L)
  list(z=z,smooth_expression=sm,modules=modules,gaps=gaps,grid_pt=grid_pt,
       bin_info=bin_info,binned_expression=mat,quality=quality,
       averaging=if(has_reps)'equal weight of available biological-replicate means' else 'pooled-cell means')
}

main <- function() {
  dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)
  run_dir <- file.path(OUT_DIR, format(Sys.time(), '%Y%m%d_%H%M%S'))
  dir.create(run_dir, recursive=TRUE, showWarnings=FALSE)
  cat('Output: ', run_dir, '\n', sep='')
  logcon <- file(file.path(run_dir, 'run.log'), open='wt', encoding='UTF-8')
  sink(logcon, split=TRUE)
  on.exit({sink(); close(logcon)}, add=TRUE)
  on.exit(writeLines(capture.output(sessionInfo()), file.path(run_dir,'sessionInfo.txt')), add=TRUE)
  saveRDS(as.list(.GlobalEnv)[intersect(ls(.GlobalEnv), c(
    'INPUT_RDS','ASSAY','TREATMENT_COL','CELLTYPE_COL','ENDODERMIS_VALUES','SAMPLE_COL',
    'ROOT_NODES','ROOT_CELL_IDS','ROOT_STATE_COL','ROOT_STATE_VALUES','ROOT_BIOLOGICALLY_VALIDATED',
    'TARGET_COMPONENT','STAGE_REFERENCE','SEED','SEURAT_REDUCTION','PREFERRED_DB_VERSION',
    'EXTRA_COMPONENT_ROOTS','MAX_HEATMAP_GENES','SMOOTH_DF','MIN_PATH_CELLS'))],
    file.path(run_dir,'analysis_parameters.rds'))
  required <- c('Seurat','SeuratObject','monocle3','Matrix','ggplot2','igraph',
                'SingleCellExperiment','SummarizedExperiment','S4Vectors')
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
  if(length(missing)) stop('请先安装依赖：',paste(missing,collapse=', '))
  suppressPackageStartupMessages(library(ggplot2))
  set.seed(SEED)
  groups <- c('CK','ABA','eBL','IAA','MeJA','SA')
  group_colors <- c(CK='#555555',ABA='#D55E00',eBL='#009E73',IAA='#0072B2',MeJA='#CC79A7',SA='#E69F00')
  csv <- function(x, name) write.csv(x,file.path(run_dir,paste0(name,'.csv')),
                                     row.names=FALSE,na='',fileEncoding='UTF-8')
  save_plot <- function(p,name,w=10,h=7) {
    ggsave(file.path(run_dir,paste0(name,'.pdf')),p,width=w,height=h,limitsize=FALSE)
    ggsave(file.path(run_dir,paste0(name,'.png')),p,width=w,height=h,dpi=300,limitsize=FALSE)
  }
  placeholder <- function(title,reason,name) {
    p <- ggplot()+annotate('text',x=0,y=0,label=paste(strwrap(reason,70),collapse='\n'),size=4)+
      xlim(-1,1)+ylim(-1,1)+theme_void()+labs(title=title)+
      theme(plot.title=element_text(hjust=0.5,face='bold',size=16),plot.margin=margin(20,20,20,20))
    save_plot(p,name,10,5)
  }
  clean <- function(x) tolower(gsub('[[:space:]_-]+','',trimws(as.character(x))))
  treatment_map <- function(x) {
    lookup <- c(control='CK',ck='CK',aba='ABA',ebl='eBL',iaa='IAA',meja='MeJA',sa='SA')
    unname(lookup[clean(x)])
  }
  is_endo <- function(x) {
    if(!is.null(ENDODERMIS_VALUES)) return(as.character(x) %in% ENDODERMIS_VALUES)
    !is.na(x) & grepl('endodermis|endodermal|内皮层',as.character(x),ignore.case=TRUE)
  }
  
  # ==================== 2. 读RDS、检查metadata、筛选 ====================
  if(!file.exists(INPUT_RDS)) stop('找不到输入RDS：',INPUT_RDS)
  cat('Reading RDS: ',INPUT_RDS,'\n',sep='')
  obj <- readRDS(INPUT_RDS)
  if(!inherits(obj,'Seurat')) stop('输入RDS不是Seurat对象。')
  md <- obj[[]]
  stopifnot(!anyDuplicated(rownames(md)),identical(rownames(md),colnames(obj)))
  text_cols <- names(md)[vapply(md,function(x)is.character(x)||is.factor(x),logical(1))]
  inventory <- do.call(rbind,lapply(names(md),function(nm) data.frame(
    column=nm,class=paste(class(md[[nm]]),collapse='/'),n_unique=length(unique(md[[nm]])),
    examples=paste(head(unique(as.character(md[[nm]])),12),collapse=' | '))))
  csv(inventory,'00_metadata_columns')
  choose_col <- function(explicit,candidates,what) {
    if(!is.null(explicit)) {
      if(!explicit %in% names(md)) stop(what,'列不存在：',explicit)
      return(explicit)
    }
    if(length(candidates)==1) return(candidates)
    if(length(candidates)>1) {
      vectors <- lapply(candidates,function(z) as.character(md[[z]]))
      if(all(vapply(vectors,function(z)identical(z,vectors[[1]]),logical(1)))) return(candidates[1])
    }
    stop(what,'列无法唯一确定；候选：',paste(candidates,collapse=', '),
         '。查看00_metadata_columns.csv并填写参数后重跑。')
  }
  tc <- choose_col(TREATMENT_COL,text_cols[vapply(md[text_cols],function(x)
    all(groups %in% treatment_map(x)),logical(1))],'TREATMENT_COL')
  tr <- treatment_map(md[[tc]])
  celltype_audit <- audit_celltype_columns(md,tr,groups,is_endo)
  csv(celltype_audit,'00b_celltype_candidate_audit')
  print(celltype_audit)
  ec <- choose_col(CELLTYPE_COL,text_cols[vapply(md[text_cols],function(x)any(is_endo(x)),logical(1))],
                   'CELLTYPE_COL')
  celltype_audit$selected <- celltype_audit$candidate_column==ec
  csv(celltype_audit,'00b_celltype_candidate_audit')
  writeLines(c(paste('Selected annotation column:',ec),
               paste('Endodermis values:',if(is.null(ENDODERMIS_VALUES))'pattern matching' else paste(ENDODERMIS_VALUES,collapse=', ')),
               'No automatic priority for new_cell_type or merged_cell_type; no mixing of annotation columns.'),
             file.path(run_dir,'00c_annotation_selection.txt'))
  cat('Treatment column: ',tc,'; Cell type column: ',ec,'\n',sep='')
  audit <- as.data.frame(table(metadata_value=as.character(md[[tc]]),useNA='ifany'),stringsAsFactors=FALSE)
  names(audit)[2] <- 'n_cells'; audit$treatment <- treatment_map(audit$metadata_value)
  audit$status <- ifelse(is.na(audit$treatment),'excluded_not_selected_or_missing','retained')
  csv(audit,'01_all_cell_treatment_audit'); print(audit)
  cat('Selected ',sum(!is.na(tr)),' cells in six treatment groups; excluded ',sum(is.na(tr)),'.\n',sep='')
  csv(as.data.frame(table(annotation=as.character(md[[ec]]),useNA='ifany')),'02_annotation_counts')
  # 保留完整注释与局部新注释的差异，方便判断是否存在后续需要人工核对的重注释。
  for(other in setdiff(intersect(c('cell_type','new_cell_type','merged_cell_type'),names(md)),ec)) {
    cross <- as.data.frame(table(selected_annotation=as.character(md[[ec]])[!is.na(tr)],
                                 alternative_annotation=as.character(md[[other]])[!is.na(tr)],useNA='ifany'))
    csv(cross,paste0('02b_annotation_comparison_',other))
  }
  keep <- !is.na(tr) & is_endo(md[[ec]])
  endo_counts <- as.data.frame(table(treatment=factor(tr[keep],levels=groups)))
  names(endo_counts)[2] <- 'n_cells'
  csv(endo_counts,'03_endodermis_counts_before_qc');print(endo_counts)
  cat('Selected endodermis across six groups: ',sum(keep),' (annotation=',ec,').\n',sep='')
  empty <- as.character(endo_counts$treatment[endo_counts$n_cells==0])
  if(length(empty)) stop('注释列 ',ec,' 中以下处理没有匹配的Endodermis：',paste(empty,collapse=', '),
                         '。已输出00b列对照表和03六组计数表；不删除空组或伪造细胞，请核对真实注释。')
  if(sum(keep)<100) stop('endodermis不足100个细胞；检查注释列及标签。')
  selected_annotation_audit <- data.frame(cell=rownames(md)[keep],treatment=tr[keep],
                                          annotation_source=ec,selected_annotation=as.character(md[[ec]])[keep],stringsAsFactors=FALSE)
  for(other in intersect(c('cell_type','new_cell_type','merged_cell_type'),names(md)))
    selected_annotation_audit[[other]] <- as.character(md[[other]])[keep]
  csv(selected_annotation_audit,'03b_selected_endodermis_annotation_audit')
  md <- md[keep,,drop=FALSE]
  md$analysis_treatment <- factor(tr[keep],levels=groups)
  md$original_endodermis_label <- as.character(md[[ec]])
  sample_candidates <- inventory[grepl('sample|replicate|rep$|orig.ident',inventory$column,ignore.case=TRUE),]
  csv(sample_candidates,'04_sample_column_candidates')
  sample_col <- SAMPLE_COL
  if(is.null(sample_col)) {
    viable <- sample_candidates$column[vapply(sample_candidates$column,function(nm) {
      x <- as.character(md[[nm]])
      !anyNA(x) && all(nzchar(trimws(x))) &&
        all(vapply(split(x,md$analysis_treatment),function(z)length(unique(z))==EXPECTED_REPLICATES,logical(1)))
    },logical(1))]
    # 不同列只有命名不同但同一批细胞分组相同，也视作相同样本划分。
    canonical_groups <- lapply(viable,function(nm) {
      x <- paste(md$analysis_treatment,md[[nm]],sep='::');match(x,unique(x))
    })
    if(length(viable)==1 || (length(viable)>1 &&
                             all(vapply(canonical_groups,function(x)identical(x,canonical_groups[[1]]),logical(1)))))
      sample_col <- viable[1]
  }
  if(!is.null(sample_col)) {
    if(!sample_col %in% names(md)) stop('SAMPLE_COL不存在；请查看00_metadata_columns.csv。')
    if(anyNA(md[[sample_col]]) || any(!nzchar(trimws(as.character(md[[sample_col]])))))
      stop('选中endodermis有缺失样本ID；请修正SAMPLE_COL。')
    md$analysis_sample <- paste(md$analysis_treatment,as.character(md[[sample_col]]),sep='::')
    writeLines(paste('Biological sample column:',sample_col),file.path(run_dir,'04_selected_sample_column.txt'))
    cat('Biological sample column: ',sample_col,'; please verify sample assignments.\n',sep='')
  } else {
    md$analysis_sample <- NA_character_
    cat('SAMPLE_COL未指定：会生成合并组描述图；设置真实生物学样本列后重跑获得两个重复的结果。\n')
  }
  # 先取子集再合并counts层，避免对全部20万细胞创建稠密矩阵。
  if(!ASSAY %in% names(obj@assays)) stop('找不到assay：',ASSAY)
  if(!SEURAT_REDUCTION %in% names(obj@reductions)) stop('缺少指定的原始UMAP：',SEURAT_REDUCTION,
                                                        '；可用reduction：',paste(names(obj@reductions),collapse=', '),'。本版不会重新降维。')
  original_umap <- SeuratObject::Embeddings(obj,SEURAT_REDUCTION)
  if(ncol(original_umap)<2L || !all(rownames(md) %in% rownames(original_umap)))
    stop('指定UMAP缺少选中细胞或不足二维。')
  original_umap <- original_umap[rownames(md),1:2,drop=FALSE]
  if(any(!is.finite(original_umap))) stop('原UMAP有非有限坐标，请先检查输入对象。')
  obj <- subset(obj,cells=rownames(md))
  if(inherits(obj[[ASSAY]],'Assay5')) {
    layers <- SeuratObject::Layers(obj[[ASSAY]])
    count_layers <- grep('^counts($|[.])',layers,value=TRUE)
    if(!length(count_layers)) stop('RNA assay没有原始counts层。')
    if(length(count_layers)>1 || count_layers!='counts')
      obj[[ASSAY]] <- SeuratObject::JoinLayers(obj[[ASSAY]],layers=count_layers,new='counts')
    counts <- SeuratObject::LayerData(obj,assay=ASSAY,layer='counts')
  } else counts <- SeuratObject::GetAssayData(obj,assay=ASSAY,slot='counts')
  if(!all(rownames(md) %in% colnames(counts))) stop('RNA counts缺少选中细胞。')
  counts <- counts[,rownames(md),drop=FALSE]
  if(!inherits(counts,'sparseMatrix')) counts <- methods::as(counts,'dgCMatrix')
  if(anyDuplicated(rownames(counts))) stop('counts有重复基因名。')
  if(any(!is.finite(counts@x)) || any(counts@x<0) || any(abs(counts@x-round(counts@x))>1e-6))
    stop('counts不是非负整数原始计数；请检查ASSAY。')
  gene_keep <- Matrix::rowSums(counts>0)>=MIN_GENE_CELLS
  counts <- counts[gene_keep,,drop=FALSE]
  cell_keep <- Matrix::colSums(counts)>0
  csv(data.frame(cell=colnames(counts),retained=cell_keep),'05_post_gene_filter_cell_qc')
  counts <- counts[,cell_keep,drop=FALSE]; md <- md[colnames(counts),,drop=FALSE]
  if(nrow(counts)<100 || ncol(counts)<100) stop('过滤后基因/细胞数过低。')
  if(any(table(md$analysis_treatment)==0)) stop('过滤后有处理组为空。')
  old_umap <- original_umap[colnames(counts),,drop=FALSE]
  csv(data.frame(cell=rownames(old_umap),UMAP_1=old_umap[,1],UMAP_2=old_umap[,2]),'05b_original_UMAP_coordinates')
  rm(original_umap)
  rm(obj); gc()
  
  # ========================== 3. 六组共同轨迹 ==========================
  gene_md <- data.frame(gene_short_name=rownames(counts),row.names=rownames(counts))
  cds <- monocle3::new_cell_data_set(counts,cell_metadata=md,gene_metadata=gene_md)
  cds <- monocle3::estimate_size_factors(cds)
  SingleCellExperiment::reducedDim(cds,'UMAP') <- old_umap
  cds <- monocle3::cluster_cells(cds,reduction_method='UMAP',random_seed=SEED)
  cds <- monocle3::learn_graph(cds,use_partition=TRUE,close_loop=FALSE)
  preserved <- SingleCellExperiment::reducedDim(cds,'UMAP')[colnames(counts),1:2,drop=FALSE]
  delta <- max(abs(unname(preserved)-unname(old_umap)))
  csv(data.frame(reduction=SEURAT_REDUCTION,max_coordinate_change=delta,
                 n_cells=nrow(preserved),preserved_exactly=delta==0),'05c_UMAP_preservation_check')
  if(delta!=0) stop('UMAP坐标发生变化；中止以防输出改变形状的图。')
  saveRDS(cds,file.path(run_dir,'cds_before_root.rds'))
  g <- monocle3::principal_graph(cds)[['UMAP']]
  vn <- igraph::V(g)$name
  # Monocle3目前将主图坐标与最近节点保存在principal_graph_aux中。
  aux <- cds@principal_graph_aux[['UMAP']]
  xy <- t(aux$dp_mst)
  if(is.null(rownames(xy)) || !all(vn %in% rownames(xy))) stop('主图坐标与节点ID不匹配。')
  nodes <- data.frame(node=vn,x=xy[vn,1],y=xy[vn,2],degree=as.integer(igraph::degree(g)[vn]))
  membership <- igraph::components(g)$membership
  nodes$component <- unname(membership[nodes$node])
  cv <- aux$pr_graph_cell_proj_closest_vertex
  raw_vertex <- as.character(cv[colnames(cds),1])
  if(all(raw_vertex %in% vn)) nearest <- raw_vertex else {
    idx <- suppressWarnings(as.integer(raw_vertex))
    if(anyNA(idx)||any(idx<1|idx>length(vn))) stop('无法解析最近主图节点索引。')
    nearest <- vn[idx]
  }
  names(nearest) <- colnames(cds)
  um <- SingleCellExperiment::reducedDim(cds,'UMAP')
  d <- data.frame(cell=colnames(cds),x=um[,1],y=um[,2],treatment=md$analysis_treatment,
                  sample=md$analysis_sample,annotation=md$original_endodermis_label,
                  cluster=as.character(monocle3::clusters(cds)),
                  partition=as.character(monocle3::partitions(cds)),nearest_node=nearest,
                  component=unname(membership[nearest]),row.names=NULL)
  nodes$n_cells <- as.integer(table(factor(nearest,levels=nodes$node)))
  edges <- igraph::as_data_frame(g,what='edges')[,c('from','to')]
  edges$x <- nodes$x[match(edges$from,nodes$node)]; edges$y <- nodes$y[match(edges$from,nodes$node)]
  edges$xend <- nodes$x[match(edges$to,nodes$node)]; edges$yend <- nodes$y[match(edges$to,nodes$node)]
  graph_layer <- function() geom_segment(data=edges,aes(x=x,y=y,xend=xend,yend=yend),
                                         inherit.aes=FALSE,color='black',linewidth=0.35)
  base_theme <- theme_classic(base_size=12)+theme(strip.background=element_blank(),
                                                  strip.text=element_text(face='bold'),legend.position='right',
                                                  plot.title=element_text(size=14,face='bold'),plot.subtitle=element_text(size=10),
                                                  plot.caption=element_text(size=8),plot.margin=margin(12,16,12,12))
  diag <- ggplot(d,aes(x,y,color=treatment))+geom_point(size=0.4,alpha=0.65)+graph_layer()+
    geom_text(data=nodes,aes(x,y,label=node),inherit.aes=FALSE,size=1.8,check_overlap=TRUE)+
    scale_color_manual(values=group_colors)+coord_equal()+base_theme+labs(x='UMAP 1',y='UMAP 2',
                                                                          title='Endodermis: pooled graph and principal-node IDs')
  save_plot(diag,'06_root_diagnostic_nodes',12,9)
  save_plot(ggplot(d,aes(x,y,color=factor(component)))+geom_point(size=0.4)+graph_layer()+
              coord_equal()+base_theme+labs(color='Component',x='UMAP 1',y='UMAP 2'),'07_components')
  csv(nodes,'08_principal_nodes'); csv(edges,'09_graph_edges')
  csv(as.data.frame(table(treatment=d$treatment,component=d$component)),'10_component_counts')
  
  # ========================== 4. 起点与拟时序 ==========================
  target <- if(is.null(TARGET_COMPONENT)) as.integer(names(which.max(table(d$component)))) else
    as.integer(TARGET_COMPONENT)
  if(!target %in% d$component) stop('TARGET_COMPONENT不存在。')
  root_source <- 'exploratory_topological_endpoint'
  chosen <- ROOT_NODES
  candidate_cells <- ROOT_CELL_IDS
  if(!length(chosen) && !length(candidate_cells) && !is.null(ROOT_STATE_COL)) {
    if(!ROOT_STATE_COL %in% names(md)) stop('ROOT_STATE_COL不存在。')
    candidate_cells <- rownames(md)[as.character(md[[ROOT_STATE_COL]]) %in% ROOT_STATE_VALUES]
    if(!length(candidate_cells)) stop('ROOT_STATE_VALUES没有匹配到细胞。')
    root_source <- 'user_early_state_annotation'
  }
  if(length(chosen)) {
    if(!all(chosen %in% vn)) stop('ROOT_NODES含不存在的节点；查看08_principal_nodes.csv。')
    comp <- unique(unname(membership[chosen]))
    if(length(comp)!=1) stop('请每次分析一个有根连通分量；ROOT_NODES跨越多个分量。')
    target <- comp; root_source <- 'user_principal_nodes'
  } else if(length(candidate_cells)) {
    if(!all(candidate_cells %in% d$cell)) stop('ROOT_CELL_IDS含不在当前子集的细胞。')
    comp <- unique(d$component[match(candidate_cells,d$cell)])
    if(length(comp)!=1) stop('早期候选细胞跨连通分量；请限定到同一分量。')
    target <- comp
    chosen <- names(which.max(table(nearest[candidate_cells])))
    if(root_source!='user_early_state_annotation') root_source <- 'user_early_cells'
  } else {
    pool <- nodes[nodes$component==target & nodes$degree==1,,drop=FALSE]
    if(!nrow(pool)) pool <- nodes[nodes$component==target,,drop=FALSE]
    # 仅为让初次运行能完整展示；该端点没有早期/成熟生物学含义。
    pool <- pool[order(-pool$n_cells,pool$node),,drop=FALSE]; chosen <- pool$node[1]
  }
  validated <- ROOT_BIOLOGICALLY_VALIDATED && root_source!='exploratory_topological_endpoint'
  orientation <- if(validated) 'User-validated root' else 'Exploratory root; direction unvalidated'
  cat('Root: ',paste(chosen,collapse=','),' | ',root_source,' | component=',target,'\n',sep='')
  cds <- monocle3::order_cells(cds,root_pr_nodes=chosen)
  pt <- monocle3::pseudotime(cds)
  d$pseudotime <- as.numeric(pt[d$cell])
  d$finite_pseudotime <- is.finite(d$pseudotime) & d$component==target
  d$pseudotime_plot <- ifelse(d$finite_pseudotime,d$pseudotime,NA_real_)
  if(sum(d$finite_pseudotime)<30) stop('有根连通分量的有限拟时序细胞过少。')
  # 共享阈值；禁止每组各自三等分（那会人为制造每组相同的比例）。
  ref <- d$pseudotime[d$finite_pseudotime & d$treatment==STAGE_REFERENCE]
  if(length(ref)<30) stop('参考组在有根分量不足30个细胞；检查起点/连通分量。')
  cuts <- unname(quantile(ref,c(1/3,2/3),type=7))
  if(cuts[1]>=cuts[2]) stop('参考组拟时序分位数重合；无法合理划分三段。')
  stage_labels <- if(validated) c('Early_PT','Middle_PT','Late_PT') else c('P1','P2','P3')
  d$pt_interval <- 'Unrooted'
  d$pt_interval[d$finite_pseudotime] <- as.character(cut(d$pseudotime[d$finite_pseudotime],
                                                         breaks=c(-Inf,cuts,Inf),labels=stage_labels,include.lowest=TRUE))
  d$pt_interval <- factor(d$pt_interval,levels=c(stage_labels,'Unrooted'))
  csv(data.frame(root_node=chosen,source=root_source,biologically_validated=validated,
                 component=target,reference=STAGE_REFERENCE,cut_1=cuts[1],cut_2=cuts[2]),'11_root_and_interval_definition')
  
  # 额外局部定根：为了分别展示各处理在断开分量内的动态，不改变主图的共同根拟时序。
  if(length(chosen)!=1L) stop('本版路径热图要求主分量恰好一个起点；请将ROOT_NODES限定为一个节点。')
  if(!all(EXTRA_COMPONENT_ROOTS %in% vn)) stop('EXTRA_COMPONENT_ROOTS中有不存在的节点。')
  if(anyDuplicated(unname(membership[EXTRA_COMPONENT_ROOTS]))) stop('其他每个分量只能指定一个起点。')
  if(any(unname(membership[EXTRA_COMPONENT_ROOTS])==target)) stop('主分量起点请用ROOT_NODES指定。')
  local_roots <- do.call(rbind,lapply(sort(unique(d$component)),function(cc) {
    if(cc==target) return(data.frame(component=cc,root_node=chosen,source=root_source,
                                     biologically_validated=validated))
    user_root <- EXTRA_COMPONENT_ROOTS[unname(membership[EXTRA_COMPONENT_ROOTS])==cc]
    if(length(user_root)) return(data.frame(component=cc,root_node=user_root,
                                            source='user_extra_component_root',biologically_validated=FALSE))
    pool <- nodes[nodes$component==cc & nodes$degree==1,,drop=FALSE]
    if(!nrow(pool)) pool <- nodes[nodes$component==cc,,drop=FALSE]
    pool <- pool[order(-pool$n_cells,pool$node),,drop=FALSE]
    data.frame(component=cc,root_node=pool$node[1],source='exploratory_local_endpoint',
               biologically_validated=FALSE)
  }))
  cds_local <- monocle3::order_cells(cds,root_pr_nodes=local_roots$root_node)
  local_pt <- monocle3::pseudotime(cds_local)
  d$local_pseudotime <- as.numeric(local_pt[d$cell])
  d$local_root <- local_roots$root_node[match(d$component,local_roots$component)]
  csv(local_roots,'11b_roots_for_separate_component_heatmaps')
  saveRDS(cds_local,file.path(run_dir,'endodermis_local_roots_separate_components.rds'))
  rm(cds_local);gc()
  
  # 图分支段：去掉度数>=3的分叉节点后，剩余连通段；分叉节点单独标注。
  # 这是确定性拓扑分段，不等价于Monocle2 State，更不能直接命名为细胞命运。
  junction <- nodes$node[nodes$degree>=3]
  segment_map <- setNames(rep(NA_character_,length(vn)),vn)
  if(length(junction)<length(vn)) {
    sg <- igraph::delete_vertices(g,junction)
    sm <- igraph::components(sg)$membership
    segment_map[names(sm)] <- sprintf('Segment_%03d',as.integer(sm))
  }
  segment_map[junction] <- paste0('Junction_',junction)
  d$graph_segment <- unname(segment_map[d$nearest_node])
  csv(data.frame(node=vn,graph_segment=unname(segment_map[vn])),'12_node_segment_map')
  csv(d,'13_cell_metadata_pseudotime_states')
  SummarizedExperiment::colData(cds)$analysis_pseudotime <- d$pseudotime
  SummarizedExperiment::colData(cds)$pt_interval <- d$pt_interval
  SummarizedExperiment::colData(cds)$graph_segment <- d$graph_segment
  SummarizedExperiment::colData(cds)$local_pseudotime <- d$local_pseudotime
  SummarizedExperiment::colData(cds)$local_root <- d$local_root
  saveRDS(cds,file.path(run_dir,'endodermis_monocle3_ordered.rds'))
  
  # ========================== 5. 图A与状态比例 ==========================
  plot_d <- d
  if(is.finite(PLOT_MAX_PER_GROUP)) {
    set.seed(SEED)
    ix <- unlist(lapply(split(seq_len(nrow(d)),d$treatment),function(ii)
      if(length(ii)>PLOT_MAX_PER_GROUP) sample(ii,PLOT_MAX_PER_GROUP) else ii),use.names=FALSE)
    plot_d <- d[ix,,drop=FALSE]
  }
  pt_limits <- range(d$pseudotime[d$finite_pseudotime])
  # 每个面板用浅灰细胞显示完整参考轮廓；只叠加该处理实际占据的分量的轨迹。
  background_d <- d[,c('x','y')]
  edges$component <- unname(membership[edges$from])
  facet_edges <- do.call(rbind,lapply(groups,function(tg) {
    ee <- edges[edges$component %in% unique(d$component[d$treatment==tg]),,drop=FALSE]
    ee$treatment <- factor(rep(tg,nrow(ee)),levels=groups);ee
  }))
  base_a <- ggplot(plot_d,aes(x,y))+
    geom_point(data=background_d,aes(x,y),inherit.aes=FALSE,color='grey93',size=POINT_SIZE)+
    geom_segment(data=facet_edges,aes(x=x,y=y,xend=xend,yend=yend),inherit.aes=FALSE,
                 color='grey35',linewidth=0.3)+
    geom_point(aes(color=pseudotime_plot),size=POINT_SIZE,alpha=0.9)+
    scale_color_gradient(low='#123047',high='#55B5F4',limits=pt_limits,na.value='grey80',name='Pseudotime')+
    geom_point(data=nodes[nodes$node %in% chosen,],aes(x,y),inherit.aes=FALSE,
               shape=23,fill='#FFD700',color='black',size=2.4)+
    coord_equal()+base_theme+labs(x='Original UMAP 1',y='Original UMAP 2',title='Endodermis | original UMAP',
                                  subtitle=paste(orientation,'| primary component',target),
                                  caption='Gold: root. Grey: outside the rooted component. Pale background: pooled reference cells.')
  save_plot(base_a+facet_wrap(~treatment,ncol=3,drop=FALSE),'A_pseudotime_6groups',14,9)
  save_plot(base_a+facet_wrap(~treatment,ncol=1,drop=FALSE),'A_pseudotime_vertical',7,22)
  save_plot(ggplot(d,aes(x,y,color=treatment))+geom_point(size=POINT_SIZE)+
              scale_color_manual(values=group_colors)+coord_equal()+base_theme+
              labs(title='Endodermis | unchanged Seurat UMAP',x='Original UMAP 1',y='Original UMAP 2'),
            'A0_original_UMAP_by_treatment',10,8)
  for(cc in local_roots$component) {
    dc <- d[d$component==cc & is.finite(d$local_pseudotime),,drop=FALSE]
    if(nrow(dc)<10) next
    ee <- edges[edges$component==cc,,drop=FALSE]
    pp <- ggplot(dc,aes(x,y,color=local_pseudotime))+geom_point(size=POINT_SIZE)+
      geom_segment(data=ee,aes(x=x,y=y,xend=xend,yend=yend),inherit.aes=FALSE,
                   color='grey35',linewidth=0.3)+
      scale_color_gradient(low='#123047',high='#55B5F4',name='Local PT')+
      facet_wrap(~treatment,ncol=3,drop=FALSE)+coord_equal()+base_theme+
      labs(title=paste('Endodermis | original UMAP | component',cc),
           subtitle='Local root; compare groups within this component only',
           x='Original UMAP 1',y='Original UMAP 2')
    save_plot(pp,paste0('A_local_component_',cc,'_6groups'),13,8)
  }
  save_plot(ggplot(d,aes(x,y,color=pt_interval))+graph_layer()+geom_point(size=POINT_SIZE)+
              facet_wrap(~treatment,ncol=3,drop=FALSE)+coord_equal()+base_theme+
              labs(x='UMAP 1',y='UMAP 2',color='PT interval',subtitle=orientation),'A_interval_6groups',14,9)
  save_plot(ggplot(d,aes(x,y,color=graph_segment))+geom_point(size=POINT_SIZE)+graph_layer()+
              coord_equal()+base_theme+labs(x='UMAP 1',y='UMAP 2',color='Graph segment'),'14_graph_segments',12,8)
  tab <- as.data.frame(table(treatment=d$treatment,interval=d$pt_interval));names(tab)[3]<-'n_cells'
  tab$proportion <- ave(tab$n_cells,tab$treatment,FUN=function(v)v/sum(v))
  csv(tab,'15_interval_proportion_all_endodermis')
  save_plot(ggplot(tab,aes(treatment,proportion,fill=interval))+geom_col()+base_theme+
              labs(y='Proportion of all selected endodermis',x=NULL,subtitle=orientation),'16_interval_proportion')
  segtab <- as.data.frame(table(treatment=d$treatment,segment=d$graph_segment));names(segtab)[3]<-'n_cells'
  segtab$proportion <- ave(segtab$n_cells,segtab$treatment,FUN=function(v)v/sum(v))
  csv(segtab,'17_graph_segment_proportion')
  finite_d <- d[d$finite_pseudotime,,drop=FALSE]
  save_plot(ggplot(finite_d,aes(treatment,pseudotime,fill=treatment))+
              geom_violin(scale='width',trim=TRUE)+geom_boxplot(width=0.12,outlier.shape=NA,fill='white')+
              scale_fill_manual(values=group_colors)+base_theme+theme(legend.position='none')+
              labs(x=NULL,y='Pseudotime (rooted component only)',subtitle=orientation),'18_pseudotime_distribution')
  coverage <- as.data.frame(table(treatment=d$treatment,finite=d$finite_pseudotime));names(coverage)[3]<-'n_cells'
  csv(coverage,'19_rooted_coverage_by_treatment')
  
  # 每个生物学重复一个点；不对细胞执行Wilcoxon/卡方并充当生物学重复。
  if(!is.null(sample_col)) {
    sample_info <- unique(d[,c('sample','treatment')])
    reps <- as.data.frame(table(treatment=sample_info$treatment));names(reps)[2]<-'n_replicates'
    csv(reps,'20_replicate_counts'); print(reps)
    if(any(reps$n_replicates!=EXPECTED_REPLICATES)) warning('样本列识别的重复数不全为2；请核对SAMPLE_COL。')
    st <- as.data.frame(table(sample=d$sample,interval=d$pt_interval)); names(st)[3]<-'n_cells'
    st$treatment <- sample_info$treatment[match(st$sample,sample_info$sample)]
    st$proportion_all <- ave(st$n_cells,st$sample,FUN=function(v)v/sum(v))
    finite_tot <- tapply(as.integer(d$finite_pseudotime),d$sample,sum)
    st$proportion_rooted <- ifelse(st$interval=='Unrooted',NA_real_,
                                   st$n_cells/as.numeric(finite_tot[as.character(st$sample)]))
    st$proportion_rooted[!is.finite(st$proportion_rooted)] <- NA_real_
    csv(st,'21_per_replicate_interval_proportions')
    sample_summary <- do.call(rbind,lapply(split(d,d$sample),function(z) {
      p <- z$pseudotime[z$finite_pseudotime]
      data.frame(sample=z$sample[1],treatment=as.character(z$treatment[1]),n_cells=nrow(z),
                 n_rooted=length(p),rooted_fraction=length(p)/nrow(z),
                 median_pseudotime=if(length(p))median(p) else NA_real_)
    }))
    sample_summary$treatment <- factor(sample_summary$treatment,levels=groups)
    csv(sample_summary,'22_per_replicate_pseudotime_summary')
    save_plot(ggplot(st,aes(treatment,proportion_all,color=treatment))+
                geom_point(position=position_jitter(width=0.08,height=0,seed=SEED),size=2.5)+
                facet_wrap(~interval,nrow=1)+scale_color_manual(values=group_colors)+base_theme+
                theme(legend.position='none')+labs(x=NULL,y='Proportion per biological replicate'),
              '23_replicate_interval_proportions',14,5)
    save_plot(ggplot(sample_summary,aes(treatment,median_pseudotime,color=treatment))+
                geom_point(position=position_jitter(width=0.08,height=0,seed=SEED),size=3)+
                scale_color_manual(values=group_colors)+base_theme+labs(x=NULL,y='Median pseudotime per replicate'),
              '24_replicate_median_pseudotime')
  }
  writeLines(c('P1/P2/P3 (or Early_PT/Middle_PT/Late_PT) are shared pseudotime intervals, not validated cell fates.',
               'Original Seurat UMAP coordinates are copied exactly; coordinate changes are audited.',
               'Heatmaps and GO are computed separately for each treatment and each root-to-leaf path.',
               'Heatmaps smooth replicate-balanced quantile-bin means, hierarchically cluster rows, and retain chronological columns.',
               'Disconnected components use separate local roots for heatmaps; their pseudotimes cannot be compared directly.',
               'Module IDs are local to each group/path; similarly named modules do not imply the same gene set.',
               'Graph segments are topological segments, not automatically biological differentiation states.',
               'No cell-level p-values are used as biological-replicate evidence.',
               'With two replicates per treatment, inspect both replicates; pooled-cell plots are descriptive.',
               'No treatment-specific trajectory is re-learned. Unrooted cells remain in coverage and proportion tables.',
               'Pseudotime is not elapsed time and does not measure a differentiation rate.',orientation),
             file.path(run_dir,'INTERPRETATION.txt'))
  
  # ==================== 6. 每组独立的路径热图与GO ====================
  # 共同图只用于固定几何和定义路径；热图表达值、筛选、聚类、模块均在各处理内计算。
  if(RUN_GENE_ANALYSIS) tryCatch({
    if(!requireNamespace('pheatmap',quietly=TRUE)) stop('扩展需要pheatmap。')
    for(tg in groups) {
      dir.create(file.path(run_dir,'B_heatmaps',tg),recursive=TRUE,showWarnings=FALSE)
      dir.create(file.path(run_dir,'C_GO',tg),recursive=TRUE,showWarnings=FALSE)
    }
    cat('Running graph_test on the pooled, unchanged UMAP graph...\n')
    gt <- monocle3::graph_test(cds,neighbor_graph='principal_graph',cores=CORES)
    gt$feature_id <- rownames(gt);csv(gt,'25_graph_test_pooled_candidate_screen')
    tested <- gt$feature_id[is.finite(gt$q_value)&is.finite(gt$morans_I)]
    sig <- gt[is.finite(gt$q_value)&is.finite(gt$morans_I)&
                gt$q_value<GRAPH_Q_CUTOFF&gt$morans_I>GRAPH_MORAN_CUTOFF,,drop=FALSE]
    sig<-sig[order(-sig$morans_I,sig$q_value),,drop=FALSE];csv(sig,'26_pooled_graph_candidate_genes')
    # 此候选筛选不代表某激素特异DE；后续模块仅描述该处理内的拟时表达模式。
    sf <- monocle3::size_factors(cds)
    if(any(!is.finite(sf)|sf<=0)) stop('无效size factors。')
    log_expr <- function(ids,ix) log1p(as.matrix(counts[ids,ix,drop=FALSE]) /
                                         rep(sf[ix],each=length(ids)))
    paths <- list();path_info <- list();path_members <- list()
    for(cc in local_roots$component) {
      root <- local_roots$root_node[local_roots$component==cc]
      leaves <- nodes$node[nodes$component==cc & nodes$degree==1 & nodes$node!=root]
      for(j in seq_along(leaves)) {
        pnodes <- igraph::as_ids(igraph::shortest_paths(g,from=root,to=leaves[j],
                                                        mode='all',weights=NA,output='vpath')$vpath[[1]])
        if(length(pnodes)<3L) next
        key<-sprintf('C%02d_P%02d',cc,j)
        paths[[key]]<-pnodes
        path_info[[key]]<-data.frame(path=key,component=cc,root_node=root,leaf_node=leaves[j],
                                     orientation=local_roots$source[local_roots$component==cc],n_principal_nodes=length(pnodes))
        path_members[[key]]<-data.frame(path=key,node=pnodes)
      }
    }
    if(!length(paths)) stop('未得到具有至少3个主图节点的起点到末端路径。')
    csv(do.call(rbind,path_info),'26b_root_to_leaf_paths')
    csv(do.call(rbind,path_members),'26c_path_node_membership')
    # 共用主干的细胞可出现在多条路径，但不会跨处理混用；输出成员表便于检查。
    module_tables<-list();jobs<-list();manifest<-list();group_status<-list()
    save_heatmap <- function(result,tg,path,ncell,file_stem) {
      z<-result$z;z[]<-pmax(-2.5,pmin(2.5,z))
      colnames(z)<-sprintf('pt_%03d',seq_len(ncol(z)))
      mods<-unique(result$modules)
      module_cols<-setNames(colorRampPalette(c('#E69F00','#56B4E9','#009E73','#CC79A7'))(length(mods)),mods)
      annrow<-data.frame(Module=factor(result$modules,levels=mods),row.names=rownames(z))
      anncol<-data.frame(Treatment=factor(rep(tg,ncol(z)),levels=groups),
                         Pseudotime=result$grid_pt,row.names=colnames(z))
      labels<-rep('',ncol(z));ticks<-unique(round(seq(1,ncol(z),length.out=5)))
      labels[ticks]<-format(round(result$grid_pt[ticks],2),trim=TRUE)
      ph<-pheatmap::pheatmap(z,cluster_rows=FALSE,cluster_cols=FALSE,
                             gaps_row=result$gaps,annotation_row=annrow,annotation_col=anncol,
                             annotation_colors=list(Treatment=group_colors,Module=module_cols,
                                                    Pseudotime=c('#E5F5F9','#2C7FB8')),
                             labels_col=labels,angle_col='0',show_colnames=TRUE,show_rownames=FALSE,
                             annotation_names_col=TRUE,annotation_names_row=FALSE,border_color=NA,
                             color=colorRampPalette(c('#2166AC','#F7F7F7','#B2182B'))(100),
                             breaks=seq(-2.5,2.5,length.out=101),legend_breaks=c(-2,-1,0,1,2),
                             main=paste0(tg,' | ',path,'\n',nrow(z),' genes; ',ncell,' cells | smooth row Z-score'),
                             fontsize=10,fontsize_col=9,silent=TRUE)
      grDevices::pdf(file.path(run_dir,paste0(file_stem,'.pdf')),width=9,height=9)
      tryCatch({grid::grid.newpage();grid::grid.draw(ph$gtable)},finally=grDevices::dev.off())
      grDevices::png(file.path(run_dir,paste0(file_stem,'.png')),width=2700,height=2700,res=300)
      tryCatch({grid::grid.newpage();grid::grid.draw(ph$gtable)},finally=grDevices::dev.off())
    }
    for(tg in groups) {
      cat('Heatmaps for ',tg,'...\n',sep='')
      np<-vapply(paths,function(pn)sum(d$treatment==tg & d$nearest_node %in% pn &
                                         is.finite(d$local_pseudotime)),integer(1))
      path_order<-names(sort(np,decreasing=TRUE));done<-0L
      for(path in path_order) {
        meta<-path_info[[path]];ncell<-unname(np[path]);job<-paste(tg,path,sep='__')
        report<-data.frame(job=job,treatment=tg,component=meta$component,path=path,n_cells=ncell,
                           status='',detail='',n_heatmap_genes=0L,stringsAsFactors=FALSE)
        if(ncell<MIN_PATH_CELLS) {
          report$status<-'skipped';report$detail<-'too_few_cells';manifest[[job]]<-report;next
        }
        if(done>=MAX_PATHS_PER_GROUP) {
          report$status<-'skipped';report$detail<-'configured_path_limit';manifest[[job]]<-report;next
        }
        outcome<-tryCatch({
          ix<-which(d$treatment==tg & d$nearest_node %in% paths[[path]] & is.finite(d$local_pseudotime))
          det<-Matrix::rowSums(counts[,ix,drop=FALSE]>0)
          detectable<-rownames(counts)[det>=max(10L,ceiling(MIN_DETECTED_FRACTION*length(ix)))]
          universe_features<-intersect(tested,detectable)
          candidates<-intersect(sig$feature_id,universe_features)
          if(length(candidates)<6L) stop('Too few expressed pooled-graph candidates in this treatment/path.')
          # 最多1500个候选进入平滑，之后按本组的动态质量选至MAX_HEATMAP_GENES。
          candidates<-head(candidates,1500L)
          ex<-log_expr(candidates,ix)
          rr<-make_smooth_heatmap(ex,d$local_pseudotime[ix],d$sample[ix])
          stem<-file.path('B_heatmaps',tg,paste0('B_heatmap_',tg,'_',path))
          save_heatmap(rr,tg,path,length(ix),stem)
          csv(data.frame(treatment=tg,path=path,rr$quality),paste0(stem,'_gene_selection'))
          csv(data.frame(treatment=tg,path=path,rr$bin_info),paste0(stem,'_bins'))
          csv(data.frame(treatment=tg,path=path,feature_id=rownames(rr$binned_expression),
                         rr$binned_expression,check.names=FALSE),paste0(stem,'_raw_binned_expression'))
          csv(data.frame(treatment=tg,path=path,feature_id=rownames(rr$z),
                         module=unname(rr$modules),rr$z,check.names=FALSE),paste0(stem,'_rowZ_matrix'))
          csv(data.frame(treatment=tg,path=path,feature_id=rownames(rr$smooth_expression),
                         rr$smooth_expression,check.names=FALSE),paste0(stem,'_smoothed_expression'))
          csv(data.frame(treatment=tg,path=path,column=seq_along(rr$grid_pt),
                         local_pseudotime=rr$grid_pt),paste0(stem,'_column_pseudotime'))
          csv(data.frame(treatment=tg,path=path,cell=d$cell[ix],sample=d$sample[ix],
                         local_pseudotime=d$local_pseudotime[ix]),paste0(stem,'_cells'))
          writeLines(c(paste('Treatment:',tg),paste('Path:',path),paste('Averaging:',rr$averaging),
                       'Rows: hierarchical clustering, modules ordered by mean peak time, within-module clustering.',
                       'Columns: increasing local pseudotime; spline display grid is not independent observations.',
                       'Different groups may retain different genes; M1/M2/M3 are local module IDs, not shared gene sets.',
                       'No comparison across disconnected components; no claim of treatment-specific differential expression.'),
                     file.path(run_dir,paste0(stem,'_notes.txt')))
          module_tables[[job]]<-data.frame(job=job,treatment=tg,component=meta$component,path=path,
                                           root_node=meta$root_node,leaf_node=meta$leaf_node,module=unname(rr$modules),feature_id=names(rr$modules))
          jobs[[job]]<-list(treatment=tg,path=path,component=meta$component,universe=universe_features,
                            n_cells=length(ix),stem=stem)
          # 附带同组、同路径的代表基因平滑曲线，与热图采用完全一致的数值。
          markers<-if(length(MARKER_GENES))intersect(MARKER_GENES,rownames(rr$smooth_expression)) else
            head(rownames(rr$smooth_expression),6L)
          if(length(markers)) {
            td<-do.call(rbind,lapply(markers,function(gene)data.frame(gene=gene,
                                                                      pseudotime=rr$grid_pt,expression=as.numeric(rr$smooth_expression[gene,]))))
            p<-ggplot(td,aes(pseudotime,expression,color=gene))+geom_line(linewidth=0.7)+base_theme+
              labs(title=paste(tg,'|',path,'| representative gene trends'),
                   x='Local pseudotime',y='Smoothed log1p normalized expression')
            save_plot(p,paste0(stem,'_gene_trends'),10,6)
          }
          nrow(rr$z)
        },error=function(e)e)
        if(inherits(outcome,'error')) {
          report$status<-'skipped_or_failed';report$detail<-conditionMessage(outcome)
        } else {done<-done+1L;report$status<-'completed';report$n_heatmap_genes<-outcome}
        manifest[[job]]<-report
      }
      group_status[[tg]]<-data.frame(treatment=tg,n_completed_heatmaps=done)
      if(done==0L) {
        reason<-'No supported dynamic heatmap for this group. See heatmap_manifest.csv for cell counts and filtering reasons.'
        placeholder(paste(tg,'| pseudotime heatmap'),reason,file.path('B_heatmaps',tg,paste0('B_heatmap_',tg,'_NO_RESULT')))
      }
    }
    csv(do.call(rbind,manifest),'heatmap_manifest');csv(do.call(rbind,group_status),'heatmap_group_summary')
    if(!length(module_tables)) stop('没有通过数据支持阈值的处理/路径热图；各组已输出原因图。')
    modules<-do.call(rbind,module_tables);csv(modules,'29_group_path_gene_modules')
    
    if(RUN_GO) {
      for(pkg in c('AnnotationDbi','clusterProfiler','GO.db'))
        if(!requireNamespace(pkg,quietly=TRUE)) stop('GO扩展需要：',pkg)
      tarball <- ORGDB_TARBALL
      if(is.null(tarball)) {
        # 只扫描指定目录，不递归扫描整个home。
        search_dirs <- ORGDB_SEARCH_DIR
        candidates <- unlist(lapply(search_dirs,function(z)list.files(z,
                                                                      pattern='^org[.]G[.]max[.]eg[.]db.*[.]tar[.]gz$',full.names=TRUE)),use.names=FALSE)
        if(length(candidates)>1) stop('发现多个OrgDb安装包；请用ORGDB_TARBALL指定所需文件。')
        if(length(candidates)==1) tarball <- candidates[1]
      }
      # 使用本次输出目录的私有R库，避免覆盖服务器原有注释包。
      if(!is.null(tarball)) {
        if(!file.exists(tarball)) stop('注释包不存在：',tarball)
        lib <- file.path(run_dir,'R_library');dir.create(lib)
        if(ORGDB_PACKAGE %in% loadedNamespaces())
          stop('注释包已在当前R进程加载；请用独立Rscript进程运行，避免加载旧包。')
        install.packages(tarball,repos=NULL,type='source',lib=lib)
        .libPaths(c(lib,.libPaths()))
      }
      if(!requireNamespace(ORGDB_PACKAGE,quietly=TRUE)) stop('找不到org.G.max.eg.db及其本地安装包。')
      db <- get(ORGDB_PACKAGE,envir=asNamespace(ORGDB_PACKAGE))
      kt <- AnnotationDbi::keytypes(db)
      if(!KEYTYPE_OVERRIDE %in% kt) stop('该OrgDb不支持GID；不自动改用ENTREZID。')
      dbmeta <- AnnotationDbi::metadata(db);csv(dbmeta,'32_orgdb_metadata')
      dbkeys <- AnnotationDbi::keys(db,keytype=KEYTYPE_OVERRIDE)
      csv(data.frame(GID=head(dbkeys,100)),'33_orgdb_GID_examples')
      db_id_info <- parse_soybean_id(dbkeys)
      csv(db_id_info,'33b_orgdb_GID_cleaning')
      db_versions <- as.data.frame(table(version=ifelse(is.na(db_id_info$version),'unversioned',db_id_info$version)))
      csv(db_versions,'33c_orgdb_version_counts')
      provenance <- c(paste('Requested package:',ORGDB_PACKAGE),paste('Keytype:',KEYTYPE_OVERRIDE),
                      paste('Package version:',as.character(utils::packageVersion(ORGDB_PACKAGE))),
                      paste('Package path:',find.package(ORGDB_PACKAGE)),paste('Tarball:',if(is.null(tarball))'installed package' else tarball),
                      paste('Preferred ID version marker:',PREFERRED_DB_VERSION),
                      paste('GIDs containing preferred version marker:',sum(db_id_info$version==PREFERRED_DB_VERSION,na.rm=TRUE)),
                      'Package version and genome annotation version are different concepts.',
                      'IDs are cleaned for matching; original database GIDs are passed to enrichGO.',
                      'Version-marker cleanup is not a gene-coordinate liftover or an assembly conversion.')
      writeLines(provenance,file.path(run_dir,'34_annotation_provenance.txt'));cat(paste(provenance,collapse='\n'),'\n')
      idmap <- match_soybean_gid(rownames(cds),dbkeys,PREFERRED_DB_VERSION)
      if(!is.null(GENE_ID_MAP_TSV)) {
        mp <- read.delim(GENE_ID_MAP_TSV,stringsAsFactors=FALSE,check.names=FALSE)
        if(!all(c('feature_id','GID') %in% names(mp))) stop('映射表必须有feature_id、GID两列。')
        if(anyDuplicated(mp$feature_id)||anyDuplicated(mp$GID)) stop('请提供经过确认的一对一feature_id-GID映射。')
        idmap$GID <- mp$GID[match(idmap$feature_id,mp$feature_id)]
        idmap$match_status <- 'user_supplied_mapping'
        idmap$db_version <- parse_soybean_id(idmap$GID)$version
      }
      idmap$matched <- !is.na(idmap$GID)&idmap$GID %in% dbkeys
      csv(idmap,'35_feature_to_GID_audit')
      csv(as.data.frame(table(status=idmap$match_status)),'35b_ID_match_status_counts')
      gid_n <- as.data.frame(table(GID=idmap$GID[idmap$matched]));names(gid_n)[2]<-'n_expression_features'
      csv(gid_n[gid_n$n_expression_features>1,,drop=FALSE],'35c_multiple_features_per_GID')
      match_rate <- mean(idmap$matched);cat('GID match rate: ',round(100*match_rate,2),'%\n',sep='')
      if(match_rate<MIN_DB_MATCH_RATE) stop('GID匹配率过低；检查v4基因ID或提供正式映射表。')
      go_results<-list();go_audit<-list()
      for(job in names(jobs)) {
        info<-jobs[[job]];tg<-info$treatment;path<-info$path
        background<-unique(idmap$GID[idmap$matched & idmap$feature_id %in% info$universe])
        csv(data.frame(treatment=tg,path=path,GID=background),
            file.path('C_GO',tg,paste0('GO_',tg,'_',path,'_input_universe')))
        job_modules<-modules[modules$job==job,,drop=FALSE]
        for(mm in unique(job_modules$module)) {
          key<-paste(job,mm,sep='__');ids<-job_modules$feature_id[job_modules$module==mm]
          genes<-unique(idmap$GID[idmap$matched & idmap$feature_id %in% ids]);genes<-intersect(genes,background)
          audit<-data.frame(job=job,treatment=tg,component=info$component,path=path,module=mm,
                            n_features=length(ids),n_GID=length(genes),n_input_background=length(background),status='')
          csv(data.frame(treatment=tg,path=path,module=mm,GID=genes),
              file.path('C_GO',tg,paste0('GO_',tg,'_',path,'_',mm,'_input_genes')))
          if(length(genes)<5L||length(background)<10L) {
            audit$status<-'insufficient_mapped_genes';go_audit[[key]]<-audit;next
          }
          eg<-tryCatch(clusterProfiler::enrichGO(gene=genes,universe=background,OrgDb=db,
                                                 keyType=KEYTYPE_OVERRIDE,ont='BP',pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1,
                                                 minGSSize=10,maxGSSize=500,readable=FALSE),error=function(e)e)
          if(inherits(eg,'error')) {
            audit$status<-paste('failed:',conditionMessage(eg));go_audit[[key]]<-audit;next
          }
          rr<-if(is.null(eg))data.frame() else as.data.frame(eg)
          audit$status<-if(nrow(rr))'tested' else 'no_testable_GO_terms'
          go_audit[[key]]<-audit
          if(nrow(rr)) {
            rr$job<-job;rr$treatment<-tg;rr$component<-info$component;rr$path<-path;rr$module<-mm
            go_results[[key]]<-rr
            if(methods::is(eg,'enrichResult')) csv(data.frame(GID=eg@universe),
                                                   file.path('C_GO',tg,paste0('GO_',tg,'_',path,'_',mm,'_effective_universe')))
          }
        }
      }
      if(length(go_audit)) csv(do.call(rbind,go_audit),'37_GO_group_path_module_audit')
      go_all<-if(length(go_results))do.call(rbind,go_results) else data.frame()
      if(nrow(go_all)) {
        go_all$FDR_across_all_group_modules<-p.adjust(go_all$pvalue,method='BH')
        csv(go_all,'38_GO_BP_all_groups_all_tested')
      }
      for(tg in groups) {
        group_go<-if(nrow(go_all))go_all[go_all$treatment==tg,,drop=FALSE] else data.frame()
        csv(group_go,file.path('C_GO',tg,paste0('GO_BP_',tg,'_all_results')))
        tg_jobs<-names(jobs)[vapply(jobs,function(x)x$treatment==tg,logical(1))]
        if(!length(tg_jobs)) {
          placeholder(paste(tg,'| GO BP'), 'No supported heatmap module for this group; GO was not tested.',
                      file.path('C_GO',tg,paste0('GO_BP_',tg,'_NO_MODULE')));next
        }
        for(job in tg_jobs) {
          path<-jobs[[job]]$path
          jj<-if(nrow(group_go))group_go[group_go$job==job,,drop=FALSE] else data.frame()
          stem<-file.path('C_GO',tg,paste0('GO_BP_',tg,'_',path))
          csv(jj,paste0(stem,'_all_results'))
          hit<-if(nrow(jj))jj[is.finite(jj$FDR_across_all_group_modules)&
                                jj$FDR_across_all_group_modules<GO_FDR_CUTOFF,,drop=FALSE] else data.frame()
          csv(hit,paste0(stem,'_significant'))
          expected_modules<-unique(modules$module[modules$job==job])
          passing_modules<-if(nrow(hit))unique(hit$module) else character(0)
          for(mm in setdiff(expected_modules,passing_modules)) {
            placeholder(paste(tg,'|',path,'|',mm,'| GO BP'),
                        'No significant GO BP terms for this module. Check the mapping and GO audit tables for details.',
                        paste0(stem,'_',mm,'_NO_SIGNIFICANT'))
          }
          if(!nrow(hit)) {
            placeholder(paste(tg,'|',path,'| GO BP'),
                        'No GO BP terms passed the global FDR threshold. See the result and mapping audit tables.',
                        paste0(stem,'_NO_SIGNIFICANT'));next
          }
          top<-do.call(rbind,lapply(split(hit,hit$module),function(z)
            head(z[order(z$FDR_across_all_group_modules,-z$Count),,drop=FALSE],GO_TOP_TERMS_PER_MODULE)))
          top$term_key<-paste(top$module,top$ID,sep='::')
          top$term_key<-factor(top$term_key,levels=rev(unique(top$term_key)))
          label_map<-setNames(vapply(top$Description,function(s)paste(strwrap(s,48),collapse='\n'),character(1)),
                              as.character(top$term_key))
          make_go_plot<-function(z,title) ggplot(z,aes(Count,term_key,fill=FDR_across_all_group_modules))+
            geom_col(width=0.7)+scale_y_discrete(labels=label_map)+
            scale_fill_gradient(low='#D6604D',high='#4393C3',name='Global FDR')+
            base_theme+theme(axis.text.y=element_text(size=10))+
            labs(title=title,subtitle='GO BP | soybean GID | treatment-specific descriptive modules',
                 x='Mapped gene count',y=NULL)
          p<-make_go_plot(top,paste(tg,'|',path,'| GO BP'))+
            facet_wrap(~module,ncol=1,scales='free_y')
          save_plot(p,stem,12,max(6,0.42*nrow(top)+1.3*length(unique(top$module))))
          for(mm in unique(top$module)) {
            zz<-top[top$module==mm,,drop=FALSE]
            save_plot(make_go_plot(zz,paste(tg,'|',path,'|',mm,'| GO BP')),
                      paste0(stem,'_',mm),11,max(4,0.5*nrow(zz)+2))
          }
        }
      }
    }
  },error=function(e) {
    msg<-paste('Core figures saved. Group heatmap/GO extension failed:',conditionMessage(e))
    writeLines(msg,file.path(run_dir,'EXTENSION_FAILED.txt'));cat(msg,'\n')
  })
  cat('Finished. Output: ',run_dir,'\n',sep='')
  invisible(run_dir)
}

main()

##GO
# -*- coding: UTF-8 -*-
# v4：从已完成的v3结果续跑GO并重绘。独立完整脚本，无需先source旧版。
# 运行：source('soybean_endodermis_resume_GO_replot_v4.txt', encoding='UTF-8')
# 或：Rscript soybean_endodermis_resume_GO_replot_v4.txt
# 本次GO空目录的已知原因：v3在发现OrgDb已加载时主动stop，尚未执行富集。
# v4优先复用已加载的同名OrgDb，不再因已加载而报错，不强制重新安装。
# 不读原始Seurat对象，不重新降维/聚类/建图/定根/graph_test/平滑/聚类模块。
# 原输出目录不变，无新时间戳目录；重绘覆盖同名图，已保存的CDS/CSV分析输入不改。

RUN_DIR <- '/data01/liuxy/Project/CK2_added_v3_Results/Endodermis_monocle3_originalUMAP_grouped_v3/20261119_163425'
REDRAW_TRAJECTORY <- TRUE
REDRAW_HEATMAPS <- TRUE
REDRAW_SUMMARIES <- TRUE
RESUME_GO <- TRUE
ORGDB_PACKAGE <- 'org.G.max.eg.db'
KEYTYPE_OVERRIDE <- 'GID'
PREFERRED_DB_VERSION <- 'v4'
ORGDB_SEARCH_DIR <- '/home/liuxy'
ORGDB_TARBALL <- NULL
GENE_ID_MAP_TSV <- NULL                 # 可选：feature_id、GID两列正式一对一映射
MIN_DB_MATCH_RATE <- 0.50
MIN_DETECTED_FRACTION <- 0.10           # 与v3一致，用于恢复各路径的GO检测背景
GO_MIN_GENES <- 5L
GO_MIN_SET_SIZE <- 10L
GO_MAX_SET_SIZE <- 500L
GO_FDR_CUTOFF <- 0.05                   # 不为出图而放宽显著性阈值
GO_TOP_TERMS <- 6L
PT_PALETTE <- 'plasma'                 # 紫→红橙→黄，替换不明显的单一蓝色
PT_COLOR_TRANSFORM <- 'sqrt'           # 只改变颜色映射，拟时序数值不改；图例标原始值
# 所有处理使用同一全局色标，绝不按每组各自归一化。另存linear线性色标版本作对照。
POINT_SIZE <- 0.55
BACKGROUND_ALPHA <- 0.30
SEED <- 20261119L

# __ID_HELPERS__

# 可测试的包加载分支：已加载 > 已安装 > 从本地tar.gz安装。
orgdb_load_route <- function(is_loaded,is_available) {
  if(is_loaded) 'reuse_loaded_namespace' else if(is_available) 'load_installed_package' else 'install_local_package'
}

# 从已保存的counts与graph_test恢复相同的检测背景；不重新做graph_test。
restore_universe <- function(count_matrix,cell_ids,tested_features,min_fraction=0.10) {
  if(!length(cell_ids)||anyDuplicated(cell_ids)) stop('路径细胞列表为空或含重复barcode。')
  if(!all(cell_ids %in% colnames(count_matrix))) stop('保存的路径细胞不全在CDS中，拒绝混用运行结果。')
  cc <- count_matrix[,cell_ids,drop=FALSE]
  n_detect <- if(inherits(cc,'Matrix')) Matrix::rowSums(cc>0) else rowSums(cc>0)
  detectable <- rownames(count_matrix)[n_detect>=max(10L,ceiling(min_fraction*length(cell_ids)))]
  intersect(tested_features,detectable)
}

main_resume <- function() {
  if(!dir.exists(RUN_DIR)) stop('找不到已完成的输出目录：',RUN_DIR)
  run_dir <- normalizePath(RUN_DIR,winslash='/',mustWork=TRUE)
  logcon <- file(file.path(run_dir,'resume_replot_v4.log'),'wt',encoding='UTF-8')
  sink(logcon,split=TRUE);on.exit({sink();close(logcon)},add=TRUE)
  on.exit(writeLines(capture.output(sessionInfo()),file.path(run_dir,'sessionInfo_v4.txt')),add=TRUE)
  cat('Using saved run, same output directory: ',run_dir,'\n',sep='')
  need <- function(pkgs) {
    bad <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
    if(length(bad)) stop('缺少依赖：',paste(bad,collapse=', '))
  }
  need('ggplot2');suppressPackageStartupMessages(library(ggplot2))
  set.seed(SEED)
  groups <- c('CK','ABA','eBL','IAA','MeJA','SA')
  group_colors <- c(CK='#4D4D4D',ABA='#D55E00',eBL='#009E73',IAA='#0072B2',MeJA='#CC79A7',SA='#E69F00')
  read_table <- function(name,required=TRUE) {
    f <- file.path(run_dir,name)
    if(!file.exists(f)) {
      if(required) stop('缺少已保存文件：',f)
      return(NULL)
    }
    read.csv(f,check.names=FALSE,stringsAsFactors=FALSE,na.strings=c('','NA'),fileEncoding='UTF-8')
  }
  csv <- function(x,name) {
    path <- file.path(run_dir,paste0(name,'.csv'))
    dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
    write.csv(x,path,row.names=FALSE,na='',fileEncoding='UTF-8')
  }
  theme_pub <- theme_classic(base_size=12)+theme(
    plot.title=element_text(size=15,face='bold'),plot.subtitle=element_text(size=10),
    plot.caption=element_text(size=8),strip.background=element_blank(),
    strip.text=element_text(size=13,face='bold'),plot.margin=margin(12,16,12,12))
  save_plot <- function(p,stem,w=11,h=7) {
    path <- file.path(run_dir,stem);dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
    ggsave(paste0(path,'.pdf'),plot=p,width=w,height=h,limitsize=FALSE,bg='white')
    ggsave(paste0(path,'.png'),plot=p,width=w,height=h,dpi=300,limitsize=FALSE,bg='white')
  }
  placeholder <- function(title,reason,stem) {
    p <- ggplot()+annotate('text',x=0,y=0,label=paste(strwrap(reason,78),collapse='\n'),size=4)+
      xlim(-1,1)+ylim(-1,1)+theme_void()+labs(title=title)+
      theme(plot.title=element_text(face='bold',hjust=0.5,size=15))
    save_plot(p,stem,11,5)
  }
  stages <- list()
  run_stage <- function(name,fun) {
    cat('\n--- ',name,' ---\n',sep='')
    out <- tryCatch({fun();list(status='completed',detail='')},error=function(e)
      list(status='failed',detail=conditionMessage(e)))
    stages[[name]] <<- data.frame(stage=name,status=out$status,detail=out$detail)
    csv(do.call(rbind,stages),'resume_replot_v4_status')
    cat(name,': ',out$status,' ',out$detail,'\n',sep='')
    invisible(out)
  }
  d <- read_table('13_cell_metadata_pseudotime_states.csv')
  for(nm in c('cell','x','y','treatment','component','pseudotime','local_pseudotime'))
    if(!nm %in% names(d)) stop('保存的细胞表缺少列：',nm)
  d$treatment <- factor(d$treatment,levels=groups)
  if(anyNA(d$treatment)||anyDuplicated(d$cell)) stop('细胞表有未知处理或重复barcode。')
  edges <- read_table('09_graph_edges.csv')
  nodes <- read_table('08_principal_nodes.csv')
  roots <- read_table('11_root_and_interval_definition.csv')
  local_roots <- read_table('11b_roots_for_separate_component_heatmaps.csv',required=FALSE)
  target <- unique(roots$component)
  if(length(target)!=1L) stop('主分量定义不唯一。')
  d$pt_plot <- ifelse(d$component==target & is.finite(d$pseudotime),d$pseudotime,NA_real_)
  edges$component <- nodes$component[match(edges$from,nodes$node)]
  if(anyNA(edges$component)) stop('保存的边表与节点表不匹配。')
  color_scale <- function(limits,transform=PT_COLOR_TRANSFORM,label='Pseudotime') {
    if(any(!is.finite(limits))||diff(limits)<=0) stop('拟时序颜色范围无效。')
    scale_color_viridis_c(option=PT_PALETTE,begin=0,end=0.97,limits=limits,trans=transform,
                          name=label,na.value='#BDBDBD',guide=guide_colorbar(barheight=grid::unit(42,'mm'),
                                                                             barwidth=grid::unit(4,'mm')))
  }
  
  # ==================== 1. 仅用已保存坐标与拟时序重绘 ====================
  if(REDRAW_TRAJECTORY) run_stage('trajectory_replot',function() {
    finite <- d$pt_plot[is.finite(d$pt_plot)]
    if(!length(finite)||min(finite)<0) stop('主拟时序为空或存在负值。')
    limits <- range(finite)
    csv(data.frame(palette=PT_PALETTE,transform=PT_COLOR_TRANSFORM,
                   min_pseudotime=limits[1],max_pseudotime=limits[2],shared_across_groups=TRUE),'v4_pseudotime_color_settings')
    facet_edges <- do.call(rbind,lapply(groups,function(tg) {
      ee<-edges[edges$component %in% unique(d$component[d$treatment==tg]),,drop=FALSE]
      ee$treatment<-factor(rep(tg,nrow(ee)),levels=groups);ee
    }))
    make_a <- function(transform) ggplot(d,aes(x,y))+
      geom_point(data=d[,c('x','y')],aes(x,y),inherit.aes=FALSE,
                 color='grey75',alpha=BACKGROUND_ALPHA,size=0.35)+
      geom_segment(data=facet_edges,aes(x=x,y=y,xend=xend,yend=yend),inherit.aes=FALSE,
                   color='grey30',linewidth=0.35)+
      geom_point(aes(color=pt_plot),size=POINT_SIZE,alpha=0.95)+color_scale(limits,transform)+
      geom_point(data=nodes[nodes$node %in% roots$root_node,,drop=FALSE],aes(x,y),
                 inherit.aes=FALSE,shape=23,fill='white',color='black',size=3,stroke=0.9)+
      coord_equal()+theme_pub+labs(title='Endodermis | original UMAP',
                                   subtitle=paste0('Saved trajectory | primary component ',target,' | ',transform,' color scale'),
                                   x='Original UMAP 1',y='Original UMAP 2',
                                   caption='Same color scale for all groups; legend shows original pseudotime. Root direction remains exploratory unless validated.')
    p<-make_a(PT_COLOR_TRANSFORM)
    save_plot(p+facet_wrap(~treatment,ncol=3,drop=FALSE),'A_pseudotime_6groups',14,9)
    save_plot(p+facet_wrap(~treatment,ncol=1,drop=FALSE),'A_pseudotime_vertical',8,22)
    save_plot(make_a('identity')+facet_wrap(~treatment,ncol=3,drop=FALSE),
              'A_pseudotime_6groups_plasma_linear',14,9)
    # 每个处理单独保存，所有图仍使用同一坐标范围与同一色标。
    for(tg in groups) {
      z<-d[d$treatment==tg,,drop=FALSE]
      ee<-edges[edges$component %in% unique(z$component),,drop=FALSE]
      one<-ggplot(z,aes(x,y,color=pt_plot))+
        geom_point(data=d[,c('x','y')],aes(x,y),inherit.aes=FALSE,
                   color='grey75',alpha=BACKGROUND_ALPHA,size=0.35)+
        geom_segment(data=ee,aes(x=x,y=y,xend=xend,yend=yend),inherit.aes=FALSE,
                     color='grey30',linewidth=0.35)+geom_point(size=POINT_SIZE,alpha=0.95)+
        color_scale(limits)+coord_equal(xlim=range(d$x),ylim=range(d$y))+theme_pub+
        labs(title=paste(tg,'| Endodermis'),subtitle=paste(PT_COLOR_TRANSFORM,'color scale; saved pseudotime'),
             x='Original UMAP 1',y='Original UMAP 2')
      save_plot(one,paste0('A_pseudotime_',tg),8,7)
    }
    if(!is.null(local_roots)) for(cc in local_roots$component) {
      z<-d[d$component==cc & is.finite(d$local_pseudotime),,drop=FALSE]
      if(nrow(z)<10L||diff(range(z$local_pseudotime))<=0) next
      ee<-edges[edges$component==cc,,drop=FALSE]
      p<-ggplot(z,aes(x,y,color=local_pseudotime))+
        geom_segment(data=ee,aes(x=x,y=y,xend=xend,yend=yend),inherit.aes=FALSE,
                     color='grey30',linewidth=0.35)+geom_point(size=POINT_SIZE)+
        color_scale(range(z$local_pseudotime),label='Local PT')+
        facet_wrap(~treatment,ncol=3,drop=FALSE)+coord_equal()+theme_pub+
        labs(title=paste('Endodermis | component',cc),
             subtitle='Saved local pseudotime; compare treatments within this component only',
             x='Original UMAP 1',y='Original UMAP 2')
      save_plot(p,paste0('A_local_component_',cc,'_6groups'),13,8)
    }
  })
  
  # ==================== 2. 从Z-score矩阵重绘，不再平滑或聚类 ====================
  modules <- read_table('29_group_path_gene_modules.csv',required=FALSE)
  if(REDRAW_HEATMAPS) run_stage('heatmap_replot',function() {
    need('pheatmap')
    files<-list.files(file.path(run_dir,'B_heatmaps'),pattern='_rowZ_matrix[.]csv$',recursive=TRUE,full.names=TRUE)
    if(!length(files)) stop('没有保存的_rowZ_matrix.csv，请查看旧版heatmap_manifest.csv。')
    checks<-list()
    for(f in files) {
      status<-tryCatch({
        ztab<-read.csv(f,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding='UTF-8')
        tg<-unique(ztab$treatment);path<-unique(ztab$path)
        if(length(tg)!=1L||!tg %in% groups||length(path)!=1L) stop('矩阵组名/路径不唯一。')
        numeric_cols<-setdiff(names(ztab),c('treatment','path','feature_id','module'))
        z<-as.matrix(ztab[,numeric_cols,drop=FALSE]);storage.mode(z)<-'double'
        rownames(z)<-ztab$feature_id
        if(any(!is.finite(z))) stop('Z-score矩阵有非有限值。')
        stem<-sub('_rowZ_matrix[.]csv$','',f)
        axis<-read.csv(paste0(stem,'_column_pseudotime.csv'),stringsAsFactors=FALSE)
        if(nrow(axis)!=ncol(z)) stop('保存的热图列和拟时序轴长度不一致。')
        if(any(diff(axis$local_pseudotime)<=0)) stop('热图列轴不是递增拟时序。')
        mods<-as.character(ztab$module);lev<-unique(mods)
        gaps<-head(cumsum(rle(mods)$lengths),-1L)
        # 行与列直接沿用保存的聚类顺序；只更新注释与色板。
        z[]<-pmax(-2.5,pmin(2.5,z));colnames(z)<-sprintf('PT%03d',seq_len(ncol(z)))
        annrow<-data.frame(Module=factor(mods,levels=lev),row.names=rownames(z))
        anncol<-data.frame(Treatment=factor(rep(tg,ncol(z)),levels=groups),
                           Pseudotime=axis$local_pseudotime,row.names=colnames(z))
        tick<-unique(round(seq(1,ncol(z),length.out=5)));labels<-rep('',ncol(z))
        labels[tick]<-format(round(axis$local_pseudotime[tick],2),trim=TRUE)
        ph<-pheatmap::pheatmap(z,cluster_rows=FALSE,cluster_cols=FALSE,gaps_row=gaps,
                               annotation_row=annrow,annotation_col=anncol,
                               annotation_colors=list(Treatment=group_colors,
                                                      Module=setNames(colorRampPalette(c('#E69F00','#009E73','#CC79A7'))(length(lev)),lev),
                                                      Pseudotime=c('#0D0887','#CC4778','#F0F921')),
                               color=colorRampPalette(c('#2166AC','#F7F7F7','#B2182B'))(100),
                               breaks=seq(-2.5,2.5,length.out=101),legend_breaks=c(-2,-1,0,1,2),
                               labels_col=labels,angle_col='0',show_rownames=FALSE,show_colnames=TRUE,
                               annotation_names_row=FALSE,border_color=NA,fontsize=10,
                               main=paste0(tg,' | ',path,'\nSaved clustered trends | row Z-score'),silent=TRUE)
        grDevices::pdf(paste0(stem,'.pdf'),width=9,height=9)
        tryCatch({grid::grid.newpage();grid::grid.draw(ph$gtable)},finally=grDevices::dev.off())
        grDevices::png(paste0(stem,'.png'),width=2700,height=2700,res=300,
                       type=if(capabilities('cairo'))'cairo' else getOption('bitmapType','Xlib'))
        tryCatch({grid::grid.newpage();grid::grid.draw(ph$gtable)},finally=grDevices::dev.off())
        'completed'
      },error=function(e)paste('failed:',conditionMessage(e)))
      checks[[f]]<-data.frame(matrix_file=basename(f),status=status)
    }
    csv(do.call(rbind,checks),'v4_heatmap_redraw_audit')
    if(any(vapply(checks,function(x)x$status!='completed',logical(1))))
      stop('部分缓存热图重绘失败，详见v4_heatmap_redraw_audit.csv；GO仍将继续。')
  })
  
  if(REDRAW_SUMMARIES) run_stage('replicate_summary_replot',function() {
    st<-read_table('21_per_replicate_interval_proportions.csv',FALSE)
    if(!is.null(st)) {
      st$treatment<-factor(st$treatment,levels=groups)
      p<-ggplot(st,aes(treatment,proportion_all,color=treatment))+
        geom_point(position=position_jitter(width=.08,height=0,seed=SEED),size=3)+
        facet_wrap(~interval,nrow=1)+scale_color_manual(values=group_colors)+theme_pub+
        theme(legend.position='none')+labs(x=NULL,y='Proportion per biological replicate')
      save_plot(p,'23_replicate_interval_proportions',14,5)
    }
    ss<-read_table('22_per_replicate_pseudotime_summary.csv',FALSE)
    if(!is.null(ss)) {
      ss$treatment<-factor(ss$treatment,levels=groups)
      p<-ggplot(ss,aes(treatment,median_pseudotime,color=treatment))+
        geom_point(position=position_jitter(width=.08,height=0,seed=SEED),size=3)+
        scale_color_manual(values=group_colors)+theme_pub+theme(legend.position='none')+
        labs(x=NULL,y='Saved median pseudotime per biological replicate')
      save_plot(p,'24_replicate_median_pseudotime')
    }
  })
  
  # ==================== 3. 修复OrgDb加载，补做此前未执行的GO ====================
  if(RESUME_GO) run_stage('GO_resume',function() {
    need(c('AnnotationDbi','clusterProfiler','GO.db','monocle3','SummarizedExperiment','Matrix'))
    if(is.null(modules)||!nrow(modules)) stop('没有保存的29_group_path_gene_modules.csv，不能凭空生成GO模块。')
    for(nm in c('job','treatment','path','component','module','feature_id'))
      if(!nm %in% names(modules)) stop('模块表缺少：',nm)
    loaded<-ORGDB_PACKAGE %in% loadedNamespaces()
    lib<-file.path(run_dir,'R_library')
    if(!loaded && dir.exists(lib)) .libPaths(unique(c(lib,.libPaths())))
    available<-if(loaded)TRUE else requireNamespace(ORGDB_PACKAGE,quietly=TRUE)
    route<-orgdb_load_route(loaded,available)
    cat('OrgDb route: ',route,'\n',sep='')
    if(route=='install_local_package') {
      tb<-ORGDB_TARBALL
      if(is.null(tb)) {
        candidates<-list.files(ORGDB_SEARCH_DIR,pattern='^org[.]G[.]max[.]eg[.]db.*[.]tar[.]gz$',full.names=TRUE)
        if(length(candidates)!=1L) stop('OrgDb未安装且安装包不唯一/不存在；请填写ORGDB_TARBALL。')
        tb<-candidates[1]
      }
      if(!file.exists(tb)) stop('安装包不存在：',tb)
      dir.create(lib,showWarnings=FALSE,recursive=TRUE)
      install.packages(tb,repos=NULL,type='source',lib=lib)
      .libPaths(unique(c(lib,.libPaths())))
      if(!requireNamespace(ORGDB_PACKAGE,quietly=TRUE)) stop('本地OrgDb安装后仍无法加载。')
    }
    # 已加载包直接复用，其命名空间路径才是实际生效的包路径。
    ns<-asNamespace(ORGDB_PACKAGE);db<-get(ORGDB_PACKAGE,envir=ns)
    if(!methods::is(db,'OrgDb')) stop('同名对象不是OrgDb。')
    if(!KEYTYPE_OVERRIDE %in% AnnotationDbi::keytypes(db)) stop('该OrgDb没有GID键类型。')
    dbkeys<-AnnotationDbi::keys(db,keytype=KEYTYPE_OVERRIDE)
    csv(AnnotationDbi::metadata(db),'32_orgdb_metadata')
    writeLines(c(paste('Load route:',route),paste('Actual namespace path:',getNamespaceInfo(ns,'path')),
                 paste('Loaded package version:',as.character(getNamespaceVersion(ns))),
                 paste('Keytype:',KEYTYPE_OVERRIDE),paste('Preferred ID marker:',PREFERRED_DB_VERSION)),
               file.path(run_dir,'34_annotation_provenance_v4.txt'))
    cds_path<-file.path(run_dir,'endodermis_monocle3_ordered.rds')
    if(!file.exists(cds_path)) stop('缺少已保存的ordered CDS，无法恢复GO背景。')
    cat('Loading saved CDS only: ',cds_path,'\n',sep='')
    cds<-readRDS(cds_path)
    counts<-SummarizedExperiment::assay(cds,'counts')
    if(!setequal(colnames(counts),d$cell)) stop('保存的CDS与细胞结果表来自不同运行。')
    gt<-read_table('25_graph_test_pooled_candidate_screen.csv')
    if(!all(c('feature_id','q_value','morans_I') %in% names(gt))) stop('保存的graph_test表缺少关键列。')
    tested<-gt$feature_id[is.finite(gt$q_value)&is.finite(gt$morans_I)]
    idmap<-match_soybean_gid(rownames(counts),dbkeys,PREFERRED_DB_VERSION)
    if(!is.null(GENE_ID_MAP_TSV)) {
      mp<-read.delim(GENE_ID_MAP_TSV,stringsAsFactors=FALSE,check.names=FALSE)
      if(!all(c('feature_id','GID') %in% names(mp))||anyDuplicated(mp$feature_id)||anyDuplicated(mp$GID))
        stop('映射表必须是feature_id/GID两列的一对一正式映射。')
      idmap$GID<-mp$GID[match(idmap$feature_id,mp$feature_id)];idmap$match_status<-'user_mapping'
    }
    idmap$matched<-!is.na(idmap$GID)&idmap$GID %in% dbkeys
    csv(idmap,'35_feature_to_GID_audit');csv(as.data.frame(table(status=idmap$match_status)),'35b_ID_match_status_counts')
    match_rate<-mean(idmap$matched);cat('GID match rate: ',round(match_rate*100,2),'%\n',sep='')
    if(match_rate<MIN_DB_MATCH_RATE) stop('GID匹配率过低；已保存匹配审计，不能可靠进行GO。')
    map_ids<-function(ids)unique(idmap$GID[idmap$matched & idmap$feature_id %in% ids])
    results<-list();audit_rows<-list();jobs<-unique(modules[,c('job','treatment','path','component')])
    if(anyDuplicated(jobs$job)) stop('同一个job对应多个处理/路径，拒绝混用模块。')
    path_nodes<-read_table('26c_path_node_membership.csv',FALSE)
    for(i in seq_len(nrow(jobs))) {
      info<-jobs[i,,drop=FALSE];job<-info$job;tg<-info$treatment;pa<-info$path
      if(!tg %in% groups) stop('模块表包含未知处理组：',tg)
      cat('GO: ',tg,' | ',pa,'\n',sep='')
      stem<-file.path('B_heatmaps',tg,paste0('B_heatmap_',tg,'_',pa))
      cells<-read_table(paste0(stem,'_cells.csv'),FALSE)
      if(!is.null(cells)) {
        if(any(cells$treatment!=tg)||any(cells$path!=pa)) stop('路径细胞表标签不匹配：',job)
        cell_ids<-cells$cell
      } else {
        if(is.null(path_nodes)) stop('缺少路径细胞表及路径节点表：',job)
        pn<-path_nodes$node[path_nodes$path==pa]
        cell_ids<-d$cell[d$treatment==tg & d$nearest_node %in% pn & is.finite(d$local_pseudotime)]
      }
      if(!all(d$treatment[match(cell_ids,d$cell)]==tg)) stop('路径细胞混入其他处理：',job)
      universe_features<-restore_universe(counts,cell_ids,tested,MIN_DETECTED_FRACTION)
      background<-map_ids(universe_features)
      # 用长度一致的向量写空表，避免零匹配时data.frame报“行数不一致”。
      csv(data.frame(GID=background),file.path('C_GO',tg,paste0('GO_',tg,'_',pa,'_input_universe')))
      job_modules<-modules[modules$job==job,,drop=FALSE]
      for(mm in unique(job_modules$module)) {
        key<-paste(job,mm,sep='__');ids<-job_modules$feature_id[job_modules$module==mm]
        genes<-intersect(map_ids(ids),background)
        ar<-data.frame(job=job,treatment=tg,component=info$component,path=pa,module=mm,
                       n_features=length(ids),n_GID=length(genes),n_input_background=length(background),
                       n_tested_terms=0L,n_significant_global=0L,status='',detail='',stringsAsFactors=FALSE)
        csv(data.frame(GID=genes),file.path('C_GO',tg,paste0('GO_',tg,'_',pa,'_',mm,'_input_genes')))
        if(length(genes)<GO_MIN_GENES||length(background)<10L) {
          ar$status<-'insufficient_mapped_genes';audit_rows[[key]]<-ar;next
        }
        eg<-tryCatch(clusterProfiler::enrichGO(gene=genes,universe=background,OrgDb=db,
                                               keyType=KEYTYPE_OVERRIDE,ont='BP',pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1,
                                               minGSSize=GO_MIN_SET_SIZE,maxGSSize=GO_MAX_SET_SIZE,readable=FALSE),error=function(e)e)
        if(inherits(eg,'error')) {
          ar$status<-'enrichment_error';ar$detail<-conditionMessage(eg);audit_rows[[key]]<-ar;next
        }
        # 保存完整检验表，避免as.data.frame按对象内阈值再次隐藏行。
        rr<-if(is.null(eg))data.frame() else eg@result
        ar$status<-if(nrow(rr))'tested' else 'no_testable_GO_terms';ar$n_tested_terms<-nrow(rr)
        audit_rows[[key]]<-ar
        if(nrow(rr)) {
          rr$job<-job;rr$treatment<-tg;rr$component<-info$component;rr$path<-pa;rr$module<-mm
          results[[key]]<-rr
          saveRDS(eg,file.path(run_dir,'C_GO',tg,paste0('GO_',tg,'_',pa,'_',mm,'_enrichResult.rds')))
          csv(data.frame(GID=eg@universe),file.path('C_GO',tg,paste0('GO_',tg,'_',pa,'_',mm,'_effective_universe')))
        }
      }
      csv(do.call(rbind,audit_rows),'37_GO_group_path_module_audit')
    }
    if(length(results)) {
      allgo<-do.call(rbind,results)
      allgo$FDR_across_all_group_modules<-p.adjust(allgo$pvalue,method='BH')
      allgo$significant_global<-is.finite(allgo$FDR_across_all_group_modules)&
        allgo$FDR_across_all_group_modules<GO_FDR_CUTOFF
    } else allgo<-data.frame(ID=character(),Description=character(),pvalue=numeric(),p.adjust=numeric(),
                             Count=integer(),job=character(),treatment=character(),component=integer(),path=character(),
                             module=character(),FDR_across_all_group_modules=numeric(),significant_global=logical())
    csv(allgo,'38_GO_BP_all_groups_all_tested')
    csv(allgo[allgo$significant_global,,drop=FALSE],'39_GO_BP_all_groups_significant')
    final_audit<-do.call(rbind,audit_rows)
    for(i in seq_len(nrow(final_audit))) {
      ar<-final_audit[i,]
      n<-sum(allgo$job==ar$job & allgo$module==ar$module & allgo$significant_global)
      final_audit$n_significant_global[i]<-n
      if(final_audit$status[i]=='tested' && n==0L) final_audit$status[i]<-'tested_no_global_significance'
    }
    csv(final_audit,'37_GO_group_path_module_audit')
    draw_go<-function(z,tg,pa,title_suffix='',exploratory=FALSE) {
      top<-do.call(rbind,lapply(split(z,z$module),function(x)
        head(x[order(x$FDR_across_all_group_modules,x$pvalue,-x$Count),,drop=FALSE],GO_TOP_TERMS)))
      top$key<-paste(top$module,top$ID,sep='::');top$key<-factor(top$key,levels=rev(unique(top$key)))
      labels<-setNames(vapply(top$Description,function(s)paste(strwrap(s,48),collapse='\n'),character(1)),as.character(top$key))
      top$FDR_strength <- -log10(pmax(top$FDR_across_all_group_modules,1e-300))
      ggplot(top,aes(Count,key,fill=FDR_strength))+geom_col(width=.72)+
        scale_y_discrete(labels=labels)+scale_fill_viridis_c(option='magma',direction=-1,name='-log10(global FDR)')+
        facet_wrap(~module,ncol=1,scales='free_y')+theme_pub+
        labs(title=paste(tg,'|',pa,'| GO BP',title_suffix),
             subtitle=if(exploratory)'Top tested terms only; NOT significant at global FDR < 0.05' else
               paste('Global FDR <',GO_FDR_CUTOFF,'| soybean GID'),x='Mapped gene count',y=NULL)
    }
    for(tg in groups) {
      group_go<-allgo[allgo$treatment==tg,,drop=FALSE]
      csv(group_go,file.path('C_GO',tg,paste0('GO_BP_',tg,'_all_results')))
      jtg<-jobs[jobs$treatment==tg,,drop=FALSE]
      if(!nrow(jtg)) placeholder(paste(tg,'| GO BP'),'No saved heatmap modules for this group.',
                                 file.path('C_GO',tg,paste0('GO_BP_',tg,'_NO_MODULE')))
      for(i in seq_len(nrow(jtg))) {
        job<-jtg$job[i];pa<-jtg$path[i];stem<-file.path('C_GO',tg,paste0('GO_BP_',tg,'_',pa))
        jj<-group_go[group_go$job==job,,drop=FALSE];hit<-jj[jj$significant_global,,drop=FALSE]
        csv(jj,paste0(stem,'_all_results'));csv(hit,paste0(stem,'_significant'))
        if(nrow(hit)) {
          nn<-sum(pmin(table(hit$module),GO_TOP_TERMS))
          save_plot(draw_go(hit,tg,pa),stem,12,max(6,.48*nn+1.3*length(unique(hit$module))))
        } else {
          aa<-final_audit[final_audit$job==job,,drop=FALSE]
          reason<-paste('No significant GO BP result. Module status:',paste(unique(aa$status),collapse='; '))
          placeholder(paste(tg,'|',pa,'| GO BP'),reason,stem)
          # 若真实检验完成但不显著，可另外看明确标注“未达显著”的描述性排序图。
          if(nrow(jj)) save_plot(draw_go(jj,tg,pa,exploratory=TRUE),paste0(stem,'_top_tested_NOT_SIGNIFICANT'),
                                 12,max(6,.48*sum(pmin(table(jj$module),GO_TOP_TERMS))+1.3*length(unique(jj$module))))
        }
        for(mm in unique(modules$module[modules$job==job])) {
          hh<-hit[hit$module==mm,,drop=FALSE]
          if(nrow(hh)) save_plot(draw_go(hh,tg,pa,paste('|',mm)),paste0(stem,'_',mm),
                                 11,max(4,.5*min(nrow(hh),GO_TOP_TERMS)+2)) else {
                                   aa<-final_audit[final_audit$job==job & final_audit$module==mm,,drop=FALSE]
                                   placeholder(paste(tg,'|',pa,'|',mm,'| GO BP'),
                                               paste('No significant GO result. Status:',paste(unique(aa$status),collapse='; ')),paste0(stem,'_',mm))
                                 }
        }
      }
    }
    writeLines(c('The previous GO step stopped before enrichment because the OrgDb was already loaded.',
                 paste('v4 load route:',route),paste('Tested GO rows:',nrow(allgo)),
                 paste('Global FDR significant rows:',sum(allgo$significant_global)),
                 'Check 37_GO_group_path_module_audit.csv to distinguish errors, mapping shortages and true lack of significance.',
                 'Files ending _NO_SIGNIFICANT from v3 are historical; use the freshly written canonical GO_BP_<group>_<path> plots and v4 audit.',
                 'No significance threshold was relaxed. Descriptive top-tested plots are explicitly labelled NOT SIGNIFICANT.'),
               file.path(run_dir,'GO_resume_v4_summary.txt'))
    if(any(final_audit$status=='enrichment_error')) stop('部分模块GO仍有计算错误；其他结果已保存，详见37审计表。')
    cat('GO enrichment completed: ',nrow(allgo),' tested rows; ',sum(allgo$significant_global),' significant rows.\n',sep='')
  })
  writeLines(c(paste('Existing run directory:',run_dir),
               'Saved coordinates, pseudotimes, counts, graph-test results and gene modules were reused.',
               'No embedding, graph, root, graph_test, smoothing or clustering was recomputed.',
               paste('Primary pseudotime palette:',PT_PALETTE,'; color transform:',PT_COLOR_TRANSFORM),
               'All groups share the same scale. Sqrt transforms colors only; legend ticks remain original pseudotime values.',
               'A_pseudotime_6groups_plasma_linear is the linear-scale comparison plot.',
               'GO was newly computed because it did not run successfully in v3.',
               'Check resume_replot_v4_status.csv for actual completion of each stage.'),
             file.path(run_dir,'REPlot_v4_README.txt'))
  if(any(vapply(stages,function(x)x$status=='failed',logical(1))))
    stop('续跑有部分步骤失败；已完成图表仍保存于原目录，查看resume_replot_v4_status.csv。',call.=FALSE)
  cat('\nFinished in the SAME directory: ',run_dir,'\n',sep='')
  invisible(run_dir)
}

main_resume()



##11.21
# -*- coding: UTF-8 -*-
# 大豆整根 CK + ABA + eBL + IAA + MeJA + SA：原UMAP、候选QC起点、WGCNA程序和根解剖。
# 运行：Rscript soybean_wholeroot_hormone_QC_v8.R
# 本脚本独立运行，不需要source旧版脚本。真实结果仅来自输入数据；截图不用于生成数值。
# 同一时点的处理比较不能证明发育速度、细胞谱系或因果机制。
# QC为组织者；QC邻近initials/meristem至分化细胞为待验证的状态顺序。
# 空间径向层次 != 拟时序；CK=control，不能自动当作细胞分裂素或早期状态。
# 依赖：Seurat, SeuratObject, Matrix, monocle3, SingleCellExperiment,
# SummarizedExperiment, igraph。GO可选AnnotationDbi, org.G.max.eg.db, GO.db。
# 图片可选png/jpeg。无自动安装、无强制连接partition、无重新UMAP。
# 修改：接入已运行的 WGCNA 模块文件 Nature_network_modules.csv（gene_name/module）。

INPUT_RDS <- '/data01/liuxy/Project/CK2_added_v3_Results/seurat_obj_with_CK2_v3.rds'
RUN_DIR <- '/data01/liuxy/Project/CK2_added_v3_Results/Endodermis_monocle3_originalUMAP_grouped_v3/20261119_163425'
# 保留原输出根目录；整根与原Endodermis对象严格分开保存。
OUT_DIR <- file.path(RUN_DIR, 'WholeRoot_hormone_QC_v8')
ASSAY <- 'RNA'
UMAP_NAME <- 'umap'
TREATMENT_COL <- 'treatment'
CELLTYPE_COL <- 'cell_type'
SAMPLE_COL <- 'orig.ident' # 必须核实为生物学重复，不是lane/处理名。
GROUPS <- c('CK','ABA','eBL','IAA','MeJA','SA')
EXPECTED_REPLICATES <- 2L
SEED <- 20261121L
CORES <- 4L
REBUILD_GRAPH <- FALSE    # FALSE优先复用本v8整根缓存；输入参数/文件改变时拒绝旧缓存。
QC_LABELS <- c('QC','Quiescent Center','Quiescent centre')
ROOT_CELL_IDS <- character(0) # 优先：经标记验证的早期细胞barcode。
ROOT_NODES <- character(0)    # 次优先：本次整根图节点，不得照抄Endodermis的Y_86。
ROOT_BIOLOGICALLY_VALIDATED <- FALSE
MIN_QC_CELLS_PER_COMPONENT <- 5L
MIN_QC_SAMPLES_PER_COMPONENT <- 2L
# 没QC/经验证root的分量保持未定向，不自动改用CK或任意Meristem节点。
MIN_GENE_CELLS <- 10L
MIN_SAMPLE_CELLS <- 20L
MIN_CELLS_PER_BIN_SAMPLE <- 5L
MIN_PATH_CELLS <- 150L
MAX_PATHS <- 12L # 超过时保留细胞数最多的12条；所有路径均在审计表列出。
N_BINS <- 24L
PNG_DPI <- 180L
MODULE_MIN_MATCHED_GENES <- 5L
MODULE_MIN_MATCH_RATE <- 0.5

# === 已运行的 WGCNA 模块输入 ===
WGCNA_MEMBERSHIP_FILE <- '/data01/liuxy/Project/muscat_results_ck_ref_paper/CK_reference/Ensemble_PI_Heat_Cold_PEG/run_05_v1_4_annotated/Nature_network_modules.csv'
WGCNA_GENE_COL   <- 'gene_name'  # Nature_network_modules.csv 的基因列
WGCNA_MODULE_COL <- 'module'     # 含 Program01..Program12 与 grey；grey 会被过滤
WGCNA_BACKGROUND_FILE <- ''      # 可选；未提供时 GO 背景用 detected whole-root genes
MODULE_GO_FILE <- ''             # 可选已计算GO表：module,ID,Description,p.adjust,Count,GeneRatio。
MARKER_FILE <- ''                # 可选真实验证marker表gene_id,cell_type,label；空时只画6个候选基因。
HISTOLOGY_MANIFEST <- ''         # CSV列file,treatment,orientation,label；可选sample,scale_bar_um。
RUN_MODULE_GO <- TRUE
ORGDB_PACKAGE <- 'org.G.max.eg.db'
GO_FDR <- 0.05
GO_MIN_SIZE <- 10L
GO_MAX_SIZE <- 500L
RUN_GRAPH_TEST <- TRUE
RUN_PSEUDOBULK_DE <- TRUE
GRAPH_Q <- 1e-5
GRAPH_I <- 0.15
MAX_DYNAMIC_GENES <- 600L
TARGET_GENES <- c('Glyma.17G152800','Glyma.06G050000','Glyma.09G139000',
                  'Glyma.06G003800','Glyma.13G094400','Glyma.12G173300')
GROUP_COLORS <- c(CK='#6A6A6A',ABA='#CC6677',eBL='#DDCC77',IAA='#44AA99',
                  MeJA='#AA4499',SA='#4477AA')
PT_COLORS <- c('#30123B','#4145AB','#22A7A2','#93D35E','#F9D057','#E34A33')
Z_COLORS <- colorRampPalette(c('#2166AC','#F7F7F7','#B2182B'))(101)

save_table <- function(x,name) write.csv(x,file.path(OUT_DIR,paste0(name,'.csv')),
                                         row.names=FALSE,na='',fileEncoding='UTF-8')
read_table <- function(path) {
  if(!nzchar(path)||!file.exists(path))stop('Input not found: ',path)
  if(grepl('\\.rds$',path,ignore.case=TRUE)) {
    z<-readRDS(path);if(!is.data.frame(z))stop('RDS must contain a data.frame: ',path);return(z)
  }
  read.table(path,header=TRUE,sep=if(grepl('\\.(tsv|txt)$',path,ignore.case=TRUE))'\t' else ',',
             quote='"',comment.char='',check.names=FALSE,stringsAsFactors=FALSE,fileEncoding='UTF-8')
}
need <- function(x,cols) {if(!all(cols%in%names(x)))stop('Required columns: ',paste(cols,collapse=', '));x}
safe_name <- function(x) gsub('[^A-Za-z0-9_.-]','_',x)
report <- function(...) {cat(format(Sys.time()),' | ',paste(...,collapse=' '),'\n',sep='')}
note <- function(name,text) {writeLines(text,file.path(OUT_DIR,paste0(name,'.txt')),useBytes=TRUE);report(text)}
draw_note <- function(title,lines) {
  plot.new();title(main=title);text(.05,.85,paste(strwrap(paste(lines,collapse=' '),70),collapse='\n'),
                                    adj=c(0,1),cex=.85)
}
save_plot <- function(name,fun,w=12,h=8) {
  pdf(file.path(OUT_DIR,paste0(name,'.pdf')),width=w,height=h,useDingbats=FALSE)
  tryCatch(fun(),finally=dev.off())
  png(file.path(OUT_DIR,paste0(name,'.png')),width=w*PNG_DPI,height=h*PNG_DPI,res=PNG_DPI,
      type=if(capabilities('cairo'))'cairo' else getOption('bitmapType','Xlib'))
  tryCatch(fun(),finally=dev.off())
}
canon_group <- function(x) {
  k<-toupper(trimws(as.character(x)));a<-c(CK='CK',CONTROL='CK',ABA='ABA',EBL='eBL',IAA='IAA',MEJA='MeJA',SA='SA')
  unname(a[k]) # 不用模糊正则把Heat/Cold/PEG分配到激素组。
}
gene_key <- function(x) {
  # 允许无版本或v4 Glyma及纯数字转录本后缀；其他版本不做隐式转换。
  z<-sub('^(gene|transcript):','',trimws(x),ignore.case=TRUE)
  ok<-grepl('^Glyma[._-]?[0-9]{2}G[0-9]{6}([._-][0-9]+)*([._-]v4)?$',z,ignore.case=TRUE)
  z[ok]<-sub('^Glyma[._-]?([0-9]{2})G([0-9]{6}).*$','Glyma.\\1G\\2',z[ok],ignore.case=TRUE)
  z
}
map_genes <- function(query,features) {
  fk<-gene_key(features);qk<-gene_key(query);ii<-split(seq_along(features),fk)
  ans<-lapply(seq_along(query),function(i){
    exact<-which(features==query[i]);hit<-if(length(exact)==1L)exact else ii[[qk[i]]]
    data.frame(gene_id=query[i],feature_id=if(length(hit)==1L)features[hit] else NA_character_,
               status=if(length(hit)==1L)'matched' else if(length(hit))'ambiguous' else 'missing')
  });do.call(rbind,ans)
}
numeric_colors <- function(v,lim,pal=PT_COLORS,na='#CFCFCF') {
  cols<-rep(na,length(v));good<-is.finite(v)
  u<-if(diff(lim)>0)(v[good]-lim[1])/diff(lim) else rep(.5,sum(good))
  cols[good]<-colorRampPalette(pal)(256)[1L+floor(255*pmax(0,pmin(1,u)))];cols
}
colorbar <- function(lim,pal,label) {
  usr<-par('usr');xr<-diff(usr[1:2]);yr<-diff(usr[3:4]);xx<-usr[1]+xr*c(.65,.94)
  bx<-seq(xx[1],xx[2],length.out=257);yy<-usr[3]+yr*c(.035,.065)
  rect(head(bx,-1),yy[1],tail(bx,-1),yy[2],col=colorRampPalette(pal)(256),border=NA)
  text(xx,yy[1]-.015*yr,format(signif(lim,3),trim=TRUE),adj=c(.5,1),cex=.55)
  text(mean(xx),yy[2]+.015*yr,label,cex=.55)
}
draw_umap <- function(d,values,title,edges=NULL,lim=NULL,pal=PT_COLORS,categorical=NULL,ix=seq_len(nrow(d))) {
  par(mar=c(3,3,3,1));plot(d$x,d$y,type='n',asp=1,xlab='Original UMAP 1',ylab='Original UMAP 2',main=title,cex.main=.9)
  points(d$x,d$y,pch=16,cex=.12,col='#E7E7E7')
  if(!is.null(edges)&&nrow(edges))segments(edges$x,edges$y,edges$xend,edges$yend,col='#646464',lwd=.5)
  if(is.null(categorical)){
    if(is.null(lim))lim<-range(values[is.finite(values)])
    if(!all(is.finite(lim))){mtext('No finite values',side=3);return(invisible(NULL))}
    co<-numeric_colors(values,lim,pal);o<-ix[order(values[ix],na.last=FALSE)]
    points(d$x[o],d$y[o],col=co[o],pch=16,cex=.24);colorbar(lim,pal,'Shared scale')
  } else {
    points(d$x[ix],d$y[ix],col=categorical[values[ix]],pch=16,cex=.22)
    for(z in intersect(names(categorical),unique(values[ix]))) {
      j<-ix[values[ix]==z];text(median(d$x[j]),median(d$y[j]),z,cex=.55,font=2)
    }
  }
}
type_palette <- function(types) {
  base<-c(QC='#D73027',Meristem='#8DAA50',Endodermis='#AA4499',Epidermis='#4477AA',
          Cortex='#66C2A5',Pericycle='#CC6677',Phloem='#DDCC77',Xylem='#8856A7',Procambium='#44AA99',
          'Root Cap'='#999999','Dividing Cell'='#E69F00')
  out<-setNames(hcl.colors(length(types),'Dark 3'),types);j<-intersect(types,names(base));out[j]<-base[j];out
}
log_expression <- function(counts,sf,genes) {
  ans<-counts[genes,,drop=FALSE]%*%Matrix::Diagonal(x=1/sf)
  ans@x<-log1p(ans@x);ans
}
mean_by_sample <- function(values,d,keys=c('treatment','cell_type')) {
  # 单个score，每个生物学重复先求均值；不把细胞当独立重复。
  ii<-split(seq_len(nrow(d)),interaction(d[,c(keys,'sample'),drop=FALSE],drop=TRUE,lex.order=TRUE))
  rows<-lapply(ii,function(j){
    good<-j[is.finite(values[j])];if(!length(good))return(NULL)
    cbind(d[j[1],c(keys,'sample'),drop=FALSE],n=length(good),value=mean(values[good]))
  });z<-do.call(rbind,rows);rownames(z)<-NULL;z
}
make_templates <- function() {
  td<-file.path(OUT_DIR,'input_templates');dir.create(td,showWarnings=FALSE)
  ww<-function(x,f)if(!file.exists(file.path(td,f)))write.csv(x,file.path(td,f),row.names=FALSE)
  ww(data.frame(gene_id=character(),module=character(),kME=numeric()),'wgcna_gene_module.csv')
  ww(data.frame(gene_id=character()),'wgcna_network_background.csv')
  ww(data.frame(module=character(),ID=character(),Description=character(),p.adjust=numeric(),
                Count=integer(),GeneRatio=character()),'module_GO.csv')
  ww(data.frame(file=character(),treatment=character(),orientation=character(),label=character(),
                sample=character(),scale_bar_um=numeric()),'histology_manifest.csv')
  ww(data.frame(gene_id=character(),cell_type=character(),label=character()),'validated_markers.csv')
}

load_or_build <- function() {
  pk<-c('Seurat','SeuratObject','monocle3','Matrix','SingleCellExperiment','SummarizedExperiment','igraph')
  missing<-pk[!vapply(pk,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing))stop('Install required packages first: ',paste(missing,collapse=', '))
  if(!file.exists(INPUT_RDS))stop('Seurat RDS not found: ',INPUT_RDS)
  signature<-list(input=normalizePath(INPUT_RDS),size=file.info(INPUT_RDS)$size,mtime=file.info(INPUT_RDS)$mtime,
                  assay=ASSAY,umap=UMAP_NAME,cols=c(TREATMENT_COL,CELLTYPE_COL,SAMPLE_COL),groups=GROUPS,
                  min_gene_cells=MIN_GENE_CELLS,seed=SEED,monocle3=as.character(packageVersion('monocle3')),schema=8L)
  cache<-file.path(OUT_DIR,'v8_wholeroot_unordered_cache.rds')
  if(file.exists(cache)&&!REBUILD_GRAPH){
    z<-readRDS(cache)
    if(!identical(z$signature,signature))stop('Cached graph/input settings differ. Set REBUILD_GRAPH=TRUE after checking inputs.')
    report('Reusing v8 whole-root graph cache.');return(z)
  }
  report('Reading full Seurat; six groups jointly fitted in unchanged UMAP.')
  obj<-readRDS(INPUT_RDS);md<-obj[[]];need(md,c(TREATMENT_COL,CELLTYPE_COL,SAMPLE_COL))
  tr<-canon_group(md[[TREATMENT_COL]])
  save_table(as.data.frame(table(original=md[[TREATMENT_COL]],selected=ifelse(is.na(tr),'EXCLUDED',tr))),'01_treatment_mapping')
  keep<-!is.na(tr);md<-md[keep,,drop=FALSE];tr<-tr[keep]
  if(anyNA(md[[CELLTYPE_COL]])||any(!nzchar(trimws(md[[CELLTYPE_COL]]))))stop('Missing cell_type annotations in selected cells.')
  if(anyNA(md[[SAMPLE_COL]])||any(!nzchar(trimws(md[[SAMPLE_COL]]))))stop('Missing biological sample IDs.')
  md$v8_treatment<-tr;md$v8_type<-as.character(md[[CELLTYPE_COL]])
  md$v8_sample<-paste(tr,as.character(md[[SAMPLE_COL]]),sep='::')
  audit<-unique(data.frame(treatment=tr,sample=md$v8_sample,original_sample=as.character(md[[SAMPLE_COL]])))
  save_table(audit,'02_biological_sample_mapping')
  nr<-table(factor(audit$treatment,levels=GROUPS))
  if(any(nr!=EXPECTED_REPLICATES))stop('Expected 2 biological replicates per group. Check 02_biological_sample_mapping.csv and SAMPLE_COL.')
  if(!UMAP_NAME%in%names(obj@reductions))stop('Original UMAP reduction absent: ',UMAP_NAME)
  xy<-SeuratObject::Embeddings(obj,UMAP_NAME)[rownames(md),1:2,drop=FALSE]
  if(any(!is.finite(xy)))stop('Invalid original coordinates.')
  obj<-subset(obj,cells=rownames(md))
  if(inherits(obj[[ASSAY]],'Assay5')){
    ll<-grep('^counts($|[.])',SeuratObject::Layers(obj[[ASSAY]]),value=TRUE)
    if(!length(ll))stop('No raw counts layers.')
    if(length(ll)!=1L||ll!='counts')obj[[ASSAY]]<-SeuratObject::JoinLayers(obj[[ASSAY]],layers=ll,new='counts')
    ct<-SeuratObject::LayerData(obj,assay=ASSAY,layer='counts')
  }else ct<-SeuratObject::GetAssayData(obj,assay=ASSAY,slot='counts')
  ct<-methods::as(ct[,rownames(md),drop=FALSE],'dgCMatrix');rm(obj);gc()
  if(any(ct@x<0)|any(!is.finite(ct@x))|any(abs(ct@x-round(ct@x))>1e-6))stop('Expected raw nonnegative integer RNA counts.')
  if(anyDuplicated(rownames(ct)))stop('Duplicated gene IDs.')
  # 保留所有gene以供候选基因绘图；检测阈值仅用于graph_test/GO背景。
  keepcell<-Matrix::colSums(ct)>0;save_table(data.frame(cell=colnames(ct),retained=keepcell),'03_cell_filter')
  ct<-ct[,keepcell,drop=FALSE];md<-md[colnames(ct),,drop=FALSE];xy<-xy[colnames(ct),,drop=FALSE]
  cds<-monocle3::new_cell_data_set(ct,cell_metadata=md,
                                   gene_metadata=data.frame(gene_short_name=rownames(ct),row.names=rownames(ct)))
  cds<-monocle3::estimate_size_factors(cds)
  SingleCellExperiment::reducedDim(cds,'UMAP')<-xy
  cds<-monocle3::cluster_cells(cds,reduction_method='UMAP',random_seed=SEED)
  cds<-monocle3::learn_graph(cds,use_partition=TRUE,close_loop=FALSE)
  delta<-max(abs(SingleCellExperiment::reducedDim(cds,'UMAP')-xy))
  if(delta!=0)stop('Original UMAP changed unexpectedly.')
  save_table(data.frame(n_cells=ncol(ct),max_coordinate_change=delta),'04_UMAP_preservation')
  z<-list(cds=cds,signature=signature);saveRDS(z,cache);z
}

root_graph <- function(cds) {
  g<-monocle3::principal_graph(cds)[['UMAP']];vn<-igraph::V(g)$name
  aux<-cds@principal_graph_aux[['UMAP']];xy<-t(aux$dp_mst)
  if(!all(vn%in%rownames(xy)))stop('Principal node coordinates do not match graph.')
  comp<-igraph::components(g)$membership
  nodes<-data.frame(node=vn,x=xy[vn,1],y=xy[vn,2],component=unname(comp[vn]),degree=as.integer(igraph::degree(g)[vn]))
  raw<-as.character(aux$pr_graph_cell_proj_closest_vertex[colnames(cds),1])
  if(all(raw%in%vn))near<-raw else {
    k<-suppressWarnings(as.integer(raw));if(anyNA(k)||any(k<1|k>length(vn)))stop('Invalid nearest nodes.');near<-vn[k]
  }
  um<-SingleCellExperiment::reducedDim(cds,'UMAP');md<-as.data.frame(SummarizedExperiment::colData(cds))
  d<-data.frame(cell=colnames(cds),x=um[,1],y=um[,2],treatment=md$v8_treatment,
                sample=md$v8_sample,cell_type=md$v8_type,node=near,component=unname(comp[near]))
  ed<-igraph::as_data_frame(g,what='edges')[,c('from','to')]
  ed$x<-nodes$x[match(ed$from,vn)];ed$y<-nodes$y[match(ed$from,vn)]
  ed$xend<-nodes$x[match(ed$to,vn)];ed$yend<-nodes$y[match(ed$to,vn)]
  save_table(nodes,'05_principal_nodes');save_table(ed,'06_graph_edges')
  roots<-character();sources<-character()
  if(length(ROOT_CELL_IDS)) {
    if(!all(ROOT_CELL_IDS%in%d$cell))stop('ROOT_CELL_IDS absent from selected data.')
    cand<-d[d$cell%in%ROOT_CELL_IDS,,drop=FALSE];src<-'user_root_cells'
  } else {cand<-d[tolower(d$cell_type)%in%tolower(QC_LABELS),,drop=FALSE];src<-'QC_annotation_candidate'}
  audit<-lapply(sort(unique(d$component)),function(cc){
    u<-cand[cand$component==cc,,drop=FALSE];root<-NA_character_;reason<-'no_supported_root'
    manual<-ROOT_NODES[ROOT_NODES%in%nodes$node[nodes$component==cc]]
    if(length(manual)>1)stop('Use one common root per component, not multiple competing roots.')
    if(length(manual)==1){root<-manual;reason<-'user_root_node'} else if(nrow(u)&&
                                                                         (length(ROOT_CELL_IDS)|| (nrow(u)>=MIN_QC_CELLS_PER_COMPONENT&&length(unique(u$sample))>=MIN_QC_SAMPLES_PER_COMPONENT))){
      # 每个样本内节点QC比例，再跨样本平均：避免单个大样本主导选根。
      tt<-table(factor(u$node,levels=nodes$node[nodes$component==cc]),u$sample)
      score<-rowMeans(sweep(tt,2,colSums(tt),'/'));root<-names(which.max(score));reason<-src
    }
    if(!is.na(root)){roots<<-c(roots,root);sources<<-c(sources,reason)}
    data.frame(component=cc,n_cells=sum(d$component==cc),candidate_cells=nrow(u),
               candidate_samples=length(unique(u$sample)),root_node=root,source=reason,biologically_validated=ROOT_BIOLOGICALLY_VALIDATED)
  });audit<-do.call(rbind,audit)
  if(length(setdiff(ROOT_NODES,vn)))stop('ROOT_NODES not in this graph: ',paste(setdiff(ROOT_NODES,vn),collapse=', '))
  d$pseudotime<-NA_real_
  if(length(roots)) {
    cds<-monocle3::order_cells(cds,root_pr_nodes=roots)
    pt<-monocle3::pseudotime(cds);d$pseudotime<-as.numeric(pt[d$cell]);d$pseudotime[!is.finite(d$pseudotime)]<-NA_real_
  }
  save_table(audit,'07_root_audit');save_table(d,'08_cells_pseudotime');
  save_table(as.data.frame(table(d$treatment,d$cell_type,d$component,useNA='ifany')),'09_component_type_group_counts')
  saveRDS(cds,file.path(OUT_DIR,'v8_wholeroot_ordered.rds'))
  save_plot('A0_root_diagnostics',function(){
    draw_umap(d,as.character(d$component),'Whole-root graph components',ed,
              categorical=setNames(hcl.colors(length(unique(d$component)),'Dark 3'),sort(unique(d$component))))
    text(nodes$x,nodes$y,nodes$node,cex=.35,pos=3)
    z<-nodes[nodes$node%in%roots,];points(z$x,z$y,pch=23,bg='#F0E442',cex=1.3)
  },14,11)
  list(cds=cds,d=d,nodes=nodes,edges=ed,graph=g,roots=roots,root_audit=audit)
}

# === 读取并接入已运行的 WGCNA 模块文件 ===
load_modules <- function(features) {
  if(!nzchar(WGCNA_MEMBERSHIP_FILE)) {
    note('WGCNA_status','WGCNA membership not supplied: no module scores, module GO or module response claims are calculated. Fill WGCNA_MEMBERSHIP_FILE; screenshot labels alone are not gene sets.')
    return(list(sets=list(),membership=NULL))
  }
  if(!file.exists(WGCNA_MEMBERSHIP_FILE))stop('WGCNA membership file not found: ',WGCNA_MEMBERSHIP_FILE)
  raw<-read_table(WGCNA_MEMBERSHIP_FILE)
  
  gene_col <- WGCNA_GENE_COL
  module_col <- WGCNA_MODULE_COL
  if(!nzchar(gene_col) || !gene_col %in% names(raw)) {
    cand <- intersect(c('gene_id','gene_name','Gene','gene','GID','feature_id'), names(raw))
    if(!length(cand)) stop('WGCNA membership file lacks a gene column. Columns: ', paste(names(raw), collapse=', '))
    gene_col <- cand[1]
    report('WGCNA gene column auto-detected: ', gene_col)
  }
  if(!nzchar(module_col) || !module_col %in% names(raw)) {
    cand <- intersect(c('module','Program','color','Module','original_module'), names(raw))
    if(!length(cand)) stop('WGCNA membership file lacks a module column. Columns: ', paste(names(raw), collapse=', '))
    module_col <- cand[1]
    report('WGCNA module column auto-detected: ', module_col)
  }
  
  m <- data.frame(
    gene_id = trimws(as.character(raw[[gene_col]])),
    module  = trimws(as.character(raw[[module_col]])),
    stringsAsFactors = FALSE
  )
  m <- m[!is.na(m$gene_id) & nzchar(m$gene_id) &
           !is.na(m$module) & nzchar(m$module), , drop = FALSE]
  if(!nrow(m))stop('WGCNA membership file contains no usable gene/module rows.')
  m <- unique(m[,c('gene_id','module')])
  
  save_table(data.frame(
    File = WGCNA_MEMBERSHIP_FILE,
    MD5  = unname(tools::md5sum(WGCNA_MEMBERSHIP_FILE)),
    Gene_column = gene_col,
    Module_column = module_col,
    N_rows = nrow(m),
    N_modules = length(unique(m$module))
  ), '20_WGCNA_membership_source')
  
  mp<-map_genes(m$gene_id,features);m$feature_id<-mp$feature_id;m$status<-mp$status
  save_table(m,'20_WGCNA_mapping')
  conflicting<-split(m$module[!is.na(m$feature_id)],m$feature_id[!is.na(m$feature_id)])
  if(any(lengths(lapply(conflicting,unique))>1L))stop('One gene assigned to multiple WGCNA modules. Resolve memberships first.')
  stats<-do.call(rbind,lapply(split(m,m$module),function(z)data.frame(module=z$module[1],
                                                                      n_input=length(unique(z$gene_id)),n_matched=length(unique(na.omit(z$feature_id))),
                                                                      match_rate=mean(!is.na(z$feature_id)))))
  stats$eligible<-stats$n_matched>=MODULE_MIN_MATCHED_GENES&stats$match_rate>=MODULE_MIN_MATCH_RATE&
    !tolower(stats$module)%in%c('grey','gray','0','megrey','megray')
  save_table(stats,'21_WGCNA_mapping_summary')
  sets<-lapply(split(m$feature_id[m$module%in%stats$module[stats$eligible]],
                     m$module[m$module%in%stats$module[stats$eligible]]),function(x)unique(na.omit(x)))
  report('WGCNA modules loaded: ',length(sets),' eligible; ',nrow(m),' membership rows.')
  list(sets=sets,membership=m)
}
module_scores <- function(counts,sf,sets) {
  # 固定全体六处理参考，每个基因log1p(count/sizefactor)后做Z标准化，再取模块平均。
  # 这是描述性module expression score，不是WGCNA eigengene、AUCell、调控活性或百分比激活。
  ans<-matrix(NA_real_,ncol(counts),length(sets),dimnames=list(colnames(counts),names(sets)))
  audit<-list()
  for(m in names(sets)) {
    genes<-sets[[m]];s<-numeric(ncol(counts));n<-0L
    for(k in split(genes,ceiling(seq_along(genes)/40))) {
      e<-as.matrix(log_expression(counts,sf,k));mu<-rowMeans(e);ss<-apply(e,1,sd);ok<-is.finite(ss)&ss>1e-8
      if(any(ok)){e<-sweep(sweep(e[ok,,drop=FALSE],1,mu[ok],'-'),1,ss[ok],'/');s<-s+colSums(e);n<-n+sum(ok)}
    }
    if(n>=MODULE_MIN_MATCHED_GENES)ans[,m]<-s/n
    audit[[m]]<-data.frame(module=m,n_matched=length(genes),n_variable=n,score_available=n>=MODULE_MIN_MATCHED_GENES)
  }
  if(length(audit))save_table(do.call(rbind,audit),'22_module_score_audit')
  ans
}
module_GO <- function(sets,features,detected) {
  if(!length(sets))return(NULL)
  if(nzchar(MODULE_GO_FILE)) {
    z<-need(read_table(MODULE_GO_FILE),c('module','ID','Description','p.adjust','Count','GeneRatio'))
    z$FDR<-as.numeric(as.character(z$p.adjust));z$Count<-as.numeric(as.character(z$Count))
    if(any(!is.finite(z$FDR)|z$FDR<0|z$FDR>1))stop('Imported GO p.adjust must be numeric in [0,1].')
    z<-z[z$module%in%names(sets),,drop=FALSE];z$source<-'imported; correction scope defined by supplied table'
    save_table(z,'26_GO_imported');return(z)
  }
  if(!RUN_MODULE_GO)return(NULL)
  pp<-c('AnnotationDbi',ORGDB_PACKAGE,'GO.db')
  if(any(!vapply(pp,requireNamespace,logical(1),quietly=TRUE))){note('GO_status','GO annotation packages missing; module scores remain available.');return(NULL)}
  db<-get(ORGDB_PACKAGE,envir=asNamespace(ORGDB_PACKAGE))
  if(!'GID'%in%AnnotationDbi::keytypes(db))stop('Expected soybean OrgDb GID keytype.')
  cols<-AnnotationDbi::columns(db);gc<-if('GOALL'%in%cols)'GOALL' else 'GO';oc<-if(gc=='GOALL')'ONTOLOGYALL' else 'ONTOLOGY'
  needcols<-c(gc,oc);if(!all(needcols%in%cols))stop('OrgDb lacks GO/ontology columns.')
  universe<-features[detected]
  if(nzchar(WGCNA_BACKGROUND_FILE)) {
    bg<-need(read_table(WGCNA_BACKGROUND_FILE),'gene_id');bgmap<-map_genes(bg$gene_id,features)
    save_table(bgmap,'23_WGCNA_background_mapping');universe<-intersect(universe,na.omit(bgmap$feature_id))
    source<-'detected genes intersect original WGCNA network background'
  }else source<-'detected whole-root genes; original WGCNA network background not supplied'
  keys<-AnnotationDbi::keys(db,keytype='GID');mapping<-map_genes(universe,keys)
  save_table(mapping,'24_GO_background_mapping');note('GO_background',source)
  if(!nrow(mapping)||mean(!is.na(mapping$feature_id))<.5){note('GO_status','Less than 50% of background genes map to OrgDb. Inspect version/mapping; no enrichment claimed.');return(NULL)}
  mapped<-unique(na.omit(mapping$feature_id))
  ann<-suppressMessages(AnnotationDbi::select(db,keys=mapped,keytype='GID',columns=needcols))
  ann<-unique(ann[!is.na(ann[[gc]])&!is.na(ann[[oc]])&ann[[oc]]=='BP',c('GID',gc)])
  bg<-unique(ann$GID);terms<-split(ann$GID,ann[[gc]]);terms<-lapply(terms,unique)
  terms<-terms[lengths(terms)>=GO_MIN_SIZE&lengths(terms)<=GO_MAX_SIZE]
  if(length(bg)<20||!length(terms)){note('GO_status','Insufficient BP-annotated background or no eligible GO term sizes.');return(NULL)}
  desc<-suppressMessages(AnnotationDbi::select(GO.db::GO.db,keys=names(terms),keytype='GOID',columns='TERM'))
  nameby<-setNames(desc$TERM,desc$GOID);out<-list();aud<-list()
  for(m in names(sets)) {
    # universe输入先限制，再映射GID；独立保留映射去重。
    ss<-mapping$feature_id[mapping$gene_id%in%sets[[m]]];ss<-intersect(na.omit(ss),bg);n<-length(ss);N<-length(bg)
    aud[[m]]<-data.frame(module=m,input_genes=length(sets[[m]]),BP_mapped_genes=n,BP_background=N,tested_terms=length(terms))
    if(n<2)next
    K<-lengths(terms);k<-vapply(terms,function(t)length(intersect(t,ss)),integer(1))
    out[[m]]<-data.frame(module=m,ID=names(terms),Description=unname(nameby[names(terms)]),
                         Count=k,GeneRatio=paste0(k,'/',n),BgRatio=paste0(K,'/',N),
                         pvalue=phyper(k-1,K,N-K,n,lower.tail=FALSE),geneID=vapply(terms,function(t)paste(intersect(t,ss),collapse='/'),character(1)))
  }
  save_table(do.call(rbind,aud),'25_GO_module_audit')
  if(!length(out)){note('GO_status','No modules with >=2 mapped BP genes.');return(NULL)}
  z<-do.call(rbind,out);z$FDR<-p.adjust(z$pvalue,'BH') # 所有module x eligibleBP terms，包括零overlap。
  z$source<-paste('Global BH over all tested module-term pairs;',source)
  save_table(z,'26_GO_all_tests');save_table(z[z$FDR<GO_FDR&z$Count>=2,],'27_GO_significant')
  note('GO_status',paste('Tested pairs:',nrow(z),'; significant with Count>=2:',sum(z$FDR<GO_FDR&z$Count>=2)))
  z
}
graph_candidates <- function(cds,signature) {
  if(!RUN_GRAPH_TEST)return(character())
  file<-file.path(OUT_DIR,'graph_test_v8_cache.rds')
  detected<-Matrix::rowSums(SummarizedExperiment::assay(cds,'counts')>0)>=MIN_GENE_CELLS
  sig<-list(graph=signature,genes=rownames(cds)[detected],test='principal_graph',cores=CORES)
  if(file.exists(file)&&!REBUILD_GRAPH) {
    z<-readRDS(file);if(!identical(z$signature,sig))stop('graph_test cache differs from current graph/genes. Set REBUILD_GRAPH=TRUE.')
    gt<-z$result
  } else {
    report('Running graph_test on detected genes; graph autocorrelation, not treatment DEGs.')
    gt<-as.data.frame(monocle3::graph_test(cds[detected,],neighbor_graph='principal_graph',cores=CORES))
    gt$gene_id<-rownames(gt);saveRDS(list(signature=sig,result=gt),file)
  }
  save_table(gt,'30_graph_test_all')
  ok<-is.finite(gt$q_value)&gt$q_value<GRAPH_Q&is.finite(gt$morans_I)&gt$morans_I>GRAPH_I
  if('status'%in%names(gt))ok<-ok&!is.na(gt$status)&tolower(gt$status)=='ok'
  selected<-gt[ok,,drop=FALSE];selected<-selected[order(-selected$morans_I,selected$q_value),,drop=FALSE]
  save_table(selected,'31_graph_test_selected_all')
  head(selected$gene_id,MAX_DYNAMIC_GENES)
}
get_paths <- function(r) {
  paths<-list();aa<-list();k<-0L
  for(rt in r$roots){
    cc<-r$nodes$component[match(rt,r$nodes$node)]
    leaves<-r$nodes$node[r$nodes$component==cc&r$nodes$degree==1&r$nodes$node!=rt]
    for(leaf in leaves){
      v<-as.character(igraph::as_ids(igraph::shortest_paths(r$graph,from=rt,to=leaf,weights=NA)$vpath[[1]]))
      ix<-which(r$d$component==cc&r$d$node%in%v&is.finite(r$d$pseudotime))
      k<-k+1L;id<-sprintf('Path_%03d',k);paths[[id]]<-list(id=id,root=rt,leaf=leaf,component=cc,nodes=v,ix=ix)
      aa[[id]]<-data.frame(path=id,root=rt,leaf=leaf,component=cc,n_cells=length(ix))
    }
  }
  if(!length(aa)) {note('paths_status','No supported root-to-leaf paths. UMAP/module/anatomy panels still generated.');return(list())}
  tab<-do.call(rbind,aa);eligible<-tab$path[tab$n_cells>=MIN_PATH_CELLS]
  eligible<-head(eligible[order(-tab$n_cells[match(eligible,tab$path)])],MAX_PATHS)
  tab$plotted<-tab$path%in%eligible;save_table(tab,'32_all_root_to_leaf_paths')
  save_table(do.call(rbind,lapply(paths,function(p)data.frame(path=p$id,node=p$nodes))),'33_path_node_membership')
  paths[eligible] # shared trunk cells legitimately occur in multiple paths; never sum path counts.
}
bounded_density <- function(x,lim,bw,grid) {
  if(length(x)<2||diff(range(x))==0)return(rep(NA_real_,length(grid)))
  kd<-density(c(x,2*lim[1]-x,2*lim[2]-x),bw=bw,n=1024,from=lim[1],to=lim[2])
  y<-approx(kd$x,kd$y,xout=grid,rule=2)$y
  a<-sum(diff(grid)*(head(y,-1)+tail(y,-1))/2);if(a>0)y/a else y
}
ridge_data <- function(d) {
  d<-d[is.finite(d$pseudotime),,drop=FALSE];lim<-range(d$pseudotime)
  if(nrow(d)<10||!all(is.finite(lim))||diff(lim)<=0)return(NULL)
  grid<-seq(lim[1],lim[2],length.out=256);bw<-max(bw.nrd0(d$pseudotime),diff(lim)/40)
  rows<-lapply(split(d,d$sample),function(z){
    y<-bounded_density(z$pseudotime,lim,bw,grid)
    data.frame(treatment=z$treatment[1],sample=z$sample[1],x=grid,density=y,n=nrow(z))
  });do.call(rbind,rows)
}
draw_ridges <- function(d,rd,title='Six treatments | equal-weight replicate density') {
  if(is.null(rd)){draw_note(title,'Insufficient finite pseudotime.');return()}
  par(mar=c(3,6,3,1));lim<-range(rd$x);plot(NA,xlim=lim,ylim=c(.4,7.1),xlab='Pseudotime (same rooted component/path)',ylab='',yaxt='n',main=title,cex.main=.82)
  ys<-rev(seq_along(GROUPS));axis(2,at=ys,labels=GROUPS,las=1,tick=FALSE,cex.axis=.8)
  pk<-max(rd$density,na.rm=TRUE);if(!is.finite(pk)||pk<=0)return()
  for(i in seq_along(GROUPS)){
    tr<-GROUPS[i];z<-rd[rd$treatment==tr&is.finite(rd$density),,drop=FALSE]
    if(!nrow(z)){text(mean(lim),ys[i]+.2,'Insufficient cells',cex=.6);next}
    mm<-aggregate(density~x,z,mean);h<-mm$density/pk*.9
    polygon(c(mm$x,rev(mm$x)),c(rep(ys[i],nrow(mm)),rev(ys[i]+h)),col=adjustcolor(GROUP_COLORS[tr],.5),border=NA)
    for(s in unique(z$sample)){zz<-z[z$sample==s,];lines(zz$x,ys[i]+zz$density/pk*.9,col=GROUP_COLORS[tr],lty=3,lwd=.6)}
    lines(mm$x,ys[i]+h,col=GROUP_COLORS[tr],lwd=1.1)
    text(lim[2],ys[i]+.8,paste0('n=',sum(d$treatment==tr),'; r=',length(unique(z$sample))),adj=1,cex=.55)
  }
}
bin_replicates <- function(values,d,breaks) {
  # values: features x cells，分箱所有处理同边界；无足够细胞的重复箱为NA，不补零。
  stopifnot(ncol(values)==nrow(d));nb<-length(breaks)-1L
  bin<-cut(d$pseudotime,breaks=breaks,include.lowest=TRUE,labels=FALSE)
  samples<-unique(d$sample);a<-array(NA_real_,c(nrow(values),nb,length(samples)),dimnames=list(rownames(values),NULL,samples))
  ct<-matrix(0L,nb,length(samples),dimnames=list(NULL,samples))
  for(j in seq_along(samples))for(b in seq_len(nb)){
    ix<-which(d$sample==samples[j]&bin==b);ct[b,j]<-length(ix)
    if(length(ix)>=MIN_CELLS_PER_BIN_SAMPLE)a[,b,j]<-if(inherits(values,'Matrix'))
      Matrix::rowMeans(values[,ix,drop=FALSE]) else rowMeans(values[,ix,drop=FALSE])
  }
  list(value=a,counts=ct,samples=samples,group=d$treatment[match(samples,d$sample)],centers=(head(breaks,-1)+tail(breaks,-1))/2)
}
group_bins <- function(b) {
  a<-array(NA_real_,c(dim(b$value)[1],dim(b$value)[2],length(GROUPS)),dimnames=list(dimnames(b$value)[[1]],NULL,GROUPS))
  for(i in seq_along(GROUPS)) {
    j<-which(b$group==GROUPS[i]);if(!length(j))next
    # Require BOTH biological replicates; missing coverage must remain visible.
    if(length(j)==EXPECTED_REPLICATES)for(k in seq_len(dim(a)[2])){
      mm<-matrix(b$value[,k,j],nrow=dim(a)[1]);ok<-rowSums(is.finite(mm))==length(j)
      a[ok,k,i]<-rowMeans(mm[ok,,drop=FALSE])
    }
  };a
}
draw_matrix <- function(m,title,xlabels=NULL,ylabs=rownames(m),lim=c(-2,2),pal=Z_COLORS) {
  if(!nrow(m)||!ncol(m)){draw_note(title,'No eligible data.');return()}
  par(mar=c(4,6,3,1));plot(NA,xlim=c(.5,ncol(m)+.5),ylim=c(.5,nrow(m)+.5),
                           xlab='',ylab='',xaxt='n',yaxt='n',main=title,cex.main=.82,xaxs='i',yaxs='i')
  for(i in seq_len(nrow(m)))rect(seq_len(ncol(m))-.5,nrow(m)-i+.5,seq_len(ncol(m))+.5,nrow(m)-i+1.5,
                                 col=numeric_colors(m[i,],lim,pal,na='#D6D6D6'),border=NA)
  if(!is.null(xlabels))axis(1,at=seq_len(ncol(m)),labels=xlabels,las=2,cex.axis=.55,tick=FALSE)
  if(length(ylabs)>0&&length(ylabs)<=35)axis(2,at=rev(seq_len(nrow(m))),labels=ylabs,las=1,tick=FALSE,cex.axis=.65)
  box();mtext(paste0('Scale ',lim[1],' to ',lim[2],'; grey = insufficient coverage'),side=1,line=2.8,cex=.55)
}

pseudobulk_DE <- function(counts,d) {
  if(!RUN_PSEUDOBULK_DE)return(invisible(NULL))
  if(!requireNamespace('edgeR',quietly=TRUE)||!requireNamespace('statmod',quietly=TRUE)){
    note('DE_status','edgeR/statmod unavailable. Descriptive plots still run; no treatment DEG claims. Install via Bioconductor and rerun for replicate-level DE.');return(invisible(NULL))
  }
  out<-list();aud<-list()
  for(ty in unique(d$cell_type)){
    ix<-which(d$cell_type==ty);meta<-d[ix,,drop=FALSE];samples<-unique(meta$sample)
    nc<-table(factor(meta$sample,levels=samples));samples<-samples[nc>=MIN_SAMPLE_CELLS]
    grp<-d$treatment[match(samples,d$sample)]
    eligible<-GROUPS[vapply(GROUPS,function(g)sum(grp==g)==EXPECTED_REPLICATES,logical(1))]
    if(!'CK'%in%eligible||length(eligible)<2){aud[[ty]]<-data.frame(cell_type=ty,status='insufficient samples',n_samples=length(samples));next}
    samples<-samples[grp%in%eligible];grp<-d$treatment[match(samples,d$sample)];ix<-which(d$cell_type==ty&d$sample%in%samples)
    sm<-Matrix::sparseMatrix(i=seq_along(ix),j=match(d$sample[ix],samples),x=1,dims=c(length(ix),length(samples)))
    pb<-as.matrix(counts[,ix,drop=FALSE]%*%sm);colnames(pb)<-samples
    group<-factor(grp,levels=eligible);design<-model.matrix(~0+group);colnames(design)<-eligible
    y<-edgeR::DGEList(pb);keep<-edgeR::filterByExpr(y,design=design)
    if(sum(keep)<20){aud[[ty]]<-data.frame(cell_type=ty,status='too few expressed genes',n_samples=length(samples));next}
    y<-edgeR::calcNormFactors(y[keep,,keep.lib.sizes=FALSE]);y<-edgeR::estimateDisp(y,design,robust=TRUE)
    fit<-edgeR::glmQLFit(y,design,robust=TRUE)
    for(gr in setdiff(eligible,'CK')){
      con<-setNames(numeric(length(eligible)),eligible);con[gr]<-1;con['CK']<--1
      tt<-edgeR::topTags(edgeR::glmQLFTest(fit,contrast=con),n=Inf)$table
      tt$gene_id<-rownames(tt);tt$cell_type<-ty;tt$comparison<-paste0(gr,'_vs_CK');out[[paste(ty,gr)]]<-tt
    }
    aud[[ty]]<-data.frame(cell_type=ty,status='tested',n_samples=length(samples))
    saveRDS(list(counts=pb,design=design,samples=data.frame(sample=samples,treatment=grp)),
            file.path(OUT_DIR,paste0('pseudobulk_',safe_name(ty),'.rds')))
  }
  if(length(aud))save_table(do.call(rbind,aud),'40_pseudobulk_audit')
  if(length(out)){
    z<-do.call(rbind,out);z$FDR_global_celltype_gene_contrast<-p.adjust(z$PValue,'BH')
    save_table(z,'41_pseudobulk_DE_all');save_table(z[z$FDR_global_celltype_gene_contrast<.05,],'42_pseudobulk_DE_significant')
  }
}
prepare_markers <- function(counts,sf,d) {
  inp<-if(nzchar(MARKER_FILE))need(read_table(MARKER_FILE),c('gene_id','cell_type','label')) else
    data.frame(gene_id=TARGET_GENES,cell_type='candidate (not validated cell-type marker)',label=TARGET_GENES)
  mp<-map_genes(inp$gene_id,rownames(counts));mp$intended_cell_type<-inp$cell_type;mp$label<-inp$label
  save_table(mp,'43_marker_mapping');rows<-list()
  for(i in which(!is.na(mp$feature_id))){
    ex<-as.numeric(log_expression(counts,sf,mp$feature_id[i]));aa<-mean_by_sample(ex,d,'cell_type')
    pp<-mean_by_sample(as.numeric(counts[mp$feature_id[i],]>0),d,'cell_type')
    aa$pct<-pp$value*100;aa<-aa[aa$n>=MIN_SAMPLE_CELLS,,drop=FALSE]
    if(!nrow(aa))next
    z<-aggregate(cbind(value,pct)~cell_type,aa,mean);z$gene_id<-mp$gene_id[i];z$label<-mp$label[i]
    z$scaled<-if(sd(z$value)>1e-8)as.numeric(scale(z$value)) else 0
    rows[[as.character(i)]]<-z
  }
  z<-if(length(rows))do.call(rbind,rows) else NULL;if(!is.null(z))save_table(z,'44_marker_dotplot_data')
  list(mapping=mp,data=z)
}
draw_dot <- function(z) {
  if(is.null(z)||!nrow(z)){draw_note('B | Marker/candidate expression','No matched markers.');return()}
  genes<-unique(z$label);types<-unique(z$cell_type);par(mar=c(8,7,3,1))
  plot(NA,xlim=c(.5,length(types)+.5),ylim=c(.5,length(genes)+.5),xlab='',ylab='',xaxt='n',yaxt='n',
       main='B | Expression by cell type',cex.main=.85)
  points(match(z$cell_type,types),match(z$label,genes),pch=21,bg=numeric_colors(z$scaled,c(-2,2),Z_COLORS),
         col='#606060',lwd=.3,cex=2.8*sqrt(z$pct/100))
  axis(1,seq_along(types),types,las=2,cex.axis=.6,tick=FALSE);axis(2,seq_along(genes),genes,las=1,cex.axis=.6,tick=FALSE)
  mtext('Size: % detected; colour: gene-wise Z of replicate means',side=3,line=.2,cex=.52)
}
similarity_matrix <- function(counts,sf,d,detected) {
  # 样本x类型先求log表达均值，再每类型等权样本均值；不是论文的单细胞PCA相关矩阵。
  typ<-sort(unique(d$cell_type));features<-rownames(counts)[detected]
  parts<-split(seq_len(nrow(d)),interaction(d$cell_type,d$sample,drop=TRUE));parts<-parts[lengths(parts)>=MIN_SAMPLE_CELLS]
  if(!length(parts)){note('similarity_status','No cell type/sample has enough cells for similarity.');return(matrix(numeric(),0,0))}
  ii<-unlist(parts,use.names=FALSE);jj<-rep(seq_along(parts),lengths(parts))
  a<-Matrix::sparseMatrix(i=ii,j=jj,x=1/rep(lengths(parts),lengths(parts)),dims=c(nrow(d),length(parts)))
  pt<-vapply(parts,function(j)d$cell_type[j[1]],character(1));pb<-matrix(0,length(features),length(parts),dimnames=list(features,NULL))
  typ<-sort(unique(pt))
  for(k in split(seq_along(features),ceiling(seq_along(features)/100)))pb[k,]<-as.matrix(log_expression(counts,sf,features[k])%*%a)
  av<-sapply(typ,function(t)if(any(pt==t))rowMeans(pb[,pt==t,drop=FALSE]) else rep(NA_real_,nrow(pb)))
  ok<-rowSums(is.finite(av))==length(typ)&apply(av,1,sd)>1e-8
  ans<-cor(av[ok,,drop=FALSE],method='pearson');save_table(data.frame(cell_type=rownames(ans),ans,check.names=FALSE),'45_celltype_expression_correlation')
  ans
}
draw_GO <- function(go) {
  if(is.null(go)||!nrow(go)){draw_note('E | WGCNA GO','No usable enrichment table; see GO_status.txt or provide membership/GO data.');return()}
  z<-go[is.finite(go$FDR)&go$FDR<GO_FDR&go$Count>=2,,drop=FALSE]
  z<-do.call(rbind,lapply(split(z,z$module),function(a)head(a[order(a$FDR,-a$Count),],2)))
  if(is.null(z)||!nrow(z)){draw_note('E | WGCNA GO','No terms pass the prespecified FDR and overlap thresholds. All tested terms are saved; thresholds were not relaxed.');return()}
  z<-head(z,22);labels<-paste(z$module,z$Description,sep=' | ')
  labels<-vapply(labels,function(a)paste(strwrap(a,52),collapse='\n'),character(1))
  par(mar=c(4,17,3,1));x<- -log10(pmax(z$FDR,.Machine$double.xmin));yy<-rev(seq_len(nrow(z)))
  plot(x,yy,pch=21,bg=numeric_colors(x,range(x),c('#9ECAE1','#8B1A4A')),cex=1+sqrt(z$Count)/3,
       yaxt='n',ylab='',xlab='-log10(FDR)',main='E | Module GO BP',cex.main=.85)
  axis(2,yy,labels,las=1,tick=FALSE,cex.axis=.53)
}
module_response <- function(scores,d) {
  if(!ncol(scores))return(NULL)
  rr<-list();delta<-list()
  for(m in colnames(scores)){
    x<-mean_by_sample(scores[,m],d);if(is.null(x))next
    x<-x[x$n>=MIN_SAMPLE_CELLS,,drop=FALSE];x$module<-m;rr[[m]]<-x
    for(ty in unique(x$cell_type)){
      z<-x[x$cell_type==ty,,drop=FALSE];ck<-z$value[z$treatment=='CK']
      for(gr in GROUPS){v<-z$value[z$treatment==gr]
      delta[[paste(m,ty,gr)]]<-data.frame(module=m,cell_type=ty,treatment=gr,n_replicates=length(v),
                                          n_CK_replicates=length(ck),mean_score=if(length(v))mean(v) else NA_real_,
                                          delta_vs_CK=if(length(v)==EXPECTED_REPLICATES&&length(ck)==EXPECTED_REPLICATES)mean(v)-mean(ck) else NA_real_)
      }
    }
  }
  if(length(rr))save_table(do.call(rbind,rr),'46_module_per_sample_celltype_scores')
  z<-if(length(delta))do.call(rbind,delta) else NULL;if(!is.null(z))save_table(z,'47_module_treatment_minus_CK')
  z
}
draw_response <- function(z,whichmodule=NULL) {
  if(is.null(z)){draw_note('F | Hormone response','Module gene membership not supplied; response cannot be calculated from a screenshot.');return()}
  if(is.null(whichmodule))whichmodule<-unique(z$module)[1]
  z<-z[z$module==whichmodule,,drop=FALSE];ty<-unique(z$cell_type)
  mat<-matrix(NA_real_,length(ty),length(GROUPS),dimnames=list(ty,GROUPS))
  mat[cbind(match(z$cell_type,ty),match(z$treatment,GROUPS))]<-z$delta_vs_CK
  vv<-abs(mat[is.finite(mat)]);lim<-if(length(vv))max(.1,max(vv)) else 1
  draw_matrix(mat,paste(whichmodule,'| treatment minus CK'),GROUPS,ty,c(-lim,lim))
}

draw_cross_section <- function(title='G | Root transverse schematic',tissue_values=NULL) {
  par(mar=c(1,1,3,1));plot(NA,xlim=c(-1.6,1.6),ylim=c(-1.25,1.35),asp=1,axes=FALSE,xlab='',ylab='',main=title,cex.main=.85)
  ty<-c('Xylem','Procambium','Pericycle','Endodermis','Cortex','Epidermis');rad<-c(0,.23,.34,.43,.52,.90,1)
  cols<-type_palette(ty)
  if(!is.null(tissue_values))cols<-setNames(numeric_colors(tissue_values[ty],c(-1,1),Z_COLORS),ty)
  for(i in seq_along(ty)){
    n<-if(i==1)8 else round(18+rad[i+1]*22);th<-seq(0,2*pi,length.out=n+1)
    for(j in seq_len(n)){
      aa<-seq(th[j],th[j+1],length.out=8)
      polygon(c(rad[i]*cos(aa),rev(rad[i+1]*cos(aa))),c(rad[i]*sin(aa),rev(rad[i+1]*sin(aa))),
              col=cols[ty[i]],border='white',lwd=.6)
    }
  }
  # 原生木质部/韧皮部位置仅为抽象符号，不宣称5日龄实测细胞数/极数/层数。
  for(a in seq(0,2*pi,length.out=6)[-6])points(.28*cos(a),.28*sin(a),pch=21,bg='#DDCC77',cex=1.3)
  text(1.05,seq(.9,-.6,length.out=length(ty)),rev(ty),adj=0,cex=.5)
  text(0,-1.15,'Anatomical schematic; not histology or measured cell layers',cex=.52)
}
draw_long_section <- function() {
  par(mar=c(1,1,3,1));plot(NA,xlim=c(-1.6,2.4),ylim=c(-.65,3.5),asp=1,axes=FALSE,xlab='',ylab='',
                           main='H | Root longitudinal schematic',cex.main=.85)
  cols<-c('#4477AA','#66C2A5','#AA4499','#CC6677','#44AA99','#8856A7','#44AA99','#CC6677','#AA4499','#66C2A5','#4477AA')
  xs<-seq(-1,1,length.out=length(cols)+1);ys<-c(.15,.32,.49,.66,.83,1,1.25,1.55,1.95,2.45,3.1)
  for(i in seq_along(cols))for(j in seq_len(length(ys)-1)){
    y<-ys[j];width<-.58+.42*min(1,y/1.4)
    rect(xs[i]*width,y,xs[i+1]*width,ys[j+1],col=adjustcolor(cols[i],.75),border='white')
  }
  polygon(c(-.6,-.4,0,.4,.6),c(.15,-.25,-.45,-.25,.15),col='#BDBDBD',border='white')
  rect(-.19,.1,.19,.27,col='#D73027',border='white');text(0,.18,'QC',cex=.6)
  text(1.1,c(.1,.65,1.65,2.75),c('Stem-cell niche','Meristem','Elongation','Differentiation'),adj=0,cex=.6)
  arrows(-1.25,.35,-1.25,2.95,length=.08);text(-1.43,1.6,'Proximal maturation',srt=90,cex=.58)
  text(0,-.57,'Schematic only; no inferred radial developmental arrows',cex=.52)
}
histology_panels <- function() {
  if(!nzchar(HISTOLOGY_MANIFEST)) {
    save_plot('G_root_anatomy_schematics',function(){par(mfrow=c(1,2));draw_cross_section();draw_long_section()},12,6)
    note('histology_status','No original histology images supplied. Schematics are labelled as schematics, not measured 5-day soybean sections. No source screenshot is relabelled as your experiment.')
    return(list(cross=draw_cross_section,long=draw_long_section))
  }
  h<-need(read_table(HISTOLOGY_MANIFEST),c('file','treatment','orientation','label'))
  if(!nrow(h))stop('Histology manifest is empty.')
  h$treatment<-canon_group(h$treatment);if(anyNA(h$treatment))stop('Histology treatment must be one of six groups.')
  paths<-h$file;rel<-!grepl('^(/|[A-Za-z]:[/\\\\])',paths)
  paths[rel]<-file.path(dirname(normalizePath(HISTOLOGY_MANIFEST)),paths[rel])
  if(any(!file.exists(paths)))stop('Histology files missing: ',paste(paths[!file.exists(paths)],collapse=', '))
  save_table(h,'48_histology_inputs')
  panels<-list(cross=draw_cross_section,long=draw_long_section);assigned<-character()
  for(i in seq_len(nrow(h))) {
    ext<-tolower(tools::file_ext(paths[i]));pkg<-if(ext=='png')'png' else if(ext%in%c('jpg','jpeg'))'jpeg' else ''
    if(!nzchar(pkg))stop('Export histology to PNG/JPEG preserving original scale bar: ',paths[i])
    if(!requireNamespace(pkg,quietly=TRUE))stop('Install package ',pkg,' to read actual histology.')
    img<-if(pkg=='png')png::readPNG(paths[i]) else jpeg::readJPEG(paths[i])
    fn<-local({im<-img;j<-i;function(){
      par(mar=c(1,1,3,1));hh<-dim(im)[1];ww<-dim(im)[2]
      plot(NA,xlim=c(0,ww),ylim=c(0,hh),asp=1,axes=FALSE,xlab='',ylab='',
           main=paste(h$treatment[j],h$orientation[j],h$label[j],sep=' | '),cex.main=.9)
      rasterImage(im,0,0,ww,hh);mtext('Original supplied image; no synthetic staining or invented scale bar',side=1,cex=.6)
    }})
    save_plot(paste0('G_real_histology_',i),fn,8,8)
    orientation<-tolower(trimws(as.character(h$orientation[i])))
    slot<-if(orientation%in%c('transverse','cross','横切'))'cross' else if(orientation%in%c('longitudinal','纵切'))'long' else ''
    if(nzchar(slot)&&!slot%in%assigned){panels[[slot]]<-fn;assigned<-c(assigned,slot)}
  }
  panels # 主版图使用manifest中每方向第一张；每张真实图片另存，缺方向才用示意图。
}

smooth_bins <- function(x,y) {
  out<-y;ii<-which(is.finite(y));if(!length(ii))return(out)
  runs<-split(ii,cumsum(c(TRUE,diff(ii)!=1)))
  for(j in runs)if(length(j)>=4&&length(unique(x[j]))>=4)
    out[j]<-predict(smooth.spline(x[j],y[j],df=min(4,length(j)-1)),x[j])$y
  out # 不跨缺失区，不向端点外推；用于作图，未生成额外显著性。
}
array_long <- function(a,centers) {
  ind<-expand.grid(feature=dimnames(a)[[1]],bin=seq_len(dim(a)[2]),treatment=dimnames(a)[[3]],stringsAsFactors=FALSE)
  ind$pseudotime<-centers[ind$bin];ind$value<-as.vector(a);ind
}
path_analysis <- function(p,r,counts,sf,scores,dynamic) {
  dd<-r$d[p$ix,,drop=FALSE];lim<-range(dd$pseudotime)
  if(length(unique(dd$pseudotime))<N_BINS/2||diff(lim)<=0)return(invisible(NULL))
  br<-seq(lim[1],lim[2],length.out=N_BINS+1L);rd<-ridge_data(dd)
  if(!is.null(rd))save_table(rd,paste0('D_',p$id,'_replicate_density'))
  save_table(dd,paste0('D_',p$id,'_cells'))
  save_plot(paste0('D_',p$id,'_sixgroup_ridges'),function()draw_ridges(dd,rd,paste(p$id,'| CK + five hormones')),9,6)
  if(ncol(scores)) {
    b<-bin_replicates(t(scores[p$ix,,drop=FALSE]),dd,br);a<-group_bins(b)
    save_table(array_long(a,b$centers),paste0('D_',p$id,'_module_bins'))
    save_table(data.frame(bin=seq_len(N_BINS),pseudotime=b$centers,b$counts,check.names=FALSE),paste0('D_',p$id,'_bin_cell_counts'))
    audit<-list()
    for(m in seq_len(ncol(scores))){
      for(gr in GROUPS){
        y<-a[m,,gr];ok<-which(is.finite(y));ck<-a[m,,'CK'];common<-which(is.finite(y)&is.finite(ck))
        audit[[paste(m,gr)]]<-data.frame(path=p$id,module=colnames(scores)[m],treatment=gr,
                                         supported_bins=length(ok),peak_pt=if(length(ok))b$centers[ok[which.max(y[ok])]] else NA_real_,
                                         peak_relative_position=if(length(ok))(b$centers[ok[which.max(y[ok])]]-lim[1])/diff(lim) else NA_real_,
                                         mean_delta_vs_CK_on_common_bins=if(length(common)>=4)mean(y[common]-ck[common]) else NA_real_,
                                         common_bins_with_CK=length(common))
      }
      mm<-m
      save_plot(paste0('D_',p$id,'_',safe_name(colnames(scores)[m]),'_program'),local({j<-mm;function(){
        yr<-range(a[j,,],na.rm=TRUE);if(!all(is.finite(yr)))yr<-c(-1,1)
        par(mfrow=c(2,3),mar=c(3,3,3,1))
        for(gr in GROUPS){
          plot(NA,xlim=lim,ylim=yr,xlab='Pseudotime',ylab='Module expression score',main=paste(colnames(scores)[j],gr),cex.main=.8)
          for(s in which(b$group==gr))lines(b$centers,b$value[j,,s],col=adjustcolor(GROUP_COLORS[gr],.5),lty=3,lwd=.8)
          lines(b$centers,smooth_bins(b$centers,a[j,,gr]),col=GROUP_COLORS[gr],lwd=2)
          if(gr!='CK')lines(b$centers,smooth_bins(b$centers,a[j,,'CK']),col=GROUP_COLORS['CK'],lwd=1,lty=2)
          abline(h=0,col='#CCCCCC',lty=3)
        }
      }}),12,7)
    }
    save_table(do.call(rbind,audit),paste0('D_',p$id,'_module_peak_and_response'))
  }
  if(!length(dynamic))return(invisible(NULL))
  ex<-log_expression(counts,sf,dynamic)[,p$ix,drop=FALSE]
  b<-bin_replicates(ex,dd,br);a<-group_bins(b);rm(ex)
  save_table(array_long(a,b$centers),paste0('D_',p$id,'_dynamic_expression_raw_bins'))
  for(i in seq_len(dim(a)[1]))for(gr in GROUPS)a[i,,gr]<-smooth_bins(b$centers,a[i,,gr])
  # 每个基因在所有有效处理x箱上共用Z标准化；各处理不能分别归一化后声称幅度可比。
  z<-a
  for(i in seq_len(dim(a)[1])){
    v<-as.vector(a[i,,]);s<-sd(v,na.rm=TRUE);mu<-mean(v,na.rm=TRUE)
    z[i,,]<-if(is.finite(s)&&s>1e-6)(a[i,,]-mu)/s else NA_real_
  }
  pooled<-apply(z,c(1,2),function(x)if(any(is.finite(x)))mean(x,na.rm=TRUE) else NA_real_)
  enough<-rowSums(is.finite(pooled))>=6
  pk<-apply(pooled,1,function(x)if(any(is.finite(x)))which.max(replace(x,!is.finite(x),-Inf)) else Inf)
  ordergenes<-which(enough)[order(pk[enough])]
  if(!length(ordergenes)) {note(paste0('D_',p$id,'_heatmap_status'),'No genes with enough replicate-supported bins.');return(invisible(NULL))}
  z<-z[ordergenes,,,drop=FALSE]
  save_table(array_long(a,b$centers),paste0('D_',p$id,'_dynamic_expression_all_bins'))
  save_table(array_long(z,b$centers),paste0('D_',p$id,'_dynamic_shared_Z'))
  save_plot(paste0('D_',p$id,'_ridges_and_dynamic_heatmaps'),function(){
    layout(rbind(c(1,1,1),c(2,3,4),c(5,6,7)),heights=c(1.2,1,1))
    draw_ridges(dd,rd,paste(p$id,'| shared trunk may occur in other paths'))
    labs<-rep('',N_BINS);ind<-unique(round(seq(1,N_BINS,length.out=5)));labs[ind]<-format(signif(b$centers[ind],2),trim=TRUE)
    for(gr in GROUPS)draw_matrix(matrix(z[,,gr],nrow=dim(z)[1],ncol=dim(z)[2]),paste(gr,'| spatially associated genes'),labs,
                                 ylabs=character(),lim=c(-2,2))
  },14,14)
}

plot_graph_summary <- function(r,scores,go) {
  for(rt in r$roots){
    cc<-r$nodes$component[match(rt,r$nodes$node)];ns<-r$nodes$node[r$nodes$component==cc]
    sub<-igraph::induced_subgraph(r$graph,ns);v<-igraph::V(sub)$name
    layout<-igraph::layout_as_tree(sub,root=match(rt,v),mode='all',flip.y=TRUE)
    # 节点dominant标签只描述附近观测细胞；不能作为谱系追踪结论。
    types<-vapply(v,function(n){ii<-which(r$d$node==n);if(!length(ii))return('');names(which.max(table(r$d$cell_type[ii])))},character(1))
    show<-igraph::degree(sub)!=2|v==rt
    lab<-ifelse(show,paste(v,types,sep='\n'),'')
    save_table(data.frame(component=cc,node=v,dominant_observed_type=types,schematic_x=layout[,1],schematic_y=layout[,2]),paste0('F_component_',cc,'_node_annotations'))
    save_plot(paste0('F_component_',cc,'_inferred_graph'),function(){
      par(mar=c(2,2,4,2));igraph::plot.igraph(sub,layout=layout,vertex.size=ifelse(v==rt,7,2.5),
                                              vertex.color=ifelse(v==rt,'#D73027','#7F7F7F'),vertex.label=lab,vertex.label.cex=.45,
                                              vertex.frame.color=NA,edge.color='#A0A0A0',edge.arrow.size=0,
                                              main=paste('Component',cc,'| observed graph topology; candidate root:',rt))
      mtext('Computational state relationships; not lineage tracing or a radial tissue-development map',side=1,cex=.65)
    },12,9)
  }
}

main_v8 <- function() {
  set.seed(SEED);dir.create(OUT_DIR,recursive=TRUE,showWarnings=FALSE)
  make_templates();report('Output root retained:',OUT_DIR)
  z<-load_or_build();r<-root_graph(z$cds);cds<-r$cds;d<-r$d
  counts<-SummarizedExperiment::assay(cds,'counts');sf<-monocle3::size_factors(cds)
  if(!is.null(names(sf)))sf<-sf[colnames(counts)]
  if(length(sf)!=ncol(counts)||any(!is.finite(sf)|sf<=0))stop('Invalid saved size factors.')
  detected<-Matrix::rowSums(counts>0)>=MIN_GENE_CELLS;tycols<-type_palette(sort(unique(d$cell_type)))
  save_table(as.data.frame(table(treatment=d$treatment,sample=d$sample,cell_type=d$cell_type)),'10_cells_per_sample_type')
  # 样本内捕获细胞比例，不能等同组织绝对数量/增殖率。
  comp<-as.data.frame(table(sample=d$sample,cell_type=d$cell_type));names(comp)[3]<-'n'
  totals<-table(d$sample);comp$proportion<-comp$n/as.numeric(totals[as.character(comp$sample)])
  comp$treatment<-d$treatment[match(comp$sample,d$sample)];save_table(comp,'11_captured_cell_proportions')
  save_plot('A1_originalUMAP_celltypes',function(){
    par(mfrow=c(2,3));for(gr in GROUPS)draw_umap(d,d$cell_type,paste('Soybean 5-day root |',gr),categorical=tycols,ix=which(d$treatment==gr))
  },15,10)
  rooted<-r$root_audit$component[!is.na(r$root_audit$root_node)]
  for(cc in rooted){
    pt<-d$pseudotime;pt[d$component!=cc]<-NA_real_;lim<-range(pt[is.finite(pt)])
    if(!all(is.finite(lim)))next
    save_table(data.frame(component=cc,pseudotime=seq(lim[1],lim[2],length.out=length(PT_COLORS)),colour=PT_COLORS),paste0('A2_component_',cc,'_PT_scale'))
    save_plot(paste0('A2_component_',cc,'_pseudotime_six_treatments'),local({val<-pt;limit<-lim;com<-cc;function(){
      par(mfrow=c(2,3));for(gr in GROUPS)draw_umap(d,val,paste(gr,'| component',com),r$edges,limit,ix=which(d$treatment==gr))
    }}),15,10)
  }
  modules<-load_modules(rownames(counts));scores<-module_scores(counts,sf,modules$sets)
  if(ncol(scores))saveRDS(list(cell=d$cell,scores=scores,definition='mean gene-wise Z of log1p(count/sizefactor)'),file.path(OUT_DIR,'v8_module_scores.rds'))
  go<-module_GO(modules$sets,rownames(counts),detected)
  markers<-prepare_markers(counts,sf,d);similarity<-similarity_matrix(counts,sf,d,detected)
  response<-module_response(scores,d)
  save_plot('B_marker_dotplot',function()draw_dot(markers$data),11,8)
  save_plot('C_celltype_expression_similarity',function()draw_matrix(similarity,'C | Replicate-averaged expression correlation',colnames(similarity),rownames(similarity),c(-1,1)),10,8)
  save_plot('E_WGCNA_GO',function()draw_GO(go),13,11)
  for(m in colnames(scores)){
    value<-scores[,m];good<-value[is.finite(value)];lim<-if(length(good))as.numeric(quantile(good,c(.01,.99))) else c(-1,1)
    save_table(data.frame(module=m,lower=lim[1],upper=lim[2],clipped=sum(value<lim[1]|value>lim[2],na.rm=TRUE)),paste0('E_',safe_name(m),'_scale'))
    save_plot(paste0('E_',safe_name(m),'_UMAP_sixgroups'),local({v<-value;ll<-lim;mm<-m;function(){
      par(mfrow=c(2,3));for(gr in GROUPS)draw_umap(d,v,paste(mm,gr,'| expression score'),lim=ll,pal=Z_COLORS,ix=which(d$treatment==gr))
    }}),14,9)
    save_plot(paste0('F_',safe_name(m),'_response_vs_CK'),local({mm<-m;function()draw_response(response,mm)}),9,8)
  }
  # 六候选基因UMAP，始终使用同一原UMAP；每个基因跨处理固定表达色标。
  mp<-map_genes(TARGET_GENES,rownames(counts));save_table(mp,'49_sixgene_mapping')
  for(i in seq_len(nrow(mp))) {
    if(is.na(mp$feature_id[i]))next
    v<-as.numeric(log_expression(counts,sf,mp$feature_id[i]));hi<-if(any(v>0))as.numeric(quantile(v[v>0],.99)) else 1
    save_plot(paste0('B_gene_',safe_name(mp$gene_id[i])),local({vv<-v;high<-hi;nm<-mp$gene_id[i];function(){
      par(mfrow=c(2,3));for(gr in GROUPS)draw_umap(d,vv,paste(nm,gr),lim=c(0,high),
                                                   pal=c('#D9D9D9','#54278F','#C51B8A','#FEC44F'),ix=which(d$treatment==gr))
    }}),14,9)
  }
  anatomy<-histology_panels();plot_graph_summary(r,scores,go)
  # 主版图先输出；耗时graph_test/DE不阻塞基础atlas和模块图。
  primary<-if(length(rooted))rooted[which.max(vapply(rooted,function(cc)sum(d$component==cc),numeric(1)))] else NA_integer_
  dd<-d[!is.na(primary)&d$component==primary&is.finite(d$pseudotime),,drop=FALSE]
  rd<-if(nrow(dd))ridge_data(dd) else NULL
  save_plot('Figure_root_hormone_overview_v8',function(){
    par(mfrow=c(3,3))
    draw_umap(d,d$cell_type,'A | Soybean root cell types',r$edges,categorical=tycols)
    pt<-d$pseudotime;pt[is.na(primary)|d$component!=primary]<-NA_real_
    if(any(is.finite(pt)))draw_umap(d,pt,paste('A | Candidate-root PT, component',primary),r$edges) else draw_note('A | Pseudotime','No supported QC root; no orientation assigned.')
    draw_dot(markers$data)
    draw_matrix(similarity,'C | Cell-type similarity',colnames(similarity),rownames(similarity),c(-1,1))
    draw_ridges(dd,rd,'D | Six-treatment pseudotime distribution')
    if(ncol(scores))draw_umap(d,scores[,1],paste('E |',colnames(scores)[1],'expression score'),pal=Z_COLORS) else draw_note('E | WGCNA programs','Membership not supplied. See input_templates.')
    draw_response(response);anatomy$cross();anatomy$long()
  },20,17)
  pseudobulk_DE(counts,d)
  dynamic<-graph_candidates(cds,z$signature);paths<-get_paths(r)
  for(p in paths){report('Processing',p$id,'cells:',length(p$ix));path_analysis(p,r,counts,sf,scores,dynamic)}
  files<-file.path(OUT_DIR,paste0('D_',names(paths),'_module_peak_and_response.csv'))
  files<-files[file.exists(files)] # 只整合当前运行选中的paths，避免读取历史遗留的其他路径。
  if(ncol(scores)&&length(files)){
    ss<-do.call(rbind,lapply(files,read.csv,stringsAsFactors=FALSE))
    ss$descriptive_peak_window<-as.character(cut(ss$peak_relative_position,c(0,1/3,2/3,1),
                                                 include.lowest=TRUE,labels=c('relative early','relative middle','relative late')))
    ss$GO_terms<-NA_character_
    if(!is.null(go)&&nrow(go))for(mm in unique(ss$module)){
      gg<-go[go$module==mm&is.finite(go$FDR)&go$FDR<GO_FDR&go$Count>=2,,drop=FALSE]
      if(nrow(gg))ss$GO_terms[ss$module==mm]<-paste(head(gg$Description[order(gg$FDR)],3),collapse='; ')
    }
    save_table(ss,'F_path_module_function_response_summary')
  }
  capture.output(sessionInfo(),file=file.path(OUT_DIR,'sessionInfo_v8.txt'))
  note('analysis_limits',c(
    'QC annotation defines candidate orientation; inspect root audit/markers and lineage connectivity before biological interpretation.',
    'Pseudotime comparisons are within a rooted component/path. Unrooted cells remain NA. Shared trunks cannot be summed over paths.',
    'Cross-sectional hormone treatment is not real elapsed developmental time; density shifts cannot establish accelerated/delayed development.',
    'Cell proportions reflect captured cells; module means reflect expression, not causal regulation or WGCNA eigengenes.',
    'Pseudobulk edgeR uses biological replicate raw counts. graph_test tests spatial autocorrelation, not hormone contrasts or branch-specific DE.',
    'Schema cartoons are not measured soybean sections. Real supplied images retain their own scale bars.',
    'WGCNA GO may be reused from a supplied table; computed GO uses the explicitly recorded background and global BH.',
    'User screenshot module labels are hypotheses. No gene memberships, q-values or activation direction were recovered from screenshots.'
  ))
  report('Finished:',OUT_DIR)
  invisible(list(root_audit=r$root_audit,n_modules=ncol(scores),n_paths=length(paths)))
}

if(sys.nframe()==0L) main_v8()
# 在RStudio中可 source(..., encoding='UTF-8') 后运行 main_v8()