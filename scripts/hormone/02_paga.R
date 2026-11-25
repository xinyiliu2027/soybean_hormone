INPUT_RDS <- '/data01/liuxy/Project/CK2_added_v3_Results/seurat_obj_with_CK2_v3.rds'
RUN_DIR <- '/data01/liuxy/Project/CK2_added_v3_Results/Endodermis_monocle3_originalUMAP_grouped_v3/20261119_163425'
INPUT_CACHE_DIR <- file.path(RUN_DIR,'QC_CK_PAGA_branches_v9') 
OUT_DIR <- Sys.getenv('PAGA_OUT_DIR', unset=getwd()) 
PYTHON <- '/home/liuxy/micromamba/envs/paga/bin/python' 

PAGA_LIB <- Sys.getenv('PAGA_LIB')
if(nzchar(PAGA_LIB)) .libPaths(unique(c(strsplit(PAGA_LIB, .Platform$path.sep)[[1]], .libPaths())))
ASSAY <- 'RNA'; UMAP_NAME <- 'umap'
TREATMENT_COL <- 'treatment'; CELLTYPE_COL <- 'cell_type'; SAMPLE_COL <- 'orig.ident'
BATCH_COL <- '' 
GROUPS <- c('CK','ABA','eBL','IAA','MeJA','SA','Heat','Cold','PEG')
REQUIRE_ALL_GROUPS <- TRUE
MIN_REPLICATES <- 2L
QC_LABELS <- c('QC','Quiescent Center','Quiescent centre')
ROOT_CELL <- '' 
ROOT_BIOLOGICALLY_VALIDATED <- FALSE
N_HVG <- 2000L; N_PCS <- 30L; N_NEIGHBORS <- 30L
CK_LEIDEN_RESOLUTION <- 0.5 
PAGA_THRESHOLD <- 0.10
THRESHOLD_GRID <- c(.05,.10,.20,.30)
NEIGHBOR_GRID <- c(15L,30L,50L)
MIN_ARM_STATES <- 2L; MIN_ARM_CELLS <- 100L
CK_SUBSAMPLES <- 20L 
EDGE_RESAMPLES <- 5L
EDGE_MIN_PER_STATE_SAMPLE <- 20L; EDGE_MAX_PER_STATE_SAMPLE <- 200L
MAP_CONFIDENCE <- 0.6; MAP_RADIUS_MULTIPLIER <- 1.5
MIN_SAMPLE_BRANCH_CELLS <- 30L
N_BINS <- 20L; MIN_BIN_SAMPLE_CELLS <- 5L
GENES_PER_BRANCH <- 30L; DE_FDR <- .05; DE_LOGFC <- .5
RUN_DE <- TRUE; RUN_GO <- TRUE
WGCNA_MEMBERSHIP_FILE <- '' 
ORGDB_PACKAGE <- 'org.G.max.eg.db'; GO_FDR <- .05
REBUILD_EXPORT <- FALSE; REBUILD_PAGA <- FALSE
SEED <- 20261121L; PNG_DPI <- 200L
GROUP_COLORS <- c(CK='#5E5E5E',ABA='#CC6677',eBL='#AA8800',IAA='#228833',MeJA='#AA4499',
                  SA='#4477AA',Heat='#D55E00',Cold='#56B4E9',PEG='#009E73')
BRANCH_COLORS <- c(Shared_trunk='#6A3D9A',Branch_1='#0072B2',Branch_2='#D55E00',
                   Branch_3='#009E73',Branch_4='#CC79A7',Outside_candidate_arms='#BDBDBD',Unmapped='#E3E3E3')
PT_COLORS <- c('#30123B','#4145AB','#22A7A2','#A0DA39','#F9D057','#E34A33')
Z_COLORS <- colorRampPalette(c('#5E9C45','#F7F7F7','#AA2676'))(101)

need <- function(x,cols){if(!all(cols%in%names(x)))stop('Required columns: ',paste(cols,collapse=', '));x}
save_table <- function(x,name)write.csv(x,file.path(OUT_DIR,paste0(name,'.csv')),row.names=FALSE,na='',fileEncoding='UTF-8')
read_table <- function(path)read.table(path,header=TRUE,sep=if(grepl('\\.(tsv|txt)$',path))'\t' else ',',
                                       quote='"',comment.char='',check.names=FALSE,stringsAsFactors=FALSE,fileEncoding='UTF-8')
safe_name <- function(x)gsub('[^A-Za-z0-9_.-]','_',x)
note <- function(name,txt){writeLines(txt,file.path(OUT_DIR,paste0(name,'.txt')),useBytes=TRUE);cat(paste(txt,collapse='\n'),'\n')}
save_plot <- function(name,fn,w=12,h=8){
  pdf(file.path(OUT_DIR,paste0(name,'.pdf')),width=w,height=h,useDingbats=FALSE);tryCatch(fn(),finally=dev.off())
  png(file.path(OUT_DIR,paste0(name,'.png')),width=w*PNG_DPI,height=h*PNG_DPI,res=PNG_DPI,
      type=if(capabilities('cairo'))'cairo' else getOption('bitmapType','Xlib'));tryCatch(fn(),finally=dev.off())
}
colors_num <- function(v,lim,pal=Z_COLORS){
  co<-rep('#CFCFCF',length(v));ok<-is.finite(v);u<-if(diff(lim)>0)(v[ok]-lim[1])/diff(lim) else rep(.5,sum(ok))
  co[ok]<-colorRampPalette(pal)(256)[1L+floor(255*pmin(1,pmax(0,u)))];co
}
message_panel <- function(title,txt){par(mar=c(2,2,3,1));plot.new();title(main=title,cex.main=1);text(.03,.85,paste(strwrap(txt,75),collapse='\n'),adj=c(0,1),cex=.8)}
PROGRESS_FILE <- file.path(OUT_DIR,'progress_v10.log')
plog <- function(...){
  line <- paste0(format(Sys.time(),'%H:%M:%S  '), paste0(...,collapse=''))
  cat(line,'\n')
  con <- tryCatch(file(PROGRESS_FILE,open='at'),error=function(e)NULL)
  if(!is.null(con)){writeLines(line,con);close(con)}
  invisible(line)
}
canon_group <- function(x){
  raw <- trimws(as.character(x))
  k <- toupper(gsub('[^A-Za-z0-9]','',raw))
  eq <- c(CK='CK',CONTROL='CK',CONT='CK',MOCK='CK',VECTOR='CK',NONTREAT='CK',
          ABA='ABA',ABSCISICACID='ABA',
          EBL='eBL',EPIBRASSINOLIDE='eBL',BRASSINOLIDE='eBL',
          IAA='IAA',AUXIN='IAA',
          MEJA='MeJA',METHYLJASMONATE='MeJA',
          SA='SA',SALICYLICACID='SA',
          HEAT='Heat',HIGHTEMP='Heat',HIGHTEMPERATURE='Heat',HEATSTRESS='Heat',
          HIGHTEMPERATURESTRESS='Heat',HIGHTEMPERATURETREATMENT='Heat',
          COLD='Cold',LOWTEMP='Cold',LOWTEMPERATURE='Cold',COLDSTRESS='Cold',
          CHILLING='Cold',LOWTEMPERATURETREATMENT='Cold',LOWTEMPERATURESTRESS='Cold',
          PEG='PEG',PEG6000='PEG',PEG8000='PEG',PEG4000='PEG',POLYETHYLENEGLYCOL='PEG')
  out <- unname(eq[k])
  i <- is.na(out); out[i & grepl('HIGH',k) & grepl('TEMP',k)] <- 'Heat'
  i <- is.na(out); out[i & grepl('LOW',k)  & grepl('TEMP',k)] <- 'Cold'
  i <- is.na(out); out[i & grepl('HEAT',k)] <- 'Heat'
  i <- is.na(out); out[i & grepl('COLD|CHILL|FROST',k)] <- 'Cold'
  i <- is.na(out); out[i & grepl('PEG',k)] <- 'PEG'
  i <- is.na(out); out[i & grepl('\u9ad8\u6e29|\u70ed\u6fc0',raw)] <- 'Heat'
  i <- is.na(out); out[i & grepl('\u4f4e\u6e29',raw)] <- 'Cold'
  i <- is.na(out); out[i & grepl('\u5e72\u65f1|\u6e17\u900f',raw)] <- 'PEG'
  out
}
gene_key <- function(x){z<-trimws(x);ok<-grepl('^Glyma[._-]?[0-9]{2}G[0-9]{6}([._-][0-9]+)*([._-]v4)?$',z,ignore.case=TRUE);z[ok]<-sub('^Glyma[._-]?([0-9]{2})G([0-9]{6}).*$','Glyma.\\1G\\2',z[ok],ignore.case=TRUE);z}
map_genes <- function(query,features){
  idx<-split(seq_along(features),gene_key(features));do.call(rbind,lapply(query,function(q){
    hit<-which(features==q);if(length(hit)!=1)hit<-idx[[gene_key(q)]]
    data.frame(gene_id=q,feature_id=if(length(hit)==1)features[hit] else NA_character_,
               status=if(length(hit)==1)'matched' else if(length(hit))'ambiguous' else 'missing')
  }))
}

export_input <- function(){
  for(p in c('Seurat','SeuratObject','Matrix','jsonlite'))if(!requireNamespace(p,quietly=TRUE))
    stop('Missing R package: ',p,' | .libPaths() = ',paste(.libPaths(),collapse=' | '),
         ' | 若包在别的库目录，用 PAGA_LIB 指过去再重跑: PAGA_LIB=<dir> setsid nohup "$RS" soybean_v9_run.R > /home/liuxy/run_v9.log 2>&1 &')
  if(!file.exists(INPUT_RDS))stop('INPUT_RDS not found.')
  sig<-list(path=normalizePath(INPUT_RDS),size=file.info(INPUT_RDS)$size,mtime=file.info(INPUT_RDS)$mtime,
            cols=c(TREATMENT_COL,CELLTYPE_COL,SAMPLE_COL,BATCH_COL),assay=ASSAY,umap=UMAP_NAME,groups=GROUPS,schema=9L)
  sigfile<-file.path(OUT_DIR,'export_signature.rds');cache<-file.path(OUT_DIR,'counts_metadata_v9.rds')
  if(file.exists(cache)&&file.exists(sigfile)&&!REBUILD_EXPORT){
    if(!identical(readRDS(sigfile),sig))stop('Input/export settings changed: set REBUILD_EXPORT=TRUE and REBUILD_PAGA=TRUE.')
    plog('Resume exported v9 counts/metadata: ',cache)
    ans<-readRDS(cache);validate_cached_export(ans)
    return(ans)
  }
  plog('readRDS: ',INPUT_RDS)
  obj<-readRDS(INPUT_RDS);md<-obj[[]];need(md,c(TREATMENT_COL,CELLTYPE_COL,SAMPLE_COL))
  plog('readRDS done: ',nrow(md),' cells x ',nrow(obj),' features')
  tr<-canon_group(md[[TREATMENT_COL]])
  aud<-as.data.frame(table(original=md[[TREATMENT_COL]],canonical=ifelse(is.na(tr),'<UNMATCHED>',tr)))
  aud<-aud[aud$Freq>0,,drop=FALSE]
  save_table(aud,'00_treatment_audit')
  plog('treatment labels: ',paste(unique(paste(aud$original,'->',aud$canonical)),collapse=' ; '))
  keep<-!is.na(tr)&tr%in%GROUPS;md<-md[keep,,drop=FALSE];tr<-tr[keep]
  if(REQUIRE_ALL_GROUPS&&length(setdiff(GROUPS,tr)))
    stop('Groups missing: ',paste(setdiff(GROUPS,tr),collapse=', '),
         ' | 未匹配的原始标签: ',paste(unique(as.character(aud$original[aud$canonical=='<UNMATCHED>'])),collapse=', '),
         ' | 映射明细见 00_treatment_audit.csv；需要补别名就改 canon_group() 的 eq 表，不要用 REQUIRE_ALL_GROUPS=FALSE 硬跳过。')
  if(anyNA(md[[CELLTYPE_COL]])||anyNA(md[[SAMPLE_COL]]))stop('Missing annotation/sample IDs.')
  if(any(!nzchar(trimws(as.character(md[[SAMPLE_COL]])))))stop('Empty sample IDs.')
  xy<-SeuratObject::Embeddings(obj,UMAP_NAME)[rownames(md),1:2,drop=FALSE]
  if(any(!is.finite(xy)))stop('Invalid/missing original UMAP.')
  d<-data.frame(cell=rownames(md),treatment=tr,sample=paste(tr,md[[SAMPLE_COL]],sep='::'),
                cell_type=as.character(md[[CELLTYPE_COL]]),x=xy[,1],y=xy[,2],row.names=NULL)
  d$batch<-if(nzchar(BATCH_COL)){need(md,BATCH_COL);as.character(md[[BATCH_COL]])} else 'not_supplied'
  samples<-unique(d[,c('treatment','sample','batch')]);save_table(samples,'01_biological_samples')
  if(any(table(samples$sample)>1))stop('Same sample spans multiple batches; check sample/batch metadata.')
  nr<-table(factor(samples$treatment,levels=intersect(GROUPS,d$treatment)))
  if(any(nr<MIN_REPLICATES))stop('Need >=2 independent biological samples per group; inspect 01_biological_samples.csv.')
  obj<-subset(obj,cells=d$cell)
  if(inherits(obj[[ASSAY]],'Assay5')){
    ll<-grep('^counts($|[.])',SeuratObject::Layers(obj[[ASSAY]]),value=TRUE);if(!length(ll))stop('No RNA counts.')
    if(length(ll)!=1||ll!='counts')obj[[ASSAY]]<-SeuratObject::JoinLayers(obj[[ASSAY]],layers=ll,new='counts')
    ct<-SeuratObject::LayerData(obj,assay=ASSAY,layer='counts')
  }else ct<-SeuratObject::GetAssayData(obj,assay=ASSAY,slot='counts')
  ct<-methods::as(ct[,d$cell,drop=FALSE],'dgCMatrix');rm(obj);gc()
  if(anyDuplicated(rownames(ct))||anyDuplicated(d$cell))stop('Duplicate gene/cell IDs.')
  if(any(!is.finite(ct@x))||any(ct@x<0)||any(abs(ct@x-round(ct@x))>1e-6))stop('Use nonnegative integer RNA counts.')
  ok<-Matrix::colSums(ct)>0;ct<-ct[,ok,drop=FALSE];d<-d[ok,,drop=FALSE]
  ct<-ct[Matrix::rowSums(ct)>0,,drop=FALSE]
  save_table(d,'input_metadata');save_table(data.frame(gene_id=rownames(ct)),'input_genes')
  Matrix::writeMM(ct,file.path(OUT_DIR,'input_counts.mtx'))
  plog('counts written: ',nrow(ct),' genes x ',ncol(ct),' cells')
  ans<-list(counts=ct,metadata=d,signature=sig);saveRDS(ans,cache);saveRDS(sig,sigfile);ans
}
run_paga <- function(z,phase=c('reference','full')){
  phase<-match.arg(phase)
  cfg<-list(out=normalizePath(OUT_DIR,winslash='/',mustWork=TRUE),seed=SEED,groups=GROUPS,qc_labels=QC_LABELS,
            root_cell=ROOT_CELL,root_validated=ROOT_BIOLOGICALLY_VALIDATED,min_qc=10L,hvgs=N_HVG,pcs=N_PCS,
            neighbors=N_NEIGHBORS,resolution=CK_LEIDEN_RESOLUTION,threshold=PAGA_THRESHOLD,
            threshold_grid=THRESHOLD_GRID,neighbor_grid=NEIGHBOR_GRID,min_arm_states=MIN_ARM_STATES,min_arm_cells=MIN_ARM_CELLS,
            resamples=CK_SUBSAMPLES,edge_resamples=EDGE_RESAMPLES,edge_min_per_state_sample=EDGE_MIN_PER_STATE_SAMPLE,
            edge_max_per_state_sample=EDGE_MAX_PER_STATE_SAMPLE,map_confidence=MAP_CONFIDENCE,
            radius_multiplier=MAP_RADIUS_MULTIPLIER,similarity_cells=600L)
  cfg$phase<-phase
  signature<-list(input=z$signature,config=cfg,engine=PYTHON_ENGINE)
  sf<-file.path(OUT_DIR,paste0('paga_',phase,'_signature_v10.rds'));final<-file.path(OUT_DIR,'cells_CK_reference_mapping.csv')
  if(file.exists(sf)&&file.exists(final)&&!REBUILD_PAGA&&!REBUILD_EXPORT){
    required<-if(phase=='reference')c('CK_nodes.csv','CK_edges_all.csv') else c('CK_branch_evidence.csv','CK_early_cell_PCs.csv')
    if(identical(readRDS(sf),signature)&&all(file.exists(file.path(OUT_DIR,required))))return(invisible(NULL))
    plog('PAGA settings/checkpoint changed; recompute this phase, keeping original counts export.')
  }
  pyfile<-file.path(OUT_DIR,'paga_v10_engine.py');conf<-file.path(OUT_DIR,paste0('paga_config_',phase,'_v10.json'))
  writeLines(PYTHON_ENGINE,pyfile,useBytes=TRUE);jsonlite::write_json(cfg,conf,auto_unbox=TRUE,pretty=TRUE)
  plog('Running actual Scanpy PAGA using ',PYTHON)
  Sys.setenv(PYTHONUNBUFFERED='1')
  pylog<-file.path(OUT_DIR,paste0('python_',phase,'_v9.log'))
  plog('Python phase ',phase,'; live log: ',pylog)
  status<-system2(PYTHON,c('-u',shQuote(pyfile),shQuote(conf)),stdout=pylog,stderr=pylog)
  plog('python engine exit status = ',status)
  if(file.exists(pylog))cat(paste(tail(readLines(pylog,warn=FALSE),60),collapse='\n'),'\n')
  if(status!=0)stop('Python ',phase,' failed. Exact error: ',pylog,
                    ' ; exception traceback: ',file.path(OUT_DIR,'python_traceback_v9.txt'),
                    ' . Existing exports and any completed reference figures are preserved.')
  saveRDS(signature,sf)
}

exact_permutation <- function(a,b){
  if(length(a)<2||length(b)<2||any(!is.finite(c(a,b))))return(NA_real_)
  v<-c(a,b);n<-length(a);obs<-abs(mean(a)-mean(b))
  if(choose(length(v),n)>20000)return(NA_real_)
  cc<-combn(seq_along(v),n);stat<-apply(cc,2,function(j)abs(mean(v[j])-mean(v[-j])))
  mean(stat>=obs-1e-12)
}
occupancy_stats <- function(d){
  lev<-unique(c('Shared_trunk',sort(unique(d$branch[grepl('^Branch_',d$branch)])),
                'Outside_candidate_arms','Unmapped'));ss<-unique(d$sample);rr<-list()
                for(s in ss){z<-d[d$sample==s,];n<-table(factor(z$branch,levels=lev));post<-sum(grepl('^Branch_',z$branch))
                rr[[s]]<-data.frame(sample=s,treatment=z$treatment[1],branch=lev,n=as.integer(n),total=nrow(z),
                                    fraction_all=as.numeric(n)/nrow(z),fraction_postbranch=if(post>0)ifelse(grepl('^Branch_',lev),as.numeric(n)/post,NA_real_) else NA_real_,
                                    unmapped_fraction=mean(z$branch=='Unmapped'))
                };z<-do.call(rbind,rr);save_table(z,'20_branch_occupancy_per_sample')
                tests<-list()
                for(b in lev)for(gr in setdiff(GROUPS,'CK'))for(metric in c('fraction_all','fraction_postbranch')){
                  aa<-z[z$branch==b&z$treatment==gr,metric];bb<-z[z$branch==b&z$treatment=='CK',metric]
                  if(!length(aa))next
                  allowed<-!nzchar(BATCH_COL)||length(unique(d$batch[d$treatment%in%c('CK',gr)]))==1
                  tests[[paste(b,gr,metric)]]<-data.frame(branch=b,treatment=gr,metric=metric,
                                                          delta=mean(aa)-mean(bb),n_treatment=length(aa),n_CK=length(bb),
                                                          pvalue=if(allowed)exact_permutation(aa,bb) else NA_real_,
                                                          test_status=if(allowed)'sample-label permutation; verify exchangeability' else 'multiple known batches: descriptive only')
                }
                tt<-do.call(rbind,tests);if(!is.null(tt)){tt$FDR<-p.adjust(tt$pvalue,'BH');save_table(tt,'21_branch_occupancy_exact_sample_tests')}
                z
}
balanced_edge_stats <- function(d=NULL){
  results<-list()
  for(gr in setdiff(GROUPS,'CK')){
    f<-file.path(OUT_DIR,paste0(gr,'_balanced_edges_replicates.csv'));if(!file.exists(f))next
    z<-read.csv(f);if(!nrow(z))next
    m<-aggregate(connectivity~sample+treatment+source+target,z,mean)
    save_table(m,paste0(gr,'_balanced_edges_sample_means'))
    rows<-lapply(split(m,paste(m$source,m$target,sep='|')),function(a){
      x<-a$connectivity[a$treatment==gr];y<-a$connectivity[a$treatment=='CK']
      allowed<-!nzchar(BATCH_COL)||(!is.null(d)&&length(unique(d$batch[d$treatment%in%c('CK',gr)]))==1)
      data.frame(treatment=gr,source=a$source[1],target=a$target[1],CK_mean=mean(y),treatment_mean=mean(x),
                 delta=mean(x)-mean(y),n_CK=length(y),n_treatment=length(x),pvalue=if(allowed)exact_permutation(x,y) else NA_real_)
    });results[[gr]]<-do.call(rbind,rows)
  }
  ans<-if(length(results))do.call(rbind,results) else NULL
  if(!is.null(ans)){ans$FDR<-p.adjust(ans$pvalue,'BH');save_table(ans,'22_balanced_edge_changes')};ans
}

expression <- function(ct,genes){
  sf<-Matrix::colSums(ct)/1e4;e<-ct[genes,,drop=FALSE]%*%Matrix::Diagonal(x=1/sf);e@x<-log1p(e@x);e
}
branch_programs <- function(ct,d){
  arms<-sort(unique(d$branch[grepl('^Branch_',d$branch)]));sets<-list();sigsets<-list();tables<-list();backgrounds<-list()
  if(length(arms)<2)return(list(sets=sets,sigsets=sigsets,background=character()))
  if(!RUN_DE||!requireNamespace('edgeR',quietly=TRUE)||!requireNamespace('statmod',quietly=TRUE)){
    note('DE_status','edgeR/statmod unavailable or RUN_DE=FALSE. Topology, occupancy and connectivity run; no gene-level DE or GO conclusions.');return(list(sets=sets,sigsets=sigsets,background=character()))
  }
  for(arm in arms){
    ix<-which(d$treatment=='CK'&d$branch%in%arms);md<-d[ix,,drop=FALSE];md$contrast_class<-ifelse(md$branch==arm,'target','other_arms')
    key<-paste(md$sample,md$contrast_class,sep='|');groups<-split(seq_along(ix),key)
    groups<-groups[lengths(groups)>=MIN_SAMPLE_BRANCH_CELLS]
    if(!length(groups))next
    sm<-do.call(rbind,lapply(groups,function(j)md[j[1],c('sample','contrast_class'),drop=FALSE]));rownames(sm)<-names(groups)
    ok<-names(which(table(sm$sample)==2));keep<-sm$sample%in%ok;groups<-groups[keep];sm<-sm[keep,,drop=FALSE]
    if(length(ok)<2)next
    pb<-do.call(cbind,lapply(groups,function(j)Matrix::rowSums(ct[,ix[j],drop=FALSE])));rownames(pb)<-rownames(ct)
    design<-model.matrix(~factor(sample)+factor(contrast_class,levels=c('other_arms','target')),sm)
    yy<-edgeR::DGEList(pb);kk<-edgeR::filterByExpr(yy,design=design);if(sum(kk)<30)next
    yy<-edgeR::calcNormFactors(yy[kk,,keep.lib.sizes=FALSE]);yy<-edgeR::estimateDisp(yy,design,robust=TRUE)
    fit<-edgeR::glmQLFit(yy,design,robust=TRUE);tt<-edgeR::topTags(edgeR::glmQLFTest(fit,coef=ncol(design)),n=Inf)$table
    tt$gene_id<-rownames(tt);tt$branch<-arm;tables[[arm]]<-tt;backgrounds[[arm]]<-rownames(tt)
  }
  if(length(tables)){
    all<-do.call(rbind,tables);all$FDR_all_branches<-p.adjust(all$PValue,'BH');save_table(all,'30_CK_branch_DE_all')
    for(arm in names(tables)){
      tt<-all[all$branch==arm,,drop=FALSE];ss<-tt[tt$FDR_all_branches<DE_FDR&tt$logFC>DE_LOGFC,,drop=FALSE]
      ss<-ss[order(ss$FDR_all_branches,-ss$logFC),,drop=FALSE]
      sigsets[[arm]]<-ss$gene_id
      if(nrow(ss)){
        sets[[arm]]<-head(ss$gene_id,GENES_PER_BRANCH);kind<-'FDR-supported positive branch-associated genes'
      }else{
        pp<-tt[tt$logFC>0,,drop=FALSE];pp<-pp[order(-pp$logFC),,drop=FALSE]
        sets[[arm]]<-head(pp$gene_id,GENES_PER_BRANCH);kind<-'Exploratory top-positive effects; no significant program'
      }
      if(length(sets[[arm]]))save_table(data.frame(branch=arm,gene_id=sets[[arm]],selection=kind),paste0('31_',arm,'_plot_genes'))
    }
  }
  list(sets=sets[lengths(sets)>0],sigsets=sigsets,background=if(length(backgrounds))Reduce(intersect,backgrounds) else character())
}
treatment_branch_DE <- function(ct,d){
  if(!RUN_DE||!requireNamespace('edgeR',quietly=TRUE)||!requireNamespace('statmod',quietly=TRUE))return(NULL)
  branches<-sort(unique(d$branch[grepl('^Branch_',d$branch)]));all<-list();audit<-list()
  for(arm in branches)for(gr in setdiff(intersect(GROUPS,d$treatment),'CK')){
    ix<-which(d$branch==arm&d$treatment%in%c('CK',gr));g<-split(ix,d$sample[ix]);g<-g[lengths(g)>=MIN_SAMPLE_BRANCH_CELLS]
    if(!length(g))next
    sm<-do.call(rbind,lapply(g,function(j)d[j[1],c('sample','treatment','batch'),drop=FALSE]));nr<-table(factor(sm$treatment,levels=c('CK',gr)))
    tag<-paste(arm,gr,sep='_');status<-'tested'
    if(any(nr<2)){audit[[tag]]<-data.frame(branch=arm,treatment=gr,status='insufficient biological replicates');next}
    sm$condition<-factor(sm$treatment,levels=c('CK',gr))
    design<-if(nzchar(BATCH_COL)&&length(unique(sm$batch))>1)model.matrix(~factor(batch)+condition,sm) else model.matrix(~condition,sm)
    if(qr(design)$rank<ncol(design)||nrow(design)<=ncol(design)){
      audit[[tag]]<-data.frame(branch=arm,treatment=gr,status='batch/treatment confounded or no residual df');next
    }
    pb<-do.call(cbind,lapply(g,function(j)Matrix::rowSums(ct[,j,drop=FALSE])));rownames(pb)<-rownames(ct)
    y<-edgeR::DGEList(pb);kk<-edgeR::filterByExpr(y,design=design);if(sum(kk)<30)next
    y<-edgeR::calcNormFactors(y[kk,,keep.lib.sizes=FALSE]);y<-edgeR::estimateDisp(y,design,robust=TRUE)
    fit<-edgeR::glmQLFit(y,design,robust=TRUE);tt<-edgeR::topTags(edgeR::glmQLFTest(fit,coef=ncol(design)),n=Inf)$table
    tt$gene_id<-rownames(tt);tt$branch<-arm;tt$treatment<-gr;all[[tag]]<-tt
    audit[[tag]]<-data.frame(branch=arm,treatment=gr,status=status)
  }
  if(length(audit))save_table(do.call(rbind,audit),'32_treatment_branch_DE_audit')
  if(length(all)){z<-do.call(rbind,all);z$FDR_global<-p.adjust(z$PValue,'BH');save_table(z,'33_treatment_branch_DE_all');return(z)}
  NULL
}
branch_interaction_DE <- function(ct,d){
  arms<-sort(unique(d$branch[grepl('^Branch_',d$branch)]))
  if(length(arms)!=2||!RUN_DE||!requireNamespace('edgeR',quietly=TRUE)||!requireNamespace('statmod',quietly=TRUE))return(NULL)
  result<-list();aud<-list()
  for(gr in setdiff(intersect(GROUPS,d$treatment),'CK')){
    ix<-which(d$treatment%in%c('CK',gr)&d$branch%in%arms);key<-paste(d$sample[ix],d$branch[ix],sep='|')
    gg<-split(ix,key);gg<-gg[lengths(gg)>=MIN_SAMPLE_BRANCH_CELLS];if(!length(gg))next
    sm<-do.call(rbind,lapply(gg,function(j)d[j[1],c('sample','treatment','branch','batch'),drop=FALSE]))
    paired<-names(which(table(sm$sample)==2));keep<-sm$sample%in%paired;sm<-sm[keep,];gg<-gg[keep]
    nr<-table(factor(unique(sm[,c('sample','treatment')])$treatment,levels=c('CK',gr)))
    if(any(nr<2)){aud[[gr]]<-data.frame(treatment=gr,status='need >=2 samples per group with both branches');next}
    if(nzchar(BATCH_COL)&&length(unique(sm$batch))>1){
      aa<-unique(sm[,c('sample','treatment','batch')]);dm<-model.matrix(~factor(batch)+factor(treatment),aa)
      if(qr(dm)$rank<ncol(dm)){aud[[gr]]<-data.frame(treatment=gr,status='batch-treatment confounded');next}
    }
    sample_design<-model.matrix(~0+factor(sample),sm);b2<-as.numeric(sm$branch==arms[2]);interaction<-b2*as.numeric(sm$treatment==gr)
    batchb<-NULL
    if(nzchar(BATCH_COL)&&length(unique(sm$batch))>1){bb<-model.matrix(~factor(batch),sm)[,-1,drop=FALSE];batchb<-bb*b2}
    design<-cbind(sample_design,branch2=b2,batchb,condition_by_branch2=interaction)
    if(qr(design)$rank<ncol(design)||nrow(design)<=ncol(design)){aud[[gr]]<-data.frame(treatment=gr,status='interaction design not estimable');next}
    pb<-do.call(cbind,lapply(gg,function(j)Matrix::rowSums(ct[,j,drop=FALSE])));rownames(pb)<-rownames(ct)
    y<-edgeR::DGEList(pb);kk<-edgeR::filterByExpr(y,design=design);if(sum(kk)<30)next
    y<-edgeR::calcNormFactors(y[kk,,keep.lib.sizes=FALSE]);y<-edgeR::estimateDisp(y,design,robust=TRUE)
    fit<-edgeR::glmQLFit(y,design,robust=TRUE);tt<-edgeR::topTags(edgeR::glmQLFTest(fit,coef=ncol(design)),n=Inf)$table
    tt$gene_id<-rownames(tt);tt$treatment<-gr;tt$contrast<-paste0('(',gr,':',arms[2],'-',arms[1],') - (CK:',arms[2],'-',arms[1],')')
    result[[gr]]<-tt;aud[[gr]]<-data.frame(treatment=gr,status='tested; branch composition may contribute')
  }
  if(length(aud))save_table(do.call(rbind,aud),'34_branch_interaction_audit')
  if(length(result)){z<-do.call(rbind,result);z$FDR_global<-p.adjust(z$PValue,'BH');save_table(z,'35_treatment_by_branch_interaction_DE');return(z)}
  NULL
}
program_scores <- function(ct,sets){
  if(!length(sets))return(matrix(numeric(),ncol(ct),0))
  ans<-matrix(NA_real_,ncol(ct),length(sets),dimnames=list(colnames(ct),names(sets)))
  for(m in names(sets)){
    genes<-intersect(sets[[m]],rownames(ct));if(!length(genes))next
    ss<-numeric(ncol(ct));n<-0L
    for(g in split(genes,ceiling(seq_along(genes)/30))){
      e<-as.matrix(expression(ct,g));mu<-rowMeans(e);sdv<-apply(e,1,sd);ok<-is.finite(sdv)&sdv>1e-8
      if(any(ok)){ss<-ss+colSums(sweep(sweep(e[ok,,drop=FALSE],1,mu[ok],'-'),1,sdv[ok],'/'));n<-n+sum(ok)}
    };if(n)ans[,m]<-ss/n
  };ans
}
branch_GO <- function(programs){
  if(!RUN_GO||!length(programs$background)||!any(lengths(programs$sigsets)>=2))return(NULL)
  pk<-c('AnnotationDbi',ORGDB_PACKAGE,'GO.db');if(any(!vapply(pk,requireNamespace,logical(1),quietly=TRUE))){note('GO_status','GO packages unavailable. No enrichment is invented.');return(NULL)}
  db<-get(ORGDB_PACKAGE,envir=asNamespace(ORGDB_PACKAGE));keys<-AnnotationDbi::keys(db,keytype='GID')
  mp<-map_genes(programs$background,keys);save_table(mp,'40_GO_mapping')
  if(mean(!is.na(mp$feature_id))<.5){note('GO_status','GO mapping <50%; check Glyma version/keytype.');return(NULL)}
  cols<-AnnotationDbi::columns(db);gc<-if('GOALL'%in%cols)'GOALL' else 'GO';oc<-if(gc=='GOALL')'ONTOLOGYALL' else 'ONTOLOGY'
  aa<-suppressMessages(AnnotationDbi::select(db,keys=unique(na.omit(mp$feature_id)),keytype='GID',columns=c(gc,oc)))
  aa<-unique(aa[!is.na(aa[[gc]])&!is.na(aa[[oc]])&aa[[oc]]=='BP',c('GID',gc)])
  bg<-unique(aa$GID);terms<-lapply(split(aa$GID,aa[[gc]]),unique);terms<-terms[lengths(terms)>=10&lengths(terms)<=500]
  if(!length(terms))return(NULL)
  nm<-suppressMessages(AnnotationDbi::select(GO.db::GO.db,keys=names(terms),keytype='GOID',columns='TERM'));labels<-setNames(nm$TERM,nm$GOID)
  all<-list()
  for(b in names(programs$sigsets)){
    genes<-intersect(na.omit(mp$feature_id[mp$gene_id%in%programs$sigsets[[b]]]),bg);if(length(genes)<2)next
    k<-vapply(terms,function(x)length(intersect(x,genes)),integer(1));K<-lengths(terms)
    all[[b]]<-data.frame(branch=b,ID=names(terms),Description=unname(labels[names(terms)]),Count=k,
                         GeneRatio=paste0(k,'/',length(genes)),pvalue=phyper(k-1,K,length(bg)-K,length(genes),lower.tail=FALSE))
  }
  if(!length(all))return(NULL)
  z<-do.call(rbind,all);z$FDR<-p.adjust(z$pvalue,'BH');save_table(z,'41_branch_GO_all_tests')
  sig<-z[z$FDR<GO_FDR&z$Count>=2,,drop=FALSE];save_table(sig,'42_branch_GO_significant');sig
}

draw_umap <- function(d,values,title,colors=NULL,edges=NULL,nodes=NULL,lim=NULL){
  if(!nrow(d)){message_panel(title,'No cells for this treatment.');return()}
  par(mar=c(3,3,3,1));plot(d$x,d$y,type='n',asp=1,xlab='Original UMAP 1',ylab='Original UMAP 2',main=title,cex.main=.95)
  if(!is.null(edges)&&nrow(edges)){
    aa<-match(edges$source,nodes$state);bb<-match(edges$target,nodes$state)
    segments(nodes$x[aa],nodes$y[aa],nodes$x[bb],nodes$y[bb],col='#777777',lwd=.5+3*edges$connectivity)
  }
  if(is.null(colors)){
    if(is.null(lim))lim<-range(values[is.finite(values)]);co<-if(all(is.finite(lim)))colors_num(values,lim,PT_COLORS) else rep('#CCCCCC',nrow(d))
    ii<-order(values,na.last=FALSE);points(d$x[ii],d$y[ii],col=co[ii],pch=16,cex=.22)
    if(all(is.finite(lim)))legend('bottomright',legend=c(paste('Early',signif(lim[1],3)),paste('Late',signif(lim[2],3))),
                                  col=PT_COLORS[c(1,length(PT_COLORS))],pch=16,bty='n',cex=.7)
  }else{
    points(d$x,d$y,col=colors[values],pch=16,cex=.22)
    legend('topleft',legend=intersect(names(colors),unique(values)),col=colors[intersect(names(colors),unique(values))],pch=16,bty='n',cex=.6)
  }
  if(!is.null(nodes)){
    points(nodes$x,nodes$y,pch=21,bg='white',cex=.8);text(nodes$x,nodes$y,nodes$state,pos=3,cex=.55)
    points(nodes$x[nodes$root],nodes$y[nodes$root],pch=23,bg='#F0E442',cex=1.5)
  }
}
draw_matrix <- function(m,title,rows=rownames(m),cols=colnames(m),lim=c(-2,2)){
  if(!is.matrix(m)||!length(m)){message_panel(title,'No eligible data.');return()}
  par(mar=c(4,6,3,1));plot(NA,xlim=c(.5,ncol(m)+.5),ylim=c(.5,nrow(m)+.5),xaxt='n',yaxt='n',xlab='',ylab='',main=title,cex.main=.9,xaxs='i',yaxs='i')
  for(i in seq_len(nrow(m)))rect(seq_len(ncol(m))-.5,nrow(m)-i+.5,seq_len(ncol(m))+.5,nrow(m)-i+1.5,
                                 col=colors_num(m[i,],lim),border=NA)
  if(length(rows)>0&&length(rows)<=35)axis(2,rev(seq_len(nrow(m))),rows,las=1,cex.axis=.65,tick=FALSE)
  if(length(cols)>0&&length(cols)<=40)axis(1,seq_len(ncol(m)),cols,las=2,cex.axis=.65,tick=FALSE)
  mtext(paste0('Colour limits: ',paste(lim,collapse=' to '),'; grey: missing'),side=1,line=2.7,cex=.6)
}
draw_model <- function(nodes,evidence,programs,go,title='QC-related candidate branch model'){
  bp<-evidence$reference_branchpoint[1];arms<-sort(unique(nodes$branch[grepl('^Branch_',nodes$branch)]))
  if(!nzchar(bp)||length(arms)<2){message_panel(title,'No supported articulation with two or more substantial downstream arms under the current criteria. The full PAGA network is retained; a Y-shaped model was not imposed.');return()}
  par(mar=c(3,2,4,1));plot(NA,xlim=c(0,1.8),ylim=c(.3,length(arms)+.9),axes=FALSE,xlab='',ylab='',main=title,cex.main=1)
  mid<-(length(arms)+1)/2;segments(.12,mid,.52,mid,lwd=4,col='#6A3D9A')
  points(.12,mid,pch=21,bg='#6A3D9A',cex=3);text(.12,mid+.35,'CK QC\n(candidate root)',cex=.8)
  points(.52,mid,pch=21,bg='black',cex=1.3);text(.52,mid+.27,bp,cex=.7)
  supported<-isTRUE(evidence$operational_two_arm_support[1])
  for(i in seq_along(arms)){
    b<-arms[i];yy<-length(arms)-i+1;co<-BRANCH_COLORS[b];if(is.na(co))co<-'#808080'
    segments(.52,mid,1,yy,lwd=3,col=co,lty=if(supported)1 else 2);points(1,yy,pch=21,bg=co,cex=2.5)
    typ<-unique(nodes$type[nodes$branch==b]);text(1.08,yy+.24,paste(b,paste(typ,collapse=' / '),sep=' | '),adj=0,cex=.72,col=co)
    gg<-if(!is.null(go))go[go$branch==b,,drop=FALSE] else NULL
    terms<-if(!is.null(gg)&&nrow(gg))head(gg$Description[order(gg$FDR)],3) else 'No FDR-supported GO terms'
    genes<-head(programs$sets[[b]],4);tx<-paste(c(terms,if(length(genes))paste(genes,collapse=', ')),collapse='\n')
    text(1.08,yy+.08,paste(strwrap(tx,53),collapse='\n'),adj=c(0,1),cex=.62)
  }
  mtext(if(supported)'Two-arm topology recovered under the stated checks; direction remains a biological hypothesis' else
    'Candidate topology only: stability checks incomplete/failed or more than two arms; dashed lines are exploratory',side=1,line=1,cex=.65)
}
draw_occupancy <- function(z){
  z<-z[grepl('^Branch_',z$branch),];if(!nrow(z)){message_panel('Branch occupancy','No candidate arms.');return()}
  arms<-unique(z$branch);par(mar=c(5,4,3,1));plot(NA,xlim=c(.5,length(GROUPS)+.5),ylim=c(0,1),xaxt='n',xlab='',ylab='Fraction among mapped downstream-arm cells',main='Per-sample branch occupancy',cex.main=.9)
  axis(1,seq_along(GROUPS),GROUPS,las=2,cex.axis=.8)
  for(i in seq_along(arms)){
    a<-z[z$branch==arms[i],];xx<-match(a$treatment,GROUPS)+(i-(length(arms)+1)/2)*.14
    points(xx,a$fraction_postbranch,pch=16,col=BRANCH_COLORS[arms[i]],cex=.8)
    mm<-aggregate(fraction_postbranch~treatment,a,mean);x<-match(mm$treatment,GROUPS)+(i-(length(arms)+1)/2)*.14
    segments(x-.06,mm$fraction_postbranch,x+.06,mm$fraction_postbranch,lwd=2,col=BRANCH_COLORS[arms[i]])
  };legend('topright',arms,col=BRANCH_COLORS[arms],pch=16,bty='n',cex=.65)
}
plot_edges_all_groups <- function(nodes,edges,changes,d){
  ref<-edges[is.finite(edges$connectivity)&edges$connectivity>=PAGA_THRESHOLD,,drop=FALSE]
  save_plot('PAGA_CK_and_all_treatments_fixed_layout',function(){
    par(mfrow=c(3,3))
    for(gr in GROUPS){
      par(mar=c(1,1,3,1));plot(nodes$x,nodes$y,type='n',asp=1,axes=FALSE,xlab='',ylab='',main=gr)
      if(gr=='CK'){
        ed<-ref;co<-rep('#555555',nrow(ed));ww<-.6+4*ed$connectivity
      }else{
        ed<-if(!is.null(changes))changes[changes$treatment==gr,,drop=FALSE] else NULL
        if(!is.null(ed)&&nrow(ed)){ed<-ed[ed$CK_mean>=PAGA_THRESHOLD|ed$treatment_mean>=PAGA_THRESHOLD,];co<-colors_num(ed$delta,c(-.5,.5),c('#2166AC','#EEEEEE','#B2182B'));ww<-.5+4*pmax(ed$CK_mean,ed$treatment_mean)}
      }
      if(!is.null(ed)&&nrow(ed)){
        ii<-match(ed$source,nodes$state);jj<-match(ed$target,nodes$state)
        segments(nodes$x[ii],nodes$y[ii],nodes$x[jj],nodes$y[jj],col=co,lwd=ww)
      }
      counts<-table(factor(d$state[d$treatment==gr],levels=nodes$state));sizes<-.6+2*sqrt(as.numeric(counts)/max(1,max(counts)))
      points(nodes$x,nodes$y,pch=21,bg=BRANCH_COLORS[nodes$branch],cex=sizes)
      text(nodes$x,nodes$y,nodes$state,pos=3,cex=.55)
      mtext(if(gr=='CK')'Full CK edges above threshold' else 'Red: stronger; blue: weaker | common-state balanced graphs',side=1,line=-1,cex=.5)
    }
  },15,13)
}

ck_program_heatmap <- function(ct,d,programs){
  allgenes<-unique(unlist(programs$sets,use.names=FALSE));arms<-names(programs$sets)
  if(!length(allgenes))return(matrix(numeric(),0,0))
  ck<-which(d$treatment=='CK'&is.finite(d$reference_PT));ex<-expression(ct,allgenes)[,ck,drop=FALSE];md<-d[ck,]
  mats<-list();binrows<-list()
  for(arm in arms){
    ix<-which(md$branch%in%c('Shared_trunk',arm));if(length(ix)<50)next
    v<-md$reference_PT[ix];lim<-range(v);if(diff(lim)<=0)next
    breaks<-seq(lim[1],lim[2],length.out=N_BINS+1);bins<-cut(v,breaks,include.lowest=TRUE,labels=FALSE)
    samples<-unique(md$sample[ix]);m<-matrix(NA_real_,length(allgenes),N_BINS,dimnames=list(allgenes,NULL))
    for(b in seq_len(N_BINS)){
      vals<-lapply(samples,function(s){jj<-ix[md$sample[ix]==s&bins==b];if(length(jj)>=MIN_BIN_SAMPLE_CELLS)Matrix::rowMeans(ex[,jj,drop=FALSE]) else NULL})
      if(length(samples)>=2&&all(lengths(vals)>0))m[,b]<-rowMeans(do.call(cbind,vals))
    }
    mats[[arm]]<-m;binrows[[arm]]<-data.frame(branch=arm,bin=seq_len(N_BINS),reference_PT=(head(breaks,-1)+tail(breaks,-1))/2)
  }
  if(!length(mats))return(matrix(numeric(),0,0))
  m<-do.call(cbind,mats);colnames(m)<-unlist(lapply(names(mats),function(a)paste0(a,'_',seq_len(N_BINS))))
  z<-m;for(i in seq_len(nrow(m))){sdv<-sd(m[i,],na.rm=TRUE);z[i,]<-if(is.finite(sdv)&&sdv>1e-8)(m[i,]-mean(m[i,],na.rm=TRUE))/sdv else NA_real_}
  save_table(data.frame(gene_id=rownames(m),m,check.names=FALSE),'50_CK_program_expression_bins')
  save_table(do.call(rbind,binrows),'51_CK_program_bin_coordinates');save_table(data.frame(gene_id=rownames(z),z,check.names=FALSE),'52_CK_program_shared_Z')
  z
}
draw_program_heatmap <- function(m){
  if(!length(m)){message_panel('CK branch-associated gene dynamics','No eligible gene program or replicate-supported bins. See DE audit; no arbitrary modules inserted.');return()}
  cols<-rep('',ncol(m));for(i in seq(1,ncol(m),by=N_BINS))cols[i]<-sub('_[0-9]+$','',colnames(m)[i])
  draw_matrix(m,'CK branch programs | shared trunk + each arm',rows=character(),cols=character())
  starts<-seq(1,ncol(m),by=N_BINS);axis(1,at=starts+(N_BINS-1)/2,
                                        labels=sub('_[0-9]+$','',colnames(m)[starts]),cex.axis=.7,tick=FALSE)
  if(ncol(m)>N_BINS)abline(v=seq(N_BINS+.5,ncol(m)-.5,by=N_BINS),col='black',lwd=2)
}
program_dot <- function(ct,d,sets){
  genes<-unique(unlist(lapply(sets,head,4),use.names=FALSE));if(!length(genes))return(NULL)
  ck<-which(d$treatment=='CK');e<-expression(ct,genes)[,ck,drop=FALSE];md<-d[ck,];rows<-list()
  for(s in unique(md$state))for(g in genes){
    ss<-split(which(md$state==s),md$sample[md$state==s]);ss<-ss[lengths(ss)>=10];if(length(ss)<2)next
    rows[[paste(s,g)]]<-data.frame(state=s,gene_id=g,mean=mean(vapply(ss,function(j)mean(as.numeric(e[g,j])),numeric(1))),
                                   pct=mean(vapply(ss,function(j)mean(as.numeric(e[g,j])>0),numeric(1)))*100)
  };z<-if(length(rows))do.call(rbind,rows) else NULL
  if(!is.null(z)){
    z$scaled<-ave(z$mean,z$gene_id,FUN=function(v){sdv<-sd(v);if(is.finite(sdv)&&sdv>1e-8)(v-mean(v))/sdv else rep(0,length(v))})
    save_table(z,'53_CK_state_marker_dotplot')
  };z
}
draw_program_dot <- function(z){
  if(is.null(z)){message_panel('Branch-associated gene expression','Gene programs require sample-level CK branch comparison.');return()}
  gs<-unique(z$gene_id);ss<-sort(unique(z$state));par(mar=c(5,9,3,1))
  plot(NA,xlim=c(.5,length(ss)+.5),ylim=c(.5,length(gs)+.5),xaxt='n',yaxt='n',xlab='',ylab='',main='CK branch-associated genes by reference state',cex.main=.9)
  points(match(z$state,ss),match(z$gene_id,gs),pch=21,bg=colors_num(z$scaled,c(-2,2)),cex=2.2*sqrt(z$pct/100),lwd=.4)
  axis(1,seq_along(ss),ss,las=2,cex.axis=.6,tick=FALSE);axis(2,seq_along(gs),gs,las=1,cex.axis=.6,tick=FALSE)
  mtext('Point area: % detected; colour: gene-wise Z across CK states',side=3,line=.2,cex=.55)
}
early_similarity <- function(){
  p<-read.csv(file.path(OUT_DIR,'CK_early_cell_PCs.csv'),check.names=FALSE);md<-read.csv(file.path(OUT_DIR,'CK_early_cells.csv'))
  if(nrow(p)<4)return(list(matrix=matrix(numeric(),0,0),metadata=md))
  raw<-cor(t(as.matrix(p[,-1,drop=FALSE])),use='pairwise.complete.obs')
  z<-t(apply(raw,1,function(v){s<-sd(v,na.rm=TRUE);if(is.finite(s)&&s>0)(v-mean(v,na.rm=TRUE))/s else rep(0,length(v))}))
  saveRDS(list(cell=p$cell,Pearson=raw,row_Z=z,metadata=md),file.path(OUT_DIR,'54_CK_early_cell_PC_similarity.rds'))
  list(matrix=z,metadata=md)
}
draw_early_similarity <- function(sim){
  draw_matrix(sim$matrix,'Early CK cells | row-Z of PC-profile Pearson r',rows=character(),cols=character())
  if(length(sim$matrix)){
    br<-sim$metadata$branch;nn<-cumsum(rle(br)$lengths);nn<-head(nn,-1)
    if(length(nn)){abline(v=nn+.5,col='black',lwd=1);abline(h=nrow(sim$matrix)-nn+.5,col='black',lwd=1)}
  }
}
program_response_curves <- function(scores,d){
  if(!ncol(scores))return()
  branches<-sort(unique(d$branch[grepl('^Branch_',d$branch)]));out<-list()
  for(b in branches){
    ck<-which(d$treatment=='CK'&d$branch%in%c('Shared_trunk',b)&is.finite(d$reference_PT));if(length(ck)<50)next
    lim<-range(d$reference_PT[ck]);if(diff(lim)<=0)next
    breaks<-seq(lim[1],lim[2],length.out=N_BINS+1);centers<-(head(breaks,-1)+tail(breaks,-1))/2
    selected<-which(d$branch%in%c('Shared_trunk',b)&is.finite(d$reference_PT));bins<-cut(d$reference_PT[selected],breaks,include.lowest=TRUE,labels=FALSE)
    for(m in colnames(scores)){
      rr<-list()
      for(s in unique(d$sample[selected]))for(j in seq_len(N_BINS)){
        ix<-selected[d$sample[selected]==s&!is.na(bins)&bins==j];valid<-length(ix)>=MIN_BIN_SAMPLE_CELLS
        rr[[paste(s,j)]]<-data.frame(branch=b,program=m,sample=s,treatment=d$treatment[match(s,d$sample)],
                                     bin=j,reference_PT=centers[j],n=length(ix),score=if(valid)mean(scores[ix,m]) else NA_real_)
      };z<-do.call(rbind,rr);out[[paste(b,m)]]<-z
      save_plot(paste0('E_',b,'_',safe_name(m),'_response_curves'),local({zz<-z;mm<-m;bb<-b;function(){
        par(mfrow=c(3,3),mar=c(3,3,3,1));yr<-range(zz$score,na.rm=TRUE);if(!all(is.finite(yr)))yr<-c(-1,1)
        for(gr in GROUPS){
          plot(NA,xlim=lim,ylim=yr,xlab='CK reference position',ylab='Program expression score',main=paste(gr,bb,mm),cex.main=.75)
          for(s in unique(zz$sample[zz$treatment==gr])){a<-zz[zz$sample==s,];lines(a$reference_PT,a$score,col=adjustcolor(GROUP_COLORS[gr],.6),lty=3)}
          for(g in unique(c('CK',gr))){
            a<-zz[zz$treatment==g,];vals<-vapply(seq_len(N_BINS),function(j){v<-a$score[a$bin==j];if(length(v)>=2&&all(is.finite(v)))mean(v) else NA_real_},numeric(1))
            lines(centers,vals,col=GROUP_COLORS[g],lwd=if(g==gr)2 else 1,lty=if(g==gr)1 else 2)
          }
        }
      }}),13,10)
    }
  };if(length(out))save_table(do.call(rbind,out),'55_program_response_bins_per_sample')
}
PAGA_PY_PKGS <- c('scanpy','anndata','scipy','numpy','pandas','sklearn',
                  'igraph','leidenalg','h5py','pynndescent')
find_paga_python <- function(){
  cand <- c(Sys.getenv('PAGA_PYTHON'),
            get0('PYTHON',envir=globalenv(),ifnot=''),
            if(nzchar(Sys.getenv('CONDA_PREFIX')))file.path(Sys.getenv('CONDA_PREFIX'),'bin','python'),
            '/home/liuxy/micromamba/envs/paga/bin/python',
            '/home/liuxy/miniconda3/envs/paga/bin/python',
            '/opt/conda/envs/paga/bin/python',
            unname(Sys.which(c('python3','python')))
  )
  cand <- unique(cand[nzchar(cand)])
  code <- paste0('import ',paste(PAGA_PY_PKGS,collapse=', '),
                 '; import sys; from importlib.metadata import version; ',
                 'print(sys.version.split()[0], version("scanpy"), version("anndata"), ',
                 'version("numpy"), version("pandas"), sep=" | ")')
  for(p in cand){
    exe <- if(grepl('/',p))p else unname(Sys.which(p))
    if(!nzchar(exe)||!file.exists(exe))next
    out <- suppressWarnings(tryCatch(system2(exe,c('-c',shQuote(code)),stdout=TRUE,stderr=TRUE),
                                     error=function(e)structure('',class='try-error')))
    ver <- if(inherits(out,'try-error'))'' else trimws(utils::tail(out[nzchar(out)],1))
    if(inherits(out,'try-error')||!is.null(attr(out,'status'))||!grepl('\\|',ver)){
      cat('skip:',exe,'  ->',if(nzchar(ver))ver else 'cannot import all packages','\n');next
    }
    cat('PAGA Python ->',exe,'\n  python | scanpy | anndata | numpy | pandas =',ver,'\n')
    return(exe)
  }
  stop('没有找到装齐 scanpy/leidenalg/igraph/h5py/pynndescent 的Python。',
       '先用 conda/micromamba 建环境(见脚本头部说明)，再把 PYTHON 设为该环境的绝对路径。',
       '系统python是',Sys.which('python3'),'时通常版本过低。')
}

main_v10 <- function(){
  set.seed(SEED);dir.create(OUT_DIR,recursive=TRUE,showWarnings=FALSE)
  plog('main_v10 fix v10.2 (UMAP resume) start | R ',R.version.string)
  prepare_existing_export_v10()
  PYTHON <<- find_paga_python()
  plog('Using PYTHON = ',PYTHON)
  plog('stage 1/2: export Seurat counts and metadata')
  z<-export_input();validate_cached_export(z)
  save_baseline_umap(z$metadata)
  run_paga(z,'reference')
  render_reference_v10()
  run_paga(z,'full')
  plog('stage 1/2 done: Python PAGA and CK mapping finished')
  d<-read.csv(file.path(OUT_DIR,'cells_CK_reference_mapping.csv'),stringsAsFactors=FALSE)
  if(!setequal(d$cell,colnames(z$counts)))stop('Mapping/count barcode mismatch.')
  d<-d[match(colnames(z$counts),d$cell),];ct<-z$counts
  stopifnot(max(abs(d$x-z$metadata$x))<1e-8,max(abs(d$y-z$metadata$y))<1e-8)
  save_table(data.frame(n_cells=nrow(d),max_UMAP_coordinate_change=max(abs(c(d$x-z$metadata$x,d$y-z$metadata$y)))),'02_UMAP_preservation')
  nodes<-read.csv(file.path(OUT_DIR,'CK_nodes.csv'));edges<-read.csv(file.path(OUT_DIR,'CK_edges_all.csv'))
  nodes$root<-tolower(as.character(nodes$root))=='true';nodes$branchpoint<-tolower(as.character(nodes$branchpoint))=='true'
  evidence<-read.csv(file.path(OUT_DIR,'CK_branch_evidence.csv'))
  for(nm in c('operational_two_arm_support','ck_replicates_agree','root_biologically_validated'))evidence[[nm]]<-tolower(as.character(evidence[[nm]]))=='true'
  allarms<-sort(unique(d$branch[grepl('^Branch_',d$branch)]));extra<-setdiff(allarms,names(BRANCH_COLORS))
  if(length(extra))BRANCH_COLORS<<-c(BRANCH_COLORS,setNames(hcl.colors(length(extra),'Dark 3'),extra))
  occ<-occupancy_stats(d);changes<-balanced_edge_stats(d)
  render_reference_v10()
  plot_treatment_network_v10(nodes,edges,changes,d)
  plog('Reference and treatment network figures saved; starting DE/GO')
  mapping<-aggregate(cbind(mapped=as.integer(d$state!='Unmapped'),cells=rep(1,nrow(d)))~sample+treatment,d,sum)
  mapping$mapped_fraction<-mapping$mapped/mapping$cells;save_table(mapping,'03_mapping_quality_by_sample')
  programs<-optional_stage('branch_programs',branch_programs(ct,d),list(sets=list(),sigsets=list(),background=character()))
  optional_stage('treatment_branch_DE',treatment_branch_DE(ct,d),NULL)
  optional_stage('branch_interaction_DE',branch_interaction_DE(ct,d),NULL)
  go<-optional_stage('branch_GO',branch_GO(programs),NULL)
  sets<-programs$sets
  if(nzchar(WGCNA_MEMBERSHIP_FILE)){
    w<-need(read_table(WGCNA_MEMBERSHIP_FILE),c('gene_id','module'));mp<-map_genes(as.character(w$gene_id),rownames(ct));mp$module<-w$module;save_table(mp,'56_WGCNA_mapping')
    ws<-lapply(split(mp$feature_id,mp$module),function(a)unique(na.omit(a)));ws<-ws[lengths(ws)>=5]
    names(ws)<-paste0('WGCNA_',names(ws));sets<-c(sets,ws)
  }
  scores<-program_scores(ct,sets);saveRDS(list(cell=d$cell,scores=scores),file.path(OUT_DIR,'program_scores_v9.rds'))
  heat<-ck_program_heatmap(ct,d,programs);dots<-program_dot(ct,d,programs$sets);sim<-early_similarity()
  ck<-d[d$treatment=='CK',];ee<-edges[is.finite(edges$connectivity)&edges$connectivity>=PAGA_THRESHOLD,]
  save_plot('A_CK_PAGA_originalUMAP',function(){par(mfrow=c(1,2));
    draw_umap(ck,ck$branch,'CK reference states and candidate arms',BRANCH_COLORS,ee,nodes)
    draw_umap(ck,ck$reference_PT,'CK diffusion pseudotime',edges=ee,nodes=nodes)
  },14,7)
  save_plot('B_CK_state_program_dotplot',function()draw_program_dot(dots),12,8)
  save_plot('C_CK_early_cell_similarity',function()draw_early_similarity(sim),10,9)
  save_plot('D_CK_branch_program_heatmap',function()draw_program_heatmap(heat),11,9)
  save_plot('F_QC_branch_function_model',function()draw_model(nodes,evidence,programs,go),16,9)
  save_plot('G_branch_occupancy_nine_groups',function()draw_occupancy(occ),12,7)
  plot_edges_all_groups(nodes,edges,changes,d)
  save_plot('A_all_treatments_CK_branch_projection',function(){par(mfrow=c(3,3));
    for(gr in GROUPS){a<-d[d$treatment==gr,];draw_umap(a,a$branch,paste(gr,'| CK reference arms'),BRANCH_COLORS,ee,nodes)}
  },15,13)
  for(m in colnames(scores)){
    value<-scores[,m];lim<-if(any(is.finite(value)))as.numeric(quantile(value[is.finite(value)],c(.01,.99))) else c(-1,1)
    save_plot(paste0('E_',safe_name(m),'_UMAP'),local({v<-value;mm<-m;ll<-lim;function(){par(mfrow=c(3,3));
      for(gr in GROUPS){ii<-which(d$treatment==gr);draw_umap(d[ii,],v[ii],paste(gr,mm),lim=ll)}
    }}),15,13)
  }
  program_response_curves(scores,d)
  save_plot('Figure_CK_QC_branch_evidence_v9',function(){
    layout(rbind(c(1,2,3),c(4,5,6),c(7,7,8)),heights=c(1,1,1.05))
    draw_umap(ck,ck$branch,'a | CK PAGA and candidate arms',BRANCH_COLORS,ee,nodes)
    draw_umap(ck,ck$reference_PT,'a | CK DPT (candidate QC root)',edges=ee,nodes=nodes)
    draw_program_dot(dots);draw_early_similarity(sim);draw_program_heatmap(heat)
    if(ncol(scores))draw_umap(ck,scores[d$treatment=='CK',1],paste('e |',colnames(scores)[1],'program expression')) else message_panel('e | Program expression','No eligible sample-supported program. See DE_status.')
    draw_model(nodes,evidence,programs,go,'f | CK QC-related branch/function model');draw_occupancy(occ)
  },20,18)
  capture.output(sessionInfo(),file=file.path(OUT_DIR,'R_sessionInfo.txt'))
  conclusion<-if(evidence$operational_two_arm_support[1])
    'CK在本次固定状态划分和预设检查下支持两条主要候选下游臂；尚不是QC细胞直接分化到两类组织的谱系实验证据。' else
      paste0('当前预设标准未给出稳定两分支结论。主阈值候选臂数=',evidence$n_candidate_arms[1],
             '；请查看重复、阈值/k值与重采样结果。不能为了获得Y形而修改图形冒充结果。')
  note('RESULT_interpretation',c(conclusion,
                                 '处理影响分别检查：21分支占比样本级检验；22平衡后的连接变化；33分支内pseudobulk差异；55表达曲线。',
                                 '每组只有2重复时，两组间双侧精确置换检验的最小p值为1/3；不能把单细胞或重采样次数当作生物学重复增加显著性。',
                                 '所有分支比较均为转录状态关联。处理细胞的reference_PT是映射到CK的位置，不是独立估计的处理DPT。',
                                 '映射外细胞、状态消失和共同状态覆盖率必须一起看；分支内差异表达可能混入状态/细胞类型组成变化。',
                                 'PAGA权重是连接性分数，不是方向概率或p值。候选图不做强制MST，不预设两个终点。',
                                 '若批次与处理完全混淆，无法单凭此数据归因激素/胁迫因果作用；DE相关比较会在不可估计设计时跳过。'))
  plog('Completed: ',OUT_DIR)
  writeLines('All required topology stages completed; inspect optional_stage_errors_v10.log for optional analyses.',file.path(OUT_DIR,'RUN_COMPLETED_v10.txt'))
  invisible(evidence)
}
SHOW_HYPOTHESIS_ARROWS <- TRUE
MIN_DPT_ARROW_DIFFERENCE <- .02
STATE_LABEL_FILE <- ''

prepare_existing_export_v10 <- function(){
  if(REBUILD_EXPORT)return(invisible(NULL))
  required<-c('counts_metadata_v9.rds','export_signature.rds','input_counts.mtx','input_metadata.csv','input_genes.csv')
  present<-file.exists(file.path(OUT_DIR,required))
  if(all(present))return(invisible(NULL))
  if(!dir.exists(INPUT_CACHE_DIR)||!all(file.exists(file.path(INPUT_CACHE_DIR,required)))){
    if(any(present))stop('Partial input export in OUT_DIR and no complete INPUT_CACHE_DIR; inspect before setting REBUILD_EXPORT=TRUE.')
    plog('No complete old export found; Seurat export will run.');return(invisible(NULL))
  }
  if(any(present))stop('Partial input export in current directory. Use the complete old result directory via PAGA_OUT_DIR, or resolve these files first.')
  for(nm in required){
    src<-normalizePath(file.path(INPUT_CACHE_DIR,nm),winslash='/',mustWork=TRUE);dst<-file.path(OUT_DIR,nm)
    if(nm%in%c('input_counts.mtx','counts_metadata_v9.rds')){
      ok<-suppressWarnings(file.symlink(src,dst))
      if(!ok)ok<-suppressWarnings(file.link(src,dst))
      if(!ok)stop('Could not link large cache ',nm,'. Set PAGA_OUT_DIR to INPUT_CACHE_DIR to run in place; no large copy was made.')
    }else if(!file.copy(src,dst,overwrite=FALSE))stop('Could not copy small input file: ',nm)
  }
  for(nm in c('00_treatment_audit.csv','01_biological_samples.csv'))
    if(file.exists(file.path(INPUT_CACHE_DIR,nm))&&!file.exists(file.path(OUT_DIR,nm)))
      file.copy(file.path(INPUT_CACHE_DIR,nm),file.path(OUT_DIR,nm),overwrite=FALSE)
  plog('Linked existing large exports from ',INPUT_CACHE_DIR,'; figures/output -> ',OUT_DIR)
}

validate_cached_export <- function(z){
  need(z$metadata,c('cell','treatment','sample','cell_type','x','y'))
  if(!inherits(z$counts,'sparseMatrix'))stop('Cached counts are not sparse.')
  if(!identical(as.character(z$metadata$cell),colnames(z$counts)))stop('Cached barcode ordering mismatch.')
  if(anyDuplicated(z$metadata$cell)||anyNA(z$metadata$treatment)||any(!is.finite(c(z$metadata$x,z$metadata$y))))stop('Invalid cached metadata.')
  if(REQUIRE_ALL_GROUPS&&length(setdiff(GROUPS,z$metadata$treatment)))stop('Cached data lack required treatments.')
  for(f in c('input_counts.mtx','input_metadata.csv','input_genes.csv'))
    if(!file.exists(file.path(OUT_DIR,f)))stop('Missing exported file: ',f,'; set REBUILD_EXPORT=TRUE to repair export.')
  md<-read.csv(file.path(OUT_DIR,'input_metadata.csv'),stringsAsFactors=FALSE)
  if(!identical(as.character(md$cell),as.character(z$metadata$cell))||
     !identical(as.character(md$treatment),as.character(z$metadata$treatment))||
     !identical(as.character(md$sample),as.character(z$metadata$sample)))stop('CSV and RDS metadata mismatch; re-export explicitly.')
  gg<-read.csv(file.path(OUT_DIR,'input_genes.csv'),stringsAsFactors=FALSE)
  if(!identical(as.character(gg$gene_id),rownames(z$counts)))stop('CSV and RDS gene order mismatch.')
  con<-file(file.path(OUT_DIR,'input_counts.mtx'),'rt');on.exit(close(con))
  repeat{ln<-readLines(con,n=1,warn=FALSE);if(!length(ln))stop('Empty Matrix Market file.');if(!startsWith(ln,'%'))break}
  dims<-scan(text=ln,quiet=TRUE)
  if(length(dims)!=3||!identical(as.numeric(dims[1:2]),as.numeric(dim(z$counts))))stop('Matrix Market dimensions mismatch.')
  samples<-unique(z$metadata[,c('treatment','sample')]);nr<-table(samples$treatment)
  if(any(nr<MIN_REPLICATES))stop('Fewer than required biological samples in cached metadata.')
  tab<-as.data.frame(table(treatment=z$metadata$treatment));tab<-tab[tab$Freq>0,]
  save_table(tab,'00_cached_treatment_counts_v10')
  old<-file.path(OUT_DIR,'00_treatment_audit.csv')
  if(file.exists(old)){a<-read.csv(old);if('Freq'%in%names(a))save_table(a[a$Freq>0,,drop=FALSE],'00_treatment_audit_nonzero_v10')}
  plog('Validated cached gene/cell ordering and treatment labels: ',nrow(z$counts),' genes x ',ncol(z$counts),' cells')
  invisible(TRUE)
}
optional_stage <- function(name,expr,fallback){
  plog('Start optional stage: ',name)
  tryCatch(force(expr),error=function(e){
    msg<-paste(format(Sys.time()),name,conditionMessage(e),sep=' | ')
    cat(msg,'\n',file=file.path(OUT_DIR,'optional_stage_errors_v10.log'),append=TRUE)
    plog('Optional stage failed (not a negative biological result): ',msg)
    fallback
  })
}
type_palette <- function(types){
  known<-c(QC='#D73027',Meristem='#8DAA35','Dividing Cell'='#E69F00',Endodermis='#CC2C91',
           Cortex='#8073AC',Epidermis='#3197C6','Root Cap'='#F46D43',Procambium='#289E8E',
           Pericycle='#70AD47',Phloem='#B8860B',Xylem='#247F4E')
  types<-sort(unique(as.character(types)));extra<-setdiff(types,names(known))
  if(length(extra))known<-c(known,setNames(hcl.colors(length(extra),'Dark 3'),extra))
  known[types]
}
decorate_nodes_v10 <- function(nodes,d){
  nodes$root<-tolower(as.character(nodes$root))%in%c('true','1')
  nodes$branchpoint<-tolower(as.character(nodes$branchpoint))%in%c('true','1')
  base<-type_palette(nodes$type);nodes$color<-unname(base[nodes$type])
  for(ty in unique(nodes$type)){
    ix<-which(nodes$type==ty);ix<-ix[order(nodes$reference_PT[ix],nodes$state[ix],na.last=TRUE)]
    pal<-colorRampPalette(c(base[ty],'#E8E8E8'))(max(3,length(ix)+1))
    nodes$color[ix]<-pal[seq_along(ix)]
  }
  nodes$label<-paste(nodes$type,nodes$state,sep=' | ')
  if(nzchar(STATE_LABEL_FILE)){
    lab<-need(read_table(STATE_LABEL_FILE),c('state','label'));if(anyDuplicated(lab$state))stop('Duplicate state labels.')
    ix<-match(nodes$state,lab$state);ok<-!is.na(ix)&nzchar(lab$label[ix]);nodes$label[ok]<-lab$label[ix[ok]]
  }
  nodes$n_cells<-as.numeric(table(factor(d$state[d$treatment=='CK'],levels=nodes$state)))
  nodes
}
edge_geometry_v10 <- function(nodes,edges){
  edges<-edges[is.finite(edges$connectivity)&edges$connectivity>=PAGA_THRESHOLD,,drop=FALSE]
  a<-match(edges$source,nodes$state);b<-match(edges$target,nodes$state)
  good<-!is.na(a)&!is.na(b);edges<-edges[good,,drop=FALSE];a<-a[good];b<-b[good]
  pa<-nodes$reference_PT[a];pb<-nodes$reference_PT[b]
  oriented<-is.finite(pa)&is.finite(pb)&abs(pa-pb)>=MIN_DPT_ARROW_DIFFERENCE
  swap<-oriented&pa>pb;aa<-a; a[swap]<-b[swap];b[swap]<-aa[swap]
  data.frame(source=nodes$state[a],target=nodes$state[b],connectivity=edges$connectivity,
             x0=nodes$x[a],y0=nodes$y[a],x1=nodes$x[b],y1=nodes$y[b],
             DPT_source=nodes$reference_PT[a],DPT_target=nodes$reference_PT[b],
             DPT_orientable=oriented,within_type=nodes$type[a]==nodes$type[b],
             color=nodes$color[b],stringsAsFactors=FALSE)
}
draw_curve_v10 <- function(x0,y0,x1,y1,col='#555555',lwd=1,arrow=FALSE,lty=1,bend=.12){
  dx<-x1-x0;dy<-y1-y0
  if(!all(is.finite(c(x0,y0,x1,y1)))||dx^2+dy^2<1e-12)return(invisible(NULL))
  t<-seq(.06,.94,length.out=60);cx<-(x0+x1)/2-bend*dy;cy<-(y0+y1)/2+bend*dx
  xx<-(1-t)^2*x0+2*t*(1-t)*cx+t^2*x1;yy<-(1-t)^2*y0+2*t*(1-t)*cy+t^2*y1
  lines(xx,yy,col=col,lwd=lwd,lty=lty)
  if(arrow)arrows(xx[55],yy[55],xx[60],yy[60],length=.065,angle=24,col=col,lwd=lwd)
}
label_nodes_v10 <- function(nodes,cex=.62){
  if(!nrow(nodes))return(invisible(NULL))
  u<-par('usr');xr<-diff(u[1:2]);yr<-diff(u[3:4]);lx<-nodes$x;ly<-nodes$y+.024*yr
  for(step in seq_len(100))for(i in seq_len(nrow(nodes)))if(i<nrow(nodes))for(j in (i+1):nrow(nodes)){
    dx<-(lx[i]-lx[j])/xr;dy<-(ly[i]-ly[j])/yr
    if(abs(dx)<.13&&abs(dy)<.026){
      sy<-if(dy==0)1 else sign(dy);ly[i]<-ly[i]+sy*.002*yr;ly[j]<-ly[j]-sy*.002*yr
    }
  }
  ly<-pmin(u[4]-.02*yr,pmax(u[3]+.02*yr,ly))
  segments(nodes$x,nodes$y,lx,ly,col='#AAAAAA',lwd=.4)
  text(lx,ly,nodes$label,cex=cex,col='#222222',xpd=NA)
}
draw_state_network_v10 <- function(nodes,edges,title,focus=NULL){
  eg<-edge_geometry_v10(nodes,edges)
  if(!is.null(focus)){
    own<-nodes$state[nodes$type==focus]
    ekeep<-eg$source%in%own|eg$target%in%own;eg<-eg[ekeep,,drop=FALSE]
    show<-unique(c(own,eg$source,eg$target));nodes<-nodes[nodes$state%in%show,,drop=FALSE]
  }
  if(!nrow(nodes)){message_panel(title,'No eligible states.');return()}
  par(mar=c(2.8,1,3,1),cex=1);xr<-range(nodes$x);yr<-range(nodes$y)
  rx<-max(1,diff(xr));ry<-max(1,diff(yr))
  plot(NA,xlim=xr+c(-.18,.18)*rx,ylim=yr+c(-.18,.22)*ry,asp=1,axes=FALSE,xlab='',ylab='',main=title,cex.main=.95)
  for(i in seq_len(nrow(eg))){
    e<-eg[i,];co<-if(e$within_type)e$color else '#999999'
    draw_curve_v10(e$x0,e$y0,e$x1,e$y1,col=co,lwd=.55+2*e$connectivity,
                   arrow=SHOW_HYPOTHESIS_ARROWS&&e$DPT_orientable,lty=if(e$within_type)1 else 3)
  }
  co<-nodes$color;if(!is.null(focus))co[nodes$type!=focus]<-'#E3E3E3'
  points(nodes$x,nodes$y,pch=21,bg=co,col='white',cex=2.1)
  points(nodes$x[nodes$root],nodes$y[nodes$root],pch=23,bg='#FFD92F',col='black',cex=2)
  label_nodes_v10(nodes,if(is.null(focus)).8 else .85)
  mtext(if(SHOW_HYPOTHESIS_ARROWS)'PAGA edges | arrows: unvalidated DPT hypothesis' else 'Within-tissue PAGA connectivity | direction not assigned',1,line=.6,cex=.7)
}
draw_state_umap_v10 <- function(d,nodes,edges,title,background=NULL){
  par(cex=1,mar=c(6,5,3,1),mgp=c(2.3,.65,0));ref<-if(is.null(background))d else background
  plot(ref$x,ref$y,type='n',asp=1,xlab='Original UMAP 1',ylab='Original UMAP 2',main=title,cex.main=.95)
  if(!is.null(background))points(background$x,background$y,pch=16,cex=.1,col='#E8E8E8')
  co<-nodes$color[match(d$state,nodes$state)];co[is.na(co)]<-'#D3D3D3'
  points(d$x,d$y,pch=16,cex=.22,col=co)
  eg<-edge_geometry_v10(nodes,edges)
  for(i in seq_len(nrow(eg))){e<-eg[i,]
  draw_curve_v10(e$x0,e$y0,e$x1,e$y1,col=if(e$within_type)e$color else '#777777',
                 lwd=.65+1.1*e$connectivity,arrow=SHOW_HYPOTHESIS_ARROWS&&e$DPT_orientable,lty=if(e$within_type)1 else 3)
  }
  points(nodes$x[nodes$root],nodes$y[nodes$root],pch=23,bg='#FFD92F',cex=1.4)
  text(nodes$x,nodes$y,nodes$state,pos=3,cex=.7)
  mtext('Curves: PAGA connections between state centers; not measured cell paths',1,line=4.3,cex=.68)
}
save_baseline_umap <- function(d){
  plog('Save baseline UMAP before Python analysis (not a trajectory result)')
  pal<-type_palette(d$cell_type)
  save_plot('00_INPUT_UMAP_celltypes_v10',function(){
    par(mfrow=c(1,2));draw_umap(d[d$treatment=='CK',],d$cell_type[d$treatment=='CK'],
                                'CK input annotations | trajectory not yet computed',pal)
    draw_umap(d,d$cell_type,'All treatments | original input UMAP',pal)
  },16,8)
}
render_reference_v10 <- function(){
  d<-read.csv(file.path(OUT_DIR,'cells_CK_reference_mapping.csv'),stringsAsFactors=FALSE)
  nodes<-decorate_nodes_v10(read.csv(file.path(OUT_DIR,'CK_nodes.csv')),d)
  edges<-read.csv(file.path(OUT_DIR,'CK_edges_all.csv'));ck<-d[d$treatment=='CK',]
  save_table(nodes,'60_CK_state_labels_colors_v10');save_table(edge_geometry_v10(nodes,edges),'61_displayed_PAGA_edges_DPT_orientation_v10')
  save_plot('P2_CK_state_network_v10',function()draw_state_network_v10(nodes,edges,'CK | cell-state connectivity and candidate directions'),13,11)
  save_plot('P3_CK_originalUMAP_candidate_paths_v10',function()draw_state_umap_v10(ck,nodes,edges,'CK | original UMAP and PAGA state connections'),12,10)
  save_plot('P2_P3_CK_network_and_UMAP_v10',function(){
    par(mfrow=c(1,2));draw_state_umap_v10(ck,nodes,edges,'a | CK original UMAP and candidate paths')
    draw_state_network_v10(nodes,edges,'b | CK state-node network')
  },21,10)
  types<-unique(nodes$type)
  save_plot('P2_CK_local_tissue_state_networks_v10',function(){
    par(mfrow=c(ceiling(length(types)/3),3))
    for(ty in types)draw_state_network_v10(nodes,edges,ty,ty)
  },16,max(6,4.4*ceiling(length(types)/3)))
  save_plot('P3_nine_groups_CK_reference_projection_v10',function(){
    par(mfrow=c(3,3));for(gr in GROUPS)draw_state_umap_v10(d[d$treatment==gr,],nodes,edges,
                                                           paste(gr,'| projected CK reference network'),d)
  },18,16)
  save_plot('P2_P3_state_legend_v10',function(){
    par(mar=c(1,1,2,1));plot.new();title('CK graph states and majority annotations')
    legend('center',nodes$label,col=nodes$color,pch=16,bty='n',cex=.85,ncol=max(1,ceiling(nrow(nodes)/20)))
  },max(8,4*ceiling(nrow(nodes)/20)),10)
  plog('P2/P3 reference figures saved to: ',OUT_DIR)
  invisible(list(nodes=nodes,edges=edges,data=d))
}
plot_treatment_network_v10 <- function(nodes,edges,changes,d){
  nodes<-decorate_nodes_v10(nodes,d)
  save_plot('P2_nine_groups_balanced_connectivity_changes_v10',function(){
    par(mfrow=c(3,3))
    for(gr in GROUPS){
      if(gr=='CK'){draw_state_network_v10(nodes,edges,'CK | reference');next}
      z<-if(!is.null(changes))changes[changes$treatment==gr,,drop=FALSE] else NULL
      par(mar=c(2,1,3,1));plot(nodes$x,nodes$y,type='n',asp=1,axes=FALSE,xlab='',ylab='',main=gr)
      if(!is.null(z)&&nrow(z)){
        z<-z[is.finite(z$CK_mean)&is.finite(z$treatment_mean)&pmax(z$CK_mean,z$treatment_mean)>=PAGA_THRESHOLD,,drop=FALSE]
        for(i in seq_len(nrow(z))){a<-match(z$source[i],nodes$state);b<-match(z$target[i],nodes$state)
        draw_curve_v10(nodes$x[a],nodes$y[a],nodes$x[b],nodes$y[b],
                       col=colors_num(z$delta[i],c(-.5,.5),c('#2166AC','#EEEEEE','#B2182B')),
                       lwd=.6+3*max(z$CK_mean[i],z$treatment_mean[i]),arrow=FALSE)
        }
      }else mtext('Insufficient common-state coverage: edges not estimated',1,line=-1,cex=.58)
      counts<-table(factor(d$state[d$treatment==gr],levels=nodes$state))
      points(nodes$x,nodes$y,pch=21,bg=ifelse(counts>0,nodes$color,'white'),cex=.7+1.7*sqrt(counts/max(1,max(counts))),col='#888888')
      text(nodes$x,nodes$y,nodes$state,pos=3,cex=.5)
      mtext('Red stronger / blue weaker vs balanced CK | no direction inferred',1,line=.5,cex=.55)
    }
  },17,15)
}

TISSUES_V11 <- c('Endodermis','Cortex','Epidermis','Root Cap','Pericycle','Procambium','Phloem','Xylem')
LOCAL_K_V11 <- 30L; LOCAL_RESOLUTION_V11 <- .5; LOCAL_EDGE_THRESHOLD_V11 <- .1
ROOT_OVERRIDES_V11 <- file.path(OUT_DIR,'root_overrides_v11.csv')
SHOW_HYPOTHESIS_ARROWS <- FALSE 

local_tag <- function(x)gsub('[^[:alnum:]]','_',x)
local_read <- function(tissue,name)read.csv(file.path(OUT_DIR,paste0('L11_',local_tag(tissue),'_',name,'.csv')),stringsAsFactors=FALSE)
local_write <- function(z,tissue,name)save_table(z,paste0('L11_',local_tag(tissue),'_',name))
local_nodes <- function(nodes,tissue){
  nodes$root<-tolower(as.character(nodes$root))=='true';nodes$branchpoint<-tolower(as.character(nodes$branchpoint))=='true'
  base<-type_palette(tissue)[1];pal<-colorRampPalette(c(base,'#D9D9D9'))(nrow(nodes)+2)
  nodes$color<-pal[seq_len(nrow(nodes))];nodes$label<-sub(paste0(local_tag(tissue),'_'),'',nodes$state,fixed=TRUE)
  nodes
}
local_umap <- function(d,nodes,title,background=NULL,pt=FALSE){
  par(cex=1,mar=c(3.5,3.5,3,1),mgp=c(2,.6,0))
  bg<-if(is.null(background))d else background
  plot(bg$x,bg$y,type='n',asp=1,xlab='Original UMAP 1',ylab='Original UMAP 2',main=title,cex.main=.95)
  if(!is.null(background))points(bg$x,bg$y,pch=16,cex=.18,col='#ECECEC')
  if(pt){co<-colors_num(d$CK_PT_percentile,c(0,1),PT_COLORS)} else{co<-nodes$color[match(d$state,nodes$state)];co[is.na(co)]<-'#BDBDBD'}
  points(d$x,d$y,pch=16,cex=.28,col=co)
  if(pt)legend('topleft',legend=c('CK percentile 0','CK percentile 1','Undefined'),col=c(PT_COLORS[1],tail(PT_COLORS,1),'#CFCFCF'),pch=16,bty='n',cex=.65)
  else {text(nodes$x,nodes$y,nodes$label,pos=3,cex=.75);legend('topleft',nodes$label,col=nodes$color,pch=16,bty='n',cex=.65,ncol=max(1,ceiling(nrow(nodes)/8)))}
}
local_stats <- function(d,nodes,tissue){
  lev<-c(nodes$state,'Unmapped');rows<-list()
  for(s in unique(d$sample)){
    z<-d[d$sample==s,];ct<-table(factor(z$state,levels=lev));rows[[s]]<-data.frame(sample=s,treatment=z$treatment[1],state=lev,n=as.numeric(ct),total=nrow(z),fraction=as.numeric(ct)/nrow(z))
  }
  occ<-do.call(rbind,rows);local_write(occ,tissue,'state_occupancy_by_sample');tests<-list()
  for(gr in setdiff(GROUPS,'CK'))for(st in lev){
    a<-occ[occ$treatment=='CK'&occ$state==st&occ$total>=30,];b<-occ[occ$treatment==gr&occ$state==st&occ$total>=30,]
    if(nrow(a)<2||nrow(b)<2)next
    knownbatch<-nzchar(BATCH_COL)
    tests[[paste(gr,st)]]<-data.frame(treatment=gr,state=st,delta=mean(b$fraction)-mean(a$fraction),p=if(knownbatch)NA_real_ else exact_permutation(a$fraction,b$fraction),n_CK=nrow(a),n_treatment=nrow(b))
  }
  if(length(tests)){z<-do.call(rbind,tests);z$FDR_within_tissue<-p.adjust(z$p,'BH');local_write(z,tissue,'state_occupancy_tests')}
  ee<-local_read(tissue,'balanced_edges');changes<-list()
  if(nrow(ee)){
    z<-aggregate(connectivity~comparison+treatment+sample+source+target,ee,mean)
    local_write(z,tissue,'balanced_edges_sample_means')
    for(key in unique(paste(z$comparison,z$source,z$target,sep='|'))){
      v<-z[paste(z$comparison,z$source,z$target,sep='|')==key,];a<-v$connectivity[v$treatment=='CK'];b<-v$connectivity[v$treatment!='CK']
      if(length(a)<2||length(b)<2)next
      changes[[key]]<-data.frame(treatment=v$comparison[1],source=v$source[1],target=v$target[1],CK_mean=mean(a),treatment_mean=mean(b),delta=mean(b)-mean(a),p=if(nzchar(BATCH_COL))NA_real_ else exact_permutation(a,b))
    }
  }
  ch<-if(length(changes))do.call(rbind,changes) else NULL
  if(!is.null(ch)){ch$FDR_within_tissue<-p.adjust(ch$p,'BH');local_write(ch,tissue,'edge_changes')}
  list(occupancy=occ,changes=ch)
}
local_response_plot <- function(stats,nodes,tissue){
  occ<-stats$occupancy
  save_plot(paste0('L11_',local_tag(tissue),'_treatment_state_occupancy'),function(){
    states<-unique(occ$state);par(mfrow=c(ceiling(length(states)/3),3),mar=c(4,4,3,1))
    for(st in states){
      z<-occ[occ$state==st,];plot(NA,xlim=c(.5,length(GROUPS)+.5),ylim=c(0,max(.05,max(z$fraction))*1.1),xaxt='n',xlab='',ylab='Within-tissue cell fraction',main=sub(paste0(local_tag(tissue),'_'),'',st,fixed=TRUE))
      axis(1,seq_along(GROUPS),GROUPS,las=2,cex.axis=.7)
      for(gr in GROUPS){v<-z$fraction[z$treatment==gr];xx<-match(gr,GROUPS);if(length(v)){points(xx+seq(-.1,.1,length.out=length(v)),v,pch=16,col=GROUP_COLORS[gr]);segments(xx-.16,mean(v),xx+.16,mean(v),lwd=2)}}
    }
  },14,max(5,3.4*ceiling(length(unique(occ$state))/3)))
  ch<-stats$changes
  save_plot(paste0('L11_',local_tag(tissue),'_treatment_connectivity'),function(){
    par(mfrow=c(3,3))
    for(gr in GROUPS){
      if(gr=='CK'){draw_state_network_v10(nodes,local_read(tissue,'edges'),paste(tissue,'CK'));next}
      par(cex=1,mar=c(3.5,1,3,1));plot(nodes$x,nodes$y,type='n',asp=1,axes=FALSE,xlab='',ylab='',main=gr)
      z<-if(!is.null(ch))ch[ch$treatment==gr,,drop=FALSE] else NULL
      if(!is.null(z)&&nrow(z)){
        for(i in seq_len(nrow(z))){if(max(z$CK_mean[i],z$treatment_mean[i])<PAGA_THRESHOLD)next
          a<-match(z$source[i],nodes$state);b<-match(z$target[i],nodes$state)
          draw_curve_v10(nodes$x[a],nodes$y[a],nodes$x[b],nodes$y[b],col=colors_num(z$delta[i],c(-.5,.5),c('#2166AC','#EEEEEE','#B2182B')),lwd=.6+4*max(z$CK_mean[i],z$treatment_mean[i]))
        }
      }else mtext('Insufficient balanced-state coverage',1,line=.2,cex=.7)
      points(nodes$x,nodes$y,pch=21,bg=nodes$color,cex=2);text(nodes$x,nodes$y,nodes$label,pos=3,cex=.8)
      mtext('Red stronger; blue weaker vs balanced CK',1,line=1.2,cex=.65)
    }
  },15,13)
}
render_lineages_v11 <- function(){
  status<-read.csv(file.path(OUT_DIR,'L11_run_status.csv'));tissues<-status$tissue[status$status%in%c('complete','cached')]
  if(!length(tissues))stop('No tissue completed; inspect python_lineages_v11.log and L11_*_error.txt')
  all<-list();net<-list();ed<-list()
  PAGA_THRESHOLD<<-LOCAL_EDGE_THRESHOLD_V11
  for(tissue in tissues){
    d<-local_read(tissue,'cells');n<-local_nodes(local_read(tissue,'nodes'),tissue);e<-local_read(tissue,'edges');ck<-d[d$treatment=='CK',]
    all[[tissue]]<-ck;net[[tissue]]<-n;ed[[tissue]]<-e
    local_write(n,tissue,'plot_labels_colors')
    save_plot(paste0('L11_',local_tag(tissue),'_CK_states_network_DPT'),function(){
      par(mfrow=c(1,3));local_umap(ck,n,paste(tissue,'CK states'));draw_state_network_v10(n,e,'Within-tissue PAGA')
      local_umap(ck,n,'Local DPT | CK percentile (not linear time)',pt=TRUE)
    },19,6)
    save_plot(paste0('L11_',local_tag(tissue),'_nine_groups_states'),function(){
      par(mfrow=c(3,3));for(gr in GROUPS)local_umap(d[d$treatment==gr,],n,paste(tissue,gr),background=d)
    },16,14)
    local_response_plot(local_stats(d,n,tissue),n,tissue)
  }
  save_plot('L11_paper_style_local_networks',function(){
    par(mfrow=c(ceiling(length(tissues)/3),3));for(tissue in tissues)draw_state_network_v10(net[[tissue]],ed[[tissue]],tissue)
  },16,4.5*ceiling(length(tissues)/3))
  save_plot('L11_paper_style_originalUMAP',function(){
    d<-do.call(rbind,all);n<-do.call(rbind,net);par(mar=c(4,4,3,1));plot(d$x,d$y,type='n',asp=1,xlab='Original UMAP 1',ylab='Original UMAP 2',main='CK | tissue-specific states (no imposed two-branch labels)')
    points(d$x,d$y,pch=16,cex=.25,col=n$color[match(d$state,n$state)])
    for(tissue in tissues){z<-all[[tissue]];text(median(z$x),median(z$y),tissue,cex=.9,font=2)}
  },12,11)
  note('L11_interpretation',c('These are within-tissue CK states. No published branch identities were transferred.',
                              'Grey Unmapped cells represent low-confidence CK projection. No candidate bifurcation does not make annotated states meaningless.',
                              'Local roots inferred from global QC DPT are unvalidated proxies. Inspect each root.csv or supply root_overrides_v11.csv.',
                              'CK_PT_percentile improves colour contrast; it is nonlinear and only comparable within a tissue against its CK reference.',
                              'PAGA edges are undirected. No between-tissue developmental arrows are inferred. Two local arms do not prove a whole-root QC bifurcation.',
                              'State marker rankings are descriptive, not DEG tests. Branch audits use fixed local clusters and require independent validation.',
                              'Two biological replicates/group limit exact two-sided permutation p-values to >=1/3. Resamples are averaged before sample-level tests.',
                              'Complete/cached tissues are plotted. Any FAILED tissue must be reviewed in L11_run_status.csv and its error file.'))
}
main_v11 <- function(){
  dir.create(OUT_DIR,recursive=TRUE,showWarnings=FALSE);set.seed(SEED);plog('Start v11: tissue-specific analysis')
  prepare_existing_export_v10();PYTHON<<-find_paga_python();z<-export_input();validate_cached_export(z);rm(z);gc()
  cfg<-list(out=normalizePath(OUT_DIR,winslash='/',mustWork=TRUE),tissues=TISSUES_V11,groups=GROUPS,seed=SEED,
            neighbors=LOCAL_K_V11,resolution=LOCAL_RESOLUTION_V11,threshold=LOCAL_EDGE_THRESHOLD_V11,thresholds=sort(unique(c(.05,.1,.2,LOCAL_EDGE_THRESHOLD_V11))),edge_resamples=3L,root_file=ROOT_OVERRIDES_V11)
  jsonlite::write_json(cfg,file.path(OUT_DIR,'lineages_v11_config.json'),auto_unbox=TRUE,pretty=TRUE)
  writeLines(PYTHON_ENGINE,file.path(OUT_DIR,'lineages_v11_engine.py'),useBytes=TRUE)
  log<-file.path(OUT_DIR,'python_lineages_v11.log');plog('Run tissue-local Python engine; live log: ',log)
  code<-system2(PYTHON,c('-u',shQuote(file.path(OUT_DIR,'lineages_v11_engine.py')),shQuote(file.path(OUT_DIR,'lineages_v11_config.json'))),stdout=log,stderr=log)
  if(code!=0)stop('Python failed: ',log)
  render_lineages_v11();plog('v11 plots/tables saved: ',OUT_DIR)
}

PYTHON_ENGINE <- r"---(
"""Actual Scanpy PAGA, CK-trained reference; embedded in the standalone R script.
UMAP coordinates are never used to infer edges or replaced. No fixed branch count.
"""
import os, sys, json, math, heapq, itertools, platform, hashlib, pickle, traceback, time
import numpy as np

# --- Patch 1: scipy>=1.13 removed csr_array.getnnz ---
import scipy.sparse as _sp_compat
if not hasattr(_sp_compat.csr_array, 'getnnz'):
    def _csr_array_getnnz(self, axis=None):
        if axis is None: return self.nnz
        elif axis == 1: return np.diff(self.indptr)
        elif axis == 0: return _sp_compat.csc_matrix(self).getnnz(axis=0)
        else: raise ValueError("axis must be None, 0, or 1")
    _sp_compat.csr_array.getnnz = _csr_array_getnnz
# --- End patch 1 ---

def _patch_umap_knn_contiguity():
    """Install in EVERY Python process, including checkpoint-resume/full phase.
    Only memory layout and distance dtype are normalized; graph settings unchanged.
    """
    import inspect
    import functools
    import umap.umap_ as _umap_mod
    current = _umap_mod.fuzzy_simplicial_set
    if getattr(current, '_paga_knn_contiguity', False):
        patched = current
        original = current.__wrapped__
    else:
        original = current
        signature = inspect.signature(original)
        @functools.wraps(original)
        def patched(*args, **kwargs):
            bound = signature.bind_partial(*args, **kwargs)
            for name in ('knn_dists', 'knn_indices'):
                value = bound.arguments.get(name)
                if value is not None:
                    bound.arguments[name] = np.ascontiguousarray(
                        value, dtype=np.float32 if name == 'knn_dists' else None)
            return original(*bound.args, **bound.kwargs)
        patched._paga_knn_contiguity = True
        _umap_mod.fuzzy_simplicial_set = patched
    # Also update an already imported Scanpy alias, if this version has one.
    connectivity = sys.modules.get('scanpy.neighbors._connectivity')
    if connectivity is not None and getattr(connectivity, 'fuzzy_simplicial_set', None) is original:
        connectivity.fuzzy_simplicial_set = patched
    return True

def components(nodes, edges, banned=None):
    adj = {v: set() for v in nodes if v != banned}
    for a, b, w in edges:
        if a in adj and b in adj:
            adj[a].add(b); adj[b].add(a)
    result = []
    while adj:
        todo = [next(iter(adj))]; seen = set()
        while todo:
            v = todo.pop()
            if v in seen: continue
            seen.add(v); todo.extend(adj[v] - seen)
        result.append(seen)
        for v in seen: del adj[v]
    return result

def distances(root, nodes, edges):
    adj = {v: [] for v in nodes}
    for a,b,w in edges:
        cost = -math.log(max(min(w, 1.0), 1e-12)) + .001
        adj[a].append((b,cost)); adj[b].append((a,cost))
    out = {root: 0.0}; heap = [(0.0,root)]
    while heap:
        dv,v = heapq.heappop(heap)
        if dv > out[v]: continue
        for u,c in adj[v]:
            nd = dv+c
            if nd < out.get(u,float('inf')):
                out[u] = nd; heapq.heappush(heap,(nd,u))
    return out

def branch_candidates(states, matrix, root, sizes, state_samples, threshold,
                      min_states=2, min_cells=100, min_samples=2):
    edges = [(a,b,float(matrix[i,j])) for i,a in enumerate(states)
             for j,b in enumerate(states) if i<j and np.isfinite(matrix[i,j])
             and matrix[i,j]>=threshold]
    reachable = next((c for c in components(states,edges) if root in c),set())
    ed = [e for e in edges if e[0] in reachable and e[1] in reachable]
    dist = distances(root,reachable,ed) if root in reachable else {}
    result = []
    for cut in sorted(reachable):
        cs = components(reachable,ed,banned=cut)
        upstream = set() if cut==root else next((c for c in cs if root in c),set())
        arms = [c for c in cs if not (root in c)]
        valid=[];small=[]
        for arm in arms:
            n=sum(sizes.get(s,0) for s in arm)
            support=set().union(*(state_samples.get(s,set()) for s in arm))
            if len(arm)>=min_states and n>=min_cells and len(support)>=min_samples: valid.append(arm)
            else: small.append(arm)
        if len(valid)>=2:
            valid.sort(key=lambda arm:(-sum(sizes.get(s,0) for s in arm),sorted(arm)))
            result.append(dict(cut=cut,arms=valid,small=small,upstream=upstream|{cut},
                root_distance=dist.get(cut,float('inf')),n_arms=len(valid),threshold=threshold))
    return sorted(result,key=lambda z:(z['root_distance'],z['cut']))

def match_arms(reference, candidate, minimum=.70):
    if candidate is None or len(reference)!=len(candidate): return False,0.0
    if len(reference)>7: return False,0.0
    best=0.0
    for order in itertools.permutations(candidate):
        score=min(len(a&b)/max(1,len(a|b)) for a,b in zip(reference,order))
        best=max(best,score)
    return best>=minimum,best

def Path_bytes():
    with open(__file__,'rb') as handle: return handle.read()

def qc_mask(series, labels):
    # Convert Categorical BEFORE fillna: '' need not be a source category.
    # Nullable StringDtype preserves real missing values until they are filled.
    if isinstance(labels, str):
        labels = [labels]
    targets = [str(x).strip().casefold() for x in labels if str(x).strip()]
    text = series.astype('string').fillna('').str.strip().str.casefold()
    return text.isin(targets).to_numpy(dtype=bool)

def test_qc_mask_types():
    import pandas as pd
    values = ['QC', 'Meristem', ' qc ', None, 'Quiescent Center', '']
    expected = np.array([True, False, True, False, True, False], dtype=bool)
    cases = [pd.Series(values, dtype=object),
             pd.Series(values, dtype='string'),
             pd.Series(pd.Categorical(values)),
             pd.Series(pd.Categorical(values, ordered=True))]
    for v in cases:
        before = v.copy(deep=True)
        got = qc_mask(v, ['QC', 'Quiescent Center'])
        if not np.array_equal(got, expected):
            raise AssertionError('QC mask type regression failed: ' + str(v.dtype))
        pd.testing.assert_series_equal(v, before)
    # Reproduces the server case: missing value, without '' as a category.
    cat = pd.Series(pd.Categorical(['QC', None, 'Meristem'], categories=['QC', 'Meristem']))
    if qc_mask(cat, 'QC').tolist() != [True, False, False]:
        raise AssertionError('Categorical without empty category failed')


def lineage_main(config_path):
    import pandas as pd
    import scipy.sparse as sp
    from scipy.io import mmread
    from scipy.sparse.csgraph import connected_components
    from sklearn.decomposition import PCA
    from sklearn.neighbors import NearestNeighbors
    import scanpy as sc
    import anndata as ad
    from importlib.metadata import version
    c=json.load(open(config_path,encoding='utf-8'));out=c['out'];seed=int(c['seed'])
    def log(*s):print(time.strftime('%H:%M:%S'),*s,flush=True)
    _patch_umap_knn_contiguity();test_qc_mask_types()
    meta=pd.read_csv(os.path.join(out,'input_metadata.csv'));genes=pd.read_csv(os.path.join(out,'input_genes.csv')).gene_id.astype(str).to_numpy()
    priorfile=os.path.join(out,'cells_CK_reference_mapping.csv')
    prior=pd.read_csv(priorfile,usecols=['cell','reference_PT']) if os.path.isfile(priorfile) else None
    meta['prior_PT']=meta.cell.map(prior.set_index('cell').reference_PT) if prior is not None else np.nan
    rootfile=c['root_file'];roots=pd.read_csv(rootfile,keep_default_na=False) if rootfile and os.path.isfile(rootfile) else pd.DataFrame()
    versions={p:version(p) for p in ['scanpy','anndata','numpy','pandas','scipy','umap-learn','leidenalg']}
    fingerprint=[(f,os.stat(os.path.join(out,f)).st_size,os.stat(os.path.join(out,f)).st_mtime_ns) for f in ['input_counts.mtx','input_metadata.csv','input_genes.csv']]
    key=hashlib.sha256(Path_bytes()+json.dumps([c,versions,fingerprint,roots.to_dict(),None if prior is None else hashlib.sha256(prior.to_csv(index=False).encode()).hexdigest()],sort_keys=True).encode()).hexdigest()
    count=None;status=[]
    for tissue in c['tissues']:
      tag=''.join(x if x.isalnum() else '_' for x in tissue);prefix='L11_'+tag+'_'
      def csv(df,name):df.to_csv(os.path.join(out,prefix+name+'.csv'),index=False)
      stamp=os.path.join(out,prefix+'complete.sha256')
      if os.path.isfile(stamp) and open(stamp).read().strip()==key:
          status.append(dict(tissue=tissue,status='cached'));log(tissue,'cached');continue
      try:
        rows=np.flatnonzero(meta.cell_type.astype(str).str.casefold().to_numpy()==tissue.casefold())
        md=meta.iloc[rows].copy().reset_index(drop=True);ck=np.flatnonzero(md.treatment.to_numpy()=='CK')
        if len(ck)<100 or md.iloc[ck]['sample'].nunique()<2:raise ValueError('Need >=100 CK cells and >=2 CK biological samples')
        if count is None:
            log('Reading existing sparse count matrix once');count=mmread(os.path.join(out,'input_counts.mtx')).T.tocsr().astype(np.float32)
            if count.shape!=(len(meta),len(genes)):raise ValueError('Count dimensions mismatch')
        log(tissue,'CK-local HVG/PCA and graph; CK cells',len(ck))
        raw=count[rows];lib=np.asarray(raw.sum(axis=1)).ravel();norm=(sp.diags(1e4/lib)@raw).tocsr();norm.data=np.log1p(norm.data)
        detect=np.asarray((raw[ck]>0).sum(axis=0)).ravel()>=max(5,int(.01*len(ck)))
        a=ad.AnnData(norm[ck][:,detect],var=pd.DataFrame(index=genes[detect]))
        sc.pp.highly_variable_genes(a,flavor='seurat',n_top_genes=min(2000,a.n_vars))
        use=np.flatnonzero(detect)[a.var.highly_variable.to_numpy()];npc=min(30,len(use)-1,len(ck)-2)
        if npc<5:raise ValueError('Too few variable genes')
        pca=PCA(n_components=npc,svd_solver='randomized',random_state=seed);fit=pca.fit_transform(norm[ck][:,use].toarray())
        pcs=np.empty((len(md),npc),dtype=np.float32)
        for j in range(0,len(md),2048):pcs[j:j+2048]=pca.transform(norm[j:j+2048][:,use].toarray())
        pcs[ck]=fit
        a=ad.AnnData(sp.csr_matrix((len(ck),0)),obs=md.iloc[ck].set_index('cell').copy());a.obsm['X_pca']=np.ascontiguousarray(fit,dtype=np.float32)
        k=min(int(c['neighbors']),len(ck)-1)
        sc.pp.neighbors(a,n_neighbors=k,use_rep='X_pca',random_state=seed)
        sc.tl.leiden(a,resolution=float(c['resolution']),key_added='state',random_state=seed,flavor='leidenalg')
        labels=sorted(a.obs.state.unique(),key=int);ren={v:tag+'_S'+str(i+1).zfill(2) for i,v in enumerate(labels)}
        a.obs['state']=pd.Categorical([ren[str(x)] for x in a.obs.state]);states=list(a.obs.state.cat.categories);lab=a.obs.state.astype(str).to_numpy()
        def paga(aa,rebuild=False):
            if rebuild:sc.pp.neighbors(aa,n_neighbors=min(k,len(aa)-1),use_rep='X_pca',random_state=seed)
            aa.obs['state']=aa.obs.state.cat.remove_unused_categories();ls=list(aa.obs.state.cat.categories)
            m=np.full((len(states),len(states)),np.nan)
            if len(ls)<2:return m
            sc.tl.paga(aa,groups='state');ix=[states.index(x) for x in ls];m[np.ix_(ix,ix)]=aa.uns['paga']['connectivities'].toarray();return m
        mat=paga(a);pt=np.full(len(ck),np.nan);rootstate=None;rootcell='';rootmethod='none';validated=False
        override=roots[roots.tissue==tissue] if 'tissue' in roots else pd.DataFrame()
        if len(override)>1:raise ValueError('Multiple root override rows')
        if len(override):
            rootcell=str(override.iloc[0].root_cell)
            if rootcell not in a.obs_names:raise ValueError('Root override must be a CK cell in this tissue')
            ri=int(a.obs_names.get_loc(rootcell));rootmethod='user_selected_CK_cell';validated=str(override.iloc[0].get('validated','False')).lower()=='true'
        else:
            # Reuse global QC-root DPT only as an explicitly unvalidated local-root proxy.
            candidates=[]
            for st in states:
                z=a.obs[a.obs.state==st];v=z.prior_PT.to_numpy();v=v[np.isfinite(v)]
                if len(v)>=10 and z['sample'].nunique()>=2:candidates.append((float(np.median(v)),st))
            if candidates:
                rootstate=min(candidates)[1];ii=np.flatnonzero(lab==rootstate)
                ri=int(ii[np.argmin(np.linalg.norm(fit[ii]-fit[ii].mean(axis=0),axis=1))]);rootcell=str(a.obs_names[ri]);rootmethod='global_QC_DPT_proxy_UNVALIDATED'
        if rootcell:
            rootstate=lab[ri];_,cc=connected_components(a.obsp['connectivities'],directed=False);ok=cc==cc[ri]
            sub=a[ok].copy()
            if len(sub)>=20:
                sub.uns['iroot']=int(np.flatnonzero(np.flatnonzero(ok)==ri)[0]);nd=min(15,len(sub)-2)
                sc.tl.diffmap(sub,n_comps=nd,random_state=seed);sc.tl.dpt(sub,n_dcs=min(10,nd),n_branchings=0)
                pt[ok]=sub.obs.dpt_pseudotime.to_numpy();pt[~np.isfinite(pt)]=np.nan
        csv(pd.DataFrame([dict(tissue=tissue,root_cell=rootcell,root_state=rootstate,root_method=rootmethod,validated=validated)]),'root')
        # Map ONLY cells annotated as the same tissue; retain unmatched cells.
        nnk=min(15,len(ck)-1);nn=NearestNeighbors(n_neighbors=nnk+1).fit(fit);ds,js=nn.kneighbors(fit)
        radius={st:np.quantile(ds[lab==st,1:].mean(axis=1),.99) for st in states}
        mapped=np.full(len(md),'Unmapped',object);proj=np.full(len(md),np.nan);conf=np.zeros(len(md))
        for j in range(0,len(md),2048):
            dd,jj=nn.kneighbors(pcs[j:j+2048],n_neighbors=nnk)
            for q in range(len(dd)):
                ls=lab[jj[q]];u,n=np.unique(ls,return_counts=True);st=u[n.argmax()];frac=n.max()/nnk;conf[j+q]=frac
                if frac>=.6 and dd[q].mean()<=radius[st]*1.5:
                    mapped[j+q]=st;valid=(ls==st)&np.isfinite(pt[jj[q]])
                    if valid.sum()>=3:proj[j+q]=np.average(pt[jj[q][valid]],weights=1/(dd[q][valid]+1e-6))
        mapped[ck]=lab;proj[ck]=pt;conf[ck]=1
        finite=np.sort(pt[np.isfinite(pt)]);rank=np.full(len(md),np.nan);ok=np.isfinite(proj)
        if len(finite):rank[ok]=np.searchsorted(finite,proj[ok],side='right')/len(finite)
        md['state']=mapped;md['reference_PT']=proj;md['CK_PT_percentile']=rank;md['confidence']=conf;md['tissue']=tissue
        sizes={st:int(np.sum(lab==st)) for st in states};support={st:set(a.obs.loc[lab==st,'sample'].value_counts().loc[lambda x:x>=5].index) for st in states}
        allcand=[];chosen=None
        for threshold in c['thresholds']:
            cs=branch_candidates(states,mat,rootstate,sizes,support,float(threshold),2,50,2) if rootstate else []
            if abs(threshold-c['threshold'])<1e-9:chosen=cs[0] if cs else None
            allcand.append(dict(threshold=threshold,n_arms=cs[0]['n_arms'] if cs else 0,cut=cs[0]['cut'] if cs else '',root_method=rootmethod))
        csv(pd.DataFrame(allcand),'branch_audit')
        repaudit=[]
        for sample in a.obs['sample'].unique():
            aa=a[a.obs['sample']==sample].copy()
            if len(aa)<50:continue
            mm=paga(aa,True);sz=aa.obs.state.astype(str).value_counts().to_dict()
            su={st:{str(sample)} for st in sz}
            cs=branch_candidates(states,mm,rootstate,sz,su,c['threshold'],2,25,1) if rootstate else []
            match,jac=match_arms(chosen['arms'],cs[0]['arms'] if cs else None) if chosen else (False,np.nan)
            repaudit.append(dict(sample=sample,n_arms=cs[0]['n_arms'] if cs else 0,match_reference_arms=match,Jaccard=jac))
        csv(pd.DataFrame(repaudit),'branch_replicate_audit')
        bm={st:'Not_assigned_to_candidate_arm' for st in states}
        if chosen:
            for st in chosen['upstream']:bm[st]='Shared_trunk'
            for bi,arm in enumerate(chosen['arms'],1):
                for st in arm:bm[st]='Candidate_arm_'+str(bi)
        md['candidate_arm']=[bm.get(st,'Unmapped') for st in mapped];csv(md,'cells')
        node=[]
        for st in states:
            ii=np.flatnonzero(lab==st);z=md.iloc[ck[ii]];xy=z[['x','y']].to_numpy();center=xy.mean(axis=0);anchor=xy[np.argmin(np.linalg.norm(xy-center,axis=1))]
            node.append(dict(state=st,type=tissue,x=anchor[0],y=anchor[1],n_cells=len(z),reference_PT=float(np.nanmedian(pt[ii])) if np.any(np.isfinite(pt[ii])) else np.nan,
                             root=st==rootstate,branchpoint=bool(chosen and st==chosen['cut']),branch=bm[st]))
        csv(pd.DataFrame(node),'nodes');csv(pd.DataFrame([dict(source=states[i],target=states[j],connectivity=mat[i,j]) for i in range(len(states)) for j in range(i+1,len(states))],columns=['source','target','connectivity']),'edges')
        # Descriptive marker ranking, replicated direction; no fabricated P values.
        marker=[]
        for st in states:
            ii=np.flatnonzero(lab==st);oo=np.flatnonzero(lab!=st)
            if len(ii)<10 or len(oo)<10:continue
            exp1=np.asarray(norm[ck[ii]].mean(axis=0)).ravel();exp0=np.asarray(norm[ck[oo]].mean(axis=0)).ravel()
            delta=exp1-exp0;pct=np.asarray((raw[ck[ii]]>0).mean(axis=0)).ravel();eligible=np.flatnonzero((delta>0)&(pct>=.1))
            top=eligible[np.argsort(delta[eligible])[::-1][:20]]
            for gi in top:
                directions=[]
                for sample in a.obs['sample'].unique():
                    ins=np.flatnonzero((lab==st)&(a.obs['sample'].to_numpy()==sample));outs=np.flatnonzero((lab!=st)&(a.obs['sample'].to_numpy()==sample))
                    if len(ins)>=10 and len(outs)>=10:directions.append(float(norm[ck[ins],gi].mean()-norm[ck[outs],gi].mean())>0)
                marker.append(dict(state=st,gene=genes[gi],delta_mean_log1p=delta[gi],pct_detected=pct[gi],n_samples_testable=len(directions),n_positive_samples=sum(directions)))
        csv(pd.DataFrame(marker,columns=['state','gene','delta_mean_log1p','pct_detected','n_samples_testable','n_positive_samples']),'markers_descriptive')
        # Equal cells/state/sample for each treatment vs CK; retain biological sample identity.
        edge_rows=[];coverage=[];rng=np.random.default_rng(seed)
        for gr in c['groups']:
            if gr=='CK':continue
            sub=md[md.treatment.isin(['CK',gr])];samples=sorted(sub['sample'].unique())
            if gr not in set(sub.treatment) or sub.groupby('treatment')['sample'].nunique().min()<2:continue
            counts=pd.crosstab(sub.state,sub['sample']).reindex(index=states,columns=samples,fill_value=0)
            cap=counts.min(axis=1).clip(upper=150).astype(int);keep=list(cap[cap>=10].index)
            coverage.extend([dict(treatment=gr,state=st,n_balanced=int(cap[st]),included=st in keep) for st in states])
            if len(keep)<2:continue
            for rep in range(int(c['edge_resamples'])):
                for sample in samples:
                    ids=np.concatenate([rng.choice(np.flatnonzero((md['sample'].to_numpy()==sample)&(mapped==st)),int(cap[st]),replace=False) for st in keep])
                    obs=md.iloc[ids].set_index('cell').copy();obs['state']=pd.Categorical(obs.state)
                    aa=ad.AnnData(sp.csr_matrix((len(ids),0)),obs=obs);aa.obsm['X_pca']=np.ascontiguousarray(pcs[ids]);mm=paga(aa,True)
                    for i in range(len(states)):
                        for j in range(i+1,len(states)):
                            if np.isfinite(mm[i,j]):edge_rows.append(dict(comparison=gr,treatment=obs.treatment.iloc[0],sample=sample,resample=rep,source=states[i],target=states[j],connectivity=mm[i,j]))
            log(tissue,'balanced connectivity',gr)
        csv(pd.DataFrame(coverage,columns=['treatment','state','n_balanced','included']),'edge_coverage')
        csv(pd.DataFrame(edge_rows,columns=['comparison','treatment','sample','resample','source','target','connectivity']),'balanced_edges')
        with open(stamp,'w') as h:h.write(key)
        status.append(dict(tissue=tissue,status='complete'));log(tissue,'complete',len(states),'states')
      except Exception:
        detail=traceback.format_exc();log(detail)
        with open(os.path.join(out,prefix+'error.txt'),'w') as h:h.write(detail)
        status.append(dict(tissue=tissue,status='FAILED; inspect error.txt'))
    pd.DataFrame(status).to_csv(os.path.join(out,'L11_run_status.csv'),index=False)
    if any('FAILED' in s['status'] for s in status):log('Some tissues failed; R will render only completed tissue results')

if __name__=='__main__':
    try:lineage_main(sys.argv[1])
    except Exception:
        traceback.print_exc();sys.exit(1)

)---"

if (sys.nframe() == 0L) main_v11()
