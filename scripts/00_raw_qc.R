# Raw 10X QC entry point. Thresholds must come from reviewed sample QC plots.
# Run from repository root after filling config/samples.csv.
required <- c("Seurat", "Matrix")
for (pkg in required) if (!requireNamespace(pkg, quietly=TRUE)) stop("Missing R package: ",pkg)
manifest <- "config/samples.csv"
if (!file.exists(manifest)) stop("Create config/samples.csv from config/samples_template.csv; do not guess QC thresholds.")
samples <- read.csv(manifest, stringsAsFactors=FALSE, check.names=FALSE)
fields <- c("sample_id","treatment","tenx_dir","min_features","max_features","min_counts","max_counts")
if (!all(fields %in% names(samples)) || !nrow(samples)) stop("Sample manifest is incomplete.")
if (anyDuplicated(samples$sample_id)) stop("Duplicate sample IDs.")
source("config/hormone_config.R", encoding="UTF-8")
canonical <- function(x) tolower(gsub("[^[:alnum:]]", "", x))
lookup <- setNames(rep(names(CONDITION_MAP),lengths(CONDITION_MAP)),canonical(unlist(CONDITION_MAP,use.names=FALSE)))
samples$treatment <- unname(lookup[canonical(samples$treatment)])
if (anyNA(samples$treatment)) stop("Unrecognized treatment in sample manifest.")
for(nm in fields[4:7]) if(any(!is.finite(samples[[nm]]))) stop("Invalid QC threshold: ",nm)
if(any(samples$min_features>samples$max_features | samples$min_counts>samples$max_counts)) stop("Inverted QC thresholds.")
dir.create("results/qc",recursive=TRUE,showWarnings=FALSE)
audit <- list()
for(i in seq_len(nrow(samples))) {
 s <- samples[i,]; counts <- Seurat::Read10X(s$tenx_dir)
 if(is.list(counts)) {if(!"Gene Expression" %in% names(counts)) stop("10X Gene Expression matrix missing.");counts<-counts[["Gene Expression"]]}
 object <- Seurat::CreateSeuratObject(counts,project=s$sample_id,min.cells=0,min.features=0)
 object$sample_id <- s$sample_id;object$treatment <- s$treatment
 object <- Seurat::RenameCells(object,add.cell.id=s$sample_id)
 md <- object[[]]
 keep <- md$nFeature_RNA>=s$min_features & md$nFeature_RNA<=s$max_features & md$nCount_RNA>=s$min_counts & md$nCount_RNA<=s$max_counts
 write.csv(data.frame(cell=rownames(md),md,keep=keep),file.path("results/qc",paste0(s$sample_id,"_cell_qc.csv")),row.names=FALSE)
 if(!any(keep)) stop("No cells pass QC for ",s$sample_id)
 saveRDS(subset(object,cells=rownames(md)[keep]),file.path("results/qc",paste0(s$sample_id,"_filtered.rds")))
 audit[[i]] <- data.frame(sample_id=s$sample_id,treatment=s$treatment,input=ncol(object),retained=sum(keep))
}
write.csv(do.call(rbind,audit),"results/qc/sample_qc_summary.csv",row.names=FALSE)
# Organelle fractions and DoubletFinder require the correct soybean feature lists
# and reviewed sample-specific parameters; no mitochondrial gene prefix is assumed.
