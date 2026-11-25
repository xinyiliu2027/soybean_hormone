library(qs)

reference_renamed <- qread(Sys.getenv("SOYBEAN_REFERENCE_QS", unset = "data/reference_renamed.qs"))

gm_objects <- ls(pattern = "Gm")
print("Objects containing 'Gm':")
print(gm_objects)

rename_map <- c(
  "Gm_1031_N327F" = "eBL1",
  "Gm_2213_eBL" = "eBL2", 
  "Gm_2224_8.37temp" = "highTemp1",
  "Gm_2225_9.6Temp" = "highTemp2",
  "Gm_2226_9.6Temp" = "lowTemp1",
  "Gm_2227_9.6Temp" = "lowTemp2",
  "Gm_2214_SA" = "SA1",
  "Gm_2215_SA" = "SA2",
  "Gm_2216_MEJA" = "MEJA1",
  "Gm_2217_MEJA" = "MEJA2",
  "Gm_2218_ABA" = "ABA1",
  "Gm_2219_ABA" = "ABA2",
  "Gm_2222_PEG" = "PEG1",
  "Gm_2223_PEG" = "PEG2",
  "Gm_2220_IAA" = "IAA1",
  "Gm_2221_IAA" = "IAA2"
)

for (old_name in gm_objects) {
  
  if (old_name %in% c("Gm_2228_CK", "Gm_2229_CK")) {
    rm(list = old_name)
    cat("Removed:", old_name, "\n")
  } 
  
  else if (old_name %in% names(rename_map)) {
    new_name <- rename_map[old_name]
    assign(new_name, get(old_name))
    rm(list = old_name)
    cat("Renamed:", old_name, "->", new_name, "\n")
  }
}

print("\nObjects after processing:")
seurat_objects <- ls()[sapply(ls(), function(x) class(get(x))[1] == "Seurat")]
print(seurat_objects)

seurat_objects <- ls()[sapply(ls(), function(x) class(get(x))[1] == "Seurat")]
print("Current Seurat objects:")
print(seurat_objects)

for (obj_name in seurat_objects) {
  obj <- get(obj_name)
  
  cat("\n=== Checking", obj_name, "===\n")
  
  col_names <- colnames(obj@meta.data)
  
  doublet_cols <- grep("doublet|DF|DoubletFinder", col_names, ignore.case = TRUE, value = TRUE)
  classification_cols <- grep("classification", col_names, ignore.case = TRUE, value = TRUE)
  singlet_cols <- grep("singlet", col_names, ignore.case = TRUE, value = TRUE)
  
  if (length(doublet_cols) > 0) {
    cat("Found doublet-related columns:", paste(doublet_cols, collapse = ", "), "\n")
  }
  if (length(classification_cols) > 0) {
    cat("Found classification columns:", paste(classification_cols, collapse = ", "), "\n")
  }
  if (length(singlet_cols) > 0) {
    cat("Found singlet-related columns:", paste(singlet_cols, collapse = ", "), "\n")
  }
  
  if (length(col_names) <= 20) {
    cat("All meta.data columns:", paste(col_names, collapse = ", "), "\n")
  } else {
    cat("Total columns:", length(col_names), "\n")
    cat("First 10 columns:", paste(head(col_names, 10), collapse = ", "), "\n")
    cat("Last 10 columns:", paste(tail(col_names, 10), collapse = ", "), "\n")
  }
}

for (obj_name in seurat_objects) {
  obj <- get(obj_name)
  
  cat("\n=== Processing", obj_name, "===\n")
  cat("Original cells:", ncol(obj), "\n")
  
  df_col <- grep("^DF\\.classifications", colnames(obj@meta.data), value = TRUE)
  
  if (length(df_col) > 1L) stop("Multiple DoubletFinder classification columns in ", obj_name, "; select one explicitly.")
  if (length(df_col) > 0) {
    cat("Found doublet column:", df_col, "\n")
    print(table(obj@meta.data[[df_col]]))
    
    cells_to_keep <- rownames(obj@meta.data)[obj@meta.data[[df_col]] == "Singlet"]
    obj_filtered <- obj[, cells_to_keep]
    
    cat("After removing doublets:", ncol(obj_filtered), "cells\n")
    cat("Removed", ncol(obj) - ncol(obj_filtered), "doublets\n")
    
    assign(obj_name, obj_filtered)
    
  } else {
    cat("No doublet information found\n")
  }
}

cat("\n=== Summary after doublet removal ===\n")
for (obj_name in seurat_objects) {
  obj <- get(obj_name)
  cat(obj_name, ":", ncol(obj), "cells\n")
}

rm(list = intersect(c("obj", "obj_filtered", "df_col", "cells_to_keep"), ls()))
rm(reference_renamed)

library(Seurat)
library(dplyr)
library(patchwork)

reference_renamed <- qread(Sys.getenv("SOYBEAN_REFERENCE_QS", unset = "data/reference_renamed.qs"))

hormone.list <- list()
for (obj_name in seurat_objects) {
  hormone.list[[obj_name]] <- get(obj_name)
}

cat("Reference object info:\n")
cat("Cells:", ncol(reference_renamed), "\n")
cat("Features:", nrow(reference_renamed), "\n")
cat("Cell types:", length(unique(reference_renamed$cell_type)), "\n")

cat("\n=== Preprocessing query datasets ===\n")
for (i in 1:length(hormone.list)) {
  cat("Processing", names(hormone.list)[i], "\n")
  hormone.list[[i]] <- NormalizeData(hormone.list[[i]], verbose = FALSE) %>% 
    FindVariableFeatures(verbose = FALSE) %>% 
    SCTransform(variable.features.n = 3000,
                method = "glmGamPoi",
                verbose = FALSE)
}

cat("\n=== Preparing reference dataset ===\n")
DefaultAssay(reference_renamed) <- "RNA"
reference_renamed <- JoinLayers(reference_renamed)

reference <- reference_renamed

reference <- SCTransform(reference,
                         variable.features.n = 3000,
                         method = "glmGamPoi",
                         verbose = FALSE)

reference <- RunPCA(reference, verbose = FALSE)
reference <- RunUMAP(reference, reduction = "pca", dims = 1:50, 
                     umap.method = "umap-learn", return.model = TRUE,
                     verbose = FALSE)

DefaultAssay(reference) <- "RNA"

cat("Reference preparation completed!\n")

cat("\n=== Performing anchor integration and label transfer ===\n")

for (i in names(hormone.list)) {
  cat("\n--- Processing", i, "---\n")
  
  query <- hormone.list[[i]]
  
  DefaultAssay(query) <- "SCT"
  DefaultAssay(reference) <- "SCT"
  
  anchors <- FindTransferAnchors(
    reference = reference,
    query = query,
    reference.assay = "SCT",
    query.assay = "SCT",
    normalization.method = "SCT",
    recompute.residuals = FALSE,
    dims = 1:50,
    reference.reduction = "pca"
  )
  
  predictions <- TransferData(
    anchorset = anchors,
    refdata = reference$cell_type,  
    dims = 1:50,
    weight.reduction = "pcaproject"
  )
  
  query <- AddMetaData(query, metadata = predictions$predicted.id, col.name = "cell_type")
  
  map_query <- MapQuery(
    anchorset = anchors, 
    reference = reference, 
    query = query, 
    refdata = list(celltype = "cell_type"), 
    reference.reduction = "pca", 
    reduction.model = "umap"
  )
  
  query@reductions$umap <- map_query@reductions$ref.umap
  
  hormone.list[[i]] <- query
  
  cat("Completed annotation for", i, "\n")
}

cat("\nLabel transfer completed for all samples!\n")

clean_with_joinlayers <- function() {
  cat("\n=== Alternative: Using JoinLayers ===\n")
  
  reference_joined <- JoinLayers(reference)
  hormone_joined <- lapply(hormone.list, JoinLayers)
  names(hormone_joined) <- names(hormone.list)
  
  merged_joined <- merge(x = reference_joined, 
                         y = hormone_joined, 
                         add.cell.ids = c("ref", names(hormone_joined)))
  
  return(merged_joined)
}
