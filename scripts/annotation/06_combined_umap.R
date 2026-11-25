library(Seurat)
library(ggplot2)
library(gridExtra)
library(dplyr)

custom_colors <- c(
  "Cortex" = "#a2d2e7",           
  "Dividing Cell" = "#67a8cd",    
  "Endodermis" = "#ffc17f",       
  "Epidermis" = "#cf9f88",        
  "Lateral Root" = "#6fb3a8",     
  "Meristem" = "#b3e19b",         
  "Pericycle" = "#50aa4b",        
  "Phloem" = "#ff9d9f",           
  "QC" = "#f36569",               
  "Root Cap" = "#3581b7",         
  "Xylem" = "#cdb6da"             
)

group_by_col <- "cell_type"  

all_samples <- c(list("REF" = reference), hormone.list)

merge_pairs <- list(
  "ABA" = c("ABA1", "ABA2"),
  "eBL" = c("eBL1", "eBL2"),
  "highTemp" = c("highTemp1", "highTemp2"),
  "IAA" = c("IAA1", "IAA2"),
  "lowTemp" = c("lowTemp1", "lowTemp2"),
  "MEJA" = c("MEJA1", "MEJA2"),
  "PEG" = c("PEG1", "PEG2"),
  "SA" = c("SA1", "SA2")
)

cat("Fixing UMAP information for merged samples...\n")

for(merged_name in names(merge_pairs)) {
  if(merged_name %in% names(merged_samples)) {
    sample_pair <- merge_pairs[[merged_name]]
    
    sample1 <- all_samples[[sample_pair[1]]]
    sample2 <- all_samples[[sample_pair[2]]]
    
    if("umap" %in% names(sample1@reductions) && "umap" %in% names(sample2@reductions)) {
      
      merged_sample <- merged_samples[[merged_name]]
      
      umap1 <- sample1@reductions$umap@cell.embeddings
      umap2 <- sample2@reductions$umap@cell.embeddings
      
      rownames(umap1) <- paste0(sample_pair[1], "_", rownames(umap1))
      rownames(umap2) <- paste0(sample_pair[2], "_", rownames(umap2))
      
      combined_umap <- rbind(umap1, umap2)
      
      if(all(colnames(merged_sample) %in% rownames(combined_umap))) {
        combined_umap <- combined_umap[colnames(merged_sample), ]
        
        merged_sample[["umap"]] <- CreateDimReducObject(
          embeddings = combined_umap,
          key = "UMAP_",
          assay = DefaultAssay(merged_sample)
        )
        
        merged_samples[[merged_name]] <- merged_sample
        
        cat("Fixed UMAP for", merged_name, "\n")
      } else {
        cat("Warning: Cell names don't match for", merged_name, "\n")
      }
    }
  }
}

merged_dimplot_list <- list()

for(sample_name in names(merged_samples)) {
  sample_obj <- merged_samples[[sample_name]]
  
  if(!"umap" %in% names(sample_obj@reductions)) {
    cat("Warning: No UMAP found for", sample_name, ", skipping...\n")
    next
  }
  
  if(!group_by_col %in% colnames(sample_obj@meta.data)) {
    cat("Warning: No", group_by_col, "found for", sample_name, ", using seurat_clusters\n")
    current_group_col <- "seurat_clusters"
  } else {
    current_group_col <- group_by_col
  }
  
  p <- DimPlot(sample_obj, 
               reduction = "umap", 
               group.by = current_group_col,
               label = FALSE,  
               pt.size = 0.3,
               cols = custom_colors) +  
    ggtitle(sample_name) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 8))
  
  merged_dimplot_list[[sample_name]] <- p
  
  cat("Created dimplot for", sample_name, "\n")
}

if(length(merged_dimplot_list) > 0) {
  
  plot_order <- c("REF", "ABA", "eBL", "highTemp", 
                  "IAA", "lowTemp", "MEJA", "PEG", "SA")
  
  ordered_plots <- list()
  for(name in plot_order) {
    if(name %in% names(merged_dimplot_list)) {
      ordered_plots[[length(ordered_plots) + 1]] <- merged_dimplot_list[[name]]
    }
  }
  
  if(length(ordered_plots) > 0) {
    grid_plot <- do.call(grid.arrange, c(ordered_plots, ncol = 3))
    
    ggsave("merged_samples_dimplot.pdf", grid_plot, width = 12, height = 12, dpi = 300)
    ggsave("merged_samples_dimplot.png", grid_plot, width = 12, height = 12, dpi = 300)
    
    print(grid_plot)
    cat("Plots saved as merged_samples_dimplot.pdf and .png\n")
  } else {
    cat("No plots to display\n")
  }
} else {
  cat("No valid plots created\n")
}

cat("\n=== Merge summary ===\n")
cat("Successfully merged samples:", length(merged_samples), "\n")
for(name in names(merged_samples)) {
  n_cells <- ncol(merged_samples[[name]])
  has_umap <- "umap" %in% names(merged_samples[[name]]@reductions)
  cat(name, ":", n_cells, "cells, UMAP:", has_umap, "\n")
}
cat("Successfully transferred UMAP coordinates to merged object\n")

save(merged_samples, file = "merged_samples.RData")
cat("Merged samples saved as merged_samples.RData.\n")

saveRDS(merged_samples, file = "merged_samples.rds")
cat("Merged samples saved as merged_samples.rds.\n")
