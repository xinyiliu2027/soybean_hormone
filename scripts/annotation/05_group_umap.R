if (exists("all_samples", inherits = FALSE)) rm(all_samples)

library(Seurat)
library(ggplot2)
library(patchwork)

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

sample_names <- c("REF", "ABA1", "ABA2", "eBL1", "eBL2", "highTemp1", "highTemp2", 
                  "IAA1", "IAA2", "lowTemp1", "lowTemp2", "MEJA1", "MEJA2", 
                  "PEG1", "PEG2", "SA1", "SA2")

all_samples <- c(list(reference), hormone.list)
names(all_samples) <- sample_names

group_by_col <- "cell_type"
cat("Using", group_by_col, "for grouping\n")

dimplot_list <- list()
for(i in 1:length(all_samples)) {
  sample_name <- sample_names[i]
  
  if(!"umap" %in% names(all_samples[[i]]@reductions)) {
    cat("Warning: No UMAP found for", sample_name, ", skipping...\n")
    next
  }
  
  if(!group_by_col %in% colnames(all_samples[[i]]@meta.data)) {
    cat("Warning: No", group_by_col, "found for", sample_name, ", using seurat_clusters\n")
    current_group_col <- "seurat_clusters"
  } else {
    current_group_col <- group_by_col
  }
  
  p <- DimPlot(all_samples[[i]], 
               reduction = "umap", 
               group.by = current_group_col,
               label = TRUE, 
               label.size = 3,
               pt.size = 0.5,
               cols = custom_colors) +  
    ggtitle(sample_name) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14))
  
  dimplot_list[[i]] <- p
  
  cat("Created dimplot for", sample_name, "\n")
}

dimplot_list <- dimplot_list[!sapply(dimplot_list, is.null)]

if(length(dimplot_list) > 0) {
  
  n_plots <- length(dimplot_list)
  n_cols <- 4
  n_rows <- ceiling(n_plots / n_cols)
  
  combined_plot <- wrap_plots(dimplot_list, ncol = n_cols)
  
  ggsave("all_samples_with_ref_dimplot_custom_colors.pdf", combined_plot, 
         width = 5 * n_cols, height = 4 * n_rows)
  ggsave("all_samples_with_ref_dimplot_custom_colors.png", combined_plot, 
         width = 5 * n_cols, height = 4 * n_rows, dpi = 300)
  
  cat("Saved individual sample plots with custom colors\n")
}

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

merged_samples <- list()
merged_samples[["REF"]] <- all_samples[["REF"]]  

for(merged_name in names(merge_pairs)) {
  sample_pair <- merge_pairs[[merged_name]]
  
  if(all(sample_pair %in% names(all_samples))) {
    
    sample1 <- all_samples[[sample_pair[1]]]
    sample2 <- all_samples[[sample_pair[2]]]
    
    if("umap" %in% names(sample1@reductions) && "umap" %in% names(sample2@reductions)) {
      
      merged_sample <- merge(sample1, sample2, 
                             add.cell.ids = sample_pair,
                             project = merged_name)
      
      merged_samples[[merged_name]] <- merged_sample
      
      cat("Successfully merged", paste(sample_pair, collapse = " and "), "into", merged_name, "\n")
    } else {
      cat("Warning: UMAP not found for", paste(sample_pair, collapse = " or "), ", skipping merge\n")
    }
  } else {
    cat("Warning: Not all samples found for", merged_name, ", skipping merge\n")
  }
}

merged_dimplot_list <- list()
merged_sample_names <- names(merged_samples)

for(i in 1:length(merged_samples)) {
  sample_name <- merged_sample_names[i]
  
  if(!"umap" %in% names(merged_samples[[i]]@reductions)) {
    cat("Warning: No UMAP found for merged", sample_name, ", skipping...\n")
    next
  }
  
  if(!group_by_col %in% colnames(merged_samples[[i]]@meta.data)) {
    cat("Warning: No", group_by_col, "found for merged", sample_name, ", using seurat_clusters\n")
    current_group_col <- "seurat_clusters"
  } else {
    current_group_col <- group_by_col
  }
  
  p <- DimPlot(merged_samples[[i]], 
               reduction = "umap", 
               group.by = current_group_col,
               label = TRUE, 
               label.size = 3,
               pt.size = 0.5,
               cols = custom_colors) +  
    ggtitle(sample_name) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14))
  
  merged_dimplot_list[[i]] <- p
  
  cat("Created dimplot for merged", sample_name, "\n")
}

merged_dimplot_list <- merged_dimplot_list[!sapply(merged_dimplot_list, is.null)]

if(length(merged_dimplot_list) > 0) {
  
  n_plots <- length(merged_dimplot_list)
  n_cols <- 3  
  n_rows <- ceiling(n_plots / n_cols)
  
  merged_combined_plot <- wrap_plots(merged_dimplot_list, ncol = n_cols)
  
  ggsave("merged_samples_dimplot_custom_colors.pdf", merged_combined_plot, 
         width = 5 * n_cols, height = 4 * n_rows)
  ggsave("merged_samples_dimplot_custom_colors.png", merged_combined_plot, 
         width = 5 * n_cols, height = 4 * n_rows, dpi = 300)
  
  cat("Saved merged sample plots with custom colors\n")
}

if("REF" %in% names(all_samples)) {
  if("umap" %in% names(all_samples[["REF"]]@reductions)) {
    legend_plot <- DimPlot(all_samples[["REF"]], 
                           reduction = "umap", 
                           group.by = group_by_col,
                           label = TRUE, 
                           label.size = 3,
                           pt.size = 0.5,
                           cols = custom_colors) +
      ggtitle("Cell Type Legend") +
      theme(plot.title = element_text(hjust = 0.5, size = 16),
            legend.text = element_text(size = 12))
    
    ggsave("cell_type_legend.pdf", legend_plot, width = 10, height = 8)
    ggsave("cell_type_legend.png", legend_plot, width = 10, height = 8, dpi = 300)
    
    cat("Saved legend plot\n")
  }
}

cat("All plots completed!\n")
