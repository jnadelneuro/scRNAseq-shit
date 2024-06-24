library(Seurat)
library(Matrix)
library(doParallel)
library(foreach)
library(ggplot2)
library(viridis)

# Set up parallel processing
# num_cores <- detectCores()
# cl <- makeCluster(num_cores)
# registerDoParallel(cl)

# List all .csv.gz files in the directory
data_dir <- "C:\\Users\\jacob\\OneDrive\\Documents\\GitHub\\scRNAseq-shit\\midbrain-striatal-cortico-assembloids\\data"
save_dir <- "C:\\Users\\jacob\\OneDrive\\Documents\\GitHub\\scRNAseq-shit\\midbrain-striatal-cortico-assembloids\\"
csv_files <- list.files(data_dir, pattern = "*.csv.gz", full.names = TRUE)

# # # Function to read and process a single file
# # read_and_process <- function(file) {
# #   data <- read.csv(file, header = TRUE, row.names = 1)
# #   return(data)
# # }

# # # Read and process files in parallel
# # processed_data <- foreach(i = 1:length(csv_files), .packages = c("Matrix")) %dopar% {
# #   list(
# #     data = read_and_process(csv_files[i]),
# #     source = ifelse(grepl("standardorg", csv_files[i]), "standardorg",
# #                     ifelse(grepl("Lam", csv_files[i]), "lam",
# #                             ifelse(grepl("silk", csv_files[i]), "silk", NA)))
# #   )
# # }

# # Stop parallel processing
# stopCluster(cl)

# # # Create Seurat objects and add source information
# # seurat_objects <- lapply(processed_data, function(item) {
# #   seurat_obj <- CreateSeuratObject(counts = as(as.matrix(item$data), "sparseMatrix"))
# #   seurat_obj$source <- item$source
# #   seurat_obj[['percent.mt']] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
# #   # clean each session first
# #   seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# #   return(seurat_obj)
# # })

# # # Save each Seurat object with source information
# # for (i in seq_along(seurat_objects)) {
# #   saveRDS(seurat_objects[[i]], file = paste0(save_dir, "seurat_object_", i, "_", seurat_objects[[i]]$source[1], ".rds"))
# # }

# Load the Seurat objects back
# rds_files <- list.files(data_dir, pattern = "*.rds", full.names = TRUE)
seurat_object <- readRDS("C:\\Users\\jacob\\OneDrive\\Documents\\GitHub\\scRNAseq-shit\\midbrain-striatal-cortico-assembloids\\data\\GSE219245_merged_data_object.rds")



# Function to process each Seurat object separately
process_seurat_object <- function(seurat_object) {
  # Normalization
  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Find variable features
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
  
  # Scaling
  seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))
  
  # PCA
  seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
  
  # Add TH expression as metadata
  seurat_object$TH_expression <- FetchData(seurat_object, vars = "TH")[,1]
  
  # t-SNE/UMAP
  seurat_object <- RunTSNE(seurat_object, dims = 1:10)
  seurat_object <- RunUMAP(seurat_object, dims = 1:10)
  
  # Clustering
  seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
  
  # Custom clustering to ensure TH-expressing cells are in the same cluster
  seurat_object <- FindClusters(seurat_object, resolution = 0.5)
  clusters <- seurat_object$seurat_clusters
  
  # Force TH-expressing cells into the same cluster
  TH_cells <- WhichCells(seurat_object, expression = TH > 0)
  clusters[TH_cells] <- "TH_cluster"
  
  # Update cluster information
  seurat_object$seurat_clusters <- as.factor(clusters)
  
  return(seurat_object)
}

# Process each source separately
processed_seurat_objects <- lapply(source_seurat_objects, process_seurat_object)

# Save processed Seurat objects
for (i in seq_along(processed_seurat_objects)) {
  saveRDS(processed_seurat_objects[[i]], file = paste0(save_dir, "processed_seurat_object_", sources[i], ".rds"))
}

# Load necessary library for visualization
processed_rds_files <- list.files(save_dir, pattern = "processed_seurat_object_.*\\.rds", full.names = TRUE)
processed_seurat_objects <- lapply(processed_rds_files, readRDS)


library(ggplot2)

# Function to create and save plots for each source
create_and_save_plots <- function(seurat_object, markers, save_dir) {
  num_clusters <- length(unique(seurat_object$seurat_clusters))
  my_palette <- viridis_pal(option = "D")(num_clusters)
  
  # Function to check if a gene is present in the dataset and plot its expression across clusters
  plot_gene_expression <- function(gene, seurat_object) {
    if (gene %in% rownames(seurat_object)) {
      plot <- VlnPlot(seurat_object, features = gene, group.by = "seurat_clusters", cols = my_palette)
      return(plot)  # Return the plot
    } else {
      cat("Gene", gene, "not found in the dataset.\n")
      return(NULL)  # Return NULL if gene not found
    }
  }
  
  # Iterate through markers and plot expression of each gene across clusters
  for (gene in markers) {
    plot <- plot_gene_expression(gene, seurat_object)
    if (!is.null(plot)) {
      # Create subfolder based on the source
      source_subfolder <- paste0(save_dir, seurat_object$source[1], "/")
      if (!dir.exists(source_subfolder)) {
        dir.create(source_subfolder)
      }
      ggsave(filename = paste0(source_subfolder, gene, ".png"), plot = plot, width = 6, height = 4, units = "in", dpi = 300)
    }
  }
}

# Define marker gene sets
F0_markers <- c("HES1", "NES", "SOX2", "TOP2A", "CCNB2", "CENPF1")
F1_markers <- c("HES5", "SOX2", "VIM", "FABP7")
F2_markers <- c("SHH", "CORIN", "FOXA2", "LMX1A", "MSX1", "OTX2")
matDA_markers <- c("GIRK2", "CALB1", "DAT", "DCX", "SYT1", "STMN2", "PBX1", "NR4A2", "EN1", "DDC", "TH", "NR3C1")

markers <- c(F0_markers, F1_markers, F2_markers, matDA_markers)
markers <- unique(markers)

# Create and save plots for each processed Seurat object
for (i in seq_along(processed_seurat_objects)) {
  create_and_save_plots(processed_seurat_objects[[i]], markers, save_dir)
}

create_and_save_dimplots <- function(seurat_object, save_dir) {
  # Create subfolder based on the source
  source_subfolder <- paste0(save_dir, seurat_object$source[1], "/")
  if (!dir.exists(source_subfolder)) {
    dir.create(source_subfolder)
  }
  
  # Create DimPlot colored by cluster
  dim_plot <- DimPlot(seurat_object, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, cols = viridis::viridis(length(unique(seurat_object$seurat_clusters))))
  ggsave(filename = paste0(source_subfolder, "dimplot_by_cluster.png"), plot = dim_plot, width = 10, height = 8, units = "in", dpi = 300)
}

# Create and save DimPlots for each processed Seurat object
for (i in seq_along(processed_seurat_objects)) {
  create_and_save_dimplots(processed_seurat_objects[[i]], save_dir)
}
