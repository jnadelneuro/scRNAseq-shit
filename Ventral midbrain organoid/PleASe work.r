library(Seurat)
library(Matrix)
library(doParallel)
library(foreach)
library(ggplot2)
library(viridis)
# Set up parallel processing
num_cores <- detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# List all .csv.gz files in the directory
data_dir <- "D:\\Data_analysis\\idk\\GSE168323_RAW"
save_dir <- "D:\\Data_analysis\\idk\\"
csv_files <- list.files(data_dir, pattern = "*.csv.gz", full.names = TRUE)

# Function to read and process a single file
# read_and_process <- function(file) {
  data <- read.csv(file, header = TRUE, row.names = 1)
  # Process your data here (e.g., quality control, normalization)
  return(data)
}

# # Read and process files in parallel
processed_data <- foreach(i = 1:length(csv_files), .packages = c("Matrix")) %dopar% {
  read_and_process(csv_files[i])
}

# # Stop parallel processing
stopCluster(cl)

# # Merge processed data into a single Seurat object
seurat_objects <- lapply(processed_data, function(data) {
  CreateSeuratObject(counts = as(as.matrix(data), "sparseMatrix"))
})

for (i in seq_along(seurat_objects)) {
  saveRDS(seurat_objects[[i]], file = paste0(save_dir, "seurat_object_", i, ".rds"))
}


rds_files <- list.files(save_dir, pattern = "*.rds", full.names = TRUE)
seurat_objects <- lapply(rds_files, readRDS)
# Combine sparse matrices
combined_seurat_object <- merge(seurat_objects[[1]], y = seurat_objects[-1], add.cell.ids = paste0("Sample", seq_along(seurat_objects)))

combined_seurat_object[["percent.mt"]] <- PercentageFeatureSet(combined_seurat_object, pattern = "^MT-")
# Create a Seurat object from the combined sparse matrix
# combined_seurat_object <- CreateSeuratObject(counts = combined_sparse_matrix)
# Save the processed Seurat object
saveRDS(combined_seurat_object, file = "processed_seurat_object.rds")

combined_seurat_object <- readRDS("processed_seurat_object.rds")

# Load necessary library for visualization
library(ggplot2)

combined_seurat_object <- subset(combined_seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# # Create a histogram of percent.mt
histogram_plot <- ggplot(combined_seurat_object@meta.data, aes(x = percent.mt)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of Percent Mitochondrial Content",
       x = "Percent Mitochondrial Content",
       y = "Frequency")

# Print the histogram plot
# print(histogram_plot)

# Create a violin plot of percent.mt
violin_plot <- VlnPlot(combined_seurat_object, features = "percent.mt") +
  theme_minimal() +
  labs(title = "Violin Plot of Percent Mitochondrial Content")

# Print the violin plot
# print(violin_plot)

# Quality control
combined_seurat_object <- subset(combined_seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalization
combined_seurat_object <- NormalizeData(combined_seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
combined_seurat_object <- FindVariableFeatures(combined_seurat_object, selection.method = "vst", nfeatures = 2000)

# Scaling
combined_seurat_object <- ScaleData(combined_seurat_object, features = rownames(combined_seurat_object))


# PCA
combined_seurat_object <- RunPCA(combined_seurat_object, features = VariableFeatures(object = combined_seurat_object))

# t-SNE/UMAP
combined_seurat_object <- RunTSNE(combined_seurat_object, dims = 1:10)
# or
combined_seurat_object <- RunUMAP(combined_seurat_object, dims = 1:10)

# Clustering
combined_seurat_object <- FindNeighbors(combined_seurat_object, dims = 1:10)
combined_seurat_object <- FindClusters(combined_seurat_object, resolution = 0.5)


combined_seurat_object <- JoinLayers(combined_seurat_object)
saveRDS(combined_seurat_object, file = "processed_seurat_object1.rds")
# Find markers for each cluster

combined_seurat_object <- readRDS("processed_seurat_object1.rds")

# cluster_markers <- FindAllMarkers(combined_seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Define marker gene sets
F0_markers <- c("HES1", "NES", "SOX2", "TOP2A", "CCNB2", "CENPF1")
F1_markers <- c("HES5", "SOX2", "VIM", "FABP7")
F2_markers <- c("SHH", "CORIN", "FOXA2", "LMX1A", "MSX1", "OTX2")
matDA_markers <- c("GIRK2", "CALB1", "DAT", "DCX", "SYT1", "STMN2", "PBX1", "NR4A2", "EN1", "DDC", "TH")

markers <- c(F0_markers, F1_markers, F2_markers, matDA_markers)

markers <- unique(markers)
# Annotate cells with marker gene sets
gene_plots <- list()

num_clusters <- length(unique(combined_seurat_object$seurat_clusters))
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

# Iterate through genes_of_interest and plot expression of each gene across clusters
for (gene in markers) {
  gene_plots[[gene]] <- plot_gene_expression(gene, combined_seurat_object)
  ggsave(filename = paste0(gene, ".png"), plot = gene_plots[[gene]], width = 6, height = 4, units = "in", dpi = 300)
}

# if ("GIRK2" %in% row.names(combined_seurat_object)) {
#   print("Gene is present in the dataset.")
# } else {
#   print("Gene is not present in the dataset.")
# }

# Annotate clusters (manually or based on known markers)
# Example: combined_seurat_object$celltype <- "annotated_cell_type"
q <- 1