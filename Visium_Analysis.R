# Load Required Libraries
library(Seurat)
library(SpatialExperiment)
library(ggplot2)
library(dplyr)
library(patchwork)

# Load Spatial Transcriptomic Data
data_dir <- "./data/"
visium_data <- Load10X_Spatial(data.dir = data_dir)

# Quality Control
visium_data["percent.mt"] <- PercentageFeatureSet(visium_data, pattern = "^MT-")
visium_data <- subset(visium_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalization
visium_data <- SCTransform(visium_data, assay = "Spatial", verbose = FALSE)

# Dimensionality Reduction
visium_data <- RunPCA(visium_data, verbose = FALSE)
visium_data <- RunUMAP(visium_data, dims = 1:30)
visium_data <- RunTSNE(visium_data, dims = 1:30)

# Clustering
visium_data <- FindNeighbors(visium_data, dims = 1:30)
visium_data <- FindClusters(visium_data, resolution = 0.5)

# Marker Gene Identification
markers <- FindAllMarkers(visium_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Visualization
SpatialDimPlot(visium_data, label = TRUE) + ggtitle("Spatial Clusters")
SpatialFeaturePlot(visium_data, features = c("MKI67", "CD3D"))

# Save Results
write.csv(markers, file = "./results/marker_genes.csv")
ggsave("./results/spatial_clusters.png")
ggsave("./results/gene_expression_maps.png")
