###SINGE CELL TRANSCRIPTOMICS####
##Najneen Rejwana
getwd()
setwd("G:/Rstudio/single_trans/")



#Load packages
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)

#Install Seurat
install.packages(c("ggplot2", "dplyr", "Matrix", "cowplot", "Rcpp"))
install.packages("Seurat")
library(Seurat)

#Install Rtools 4.4
browseURL("https://cran.r-project.org/bin/windows/Rtools/")
system("g++ --version")

#Install dplyr
install.packages("dplyr", dependencies = TRUE)
library(dplyr)

#Install Matrix
install.packages("Matrix")
library(Matrix)


#Install ggplot2 
install.packages("ggplot2")
library(ggplot2)


str(seurat)     # Check the structure of your Seurat object
assay(seurat)   # List all assays in the Seurat object
names(seurat@assays$RNA@layers)   # Check layers available in the RNA assay
rna_counts <- seurat@assays$RNA@layers$counts  #access the count data
rna_data <- seurat@assays$RNA@layers$data   #access the normalized data



counts <- readMM("matrix.mtx.gz")
barcodes <- read.table("barcodes.tsv.gz", stringsAsFactors=F)[,1]
features <- read.csv("features.tsv.gz", stringsAsFactors=F, sep="\t", header=F)
rownames(counts) <- make.unique(features[,2])
colnames(counts) <- barcodes
counts


seurat <- CreateSeuratObject(counts, project="hb1")
#Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#Warning: Data is of class dgTMatrix. Coercing to dgCMatrix
#Just ignore this error
# Check the Seurat object summary
seurat



seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT[-\\.]")
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)


library(patchwork)
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2



seurat <- subset(seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5)
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, nfeatures = 3000)

top_features <- head(VariableFeatures(seurat), 20)
plot1 <- VariableFeaturePlot(seurat)
print(plot1)

plot2 <- LabelPoints(plot = plot1, points = top_features, repel = TRUE, xnudge = 0, ynudge = 0, max.overlaps = 40, width = 10, height = 8, dpi = 300)
print(plot2)



#### Data scaling
seurat <- ScaleData(seurat)

## Centering and scaling data matrix
seurat <- RunPCA(seurat, npcs = 50)
# Plotting the PCA visualization
pca_plot <- DimPlot(seurat, reduction = "pca", dims = c(1, 2))  # First and second PCA components
print(pca_plot)



##### Elbow plot to choose number of PCs
ElbowPlot(seurat, ndims = ncol(Embeddings(seurat, "pca")))



PCHeatmap(seurat, dims = 1:15, cells = 500, balanced = TRUE, ncol = 4)

# Generate PC heatmaps and arrange them in 3 columns, 2 rows
pc_heatmap <- PCHeatmap(seurat, dims = 1:6, cells = 500, balanced = TRUE)  # Showing 6 PCs
pc_heatmap + plot_layout(ncol = 3, nrow = 2)




#### Clustering based on selected PCs
seurat <- RunTSNE(seurat, dims = 1:15)
seurat <- RunUMAP(seurat, dims = 1:15)

## plots
plot1 <- TSNEPlot(seurat)
plot2 <- UMAPPlot(seurat)
plot1 + plot2




## checking and plotting canonical markers of the cell types
plot1 <- FeaturePlot(seurat, c("MKI67","NES","DCX","FOXG1","DLX2","EMX1","OTX2","LHX9","TFAP2A"),
                     ncol=3, reduction = "tsne")
plot2 <- FeaturePlot(seurat, c("MKI67","NES","DCX","FOXG1","DLX2","EMX1","OTX2","LHX9","TFAP2A"),
                     ncol=3, reduction = "umap")
plot1 / plot2



#Diffuse Large B-Cell Lymphoma (DLBCL) markers
plot1 <- FeaturePlot(seurat, c("CD19", "CD20", "CD22", "CD79a", "BCL6", "MYC", "MUM1"),
                     ncol=3, reduction = "tsne")
plot2 <- FeaturePlot(seurat, c("CD19", "CD20", "CD22", "CD79a", "BCL6", "MYC", "MUM1"),
                     ncol=3, reduction = "umap")
plot1 / plot2


#Mantle Cell Lymphoma (MCL)
plot1 <- FeaturePlot(seurat, c("CD5", "Cyclin D1", "SOX11"),
                     ncol=3, reduction = "tsne")
plot2 <- FeaturePlot(seurat, c("CD5", "Cyclin D1", "SOX11"),
                     ncol=3, reduction = "umap")
plot1 / plot2





## K- means clustering
seurat <- FindNeighbors(seurat, dims = 1:15)
seurat <- FindClusters(seurat, resolution = 1)


## tSNE & UMAP results

plot1 <- DimPlot(seurat, reduction = "tsne", label = TRUE)
plot2 <- DimPlot(seurat, reduction = "umap", label = TRUE)
plot1 + plot2



# Prior knowledge of markers (from database/ Publications). These markers here are just used as an example.
colon_markers <- c("CEA", "CA 19-9", "mSEPT9", "FIT", "ctDNA", "KRAS", "BRAF", "TP53", "MSI-H", "HSH2D", "CDX2", "EGFR", "NRAS", "dMMR", "HER2", "CMS classification")
DoHeatmap(seurat, features = ct_markers) + NoLegend()



## Identify markers for each of the clusters.
cl_markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))


# Heatmap for visualizing identified markers
DoHeatmap(seurat, features = cl_markers$gene) + NoLegend()
