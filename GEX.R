library(dplyr)
library(Seurat)
library(patchwork)

####################
### LN antigen #####
####################
# Load in the UMI matrix
data_dir <- './aggr_LN_antigen/'
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
gex.data <- Read10X(data.dir = data_dir)

# Setup Seurat object
LN.antigen <- CreateSeuratObject(counts = gex.data, project = "LN.antigen", min.cells = 3, min.features = 200)

## Data QC human MT and mouse mt
LN.antigen[["percent.mt"]] <- PercentageFeatureSet(LN.antigen, pattern = "^mt-")
# Visualize QC metrics as a violin plot
VlnPlot(LN.antigen, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(LN.antigen, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(LN.antigen, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
## filter data
LN.antigen <- subset(LN.antigen, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)

# Removing BCR-genes
LN.antigen <- LN.antigen[!grepl("^IG[HKL]V", rownames(LN.antigen), ignore.case = TRUE), ]
LN.antigen <- LN.antigen[!grepl("^IG[HKL]J", rownames(LN.antigen), ignore.case = TRUE), ]
LN.antigen <- LN.antigen[!grepl("^IG[KL]C", rownames(LN.antigen), ignore.case = TRUE), ]
LN.antigen <- LN.antigen[!grepl("^IGH[ADEGM]", rownames(LN.antigen), ignore.case = TRUE), ]

# Normalize RNA data with log normalization
LN.antigen <- NormalizeData(LN.antigen, normalization.method = "LogNormalize", scale.factor = 10000)
# Find and scale variable features
LN.antigen <- FindVariableFeatures(LN.antigen, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(LN.antigen), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(LN.antigen)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
## Scaling the data
all.genes <- rownames(LN.antigen)
LN.antigen <- ScaleData(LN.antigen, features = all.genes)

## Perform linear dimensional reduction
LN.antigen <- RunPCA(LN.antigen, features = VariableFeatures(object = LN.antigen))
print(LN.antigen[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(LN.antigen, dims = 1:2, reduction = "pca")
DimPlot(LN.antigen, reduction = "pca")
DimHeatmap(LN.antigen, dims = 1:15, cells = 500, balanced = TRUE)

## Determine the ‘dimensionality’ of the dataset
# LN.antigen <- JackStraw(LN.antigen, num.replicate = 100)
# LN.antigen <- ScoreJackStraw(LN.antigen, dims = 1:20)
# JackStrawPlot(LN.antigen, dims = 1:15)
ElbowPlot(LN.antigen)

## Cluster the cells
# select number of PCs used
LN.antigen <- FindNeighbors(LN.antigen, dims = 1:15)
# select resolution setting this parameter between 0.4-1.2 typically returns good results
LN.antigen <- FindClusters(LN.antigen, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(LN.antigen), 5)

LN.antigen <- RunUMAP(LN.antigen, dims = 1:15)
DimPlot(LN.antigen, reduction = "umap")
## save projects
saveRDS(LN.antigen, file = "./LN.antigen.rds")

### find cluster markers
# find all markers of cluster 2
# cluster2.markers <- FindMarkers(LN.antigen, ident.1 = 2, min.pct = 0.25)
# head(cluster2.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster5.markers <- FindMarkers(LN.antigen, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive ones
LN.antigen.markers <- FindAllMarkers(LN.antigen, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
LN.antigen.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
# top 10 marker heatmap
LN.antigen.markers %>% filter(cluster %in% as.character(c(0:9))) %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(LN.antigen, features = top5$gene) + NoLegend()

## assign cell type name
new.cluster.ids <- c(0:17)
names(new.cluster.ids) <- levels(LN.antigen)
LN.antigen <- RenameIdents(LN.antigen, new.cluster.ids)
DimPlot(LN.antigen, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(LN.antigen, file = "./LN.antigen.rds")
