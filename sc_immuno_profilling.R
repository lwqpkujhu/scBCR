library(dplyr)
library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)

source("load_data.R")

gex.combined <- readRDS("gex.combined.final.rds")

# Removing BCR-genes
gex.combined <- gex.combined[!grepl("^IG[HKL]V", rownames(gex.combined), ignore.case = TRUE), ]
gex.combined <- gex.combined[!grepl("^IG[HKL]J", rownames(gex.combined), ignore.case = TRUE), ]
gex.combined <- gex.combined[!grepl("^IG[KL]C", rownames(gex.combined), ignore.case = TRUE), ]
gex.combined <- gex.combined[!grepl("^IGH[ADEGM]", rownames(gex.combined), ignore.case = TRUE), ]

### annotate cell types
gex.combined <- RunAzimuth(gex.combined, reference = "pbmcref")
p1 <- DimPlot(gex.combined, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) 
p1
p2 <- DimPlot(gex.combined, group.by = "predicted.celltype.l1", label = TRUE, label.size = 3)
p2
### Idents(gex.combined) <- "predicted.celltype.l2"
# B cell marker plot
FeaturePlot(gex.combined, features = c("MS4A1", "CD19"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD79A", "CD79B"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD40", "CD27"), min.cutoff = NA)
# CD4 T cell
FeaturePlot(gex.combined, features = c("IL7R", "CCR7"), min.cutoff = NA)
# CD8 T cell
FeaturePlot(gex.combined, features = c("CD8A", "CD8B"), min.cutoff = NA)
# T cell markers
FeaturePlot(gex.combined, features = c("CD3E", "CD3G"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD247", "CD3D"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD44", "PTPRC"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD197", "S100A4"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD69", "IL2RA"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("SELL", "IL2RA"), min.cutoff = NA)
#PTPRC is CD45
# DC cell
FeaturePlot(gex.combined, features = c("FCER1A", "CD83"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("BST2", "CST3"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("ITGAX", "XCR1"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("ITGAE", "ITGAM"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD24", "CD24A"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CLEC9A", "CD209"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CADM1", "LY75"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD207", "CD80"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD86", "CD40"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CX3CR1", "CLEC9A"), min.cutoff = NA)
# Mono cell
FeaturePlot(gex.combined, features = c("CD14", "LYZ"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("FCGR3A", "MS4A7"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CSF1R", "CCR2"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("ITGA2", "SELL"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("ITGAX", "ITGAM"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("SPN", "CX3CR1"), min.cutoff = NA)
# NK cell
FeaturePlot(gex.combined, features = c("ITGAM", "NKG7"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD27", "IL2RB"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("KLRB1", "NCR1"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("ITGA2", "KLRK1"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD3E", "CD247"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD3D", "CD69"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("GZMA", "GZMB"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("GZMC", "GZMD"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("GZME", "GZMF"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("GZMG", "GZMH"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("GZMI", "GZMJ"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("GZMK", "PRF1"), min.cutoff = NA)

## feature selection
gex.combined <- NormalizeData(gex.combined, normalization.method = "LogNormalize", scale.factor = 10000)
gex.combined <- FindVariableFeatures(gex.combined, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(gex.combined)
gex.combined <- ScaleData(gex.combined, features = all.genes)
gex.combined <- RunPCA(gex.combined, features = VariableFeatures(object = gex.combined))
ElbowPlot(gex.combined)

gex.combined <- FindNeighbors(gex.combined, dims = 1:20)
gex.combined <- FindClusters(gex.combined, resolution = 0.22)
gex.combined <- RunUMAP(gex.combined, dims = 1:20)

# Visualization
DimPlot(gex.combined, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(gex.combined, reduction = "umap", split.by = "adjuvant", label = T) + NoLegend()
DimPlot(gex.combined, reduction = "umap", split.by = "day", label = T) + NoLegend()
DimPlot(gex.combined, reduction = "umap", split.by = "tissue", label = T) + NoLegend()
DimPlot(gex.combined, reduction = "umap", split.by = "sorting", label = T) + NoLegend()
df1 <- ncol(subset(x = gex.combined, subset = adjuvant == "RSV and LiNA2" & sorting == "All B Cells"))
df2 <- ncol(subset(x = gex.combined, subset = adjuvant == "RSV and LiNA2" & sorting == "Antigen Specific B Cells"))
df3 <- ncol(subset(x = gex.combined, subset = adjuvant == "RSV Only" & sorting == "All B Cells"))
df4 <- ncol(subset(x = gex.combined, subset = adjuvant == "RSV Only" & sorting == "Antigen Specific B Cells"))
df12 <- subset(x = gex.combined, subset = adjuvant == "RSV and LiNA2")
df34 <- subset(x = gex.combined, subset = adjuvant == "RSV Only")
DimPlot(df12, reduction = "umap", split.by = "sorting", label = T) + NoLegend()
DimPlot(df34, reduction = "umap", split.by = "sorting", label = T) + NoLegend()

## save and load object
saveRDS(gex.combined, file = "./gex.combined.final.rds")
gex.combined <- readRDS("gex.combined.final.rds")
Idents(gex.combined) <- "seurat_clusters"

### calculate % for each cluster
LiNA2_total <- ncol(subset(x = gex.combined, subset = adjuvant == "RSV and LiNA2"))
RSVOnly_total <- ncol(subset(x = gex.combined, subset = adjuvant == "RSV Only"))
df <- data.frame(matrix(ncol = 2, nrow = 14))
colnames(df) <- c("LiNA2", "RSVOnly")
df[15,] <- c(LiNA2_total, RSVOnly_total)
for (i in c(0:13)) {
  n1 <- ncol(subset(x = gex.combined, subset = adjuvant == "RSV and LiNA2" & seurat_clusters == as.character(i)))*100/LiNA2_total
  n2 <- ncol(subset(x = gex.combined, subset = adjuvant == "RSV Only" & seurat_clusters == as.character(i)))*100/RSVOnly_total
  df[i+1,] <- c(n1,n2)
}
df <- format(df, scientific = FALSE)

# Find cluster marker genes
DefaultAssay(gex.combined) <- "RNA"
cluster2.markers <- FindMarkers(gex.combined, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers)
# plot cluster marker genes
FeaturePlot(gex.combined, features = c("TCEA1", "PCMTD1", "MYBL1", "MYBL1"), min.cutoff = "q9")
VlnPlot(gex.combined, features = c("TCEA1", "PCMTD1", "MYBL1", "MYBL1"))


# find markers for every cluster compared to all remaining cells, report only the positive
# ones
cluster.markers <- FindAllMarkers(gex.combined, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers <- cluster.markers %>% filter(p_val_adj < 0.01)
write.table(cluster.markers, "./cluster.markers.txt", quote = F, row.names = F, sep = "\t")
cluster.markers2 <- cluster.markers %>% filter(p_val_adj < 0.01 & (avg_log2FC > 1 | avg_log2FC < -1 ))
write.table(cluster.markers2, "./cluster.markers_qc.txt", quote = F, row.names = F, sep = "\t")

## reload data to avoid conflict of umap
gex.combined <- readRDS("gex.combined.final.rds")
cluster.markers <- read.delim("cluster.markers.txt")
Idents(gex.combined) <- "seurat_clusters"
cluster.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10

DoHeatmap(gex.combined, features = top10$gene) + NoLegend()

DoHeatmap(gex.combined, features = top10$gene, cells = 1:3000, size = 4,
          angle = 90) + NoLegend()

DoHeatmap(gex.combined, features = top10$gene, cells = 1:300, size = 4,
          angle = 90) 

## find markers
cluster.markers %>%
  filter(cluster == 2) %>%
  slice_max(n = 4, order_by = avg_log2FC)
FeaturePlot(gex.combined, features = c("BASP1", "JCHAIN", "MEF2B", "RGS13"), min.cutoff = "q9")

cluster.markers %>%
  filter(cluster == 4) %>%
  slice_max(n = 4, order_by = avg_log2FC)
FeaturePlot(gex.combined, features = c("PCLAF", "MKI67", "STMN1", "MKI67"), min.cutoff = "q9")

cluster.markers %>%
  filter(cluster == 5) %>%
  slice_max(n = 4, order_by = avg_log2FC)
FeaturePlot(gex.combined, features = c("PLAC8", "LY6D", "HEXB", "PSAP"), min.cutoff = "q9")

cluster.markers %>%
  filter(cluster == 8) %>%
  slice_max(n = 4, order_by = avg_log2FC)
FeaturePlot(gex.combined, features = c("EGR1", "VIM", "GIMAP4", "CAPG"), min.cutoff = "q9")

cluster.markers %>%
  filter(cluster == 9) %>%
  slice_max(n = 4, order_by = avg_log2FC)
FeaturePlot(gex.combined, features = c("FCER1G", "TYROBP", "CST3", "LST1"), min.cutoff = "q9")


### find differential expressed gene between conditions
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
### extract cluster cells
cluster2.cells <- subset(gex.combined, idents = "2")
Idents(cluster2.cells) <- "seurat_clusters"
avg.cluster2.cells <- as.data.frame(log1p(AverageExpression(cluster2.cells, verbose = FALSE)$RNA))
avg.cluster2.cells$gene <- rownames(avg.cluster2.cells)

cluster4.cells <- subset(gex.combined, idents = "4")
Idents(cluster4.cells) <- "seurat_clusters"
avg.cluster4.cells <- as.data.frame(log1p(AverageExpression(cluster4.cells, verbose = FALSE)$RNA))
avg.cluster4.cells$gene <- rownames(avg.cluster4.cells)

### find differential expressed gene between conditions
Idents(gex.combined) <- "seurat_clusters"
gex.combined$celltype.stim <- paste(Idents(gex.combined), gex.combined$adjuvant, sep = "_")
gex.combined$celltype <- Idents(gex.combined)
Idents(gex.combined) <- "celltype.stim"
LiNA2.response.2 <- FindMarkers(gex.combined, ident.1 = "2_RSV and LiNA2", ident.2 = "2_RSV Only", verbose = FALSE)
head(LiNA2.response.2, n = 15)
LiNA2.response.4 <- FindMarkers(gex.combined, ident.1 = "4_RSV and LiNA2", ident.2 = "4_RSV Only", verbose = FALSE)
head(LiNA2.response.4, n = 15)
## plot top genes
genes.to.label.2 = rownames(LiNA2.response.2)[1:12]
genes.to.label.4 = rownames(LiNA2.response.4)[1:12]
p1 <- ggplot(avg.cluster2.cells, aes("RSV Only", "RSV and LiNA2")) + geom_point() + ggtitle("Cluster 2 Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label.2, repel = TRUE)
p2 <- ggplot(avg.cluster4.cells, aes("RSV Only", "RSV and LiNA2")) + geom_point() + ggtitle("Cluster 4 Cells")
p2 <- LabelPoints(plot = p2, points = genes.to.label.4, repel = TRUE)
p1 + p2

FeaturePlot(gex.combined, features = rownames(LiNA2.response.0)[1:3], split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

plots <- VlnPlot(gex.combined, features = rownames(LiNA2.response.0)[1:3], split.by = "orig.ident", group.by = "celltype",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)


###############
## plot metadata
###############

FeaturePlot(gex.combined, features = c("percent.mt"), max.cutoff = 20,
            cols = c("grey", "red"))

### calculate % for each cluster and for each hashtag
hash <- "Hash-tag1"
num <- 21
hash_total <- ncol(subset(x = gex.combined, subset = hash.ID == hash & day == num))
df <- data.frame(matrix(ncol = 1, nrow = 14))
colnames(df) <- hash
df[15,] <- hash_total
for (i in c(0:13)) {
  n1 <- ncol(subset(x = gex.combined, subset = seurat_clusters == as.character(i) & hash.ID == hash & day == num))*100/hash_total
  df[i+1,] <- round(n1, digits = 2)
}
df <- format(df, scientific = FALSE)


### calculate B cell percentage
Idents(gex.combined) <- "adjuvant"
LiNA2_total <- ncol(subset(x = gex.combined, subset = adjuvant == "RSV and LiNA2"))
RSVOnly_total <- ncol(subset(x = gex.combined, subset = adjuvant == "RSV Only"))
df <- data.frame(matrix(ncol = 2, nrow = 14))
colnames(df) <- c("LiNA2", "RSVOnly")
df[15,] <- c(LiNA2_total, RSVOnly_total)
for (i in c(0:13)) {
  n1 <- ncol(subset(x = gex.combined, subset = adjuvant == "RSV and LiNA2" & seurat_clusters == as.character(i)))*100/LiNA2_total
  n2 <- ncol(subset(x = gex.combined, subset = adjuvant == "RSV Only" & seurat_clusters == as.character(i)))*100/RSVOnly_total
  df[i+1,] <- c(n1,n2)
}
df <- format(df, scientific = FALSE)


### B Cell subtype makers
Idents(object = gex.combined) <- "seurat_clusters"
#check gene availability
test <- gex.combined@assays[["RNA"]]@counts@Dimnames[[1]]
test[str_detect(test, "CD20")]
# Naive B cell marker plot
FeaturePlot(gex.combined, features = c("CD19", "CD220"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("IgD", "CD20"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD22", "CD27"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD19", "CD22", "CD27"), min.cutoff = NA)
# Activated B cell marker plot
FeaturePlot(gex.combined, features = c("CD19", "CD220"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD25", "CD30"), min.cutoff = NA)
# Transitional B cell marker plot
FeaturePlot(gex.combined, features = c("CD93", "CD220"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD24", "BR3"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("IgM", "TACI"), min.cutoff = NA)
# Plasma cell marker plot
FeaturePlot(gex.combined, features = c("CD27", "CD138"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("BCMA", "CXCR4"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD126", "CD184"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD320", "BLIMP1"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("IRF4", "XBP1"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("IL6", "IGD"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CXCR4","IRF4", "XBP1"), min.cutoff = NA)
# Follicular B cell
FeaturePlot(gex.combined, features = c("CD19", "CD220"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD38", "IGD"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD23", "CD22"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CXCR5", "PAX5"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD19","CD38","CD22","CXCR5", "PAX5"), min.cutoff = NA)
# Germinal Center B cell
FeaturePlot(gex.combined, features = c("GL7", "CD95","PNA"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("IGM", "IGD"), min.cutoff = NA)
# Marginal Zone B cell
FeaturePlot(gex.combined, features = c("CD19", "CD220"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD9", "IGM"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD23", "CD22"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD21", "CD35"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD1D", "PAX5"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD19","CD9","CD22","CD1D", "PAX5"), min.cutoff = NA)
# Memory B cell
FeaturePlot(gex.combined, features = c("CD19", "CD220"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD80", "CD73"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD27", "CD38"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD273", "CD84"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD86"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD19","CD80","CD27", "CD38", "CD84", "CD86"), min.cutoff = NA)
# Regulatory B cells
FeaturePlot(gex.combined, features = c("CD19", "CD21"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("IGM", "IGD"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD40", "CD1D"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("TIM1", "PAX5"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD62L"), min.cutoff = NA)
FeaturePlot(gex.combined, features = c("CD19","CD40", "CD1D", "PAX5"), min.cutoff = NA)


### Cell-Cycle Scoring and Regression
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
gex.combined <- CellCycleScoring(gex.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# Visualize the distribution of cell cycle markers across
RidgePlot(gex.combined, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely
gex.combined <- RunPCA(gex.combined, features = c(s.genes, g2m.genes))
DimPlot(gex.combined)
DimPlot(gex.combined, reduction = "umap", split.by = "adjuvant", label = F) 
DimPlot(gex.combined, reduction = "umap", split.by = "day", label = F) 
DimPlot(gex.combined, reduction = "umap", split.by = "tissue", label = F) 
DimPlot(gex.combined, reduction = "umap", split.by = "sorting", label = F) 
DimPlot(gex.combined, reduction = "umap", split.by = "sorting", label = F) 
DimPlot(gex.combined, reduction = "umap", split.by = "seurat_clusters", label = F, ncol = 4) 
DimPlot(gex.combined, reduction = "umap", split.by = "BCR", label = F) 
DimPlot(gex.combined, reduction = "umap", split.by = "BCR2", label = F) 
