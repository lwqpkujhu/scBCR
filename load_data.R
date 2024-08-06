library(dplyr)
library(Seurat)
library(patchwork)

# Load in the UMI matrix
gex1.data <- Read10X(data.dir = "./21_Tube_1_results/21_Tube_1_results/count/sample_filtered_feature_bc_matrix/")
gex2.data <- Read10X(data.dir = "./21_Tube_2_results/21_Tube_2_results/count/sample_filtered_feature_bc_matrix/")
gex3.data <- Read10X(data.dir = "./21_Tube_3_results/21_Tube_3_results/count/sample_filtered_feature_bc_matrix/")
gex4.data <- Read10X(data.dir = "./21_Tube_4_results/21_Tube_4_results/count/sample_filtered_feature_bc_matrix/")
gex5.data <- Read10X(data.dir = "./21_Tube_5_results/21_Tube_5_results/count/sample_filtered_feature_bc_matrix/")
gex6.data <- Read10X(data.dir = "./21_Tube_6_results/21_Tube_6_results/count/sample_filtered_feature_bc_matrix/")
gex7.data <- Read10X(data.dir = "./21_Tube_7_results/21_Tube_7_results/count/sample_filtered_feature_bc_matrix/")
gex11.data <- Read10X(data.dir = "./40_Tube_1_results/40_Tube_1_results/count/sample_filtered_feature_bc_matrix/")
gex12.data <- Read10X(data.dir = "./40_Tube_2_results/40_Tube_2_results/count/sample_filtered_feature_bc_matrix/")
gex13.data <- Read10X(data.dir = "./40_Tube_3_results/40_Tube_3_results/count/sample_filtered_feature_bc_matrix/")
gex14.data <- Read10X(data.dir = "./40_Tube_4_results/40_Tube_4_results/count/sample_filtered_feature_bc_matrix/")
gex15.data <- Read10X(data.dir = "./40_Tube_5_results/40_Tube_5_results/count/sample_filtered_feature_bc_matrix/")
gex16.data <- Read10X(data.dir = "./40_Tube_6_results/40_Tube_6_results/count/sample_filtered_feature_bc_matrix/")
gex17.data <- Read10X(data.dir = "./40_Tube_7_results/40_Tube_7_results/count/sample_filtered_feature_bc_matrix/")
gex18.data <- Read10X(data.dir = "./40_Tube_8_results/40_Tube_8_results/count/sample_filtered_feature_bc_matrix/")

# Select cell barcodes detected by both RNA and HTO In the example datasets we have already
gex1.umis <- gex1.data$`Gene Expression`
gex1.htos <- gex1.data$`Antibody Capture`
joint.bcs <- intersect(colnames(gex1.umis), colnames(gex1.htos))
gex1.umis <- gex1.umis[, joint.bcs]
gex1.htos <- as.matrix(gex1.htos[, joint.bcs])
gex2.umis <- gex2.data$`Gene Expression`
gex2.htos <- gex2.data$`Antibody Capture`
joint.bcs <- intersect(colnames(gex2.umis), colnames(gex2.htos))
gex2.umis <- gex2.umis[, joint.bcs]
gex2.htos <- as.matrix(gex2.htos[, joint.bcs])
gex3.umis <- gex3.data$`Gene Expression`
gex3.htos <- gex3.data$`Antibody Capture`
joint.bcs <- intersect(colnames(gex3.umis), colnames(gex3.htos))
gex3.umis <- gex3.umis[, joint.bcs]
gex3.htos <- as.matrix(gex3.htos[, joint.bcs])
gex4.umis <- gex4.data$`Gene Expression`
gex4.htos <- gex4.data$`Antibody Capture`
joint.bcs <- intersect(colnames(gex4.umis), colnames(gex4.htos))
gex4.umis <- gex4.umis[, joint.bcs]
gex4.htos <- as.matrix(gex4.htos[, joint.bcs])
gex5.umis <- gex5.data$`Gene Expression`
gex5.htos <- gex5.data$`Antibody Capture`
joint.bcs <- intersect(colnames(gex5.umis), colnames(gex5.htos))
gex5.umis <- gex5.umis[, joint.bcs]
gex5.htos <- as.matrix(gex5.htos[, joint.bcs])
gex6.umis <- gex6.data$`Gene Expression`
gex6.htos <- gex6.data$`Antibody Capture`
joint.bcs <- intersect(colnames(gex6.umis), colnames(gex6.htos))
gex6.umis <- gex6.umis[, joint.bcs]
gex6.htos <- as.matrix(gex6.htos[, joint.bcs])
gex7.umis <- gex7.data$`Gene Expression`
gex7.htos <- gex7.data$`Antibody Capture`
joint.bcs <- intersect(colnames(gex7.umis), colnames(gex7.htos))
gex7.umis <- gex7.umis[, joint.bcs]
gex7.htos <- as.matrix(gex7.htos[, joint.bcs])
gex11.umis <- gex11.data$`Gene Expression`
gex11.htos <- gex11.data$`Antibody Capture`
joint.bcs <- intersect(colnames(gex11.umis), colnames(gex11.htos))
gex11.umis <- gex11.umis[, joint.bcs]
gex11.htos <- as.matrix(gex11.htos[, joint.bcs])
gex12.umis <- gex12.data$`Gene Expression`
gex12.htos <- gex12.data$`Antibody Capture`
joint.bcs <- intersect(colnames(gex12.umis), colnames(gex12.htos))
gex12.umis <- gex12.umis[, joint.bcs]
gex12.htos <- as.matrix(gex12.htos[, joint.bcs])
gex13.umis <- gex13.data$`Gene Expression`
gex13.htos <- gex13.data$`Antibody Capture`
joint.bcs <- intersect(colnames(gex13.umis), colnames(gex13.htos))
gex13.umis <- gex13.umis[, joint.bcs]
gex13.htos <- as.matrix(gex13.htos[, joint.bcs])
gex14.umis <- gex14.data$`Gene Expression`
gex14.htos <- gex14.data$`Antibody Capture`
joint.bcs <- intersect(colnames(gex14.umis), colnames(gex14.htos))
gex14.umis <- gex14.umis[, joint.bcs]
gex14.htos <- as.matrix(gex14.htos[, joint.bcs])
gex15.umis <- gex15.data$`Gene Expression`
gex15.htos <- gex15.data$`Antibody Capture`
joint.bcs <- intersect(colnames(gex15.umis), colnames(gex15.htos))
gex15.umis <- gex15.umis[, joint.bcs]
gex15.htos <- as.matrix(gex15.htos[, joint.bcs])
gex16.umis <- gex16.data$`Gene Expression`
gex16.htos <- gex16.data$`Antibody Capture`
joint.bcs <- intersect(colnames(gex16.umis), colnames(gex16.htos))
gex16.umis <- gex16.umis[, joint.bcs]
gex16.htos <- as.matrix(gex16.htos[, joint.bcs])
gex17.umis <- gex17.data$`Gene Expression`
gex17.htos <- gex17.data$`Antibody Capture`
joint.bcs <- intersect(colnames(gex17.umis), colnames(gex17.htos))
gex17.umis <- gex17.umis[, joint.bcs]
gex17.htos <- as.matrix(gex17.htos[, joint.bcs])
gex18.umis <- gex18.data$`Gene Expression`
gex18.htos <- gex18.data$`Antibody Capture`
joint.bcs <- intersect(colnames(gex18.umis), colnames(gex18.htos))
gex18.umis <- gex18.umis[, joint.bcs]
gex18.htos <- as.matrix(gex18.htos[, joint.bcs])

# Setup Seurat object
gex1.hashtag <- CreateSeuratObject(counts = gex1.umis, min.cells = 3, min.features = 300)
gex1.htos <- as.matrix(gex1.htos[, Cells(gex1.hashtag)])
gex1.hashtag[["HTO"]] <- CreateAssayObject(counts = gex1.htos)
gex2.hashtag <- CreateSeuratObject(counts = gex2.umis, min.cells = 3, min.features = 300)
gex2.htos <- as.matrix(gex2.htos[, Cells(gex2.hashtag)])
gex2.hashtag[["HTO"]] <- CreateAssayObject(counts = gex2.htos)
gex3.hashtag <- CreateSeuratObject(counts = gex3.umis, min.cells = 3, min.features = 300)
gex3.htos <- as.matrix(gex3.htos[, Cells(gex3.hashtag)])
gex3.hashtag[["HTO"]] <- CreateAssayObject(counts = gex3.htos)
gex4.hashtag <- CreateSeuratObject(counts = gex4.umis, min.cells = 3, min.features = 300)
gex4.htos <- as.matrix(gex4.htos[, Cells(gex4.hashtag)])
gex4.hashtag[["HTO"]] <- CreateAssayObject(counts = gex4.htos)
gex5.hashtag <- CreateSeuratObject(counts = gex5.umis, min.cells = 3, min.features = 300)
gex5.htos <- as.matrix(gex5.htos[, Cells(gex5.hashtag)])
gex5.hashtag[["HTO"]] <- CreateAssayObject(counts = gex5.htos)
gex6.hashtag <- CreateSeuratObject(counts = gex6.umis, min.cells = 3, min.features = 300)
gex6.htos <- as.matrix(gex6.htos[, Cells(gex6.hashtag)])
gex6.hashtag[["HTO"]] <- CreateAssayObject(counts = gex6.htos)
gex7.hashtag <- CreateSeuratObject(counts = gex7.umis, min.cells = 3, min.features = 300)
gex7.htos <- as.matrix(gex7.htos[, Cells(gex7.hashtag)])
gex7.hashtag[["HTO"]] <- CreateAssayObject(counts = gex7.htos)
gex11.hashtag <- CreateSeuratObject(counts = gex11.umis, min.cells = 3, min.features = 300)
gex11.htos <- as.matrix(gex11.htos[, Cells(gex11.hashtag)])
gex11.hashtag[["HTO"]] <- CreateAssayObject(counts = gex11.htos)
gex12.hashtag <- CreateSeuratObject(counts = gex12.umis, min.cells = 3, min.features = 300)
gex12.htos <- as.matrix(gex12.htos[, Cells(gex12.hashtag)])
gex12.hashtag[["HTO"]] <- CreateAssayObject(counts = gex12.htos)
gex13.hashtag <- CreateSeuratObject(counts = gex13.umis, min.cells = 3, min.features = 300)
gex13.htos <- as.matrix(gex13.htos[, Cells(gex13.hashtag)])
gex13.hashtag[["HTO"]] <- CreateAssayObject(counts = gex13.htos)
gex14.hashtag <- CreateSeuratObject(counts = gex14.umis, min.cells = 3, min.features = 300)
gex14.htos <- as.matrix(gex14.htos[, Cells(gex14.hashtag)])
gex14.hashtag[["HTO"]] <- CreateAssayObject(counts = gex14.htos)
gex15.hashtag <- CreateSeuratObject(counts = gex15.umis, min.cells = 3, min.features = 300)
gex15.htos <- as.matrix(gex15.htos[, Cells(gex15.hashtag)])
gex15.hashtag[["HTO"]] <- CreateAssayObject(counts = gex15.htos)
gex16.hashtag <- CreateSeuratObject(counts = gex16.umis, min.cells = 3, min.features = 300)
gex16.htos <- as.matrix(gex16.htos[, Cells(gex16.hashtag)])
gex16.hashtag[["HTO"]] <- CreateAssayObject(counts = gex16.htos)
gex17.hashtag <- CreateSeuratObject(counts = gex17.umis, min.cells = 3, min.features = 300)
gex17.htos <- as.matrix(gex17.htos[, Cells(gex17.hashtag)])
gex17.hashtag[["HTO"]] <- CreateAssayObject(counts = gex17.htos)
gex18.hashtag <- CreateSeuratObject(counts = gex18.umis, min.cells = 3, min.features = 300)
gex18.htos <- as.matrix(gex18.htos[, Cells(gex18.hashtag)])
gex18.hashtag[["HTO"]] <- CreateAssayObject(counts = gex18.htos)

## assign meta data
gex1.hashtag$day <- "21"
gex2.hashtag$day <- "21"
gex3.hashtag$day <- "21"
gex4.hashtag$day <- "21"
gex5.hashtag$day <- "21"
gex6.hashtag$day <- "21"
gex7.hashtag$day <- "21"
gex11.hashtag$day <- "40"
gex12.hashtag$day <- "40"
gex13.hashtag$day <- "40"
gex14.hashtag$day <- "40"
gex15.hashtag$day <- "40"
gex16.hashtag$day <- "40"
gex17.hashtag$day <- "40"
gex18.hashtag$day <- "40"

gex1.hashtag$tissue <- "LN"
gex2.hashtag$tissue <- "LN"
gex3.hashtag$tissue <- "LN"
gex4.hashtag$tissue <- "LN"
gex5.hashtag$tissue <- "Spleen"
gex6.hashtag$tissue <- "Spleen"
gex7.hashtag$tissue <- "Spleen"
gex11.hashtag$tissue <- "LN"
gex12.hashtag$tissue <- "LN"
gex13.hashtag$tissue <- "LN"
gex14.hashtag$tissue <- "LN"
gex15.hashtag$tissue <- "Spleen"
gex16.hashtag$tissue <- "Spleen"
gex17.hashtag$tissue <- "Spleen"
gex18.hashtag$tissue <- "Spleen"

gex1.hashtag$sorting <- "All B Cells"
gex2.hashtag$sorting <- "All B Cells"
gex3.hashtag$sorting <- "Antigen Specific B Cells"
gex4.hashtag$sorting <- "Antigen Specific B Cells"
gex5.hashtag$sorting <- "All B Cells"
gex6.hashtag$sorting <- "Antigen Specific B Cells"
gex7.hashtag$sorting <- "Antigen Specific B Cells"
gex11.hashtag$sorting <- "All B Cells"
gex12.hashtag$sorting <- "All B Cells"
gex13.hashtag$sorting <- "Antigen Specific B Cells"
gex14.hashtag$sorting <- "Antigen Specific B Cells"
gex15.hashtag$sorting <- "All B Cells"
gex16.hashtag$sorting <- "All B Cells"
gex17.hashtag$sorting <- "Antigen Specific B Cells"
gex18.hashtag$sorting <- "Antigen Specific B Cells"

#############
##combine####
#############

gex.combined <- merge(gex1.hashtag, y = list(gex2.hashtag, gex3.hashtag,gex4.hashtag,gex5.hashtag,gex6.hashtag,gex7.hashtag,gex11.hashtag,gex12.hashtag,gex13.hashtag,gex14.hashtag,gex15.hashtag,gex16.hashtag,gex17.hashtag,gex18.hashtag), project = "RSVAdjuvant")

## Data QC human MT and mouse mt
gex.combined[["percent.mt"]] <- PercentageFeatureSet(gex.combined, pattern = "^mt-")
# Visualize QC metrics as a violin plot
VlnPlot(gex.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(gex.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(gex.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
## filter data
gex.combined <- subset(gex.combined, subset = nFeature_RNA > 300 & nFeature_RNA < 7500 & percent.mt < 20)
gex.combined 
# Normalize RNA data with log normalization
gex.combined <- NormalizeData(gex.combined)
# Find and scale variable features
gex.combined <- FindVariableFeatures(gex.combined, selection.method = "mean.var.plot")
gex.combined <- ScaleData(gex.combined, features = VariableFeatures(gex.combined))
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
gex.combined <- NormalizeData(gex.combined, assay = "HTO", normalization.method = "CLR")
# If you have a very large dataset we suggest using kfunc = "clara". This is a k-medoid
# clustering function for large applications You can also play with additional parameters
gex.combined <- HTODemux(gex.combined, assay = "HTO", positive.quantile = 0.9999999)
# Global classification results
table(gex.combined$HTO_classification.global)
# Classification results
table(gex.combined$hash.ID)

## remove singlet and doublet
set1 <- subset(x = gex.combined, subset = hash.ID =="Hash-tag5" | hash.ID =="Hash-tag4" |hash.ID =="Hash-tag3" |hash.ID =="Hash-tag2" |hash.ID =="Hash-tag1")
table(set1$hash.ID)
set2 <- subset(x = gex.combined, subset = hash.ID =="Hash-tag6" | hash.ID =="Hash-tag7" |hash.ID =="Hash-tag8" |hash.ID =="Hash-tag9" |hash.ID =="Hash-tag10")
table(set2$hash.ID)

set1$adjuvant <- "RSV Only"
set2$adjuvant <- "RSV and LiNA2"

gex.final <- merge(set1,set2)

saveRDS(gex.final, file = "./gex.combined.final.rds")
