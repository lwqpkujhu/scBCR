library(dplyr)
library(tidyr)
library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)
library(scRepertoire)
library(ggplot2)
library(webr)
require(moonBook)
library(immunarch)
library(stringr)

####################
## scRepertoire ####
####################

bcr1.data <- read.csv("./21_Tube_1_results/21_Tube_1_results/vdj_b/filtered_contig_annotations.csv")
bcr2.data <- read.csv("./21_Tube_2_results/21_Tube_2_results/vdj_b/filtered_contig_annotations.csv")
bcr3.data <- read.csv("./21_Tube_3_results/21_Tube_3_results/vdj_b/filtered_contig_annotations.csv")
bcr4.data <- read.csv("./21_Tube_4_results/21_Tube_4_results/vdj_b/filtered_contig_annotations.csv")
bcr5.data <- read.csv("./21_Tube_5_results/21_Tube_5_results/vdj_b/filtered_contig_annotations.csv")
bcr6.data <- read.csv("./21_Tube_6_results/21_Tube_6_results/vdj_b/filtered_contig_annotations.csv")
bcr7.data <- read.csv("./21_Tube_7_results/21_Tube_7_results/vdj_b/filtered_contig_annotations.csv")
bcr11.data <- read.csv("./40_Tube_1_results/40_Tube_1_results/vdj_b/filtered_contig_annotations.csv")
bcr12.data <- read.csv("./40_Tube_2_results/40_Tube_2_results/vdj_b/filtered_contig_annotations.csv")
bcr13.data <- read.csv("./40_Tube_3_results/40_Tube_3_results/vdj_b/filtered_contig_annotations.csv")
bcr14.data <- read.csv("./40_Tube_4_results/40_Tube_4_results/vdj_b/filtered_contig_annotations.csv")
bcr15.data <- read.csv("./40_Tube_5_results/40_Tube_5_results/vdj_b/filtered_contig_annotations.csv")
bcr16.data <- read.csv("./40_Tube_6_results/40_Tube_6_results/vdj_b/filtered_contig_annotations.csv")
bcr17.data <- read.csv("./40_Tube_7_results/40_Tube_7_results/vdj_b/filtered_contig_annotations.csv")
bcr18.data <- read.csv("./40_Tube_8_results/40_Tube_8_results/vdj_b/filtered_contig_annotations.csv")

bcr1.data$barcode <- paste0(bcr1.data$barcode, "_1")
bcr2.data$barcode <- paste0(bcr2.data$barcode, "_2")
bcr3.data$barcode <- paste0(bcr3.data$barcode, "_3")
bcr4.data$barcode <- paste0(bcr4.data$barcode, "_4")
bcr5.data$barcode <- paste0(bcr5.data$barcode, "_5")
bcr6.data$barcode <- paste0(bcr6.data$barcode, "_6")
bcr7.data$barcode <- paste0(bcr7.data$barcode, "_7")
bcr11.data$barcode <- paste0(bcr11.data$barcode, "_8")
bcr12.data$barcode <- paste0(bcr12.data$barcode, "_9")
bcr13.data$barcode <- paste0(bcr13.data$barcode, "_10")
bcr14.data$barcode <- paste0(bcr14.data$barcode, "_11")
bcr15.data$barcode <- paste0(bcr15.data$barcode, "_12")
bcr16.data$barcode <- paste0(bcr16.data$barcode, "_13")
bcr17.data$barcode <- paste0(bcr17.data$barcode, "_14")
bcr18.data$barcode <- paste0(bcr18.data$barcode, "_15")

bcr.data <- bind_rows(bcr1.data, bcr2.data) %>% bind_rows(bcr3.data) %>% bind_rows(bcr4.data) %>% 
       bind_rows(bcr5.data) %>% bind_rows(bcr6.data) %>% bind_rows(bcr7.data) %>% 
       bind_rows(bcr11.data) %>% bind_rows(bcr12.data) %>% bind_rows(bcr13.data) %>%   
       bind_rows(bcr14.data) %>% bind_rows(bcr15.data) %>% bind_rows(bcr16.data) %>% 
       bind_rows(bcr17.data) %>% bind_rows(bcr18.data) 

write.csv(bcr.data, "BCR_combined_data.csv", row.names = F)

## reload data from file
bcr.data <- read.csv("BCR_combined_data.csv")
gex.combined <- readRDS("gex.combined.final.rds")

# assign BCR clones
bcr.combined <- combineBCR(bcr.data, samples = "Batch1", threshold = 0.85, removeNA = T, removeMulti = T)
clone_counts <- as.data.frame(table(bcr.combined[["Batch1"]][["CTstrict"]]))

# change barcode
newbarcodes <- stringr::str_replace(bcr.combined[["Batch1"]][["barcode"]], "^Batch1_", "")
bcr.combined[["Batch1"]][["barcode"]] <- newbarcodes

## generate dataframe
df_bcr <- as_tibble(bcr.combined[["Batch1"]])
df_gex <- FetchData(object = gex.combined, vars = c("day", "tissue", "sorting", "hash.ID", "adjuvant", "seurat_clusters"))
df_gex$barcode <- rownames(df_gex)

df_combined <- left_join(df_gex, df_bcr)

# backup and load data
write.table(df_combined, "df_combined_gex_bcr.txt", quote = F, row.names = F, sep = "\t")
df_combined <- read.delim("df_combined_gex_bcr.txt")

## clone analysis
df_clone <- df_combined %>% select("barcode", "CTstrict","day", "tissue", "sorting", "hash.ID", "adjuvant", "seurat_clusters") %>%
      filter(!is.na(CTstrict) & sorting == "Antigen Specific B Cells")

#############
## mouse 1 ##
#############
df_hash_1 <- df_clone %>% filter(hash.ID == "Hash-tag3")
df_hash_1_clone_type <- df_hash_1 %>% count(CTstrict)
tmp1 <- df_hash_1 %>% count(CTstrict, day) %>% count(CTstrict) %>% filter(n > 1)

df_hash_1_clone_type$CloneType <- "Singletons"
df_hash_1_clone_type[df_hash_1_clone_type$n>1,3] <- "Expanded"
df_hash_1_clone_type[df_hash_1_clone_type$CTstrict %in% tmp1$CTstrict,3] <- "Persisting"

df_hash_1 <- left_join(df_hash_1, df_hash_1_clone_type)

PD = df_hash_1 %>% group_by(CloneType, CTstrict, day) %>% summarise(n = n())
clonename <- PD %>% ungroup() %>% select(CTstrict) %>% distinct() 
clonename <- clonename %>% bind_cols(paste0("C", c(1:nrow(clonename)))) %>% rename(Clone = `...2`)
PD <- left_join(PD,clonename)

PD1 = PD %>% filter(day == 21)
PieDonut(PD1, aes(CloneType, Clone, count=n), title = " ", 
         showRatioThreshold = 0.05, labelposition=0, ratioByGroup = F,
         r0 = 0.4, r1 = 1.0, r2 = 1.2)
sum(PD1$n)

PD2 = PD %>% filter(day == 40)
PieDonut(PD2, aes(CloneType, Clone, count=n), title = " ", 
         showRatioThreshold = 0.05, labelposition=0, ratioByGroup = F,
         r0 = 0.4, r1 = 1.0, r2 = 1.2)
sum(PD2$n)


#############
## LiNA2  ###
#############
df_hash_1 <- df_clone %>% filter(hash.ID %in% c("Hash-tag10","Hash-tag9","Hash-tag8","Hash-tag7","Hash-tag6"))
df_hash_1_clone_type <- df_hash_1 %>% count(CTstrict)
tmp1 <- df_hash_1 %>% count(CTstrict, day) %>% count(CTstrict) %>% filter(n > 1)

df_hash_1_clone_type$CloneType <- "Singletons"
df_hash_1_clone_type[df_hash_1_clone_type$n>1,3] <- "Expanded"
df_hash_1_clone_type[df_hash_1_clone_type$CTstrict %in% tmp1$CTstrict,3] <- "Persisting"

df_hash_1 <- left_join(df_hash_1, df_hash_1_clone_type)

PD = df_hash_1 %>% group_by(CloneType, CTstrict, day) %>% summarise(n = n()) %>% 
  ungroup()  %>% group_by(CloneType, day) %>% arrange(desc(n),.by_group = T)
clonename <- PD %>% ungroup() %>% select(CTstrict) %>% distinct() 
temp2 <- paste0("C", c(1:nrow(clonename)))
temp3 <- temp2[order(temp2, decreasing = T)]
clonename <- clonename %>% bind_cols(temp3) %>% rename(Clone = `...2`)
PD <- left_join(PD,clonename)

PD1 = PD %>% filter(day == 21)
PieDonut(PD1, aes(CloneType, Clone, count=n), title = " ", 
         showRatioThreshold = 0.05, labelposition=0, ratioByGroup = F,
         r0 = 0.4, r1 = 1.0, r2 = 1.2)
sum(PD1$n)

PD2 = PD %>% filter(day == 40)
PieDonut(PD2, aes(CloneType, Clone, count=n), title = " ", 
         showRatioThreshold = 0.05, labelposition=0, ratioByGroup = F,
         r0 = 0.4, r1 = 1.0, r2 = 1.2)
sum(PD2$n)

df_LiNA2 <- df_hash_1

#############
## RSV Only  ###
#############
df_hash_1 <- df_clone %>% filter(hash.ID %in% c("Hash-tag1","Hash-tag2","Hash-tag3","Hash-tag4","Hash-tag5"))
df_hash_1_clone_type <- df_hash_1 %>% count(CTstrict)
tmp1 <- df_hash_1 %>% count(CTstrict, day) %>% count(CTstrict) %>% filter(n > 1)

df_hash_1_clone_type$CloneType <- "Singletons"
df_hash_1_clone_type[df_hash_1_clone_type$n>1,3] <- "Expanded"
df_hash_1_clone_type[df_hash_1_clone_type$CTstrict %in% tmp1$CTstrict,3] <- "Persisting"

df_hash_1 <- left_join(df_hash_1, df_hash_1_clone_type)

PD = df_hash_1 %>% group_by(CloneType, CTstrict, day) %>% summarise(n = n()) %>% 
  ungroup()  %>% group_by(CloneType, day) %>% arrange(desc(n),.by_group = T)
clonename <- PD %>% ungroup() %>% select(CTstrict) %>% distinct() 
temp2 <- paste0("C", c(1:nrow(clonename)))
temp3 <- temp2[order(temp2, decreasing = T)]
clonename <- clonename %>% bind_cols(temp3) %>% rename(Clone = `...2`)
PD <- left_join(PD,clonename)

PD1 = PD %>% filter(day == 21)
PieDonut(PD1, aes(CloneType, Clone, count=n), title = " ", 
         showRatioThreshold = 0.05, labelposition=0, ratioByGroup = F,
         r0 = 0.4, r1 = 1.0, r2 = 1.2)
sum(PD1$n)

PD2 = PD %>% filter(day == 40)
PieDonut(PD2, aes(CloneType, Clone, count=n), title = " ", 
         showRatioThreshold = 0.05, labelposition=0, ratioByGroup = F,
         r0 = 0.4, r1 = 1.0, r2 = 1.2)
sum(PD2$n)

df_RSV <- df_hash_1

### cluster vs cloneType
df_cluster <- bind_rows(df_LiNA2, df_RSV)
df_cluster_count <- df_cluster %>% count(seurat_clusters, CloneType) %>% 
  group_by(seurat_clusters) %>%  mutate(per = n*100/sum(n), total = sum(n))
df_cluster_table <- df_cluster_count %>% select(seurat_clusters, CloneType, per) %>%
     spread(CloneType, per) %>% left_join(distinct(df_cluster_count[,c(1,5)]))


df_cluster <- bind_rows(df_LiNA2, df_RSV)
df_cluster_count <- df_cluster %>% count(seurat_clusters, CloneType) %>% 
  group_by(CloneType) %>%  mutate(per = n*100/sum(n), total = sum(n))
df_cluster_table <- df_cluster_count %>% select(seurat_clusters, CloneType, per) %>%
  spread(seurat_clusters, per) %>% left_join(distinct(df_cluster_count[,c(2,5)]))


### map clone type to seruat object
df_cluster <- bind_rows(df_LiNA2, df_RSV) %>% select(barcode, CloneType)
df_combined <- left_join(df_combined, df_cluster)

gex.combined$BCR <- df_combined$CloneType
Idents(gex.combined) <- "BCR"
bcr.markers <- FindAllMarkers(gex.combined, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
write.table(bcr.markers, "BCR_CloneType_markers.txt", quote = F, row.names = F, sep = "\t")

## visualize
p1 <- DimPlot(gex.combined, reduction = "umap", split.by = "BCR")
p1
# expanded
FeaturePlot(gex.combined, features = c("BASP1", "JCHAIN", "MEF2B"), split.by = "BCR", max.cutoff = 3,
            cols = c("grey", "red"))
# Persisting
FeaturePlot(gex.combined, features = c("RGS13", "RGS1", "TXN"), split.by = "BCR", max.cutoff = 3,
            cols = c("grey", "red"))
# singleton
FeaturePlot(gex.combined, features = c("EGR1", "LY6D", "IER2"), split.by = "BCR", max.cutoff = 3,
            cols = c("grey", "red"))


## shared clones between adjuvant
#############
## D21  ###
#############
df_hash_1 <- df_clone %>% filter(day == 21)
df_hash_1_clone_type <- df_hash_1 %>% count(CTstrict)
tmp1 <- df_hash_1 %>% count(CTstrict, adjuvant) %>% count(CTstrict) %>% filter(n > 1)

df_hash_1_clone_type$CloneType <- "Singletons"
df_hash_1_clone_type[df_hash_1_clone_type$n>1,3] <- "Unique"
df_hash_1_clone_type[df_hash_1_clone_type$CTstrict %in% tmp1$CTstrict,3] <- "Shared"

df_hash_1 <- left_join(df_hash_1, df_hash_1_clone_type)

PD = df_hash_1 %>% group_by(CloneType, CTstrict, adjuvant) %>% summarise(n = n()) %>% 
     ungroup()  %>% group_by(CloneType, adjuvant) %>% arrange(desc(n),.by_group = T)
clonename <- PD %>% ungroup() %>% select(CTstrict) %>% distinct() 
temp2 <- paste0("C", c(1:nrow(clonename)))
temp3 <- temp2[order(temp2, decreasing = T)]
clonename <- clonename %>% bind_cols(temp3) %>% rename(Clone = `...2`)
PD <- left_join(PD,clonename)

PD1 = PD %>% filter(adjuvant == "RSV Only")
PieDonut(PD1, aes(CloneType, Clone, count=n), title = " ", 
         showRatioThreshold = 0.05, labelposition=0, ratioByGroup = F,
         r0 = 0.4, r1 = 1.0, r2 = 1.2)
sum(PD1$n)

PD2 = PD %>% filter(adjuvant == "RSV and LiNA2")
PieDonut(PD2, aes(CloneType, Clone, count=n), title = " ", 
         showRatioThreshold = 0.05, labelposition=0, ratioByGroup = F,
         r0 = 0.4, r1 = 1.0, r2 = 1.2)
sum(PD2$n)

df_21 <- df_hash_1


#############
## D40  ###
#############
df_hash_1 <- df_clone %>% filter(day == 40)
df_hash_1_clone_type <- df_hash_1 %>% count(CTstrict)
tmp1 <- df_hash_1 %>% count(CTstrict, adjuvant) %>% count(CTstrict) %>% filter(n > 1)

df_hash_1_clone_type$CloneType <- "Singletons"
df_hash_1_clone_type[df_hash_1_clone_type$n>1,3] <- "Unique"
df_hash_1_clone_type[df_hash_1_clone_type$CTstrict %in% tmp1$CTstrict,3] <- "Shared"

df_hash_1 <- left_join(df_hash_1, df_hash_1_clone_type)

PD = df_hash_1 %>% group_by(CloneType, CTstrict, adjuvant) %>% summarise(n = n()) %>% 
  ungroup()  %>% group_by(CloneType, adjuvant) %>% arrange(desc(n),.by_group = T)
clonename <- PD %>% ungroup() %>% select(CTstrict) %>% distinct() 
temp2 <- paste0("C", c(1:nrow(clonename)))
temp3 <- temp2[order(temp2, decreasing = T)]
clonename <- clonename %>% bind_cols(temp3) %>% rename(Clone = `...2`)
PD <- left_join(PD,clonename)

PD1 = PD %>% filter(adjuvant == "RSV Only")
PieDonut(PD1, aes(CloneType, Clone, count=n), title = " ", 
         showRatioThreshold = 0.05, labelposition=0, ratioByGroup = F,
         r0 = 0.4, r1 = 1.0, r2 = 1.2)
sum(PD1$n)

PD2 = PD %>% filter(adjuvant == "RSV and LiNA2")
PieDonut(PD2, aes(CloneType, Clone, count=n), title = " ", 
         showRatioThreshold = 0.05, labelposition=0, ratioByGroup = F,
         r0 = 0.4, r1 = 1.0, r2 = 1.2)
sum(PD2$n)

df_40 <- df_hash_1


### map clone type to seruat object
df_cluster <- bind_rows(df_21, df_40) %>% select(barcode, CloneType) %>% rename(CloneType2 = CloneType)
df_combined <- left_join(df_combined, df_cluster)

gex.combined$BCR2 <- df_combined$CloneType2
Idents(gex.combined) <- "BCR2"
bcr.markers <- FindAllMarkers(gex.combined, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
write.table(bcr.markers, "BCR_CloneType2_markers.txt", quote = F, row.names = F, sep = "\t")

## visualize
p1 <- DimPlot(gex.combined, reduction = "umap", split.by = "BCR2")
p1
# expanded
FeaturePlot(gex.combined, features = c("BASP1", "JCHAIN", "MEF2B"), split.by = "BCR", max.cutoff = 3,
            cols = c("grey", "red"))
# Persisting
FeaturePlot(gex.combined, features = c("RGS13", "RGS1", "TXN"), split.by = "BCR", max.cutoff = 3,
            cols = c("grey", "red"))
# singleton
FeaturePlot(gex.combined, features = c("EGR1", "LY6D", "IER2"), split.by = "BCR", max.cutoff = 3,
            cols = c("grey", "red"))

# backup and load data
write.table(df_combined, "df_combined_gex_bcr.txt", quote = F, row.names = F, sep = "\t")
df_combined <- read.delim("df_combined_gex_bcr.txt")

saveRDS(gex.combined , file = "./gex.combined.final.rds")
gex.combined <- readRDS("gex.combined.final.rds")
