library(alakazam)
library(shazam)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(ggpubr)
library(webr)
require(moonBook)

#### prepare data ###
db_obs1 <- read.delim("./SHM/21_Tube_3_results/heavy_parse-select.tsv")
db_obs2 <- read.delim("./SHM/21_Tube_4_results/heavy_parse-select.tsv")
db_obs3 <- read.delim("./SHM/21_Tube_6_results/heavy_parse-select.tsv")
db_obs4 <- read.delim("./SHM/21_Tube_7_results/heavy_parse-select.tsv")
db_obs5 <- read.delim("./SHM/40_Tube_3_results/heavy_parse-select.tsv")
db_obs6 <- read.delim("./SHM/40_Tube_4_results/heavy_parse-select.tsv")
db_obs7 <- read.delim("./SHM/40_Tube_7_results/heavy_parse-select.tsv")
db_obs8 <- read.delim("./SHM/40_Tube_8_results/heavy_parse-select.tsv")

db_obs1$cell_id <- paste0(db_obs1$cell_id, "_3")
db_obs2$cell_id <- paste0(db_obs2$cell_id, "_4")
db_obs3$cell_id <- paste0(db_obs3$cell_id, "_6")
db_obs4$cell_id <- paste0(db_obs4$cell_id, "_7")
db_obs5$cell_id <- paste0(db_obs5$cell_id, "_10")
db_obs6$cell_id <- paste0(db_obs6$cell_id, "_11")
db_obs7$cell_id <- paste0(db_obs7$cell_id, "_14")
db_obs8$cell_id <- paste0(db_obs8$cell_id, "_15")

SHM.data <- bind_rows(db_obs1, db_obs2) %>% bind_rows(db_obs3) %>% bind_rows(db_obs4) %>% bind_rows(db_obs5) %>% 
  bind_rows(db_obs6) %>% bind_rows(db_obs7) %>% bind_rows(db_obs8)

write.table(SHM.data, "./SHM/combined_heavy_parse-select.tsv", quote = F, row.names = F, sep = "\t")

#######
## need to run SHM_per_sample.sh first, then run the above code to combined files, then run SHM_combine_sample.sh
#######

##################################
## Immcantation Portal ####
##################################

df1 <- read.delim("./SHM/by_clones_combined_heavy_parse-select_clone-pass_germ-pass.tsv")

# Calculate R and S mutation counts
db_obs <- observedMutations(df1, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=FALSE, 
                            nproc=1)
db_obs %>% 
  select(sequence_id, starts_with("mu_count_")) %>%
  head(n=4)

# Calculate R and S mutation frequencies
db_obs <- observedMutations(db_obs, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=TRUE, 
                            nproc=1)
# Show new mutation frequency columns
db_obs %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)

# Calculate combined R and S mutation frequencies
db_obs1 <- observedMutations(df1, sequenceColumn="sequence_alignment",
                             germlineColumn="germline_alignment_d_mask",
                             regionDefinition=NULL,
                             frequency=TRUE, 
                             combine=TRUE,
                             nproc=1)
# Show new mutation frequency columns
db_obs1 %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)

## We can plot the mutation frequencies a explore differences between samples or isotypes.
g1 <- ggplot(db_obs1, aes(x=c_call, y=mu_freq, fill=c_call)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot() + ylim(0,0.2)
plot(g1)

######################################################################################################

df1 <- read.delim("./SHM/combined_heavy_parse-select_clone-pass_germ-pass.tsv")

# Calculate R and S mutation counts
db_obs <- observedMutations(df1, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=FALSE, 
                            nproc=1)
db_obs %>% 
  select(sequence_id, starts_with("mu_count_")) %>%
  head(n=4)

# Calculate R and S mutation frequencies
db_obs <- observedMutations(db_obs, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=TRUE, 
                            nproc=1)
# Show new mutation frequency columns
db_obs %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)

# Calculate combined R and S mutation frequencies
db_obs1 <- observedMutations(df1, sequenceColumn="sequence_alignment",
                             germlineColumn="germline_alignment_d_mask",
                             regionDefinition=NULL,
                             frequency=TRUE, 
                             combine=TRUE,
                             nproc=1)
# Show new mutation frequency columns
db_obs1 %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)

db_obs1[str_detect(db_obs1[,18], "IGHG"),18] <- "IGHG"
db_obs1 <- db_obs1 %>% filter(c_call != "")

## We can plot the mutation frequencies a explore differences between samples or isotypes.
g1 <- ggplot(db_obs1, aes(x=c_call, y=mu_freq, fill=c_call)) +
  theme_bw() + ggtitle("") +
  xlab("Isotype") + ylab("SHM frequency") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot() + ylim(0,0.2)
plot(g1)

g2 <- ggboxplot(db_obs1, x = "c_call", y = "mu_freq", color = "c_call", palette = c("#377EB8", "#FF7F00", "#E41A1C", "#4DAF4A", "#984EA3", "#E5C494", "#FFD92F"),
                order = c("IGHA", "IGHD","IGHE", "IGHG", "IGHM")) +
      theme_bw() + ggtitle("") + xlab("Isotype") + ylab("SHM frequency") + labs(color ='Isotype') + ylim(0,0.2) + 
      geom_hline(yintercept = mean(db_obs1$mu_freq), linetype = 2)+ # Add horizontal line at base mean
      stat_compare_means(method = "anova", label.y = 0.2)+        # Add global annova p-value
      stat_compare_means(label = "p.signif", method = "t.test", label.y = 0.15,
                         ref.group = ".all.")                      # Pairwise comparison against all
plot(g2)


#############################################
#############################################
## combined with meta data
db_obs1 <- db_obs1 %>% group_by(cell_id) %>% slice(1:1) %>% ungroup()

df_meta <- read.delim("./df_combined_gex_bcr.txt")
df_meta <- df_meta %>% select(day, tissue, sorting, adjuvant, seurat_clusters, barcode, CTstrict, hash.ID) %>%
  rename(cell_id = barcode)

df_combined <- left_join(df_meta, db_obs1)


gg_adjuvant <- ggboxplot(df_combined, x = "adjuvant", y = "mu_freq", color = "adjuvant", palette = c("#377EB8", "#FF7F00", "#E41A1C", "#4DAF4A", "#984EA3", "#E5C494", "#FFD92F"),
                order = c("RSV Only", "RSV and LiNA2")) +
  theme_bw() + ggtitle("") + xlab("Adjuvant") + ylab("SHM frequency") + labs(color ='Adjuvant') + ylim(0,0.2) + 
  stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5, label.y = 0.15)                 # Pairwise comparison
plot(gg_adjuvant)

gg_adjuvant <- ggboxplot(df_combined, x = "day", y = "mu_freq", color = "day", palette = c("#377EB8", "#FF7F00", "#E41A1C", "#4DAF4A", "#984EA3", "#E5C494", "#FFD92F"),
                         order = c("21", "40")) +
  theme_bw() + ggtitle("") + xlab("Day") + ylab("SHM frequency") + labs(color ='day') + ylim(0,0.2) + 
  stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5, label.y = 0.15)                 # Pairwise comparison
plot(gg_adjuvant)

gg_adjuvant <- ggboxplot(df_combined, x = "tissue", y = "mu_freq", color = "tissue", palette = c("#377EB8", "#FF7F00", "#E41A1C", "#4DAF4A", "#984EA3", "#E5C494", "#FFD92F"),
                         order = c("Spleen", "LN")) +
  theme_bw() + ggtitle("") + xlab("tissue") + ylab("SHM frequency") + labs(color ='tissue') + ylim(0,0.2) + 
  stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5, label.y = 0.15)                 # Pairwise comparison
plot(gg_adjuvant)

df_combined2 <- df_combined %>% filter(c_call !="" & day == 21)
gg_adjuvant <- ggboxplot(df_combined2, x = "hash.ID", y = "mu_freq", color = "hash.ID", 
                         order = c("Hash-tag1", "Hash-tag2", "Hash-tag3", "Hash-tag4", "Hash-tag5", "Hash-tag6", "Hash-tag7", "Hash-tag8", "Hash-tag9", "Hash-tag10")) +
  theme_bw() + ggtitle("") + xlab("hash.ID") + ylab("SHM frequency") + labs(color ='hash.ID') + ylim(0,0.2) + 
  geom_hline(yintercept = mean(df_combined2$mu_freq), linetype = 2) +      # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 0.2)+                # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test", label.y = 0.15,
                     ref.group = ".all.")                             # Pairwise comparison against all
plot(gg_adjuvant)


df_combined2 <- df_combined %>% filter(c_call !="")
gg_adjuvant2 <- ggboxplot(df_combined2, x = "c_call", y = "mu_freq", color = "adjuvant", palette = c("#377EB8", "#FF7F00", "#E41A1C", "#4DAF4A", "#984EA3", "#E5C494", "#FFD92F"),
                          order = c("IGHA", "IGHD","IGHE", "IGHG", "IGHM")) +
  theme_bw() + ggtitle("") + xlab("Isotype") + ylab("SHM frequency") + labs(color ='Adjuvant') + ylim(0,0.2) + 
  stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5, label.y = 0.15, aes(group = adjuvant))                 # Pairwise comparison
plot(gg_adjuvant2)


df_combined2 <- df_combined %>% filter(!(seurat_clusters %in% c(6,9,12,13)))
gg_adjuvant2 <- ggboxplot(df_combined2, x = "seurat_clusters", y = "mu_freq", color = "adjuvant", palette = c("#377EB8", "#FF7F00", "#E41A1C", "#4DAF4A", "#984EA3", "#E5C494", "#FFD92F"),
                          order = c("0", "1","2", "3", "4", "5", "6","7", "8", "9", "10", "11","12", "13")) +
  theme_bw() + ggtitle("") + xlab("Seurat Cell Clusters") + ylab("SHM frequency") + labs(color ='Adjuvant') + ylim(0,0.2) + 
  stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5, label.y = 0.15, aes(group = adjuvant))                 # Pairwise comparison
plot(gg_adjuvant2)


df_combined2 <- df_combined %>% filter(!(seurat_clusters %in% c(13))) %>% filter(mu_freq != "")
gg_adjuvant2 <- ggboxplot(df_combined2, x = "seurat_clusters", y = "mu_freq", color = "seurat_clusters",order = c("0", "1","2", "3", "4", "5", "6","7", "8", "9", "10", "11","12", "13")) +
  theme_bw() + ggtitle("") + xlab("Seurat Cell Clusters") + ylab("SHM frequency") + labs(color ='Seurat Cell Clusters') + ylim(0,0.2) + 
  geom_hline(yintercept = mean(df_combined2$mu_freq), linetype = 2) +      # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 0.2)+                # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test", label.y = 0.15,
                     ref.group = ".all.")                             # Pairwise comparison against all
plot(gg_adjuvant2)

############################################################################################
############################################################################################
############################################################################################
## clone analysis
df_clone <- read.delim("df_combined_gex_bcr.txt") %>% select("barcode", "CloneType", "CloneType2") %>% rename(cell_id = barcode )


########################
### SHM ################
########################
df_combined2 <- df_combined %>% filter(CloneType !="")
my_comparisons <- list( c("Singletons", "Expanded"), c("Expanded", "Persisting"), c("Singletons", "Persisting") )
gg_adjuvant2 <- ggboxplot(df_combined2, x = "CloneType", y = "mu_freq", color = "CloneType", palette = c("#377EB8", "#FF7F00", "#E41A1C", "#4DAF4A", "#984EA3", "#E5C494", "#FFD92F"),
                          order = c("Singletons", "Expanded", "Persisting")) +
  theme_bw() + ggtitle("") + xlab("\nClone Type") + ylab("SHM frequency") + labs(color ='Clone Type') + ylim(0,0.2) + 
  stat_compare_means(label = "p.signif", method = "t.test",  
                     comparisons = my_comparisons, label.y = c(0.15, 0.16, 0.18))  # Pairwise comparison
plot(gg_adjuvant2)

df_combined2 <- df_combined %>% filter(CloneType !="")
gg_adjuvant2 <- ggboxplot(df_combined2, x = "CloneType", y = "mu_freq", color = "adjuvant", palette = c("#377EB8", "#FF7F00", "#E41A1C", "#4DAF4A", "#984EA3", "#E5C494", "#FFD92F"),
                          order = c("Singletons", "Expanded", "Persisting")) +
  theme_bw() + ggtitle("") + xlab("\nClone Type") + ylab("SHM frequency") + labs(color ='Adjuvant') + ylim(0,0.2) + 
  stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5, label.y = 0.15, aes(group = adjuvant))                 # Pairwise comparison
plot(gg_adjuvant2)



df_combined2 <- df_combined %>% filter(CloneType2 !="")
my_comparisons <- list( c("Singletons", "Unique"), c("Unique", "Shared"), c("Singletons", "Shared") )
gg_adjuvant2 <- ggboxplot(df_combined2, x = "CloneType2", y = "mu_freq", color = "CloneType2", palette = c("#377EB8", "#FF7F00", "#E41A1C", "#4DAF4A", "#984EA3", "#E5C494", "#FFD92F"),
                          order = c("Singletons", "Unique", "Shared")) +
  theme_bw() + ggtitle("") + xlab("\nClone Type") + ylab("SHM frequency") + labs(color ='Clone Type') + ylim(0,0.2) + 
  stat_compare_means(label = "p.signif", method = "t.test",  
                     comparisons = my_comparisons, label.y = c(0.15, 0.16, 0.18))  # Pairwise comparison
plot(gg_adjuvant2)


df_combined2 <- df_combined %>% filter(CloneType2 !="")
gg_adjuvant2 <- ggboxplot(df_combined2, x = "CloneType2", y = "mu_freq", color = "adjuvant", palette = c("#377EB8", "#FF7F00", "#E41A1C", "#4DAF4A", "#984EA3", "#E5C494", "#FFD92F"),
                          order = c("Singletons", "Unique", "Shared")) +
  theme_bw() + ggtitle("") + xlab("\nClone Type") + ylab("SHM frequency") + labs(color ='Adjuvant') + ylim(0,0.2) + 
  stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5, label.y = 0.15, aes(group = adjuvant))                 # Pairwise comparison
plot(gg_adjuvant2)


######################################################################################################
## CDR3                 ##############################################################################
######################################################################################################

df_combined$CDR3_length <- nchar(df_combined$cdr3)/3

gg_adjuvant <- ggboxplot(df_combined, x = "adjuvant", y = "CDR3_length", color = "adjuvant", palette = c("#377EB8", "#FF7F00", "#E41A1C", "#4DAF4A", "#984EA3", "#E5C494", "#FFD92F"),
                         order = c("RSV Only", "RSV and LiNA2")) +
  theme_bw() + ggtitle("") + xlab("\nAdjuvant") + ylab("CDR3 AA Length") + labs(color ='Adjuvant') + ylim(0,25) + 
  stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5, label.y = 22)                 # Pairwise comparison
plot(gg_adjuvant)


df_combined2 <- df_combined %>% filter(CDR3_length > 0 & adjuvant == "RSV Only")
p <- ggplot(df_combined2, aes(x=CDR3_length) ) + theme_bw() +
  geom_histogram( binwidth=1, fill="#377EB8", color="#e9ecef", alpha=0.9) + xlim(4,21) +
  xlab("CDR3 Length") + ylab("Count")
p

df_combined2 <- df_combined %>% filter(CDR3_length > 0 & adjuvant == "RSV and LiNA2")
p <- ggplot(df_combined2, aes(x=CDR3_length) ) + theme_bw() +
  geom_histogram( binwidth=1, fill="#FF7F00", color="#e9ecef", alpha=0.9) + xlim(4,21) +
  xlab("CDR3 Length") + ylab("Count")
p

gg_adjuvant <- ggboxplot(df_combined, x = "day", y = "CDR3_length", color = "day", palette = c("#377EB8", "#FF7F00", "#E41A1C", "#4DAF4A", "#984EA3", "#E5C494", "#FFD92F"),
                         order = c("21", "40")) +
  theme_bw() + ggtitle("") + xlab("Day") + ylab("CDR3 AA Length") + labs(color ='day') + ylim(0,25) + 
  stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5, label.y = 22)                 # Pairwise comparison
plot(gg_adjuvant)

gg_adjuvant <- ggboxplot(df_combined, x = "tissue", y = "CDR3_length", color = "tissue", palette = c("#377EB8", "#FF7F00", "#E41A1C", "#4DAF4A", "#984EA3", "#E5C494", "#FFD92F"),
                         order = c("Spleen", "LN")) +
  theme_bw() + ggtitle("") + xlab("tissue") + ylab("CDR3 AA Length") + labs(color ='tissue') + ylim(0,25) + 
  stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5, label.y = 22)                 # Pairwise comparison
plot(gg_adjuvant)


### clonetype
df_combined2 <- df_combined %>% filter(CloneType !="")
my_comparisons <- list( c("Singletons", "Expanded"), c("Expanded", "Persisting"), c("Singletons", "Persisting") )
gg_adjuvant2 <- ggboxplot(df_combined2, x = "CloneType", y = "CDR3_length", color = "CloneType", palette = c("#377EB8", "#FF7F00", "#E41A1C", "#4DAF4A", "#984EA3", "#E5C494", "#FFD92F"),
                          order = c("Singletons", "Expanded", "Persisting")) +
  theme_bw() + ggtitle("") + xlab("\nClone Type") + ylab("CDR3 length") + labs(color ='Clone Type') + ylim(0,25) + 
  stat_compare_means(label = "p.signif", method = "t.test",  
                     comparisons = my_comparisons, label.y = c(21, 22, 24))  # Pairwise comparison
plot(gg_adjuvant2)

df_combined2 <- df_combined %>% filter(CloneType !="")
gg_adjuvant2 <- ggboxplot(df_combined2, x = "CloneType", y = "CDR3_length", color = "adjuvant", palette = c("#377EB8", "#FF7F00", "#E41A1C", "#4DAF4A", "#984EA3", "#E5C494", "#FFD92F"),
                          order = c("Singletons", "Expanded", "Persisting")) +
  theme_bw() + ggtitle("") + xlab("\nClone Type") + ylab("CDR3 length") + labs(color ='Adjuvant') + ylim(0,25) + 
  stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5, label.y = 22, aes(group = adjuvant))                 # Pairwise comparison
plot(gg_adjuvant2)



df_combined2 <- df_combined %>% filter(CloneType2 !="")
my_comparisons <- list( c("Singletons", "Unique"), c("Unique", "Shared"), c("Singletons", "Shared") )
gg_adjuvant2 <- ggboxplot(df_combined2, x = "CloneType2", y = "CDR3_length", color = "CloneType2", palette = c("#377EB8", "#FF7F00", "#E41A1C", "#4DAF4A", "#984EA3", "#E5C494", "#FFD92F"),
                          order = c("Singletons", "Unique", "Shared")) +
  theme_bw() + ggtitle("") + xlab("\nClone Type") + ylab("CDR3 length") + labs(color ='Clone Type') + ylim(0,25) + 
  stat_compare_means(label = "p.signif", method = "t.test",  
                     comparisons = my_comparisons, label.y = c(21, 22, 24))  # Pairwise comparison
plot(gg_adjuvant2)


df_combined2 <- df_combined %>% filter(CloneType2 !="")
gg_adjuvant2 <- ggboxplot(df_combined2, x = "CloneType2", y = "CDR3_length", color = "adjuvant", palette = c("#377EB8", "#FF7F00", "#E41A1C", "#4DAF4A", "#984EA3", "#E5C494", "#FFD92F"),
                          order = c("Singletons", "Unique", "Shared")) +
  theme_bw() + ggtitle("") + xlab("\nClone Type") + ylab("CDR3 length") + labs(color ='Adjuvant') + ylim(0,25) + 
  stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5, label.y = 22, aes(group = adjuvant))                 # Pairwise comparison
plot(gg_adjuvant2)


########################
### save file ##########
########################

write.table(df_combined, "df_combined_seurat_BCR_SHM_CloneType_CDR3.txt", quote = F, row.names = F, sep = "\t")

###############
#### seurat
###############
## reload data from file
df_combined <- read.delim("df_combined_seurat_BCR_SHM_CloneType_CDR3.txt")
gex.combined <- readRDS("gex.combined.final.rds")

### map clone type to seruat object
Idents(gex.combined) <- "day"
DimPlot(gex.combined, reduction = "umap", split.by = "hash.ID", label = F, group.by = "day",
        ncol = 5 ) 

gex.combined$SHM <- df_combined$mu_freq
Idents(gex.combined) <- "day"
gex.combined$CDR3 <- df_combined$CDR3_length
FeaturePlot(gex.combined, features = c("SHM", "CDR3"), min.cutoff = 0.001)



######################################################################################################
## old scripts as below ##############################################################################
######################################################################################################

##################################
## Immcantation Portal ####
##################################

df1 <- read.delim("./21_Tube_3_results/heavy_parse-select_clone-pass_germ-pass.tsv")

# Calculate R and S mutation counts
db_obs <- observedMutations(df1, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=FALSE, 
                            nproc=1)
db_obs %>% 
  select(sequence_id, starts_with("mu_count_")) %>%
  head(n=4)

# Calculate R and S mutation frequencies
db_obs <- observedMutations(db_obs, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=TRUE, 
                            nproc=1)
# Show new mutation frequency columns
db_obs %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)

# Calculate combined R and S mutation frequencies
db_obs1 <- observedMutations(df1, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=TRUE, 
                            combine=TRUE,
                            nproc=1)
# Show new mutation frequency columns
db_obs1 %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)

## We can plot the mutation frequencies a explore differences between samples or isotypes.
g1 <- ggplot(db_obs1, aes(x=c_call, y=mu_freq, fill=c_call)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot() + ylim(0,0.2)
plot(g1)

######################################################################################################

df2 <- read.delim("./21_Tube_4_results/heavy_parse-select_clone-pass_germ-pass.tsv")

# Calculate R and S mutation counts
db_obs <- observedMutations(df2, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=FALSE, 
                            nproc=1)
db_obs %>% 
  select(sequence_id, starts_with("mu_count_")) %>%
  head(n=4)

# Calculate R and S mutation frequencies
db_obs <- observedMutations(db_obs, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=TRUE, 
                            nproc=1)
# Show new mutation frequency columns
db_obs %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)

# Calculate combined R and S mutation frequencies
db_obs2 <- observedMutations(df2, sequenceColumn="sequence_alignment",
                             germlineColumn="germline_alignment_d_mask",
                             regionDefinition=NULL,
                             frequency=TRUE, 
                             combine=TRUE,
                             nproc=1)
# Show new mutation frequency columns
db_obs2 %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)

## We can plot the mutation frequencies a explore differences between samples or isotypes.
g2 <- ggplot(db_obs2, aes(x=c_call, y=mu_freq, fill=c_call)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot() + ylim(0,0.2)
plot(g2)


######################################################################################################

df3 <- read.delim("./21_Tube_6_results/heavy_parse-select_clone-pass_germ-pass.tsv")

# Calculate R and S mutation counts
db_obs <- observedMutations(df3, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=FALSE, 
                            nproc=1)
db_obs %>% 
  select(sequence_id, starts_with("mu_count_")) %>%
  head(n=4)

# Calculate R and S mutation frequencies
db_obs <- observedMutations(db_obs, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=TRUE, 
                            nproc=1)
# Show new mutation frequency columns
db_obs %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)

# Calculate combined R and S mutation frequencies
db_obs3 <- observedMutations(df3, sequenceColumn="sequence_alignment",
                             germlineColumn="germline_alignment_d_mask",
                             regionDefinition=NULL,
                             frequency=TRUE, 
                             combine=TRUE,
                             nproc=1)
# Show new mutation frequency columns
db_obs3 %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)

## We can plot the mutation frequencies a explore differences between samples or isotypes.
g3 <- ggplot(db_obs3, aes(x=c_call, y=mu_freq, fill=c_call)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot() + ylim(0,0.2)
plot(g3)


######################################################################################################

df4 <- read.delim("./21_Tube_7_results/heavy_parse-select_clone-pass_germ-pass.tsv")

# Calculate R and S mutation counts
db_obs <- observedMutations(df4, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=FALSE, 
                            nproc=1)
db_obs %>% 
  select(sequence_id, starts_with("mu_count_")) %>%
  head(n=4)

# Calculate R and S mutation frequencies
db_obs <- observedMutations(db_obs, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=TRUE, 
                            nproc=1)
# Show new mutation frequency columns
db_obs %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)

# Calculate combined R and S mutation frequencies
db_obs4 <- observedMutations(df4, sequenceColumn="sequence_alignment",
                             germlineColumn="germline_alignment_d_mask",
                             regionDefinition=NULL,
                             frequency=TRUE, 
                             combine=TRUE,
                             nproc=1)
# Show new mutation frequency columns
db_obs4 %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)

## We can plot the mutation frequencies a explore differences between samples or isotypes.
g4 <- ggplot(db_obs4, aes(x=c_call, y=mu_freq, fill=c_call)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot() + ylim(0,0.2)
plot(g4)

######################################################################################################

df5 <- read.delim("./40_Tube_3_results/heavy_parse-select_clone-pass_germ-pass.tsv")

# Calculate R and S mutation counts
db_obs <- observedMutations(df5, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=FALSE, 
                            nproc=1)
db_obs %>% 
  select(sequence_id, starts_with("mu_count_")) %>%
  head(n=4)

# Calculate R and S mutation frequencies
db_obs <- observedMutations(db_obs, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=TRUE, 
                            nproc=1)
# Show new mutation frequency columns
db_obs %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)

# Calculate combined R and S mutation frequencies
db_obs5 <- observedMutations(df5, sequenceColumn="sequence_alignment",
                             germlineColumn="germline_alignment_d_mask",
                             regionDefinition=NULL,
                             frequency=TRUE, 
                             combine=TRUE,
                             nproc=1)
# Show new mutation frequency columns
db_obs5 %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)

## We can plot the mutation frequencies a explore differences between samples or isotypes.
g5 <- ggplot(db_obs5, aes(x=c_call, y=mu_freq, fill=c_call)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot() + ylim(0,0.2)
plot(g5)


######################################################################################################

df6 <- read.delim("./40_Tube_4_results/heavy_parse-select_clone-pass_germ-pass.tsv")

# Calculate R and S mutation counts
db_obs <- observedMutations(df6, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=FALSE, 
                            nproc=1)
db_obs %>% 
  select(sequence_id, starts_with("mu_count_")) %>%
  head(n=4)

# Calculate R and S mutation frequencies
db_obs <- observedMutations(db_obs, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=TRUE, 
                            nproc=1)
# Show new mutation frequency columns
db_obs %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)

# Calculate combined R and S mutation frequencies
db_obs6 <- observedMutations(df6, sequenceColumn="sequence_alignment",
                             germlineColumn="germline_alignment_d_mask",
                             regionDefinition=NULL,
                             frequency=TRUE, 
                             combine=TRUE,
                             nproc=1)
# Show new mutation frequency columns
db_obs6 %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)

## We can plot the mutation frequencies a explore differences between samples or isotypes.
g6 <- ggplot(db_obs6, aes(x=c_call, y=mu_freq, fill=c_call)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot() + ylim(0,0.2)
plot(g6)


######################################################################################################

df7 <- read.delim("./40_Tube_7_results/heavy_parse-select_clone-pass_germ-pass.tsv")

# Calculate R and S mutation counts
db_obs <- observedMutations(df7, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=FALSE, 
                            nproc=1)
db_obs %>% 
  select(sequence_id, starts_with("mu_count_")) %>%
  head(n=4)

# Calculate R and S mutation frequencies
db_obs <- observedMutations(db_obs, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=TRUE, 
                            nproc=1)
# Show new mutation frequency columns
db_obs %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)

# Calculate combined R and S mutation frequencies
db_obs7 <- observedMutations(df7, sequenceColumn="sequence_alignment",
                             germlineColumn="germline_alignment_d_mask",
                             regionDefinition=NULL,
                             frequency=TRUE, 
                             combine=TRUE,
                             nproc=1)
# Show new mutation frequency columns
db_obs7 %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)

## We can plot the mutation frequencies a explore differences between samples or isotypes.
g7 <- ggplot(db_obs7, aes(x=c_call, y=mu_freq, fill=c_call)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot() + ylim(0,0.2)
plot(g7)


######################################################################################################

df8 <- read.delim("./40_Tube_8_results/heavy_parse-select_clone-pass_germ-pass.tsv")

# Calculate R and S mutation counts
db_obs <- observedMutations(df8, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=FALSE, 
                            nproc=1)
db_obs %>% 
  select(sequence_id, starts_with("mu_count_")) %>%
  head(n=4)

# Calculate R and S mutation frequencies
db_obs <- observedMutations(db_obs, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=TRUE, 
                            nproc=1)
# Show new mutation frequency columns
db_obs %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)

# Calculate combined R and S mutation frequencies
db_obs8 <- observedMutations(df8, sequenceColumn="sequence_alignment",
                             germlineColumn="germline_alignment_d_mask",
                             regionDefinition=NULL,
                             frequency=TRUE, 
                             combine=TRUE,
                             nproc=1)
# Show new mutation frequency columns
db_obs8 %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)

## We can plot the mutation frequencies a explore differences between samples or isotypes.
g8 <- ggplot(db_obs8, aes(x=c_call, y=mu_freq, fill=c_call)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot() + ylim(0,0.2)
plot(g8)


#############################################
#############################################
## combined with meta data

df_meta <- read.delim("../df_combined_gex_bcr.txt")
df_meta <- df_meta %>% select(day, tissue, sorting, adjuvant, seurat_clusters, barcode, CTstrict) %>%
           rename(cell_id = barcode)

db_obs1$cell_id <- paste0(db_obs1$cell_id, "_3")
db_obs2$cell_id <- paste0(db_obs2$cell_id, "_4")
db_obs3$cell_id <- paste0(db_obs3$cell_id, "_6")
db_obs4$cell_id <- paste0(db_obs4$cell_id, "_7")
db_obs5$cell_id <- paste0(db_obs5$cell_id, "_10")
db_obs6$cell_id <- paste0(db_obs6$cell_id, "_11")
db_obs7$cell_id <- paste0(db_obs7$cell_id, "_14")
db_obs8$cell_id <- paste0(db_obs8$cell_id, "_15")

SHM.data <- bind_rows(db_obs1, db_obs2) %>% bind_rows(db_obs3) %>% bind_rows(db_obs4) %>% bind_rows(db_obs5) %>% 
          bind_rows(db_obs6) %>% bind_rows(db_obs7) %>% bind_rows(db_obs8)


df_combined <- left_join(df_meta, SHM.data)

df_rsv <- df_combined %>% filter(adjuvant =="RSV Only") %>% filter(!is.na(mu_freq))
df_lina2 <- df_combined %>% filter(adjuvant =="RSV and LiNA2") %>% filter(!is.na(mu_freq))


## We can plot the mutation frequencies a explore differences between samples or isotypes.
g_rsv <- ggplot(df_rsv, aes(x=c_call, y=mu_freq, fill=c_call)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot() + ylim(0,0.2)
plot(g_rsv)

## We can plot the mutation frequencies a explore differences between samples or isotypes.
g_lina2 <- ggplot(df_lina2, aes(x=c_call, y=mu_freq, fill=c_call)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot() + ylim(0,0.2)
plot(g_lina2)


############################################################################################
############################################################################################
############################################################################################

df_rsv <- df_combined %>% filter(adjuvant =="RSV Only" & seurat_clusters == 2) %>% filter(!is.na(mu_freq))
df_lina2 <- df_combined %>% filter(adjuvant =="RSV and LiNA2" & seurat_clusters == 2) %>% filter(!is.na(mu_freq))


## We can plot the mutation frequencies a explore differences between samples or isotypes.
g_rsv <- ggplot(df_rsv, aes(x=c_call, y=mu_freq, fill=c_call)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot() + ylim(0,0.2)
plot(g_rsv)

## We can plot the mutation frequencies a explore differences between samples or isotypes.
g_lina2 <- ggplot(df_lina2, aes(x=c_call, y=mu_freq, fill=c_call)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot() + ylim(0,0.2)
plot(g_lina2)


############################################################################################
############################################################################################
############################################################################################

df_rsv <- df_combined %>% filter(adjuvant =="RSV Only" & seurat_clusters == 4) %>% filter(!is.na(mu_freq))
df_lina2 <- df_combined %>% filter(adjuvant =="RSV and LiNA2" & seurat_clusters == 4) %>% filter(!is.na(mu_freq))


## We can plot the mutation frequencies a explore differences between samples or isotypes.
g_rsv <- ggplot(df_rsv, aes(x=c_call, y=mu_freq, fill=c_call)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot() + ylim(0,0.2)
plot(g_rsv)

## We can plot the mutation frequencies a explore differences between samples or isotypes.
g_lina2 <- ggplot(df_lina2, aes(x=c_call, y=mu_freq, fill=c_call)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot() + ylim(0,0.2)
plot(g_lina2)


############################################################################################
############################################################################################
############################################################################################

df_rsv <- df_combined %>% filter(adjuvant =="RSV Only" & seurat_clusters == 5) %>% filter(!is.na(mu_freq))
df_lina2 <- df_combined %>% filter(adjuvant =="RSV and LiNA2" & seurat_clusters == 5) %>% filter(!is.na(mu_freq))


## We can plot the mutation frequencies a explore differences between samples or isotypes.
g_rsv <- ggplot(df_rsv, aes(x=c_call, y=mu_freq, fill=c_call)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot() + ylim(0,0.2)
plot(g_rsv)

## We can plot the mutation frequencies a explore differences between samples or isotypes.
g_lina2 <- ggplot(df_lina2, aes(x=c_call, y=mu_freq, fill=c_call)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot() + ylim(0,0.2)
plot(g_lina2)


############################################################################################
############################################################################################
############################################################################################

df_rsv <- df_combined %>% filter(adjuvant =="RSV Only" & seurat_clusters == 8) %>% filter(!is.na(mu_freq))
df_lina2 <- df_combined %>% filter(adjuvant =="RSV and LiNA2" & seurat_clusters == 8) %>% filter(!is.na(mu_freq))


## We can plot the mutation frequencies a explore differences between samples or isotypes.
g_rsv <- ggplot(df_rsv, aes(x=c_call, y=mu_freq, fill=c_call)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot() + ylim(0,0.2)
plot(g_rsv)

## We can plot the mutation frequencies a explore differences between samples or isotypes.
g_lina2 <- ggplot(df_lina2, aes(x=c_call, y=mu_freq, fill=c_call)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot() + ylim(0,0.2)
plot(g_lina2)


############################################################################################
############################################################################################
############################################################################################

df_rsv <- df_combined %>% filter(adjuvant =="RSV Only" & seurat_clusters == 9) %>% filter(!is.na(mu_freq))
df_lina2 <- df_combined %>% filter(adjuvant =="RSV and LiNA2" & seurat_clusters == 9) %>% filter(!is.na(mu_freq))


## We can plot the mutation frequencies a explore differences between samples or isotypes.
g_rsv <- ggplot(df_rsv, aes(x=c_call, y=mu_freq, fill=c_call)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot() + ylim(0,0.2)
plot(g_rsv)

## We can plot the mutation frequencies a explore differences between samples or isotypes.
g_lina2 <- ggplot(df_lina2, aes(x=c_call, y=mu_freq, fill=c_call)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Isotype") + ylab("Mutation frequency") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot() + ylim(0,0.2)
plot(g_lina2)
