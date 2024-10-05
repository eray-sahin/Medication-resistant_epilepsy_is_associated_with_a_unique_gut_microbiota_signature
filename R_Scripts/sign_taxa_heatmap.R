rm(list = ls())

library(tidyverse)
library(ggplot2)
library(vegan)
library(FSA)
library(gghalves)
library(readxl)

# setwd(PATH)

# import genus ra table
genus.ra.raw <- read.table("asv.genus.ra.tsv", 
                           header = T, row.names = 1, sep = "\t")

# import metadata, and sort based on the order of samples on genus ra table
metadata <- readxl::read_xlsx("Epi_Diet_Final_Metadata.xlsx") %>% column_to_rownames(var = "SampleID")
metadata <- metadata[match(rownames(genus.ra.raw), rownames(metadata)), ]
# factor the groups in desired order
metadata$Group <- factor(metadata$Group, levels = c("HC", "MS", "MR"))

  

#minimm of the genus rel abundance
genus_min <- min(apply(genus.ra.raw, 2, function(x) min(x[x != 0])))
# pseudo-count as half of min
genus.rel.abun.no.zero <- replace(genus.ra.raw, genus.ra.raw == 0, (genus_min/2))

# genus.rel.abun.no.zero with group info
genus.rel.abun.no.zero.group <- genus.rel.abun.no.zero %>% rownames_to_column(var = "sample.id") %>%
  left_join(., metadata %>% rownames_to_column(var = "sample.id") %>% select(sample.id, Group), by = "sample.id")  %>%
  column_to_rownames(var = "sample.id") %>%
  as.data.frame()

# clr transformation
genus.clr <- vegan::decostand(genus.rel.abun.no.zero, method = "clr", MARGIN = 2)

# sign genera only
sign.genera.maaslin2 <- read_xlsx("sign_genera_from_maaslin2.xlsx", col_names = T)
sign.genera.clr <- genus.clr %>%
  select(all_of(unique(sign.genera.maaslin2$feature)))

# correcting Eubacterium name
colnames(sign.genera.clr)[2] <- "[Eubacterium] siraeum group"

sign.genera.clr.meta <- sign.genera.clr %>% rownames_to_column(var = "sample.id") %>%
  left_join(., metadata %>% rownames_to_column(var = "sample.id"), by = "sample.id")  %>%
  column_to_rownames(var = "sample.id") %>%
  as.data.frame()


##########################
# Significance
#########################

# calculating fold changes
# MS vs HC
ms_vs_hc_sign_genera_fc <- data.frame(matrix(nrow = 1, ncol = 5))
colnames(ms_vs_hc_sign_genera_fc) <- c("Genus", "MS_Median", "HC_Median", "MS_HC_Log2_FC", "qval")
ms_vs_hc_sign_genera_fc$Genus <- sign.genera.maaslin2$feature[sign.genera.maaslin2$comparison == "MS_vs_HC"]
ms_vs_hc_sign_genera_fc$MS_Median <- median(genus.rel.abun.no.zero.group$Hungatella[genus.rel.abun.no.zero.group$Group == "MS"])
ms_vs_hc_sign_genera_fc$HC_Median <- median(genus.rel.abun.no.zero.group$Hungatella[genus.rel.abun.no.zero.group$Group == "HC"])
ms_vs_hc_sign_genera_fc$MS_HC_Log2_FC <- log2(ms_vs_hc_sign_genera_fc$MS_Median) - log2(ms_vs_hc_sign_genera_fc$HC_Median)
ms_vs_hc_sign_genera_fc$qval <- sign.genera.maaslin2$qval[sign.genera.maaslin2$feature == "Hungatella" & sign.genera.maaslin2$comparison == "MS_vs_HC"]

# MR vs HC

mr_vs_hc_sign_genera <- sign.genera.maaslin2$feature[sign.genera.maaslin2$comparison == "MR_vs_HC"]
mr_vs_hc_sign_genera_fc <- data.frame(matrix(nrow = length(mr_vs_hc_sign_genera), ncol = 5))
colnames(mr_vs_hc_sign_genera_fc) <- c("Genus", "MR_Median", "HC_Median", "MR_HC_Log2_FC", "qval")
for (i in 1:length(mr_vs_hc_sign_genera)) {
  
  mr_vs_hc_sign_genera_fc$Genus[i] <- sign.genera.maaslin2$feature[sign.genera.maaslin2$comparison == "MR_vs_HC"][i]
  mr_vs_hc_sign_genera_fc$MR_Median[i] <- median(genus.rel.abun.no.zero.group[genus.rel.abun.no.zero.group$Group == "MR", as.character(mr_vs_hc_sign_genera_fc$Genus[i])])
  mr_vs_hc_sign_genera_fc$HC_Median[i] <- median(genus.rel.abun.no.zero.group[genus.rel.abun.no.zero.group$Group == "HC", as.character(mr_vs_hc_sign_genera_fc$Genus[i])])
  mr_vs_hc_sign_genera_fc$MR_HC_Log2_FC[i] <- log2(mr_vs_hc_sign_genera_fc$MR_Median)[i] - log2(mr_vs_hc_sign_genera_fc$HC_Median)[i]
  mr_vs_hc_sign_genera_fc$qval[i] <- sign.genera.maaslin2$qval[sign.genera.maaslin2$feature == as.character(mr_vs_hc_sign_genera_fc$Genus[i]) & sign.genera.maaslin2$comparison == "MR_vs_HC"]
  
}

# MR vs MS

mr_vs_ms_sign_genera <- sign.genera.maaslin2$feature[sign.genera.maaslin2$comparison == "MR_vs_MS"]
mr_vs_ms_sign_genera_fc <- data.frame(matrix(nrow = length(mr_vs_ms_sign_genera), ncol = 5))
colnames(mr_vs_ms_sign_genera_fc) <- c("Genus", "MR_Median", "MS_Median", "MR_MS_Log2_FC", "qval")

for (i in 1:length(mr_vs_ms_sign_genera)) {
  
  mr_vs_ms_sign_genera_fc$Genus[i] <- sign.genera.maaslin2$feature[sign.genera.maaslin2$comparison == "MR_vs_MS"][i]
  mr_vs_ms_sign_genera_fc$MR_Median[i] <- median(genus.rel.abun.no.zero.group[genus.rel.abun.no.zero.group$Group == "MR", as.character(mr_vs_ms_sign_genera_fc$Genus[i])])
  mr_vs_ms_sign_genera_fc$MS_Median[i] <- median(genus.rel.abun.no.zero.group[genus.rel.abun.no.zero.group$Group == "MS", as.character(mr_vs_ms_sign_genera_fc$Genus[i])])
  mr_vs_ms_sign_genera_fc$MR_MS_Log2_FC[i] <- log2(mr_vs_ms_sign_genera_fc$MR_Median)[i] - log2(mr_vs_ms_sign_genera_fc$MS_Median)[i]
  mr_vs_ms_sign_genera_fc$qval[i] <- sign.genera.maaslin2$qval[sign.genera.maaslin2$feature == as.character(mr_vs_ms_sign_genera_fc$Genus[i]) & sign.genera.maaslin2$comparison == "MR_vs_MS"]
  
}

mr_vs_ms_sign_genera_fc$Genus <- "[Eubacterium] siraeum group"

# genus
ms_vs_hc_sign_genera_fc_heatmap <- ms_vs_hc_sign_genera_fc %>%
  select(Genus, MS_HC_Log2_FC, qval) %>%
  mutate(genus_ann = case_when(qval > 0.01 & qval < 0.05 & MS_HC_Log2_FC > 0 ~ paste0("+", round(MS_HC_Log2_FC, 2), "*"),
                               qval > 0.01 & qval < 0.05 & MS_HC_Log2_FC < 0 ~ paste0(round(MS_HC_Log2_FC, 2), "*"),
                               qval > 0.001 & qval <= 0.01 & MS_HC_Log2_FC > 0 ~ paste0("+", round(MS_HC_Log2_FC, 2), "**"),
                               qval > 0.001 & qval <= 0.01 & MS_HC_Log2_FC < 0 ~ paste0(round(MS_HC_Log2_FC, 2), "**"),
                               qval <= 0.001 & MS_HC_Log2_FC > 0 ~ paste0("+", round(MS_HC_Log2_FC, 2), "***"),
                               qval <= 0.001 & MS_HC_Log2_FC < 0 ~ paste0(round(MS_HC_Log2_FC, 2), "***")))

mr_vs_hc_sign_genera_fc_heatmap <- mr_vs_hc_sign_genera_fc %>%
  select(Genus, MR_HC_Log2_FC, qval) %>%
  mutate(genus_ann = case_when(qval > 0.01 & qval < 0.05 & MR_HC_Log2_FC > 0 ~ paste0("+", round(MR_HC_Log2_FC, 2), "*"),
                               qval > 0.01 & qval < 0.05 & MR_HC_Log2_FC < 0 ~ paste0(round(MR_HC_Log2_FC, 2), "*"),
                               qval > 0.001 & qval <= 0.01 & MR_HC_Log2_FC > 0 ~ paste0("+", round(MR_HC_Log2_FC, 2), "**"),
                               qval > 0.001 & qval <= 0.01 & MR_HC_Log2_FC < 0 ~ paste0(round(MR_HC_Log2_FC, 2), "**"),
                               qval <= 0.001 & MR_HC_Log2_FC > 0 ~ paste0("+", round(MR_HC_Log2_FC, 2), "***"),
                               qval <= 0.001 & MR_HC_Log2_FC < 0 ~ paste0(round(MR_HC_Log2_FC, 2), "***")))

mr_vs_ms_sign_genera_fc_heatmap <- mr_vs_ms_sign_genera_fc %>%
  select(Genus, MR_MS_Log2_FC, qval) %>%
  mutate(genus_ann = case_when(qval > 0.01 & qval < 0.05 & MR_MS_Log2_FC > 0 ~ paste0("+", round(MR_MS_Log2_FC, 2), "*"),
                               qval > 0.01 & qval < 0.05 & MR_MS_Log2_FC < 0 ~ paste0(round(MR_MS_Log2_FC, 2), "*"),
                               qval > 0.001 & qval <= 0.01 & MR_MS_Log2_FC > 0 ~ paste0("+", round(MR_MS_Log2_FC, 2), "**"),
                               qval > 0.001 & qval <= 0.01 & MR_MS_Log2_FC < 0 ~ paste0(round(MR_MS_Log2_FC, 2), "**"),
                               qval <= 0.001 & MR_MS_Log2_FC > 0 ~ paste0("+", round(MR_MS_Log2_FC, 2), "***"),
                               qval <= 0.001 & MR_MS_Log2_FC < 0 ~ paste0(round(MR_MS_Log2_FC, 2), "***")))


# combine
sign.2f.taxa_scaled_t_significant_fc <- as.data.frame(t(sign.genera.clr)) %>%
  rownames_to_column("taxa") %>%
  select(taxa) %>%
  as.data.frame() %>%
  left_join(ms_vs_hc_sign_genera_fc_heatmap %>% select(Genus, genus_ann), by = c("taxa" = "Genus")) %>% 
  mutate(MS_vs_HC_Log2.FC = coalesce(genus_ann), .keep = "unused") %>%
  left_join(mr_vs_hc_sign_genera_fc_heatmap %>% select(Genus, genus_ann), by = c("taxa" = "Genus")) %>% 
  mutate(MR_vs_HC_Log2.FC = coalesce(genus_ann), .keep = "unused") %>%
  left_join(mr_vs_ms_sign_genera_fc_heatmap %>% select(Genus, genus_ann), by = c("taxa" = "Genus")) %>% 
  mutate(MR_vs_MS_Log2.FC = coalesce(genus_ann), .keep = "unused")

# removing NA's
sign.2f.taxa_scaled_t_significant_fc[is.na(sign.2f.taxa_scaled_t_significant_fc)] <- ""

library(ComplexHeatmap)
set.seed(123)
m = as.matrix(t(sign.genera.clr))
rownames(m)[2] <- "[Eubacterium] siraeum group"
fa = sign.genera.clr.meta$Group
fa_col = c("HC" = "#F4B942", "MS" = "#97D8C4", "MR" = "#4059AD")
# dend1 = cluster_between_groups(m, fa)
col_fun = viridis::viridis(100, begin = 0, end = 1)   

png("Results/heatmap/sign_taxa_heatmap.png",
    units = "in", width = 12, height = 6, res = 800)

a <- Heatmap(m, 
             cluster_columns = F, 
             column_split = fa, 
             show_column_dend = F, 
             # column_title = NULL,   
             row_title = NULL, col = col_fun, show_column_names = F, 
             top_annotation = HeatmapAnnotation(Group = fa,
                                                col = list(Group = fa_col),
                                                show_legend = F,
                                                show_annotation_name = F),
             cluster_rows = F, 
             cluster_row_slices = FALSE, show_row_dend = F, 
             heatmap_legend_param = list(title = "Clr rel. abund.",
                                         title_position = "lefttop",
                                         direction = "horizontal",
                                         title_gp = gpar(fontsize = 15), 
                                         labels_gp = gpar(fontsize = 15)),
             right_annotation = rowAnnotation(
               `MS vs HC` = anno_text(sign.2f.taxa_scaled_t_significant_fc$MS_vs_HC_Log2.FC,
                                            location = 0.5, just = "center", 
                                            unit(12, "mm"), show_name = T, 
                                            gp = gpar(col = "black", fontsize = 10, border = "black", fontface = "bold")),
               `MR vs HC` = anno_text(sign.2f.taxa_scaled_t_significant_fc$MR_vs_HC_Log2.FC,
                                         location = 0.5, just = "center", 
                                         unit(12, "mm"), show_name = T,
                                         gp = gpar(col = "black", fontsize = 10, border = "black", fontface = "bold")),
               `MR vs MS` = anno_text(sign.2f.taxa_scaled_t_significant_fc$MR_vs_MS_Log2.FC,
                                           location = 0.5, just = "center", 
                                           unit(12, "mm"), show_name = T,
                                           gp = gpar(col = "black", fontsize = 10, border = "black", fontface = "bold")),
               annotation_name_side = c("top", "top", "top"),
               annotation_name_gp = gpar(fontsize = 15, fontface = "bold"),
               annotation_name_rot = c(45, 45, 45),
               annotation_name_offset =  unit(2, "mm"),
               annotation_name_align = T),
             row_names_max_width = unit(10, "cm"),
             show_heatmap_legend = T, 
             column_title_gp = grid::gpar(fontsize = 20, fontface = "bold"),
             row_names_gp = grid::gpar(fontsize = 15, fontface = "italic"))
draw(a, heatmap_legend_side="bottom")
dev.off()

