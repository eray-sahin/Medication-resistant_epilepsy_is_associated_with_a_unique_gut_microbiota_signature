rm(list = ls())

library(phyloseq)
library(tidyverse)
library(vegan)
library(ggplot2)
library(gghalves)
library(Maaslin2)
library(readxl)
library(egg)

# setwd(PATH)

#########
# Import
#########

# ASV count table
asv.counts <- read_xlsx("epi_diet_table_feature_table_arranged.xlsx", col_names = T) %>%
  column_to_rownames("ASV")

# ASV taxonomy table
asv.tax <- read_xlsx("epi_diet_taxonomy_arranged.xlsx", col_names = T, na = "") %>%
  column_to_rownames("ASV")
# Replacing NA with Unidentified
asv.tax[is.na(asv.tax)] <- "Unidentified"

# Import metadata, and sort based on the order of samples on genus ra table
metadata <- readxl::read_xlsx("Epi_Diet_Final_Metadata.xlsx") %>% column_to_rownames(var = "SampleID") %>%
  select(-SeqID)
metadata <- metadata[match(colnames(asv.counts), rownames(metadata)), ]
# factor the groups in desired order
metadata$Group <- factor(metadata$Group, levels = c("HC", "MS", "MR"))

# checking order of the samples
all(rownames(metadata) %in% colnames(asv.counts)) #TRUE
all(rownames(metadata) == colnames(asv.counts)) #TRUE

# Import rooted phylogenetic tree
tree <- read_tree("epi_diet_tree.nwk")

# creating phlyseq object
phylo <- phyloseq(otu_table(as.matrix(asv.counts), taxa_are_rows = T),
                  sample_data(metadata),
                  tax_table(as.matrix(asv.tax)))
phylo <- merge_phyloseq(phylo, tree)
ape::is.rooted(phy_tree(phylo)) #TRUE

# save(phylo, file = "phylo.RData")

# rarefication for alpha diversity analysis using the minimum feature size across the samples
# phylo.rarefied = rarefy_even_depth(phylo, rngseed=21, sample.size=min(sample_sums(phylo)), replace=F)
# save(phylo.rarefied, file = "phylo.rarefied.RData")
load("phylo.rarefied.RData")



########################
# Beta Diversity
########################

### Bray-Curtis ###

bray.curtis_dist = phyloseq::distance(phylo, method="bray")
bray.curtis_ordination = ordinate(phylo, method="PCoA", distance=bray.curtis_dist)

sampledf <- data.frame(sample_data(phylo))

# PERMANOVA
set.seed(123)
adonis2(bray.curtis_dist ~ Group + Age + Sex + Carbohydrate + Protein + Fat + Fiber, 
        data = sampledf, permutations = 9999, by = "terms")

# Bray Curtis PCoA plot

bray.curtis_df <- data.frame(bray.curtis_ordination$vectors) %>%
  rownames_to_column(var = "SampleID") %>%
  select(SampleID, Axis.1, Axis.2) %>%
  left_join(., metadata %>% rownames_to_column(var = "SampleID"), by = "SampleID") 

bray <- ggplot(bray.curtis_df, aes(x=Axis.1, y=Axis.2, color = Group)) +
  geom_point(alpha=0.9, size=3) + #alpha controls transparency and helps when points are overlapping
  stat_ellipse(geom = "polygon", alpha = 1/7, aes(fill = Group), linewidth = 1) +
  labs(x=paste0("Axis.1 (", round(bray.curtis_ordination$values$Relative_eig[1] * 100, digits = 2), "%)"), 
       y=paste0("Axis.2 (", round(bray.curtis_ordination$values$Relative_eig[2] * 100, digits = 2), "%)"),
       title = "PCoA (Bray-Curtis)") + theme(axis.line = element_line(colour = "grey"),
                                             panel.grid.major = element_blank(),
                                             panel.grid.minor = element_blank(),
                                             panel.border = element_blank(),
                                             panel.background = element_blank(),
                                             axis.title = element_text(size = 20),
                                             axis.text.x = element_text(size = 15),
                                             axis.text.y = element_text(size = 15),
                                             plot.title = element_text(size = 20),
                                             legend.text = element_text(size = 15),
                                             legend.title = element_text(size=15),
                                             strip.text = element_text(size=15),
                                             legend.position='right') +
  scale_fill_manual(values = c("#F4B942", "#97D8C4", "#4059AD")) + 
  scale_colour_manual(values = c("#F4B942", "#97D8C4", "#4059AD")) 

# Pairwise PERMANOVA
set.seed(123)
bray_pairwise_adonis <- pairwiseAdonis::pairwise.adonis2(bray.curtis_dist ~ Group + Age + Sex + Carbohydrate + Protein + Fat + Fiber, 
                                                         data = sampledf, nperm = 9999)
bray_pairwise_adonis_summary <- rbind(
  bray_pairwise_adonis$MS_vs_HC %>% 
    as.data.frame() %>% 
    rownames_to_column("Variable") %>%
    dplyr::filter(Variable == "Model") %>%
    mutate(Comparison = "MS vs HC", .before = Variable) %>%
    dplyr::select(Comparison, `F`, `Pr(>F)`),
  
  bray_pairwise_adonis$MR_vs_HC %>% 
    as.data.frame() %>% 
    rownames_to_column("Variable") %>%
    dplyr::filter(Variable == "Model") %>%
    mutate(Comparison = "MR vs HC", .before = Variable) %>%
    dplyr::select(Comparison, `F`, `Pr(>F)`) ,
  
  bray_pairwise_adonis$MR_vs_MS %>% 
    as.data.frame() %>% 
    rownames_to_column("Variable") %>%
    dplyr::filter(Variable == "Model") %>%
    mutate(Comparison = "MR vs MS", .before = Variable) %>%
    dplyr::select(Comparison, `F`, `Pr(>F)`)
)
# Summary table plot
bray_pairwise_adonis_summary_p <- ggpubr::ggtexttable(bray_pairwise_adonis_summary, 
                                                      rows = NULL, theme = ggpubr::ttheme("blank")) %>%
  ggpubr::tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)


### Jaccard ###

jaccard_dist = phyloseq::distance(phylo, method="jaccard")
jaccard_ordination = ordinate(phylo, method="PCoA", distance = jaccard_dist)

# PERMANOVA
set.seed(123)
adonis2(jaccard_dist ~ Group + Age + Sex + Carbohydrate + Protein + Fat + Fiber, 
        data = sampledf, 
        permutations = 9999, by = "terms")


# Jaccard PCoA plot
jaccard_df <- data.frame(jaccard_ordination$vectors) %>%
  rownames_to_column(var = "SampleID") %>%
  select(SampleID, Axis.1, Axis.2) %>%
  left_join(., metadata %>% rownames_to_column(var = "SampleID"), by = "SampleID") 

jaccard <- ggplot(jaccard_df, aes(x=Axis.1, y=Axis.2, color = Group)) +
  geom_point(alpha = 0.9, size = 3) + #alpha controls transparency and helps when points are overlapping
  stat_ellipse(geom = "polygon", alpha = 1/7, aes(fill = Group), linewidth = 1) +
  labs(x=paste0("Axis.1 (", round(jaccard_ordination$values$Relative_eig[1] * 100, digits = 2), "%)"), 
       y=paste0("Axis.2 (", round(jaccard_ordination$values$Relative_eig[2] * 100, digits = 2), "%)"),
       title = "PCoA (Jaccard)") + theme(axis.line = element_line(colour = "grey"),
                                             panel.grid.major = element_blank(),
                                             panel.grid.minor = element_blank(),
                                             panel.border = element_blank(),
                                             panel.background = element_blank(),
                                             axis.title = element_text(size = 20),
                                             axis.text.x = element_text(size = 15),
                                             axis.text.y = element_text(size = 15),
                                             plot.title = element_text(size = 20),
                                             legend.text = element_text(size = 15),
                                         legend.title = element_text(size=15),
                                             strip.text = element_text(size=15),
                                         legend.position='right') +
  scale_fill_manual(values = c("#F4B942", "#97D8C4", "#4059AD")) + 
  scale_colour_manual(values = c("#F4B942", "#97D8C4", "#4059AD"))

# Pairwise PERMANOVA

set.seed(123)
jaccard_pairwise_adonis <- pairwiseAdonis::pairwise.adonis2(jaccard_dist ~ Group + Age + Sex + Carbohydrate + Protein + Fat + Fiber, 
                                                         data = sampledf, nperm = 9999)
jaccard_dist_pairwise_adonis_summary <- rbind(
  jaccard_pairwise_adonis$MS_vs_HC %>% 
    as.data.frame() %>% 
    rownames_to_column("Variable") %>%
    dplyr::filter(Variable == "Model") %>%
    mutate(Comparison = "MS vs HC", .before = Variable) %>%
    dplyr::select(Comparison, `F`, `Pr(>F)`),
  
  jaccard_pairwise_adonis$MR_vs_HC %>% 
    as.data.frame() %>% 
    rownames_to_column("Variable") %>%
    dplyr::filter(Variable == "Model") %>%
    mutate(Comparison = "MR vs HC", .before = Variable) %>%
    dplyr::select(Comparison, `F`, `Pr(>F)`) ,
  
  jaccard_pairwise_adonis$MR_vs_MS %>% 
    as.data.frame() %>% 
    rownames_to_column("Variable") %>%
    dplyr::filter(Variable == "Model") %>%
    mutate(Comparison = "MR vs MS", .before = Variable) %>%
    dplyr::select(Comparison, `F`, `Pr(>F)`)
)
# Summary table plot
jaccard_dist_pairwise_adonis_summary_p <- ggpubr::ggtexttable(jaccard_dist_pairwise_adonis_summary, 
                                                      rows = NULL, theme = ggpubr::ttheme("blank")) %>%
  ggpubr::tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2) 

library(patchwork)
png("Results/beta_diversity/bray_jaccard_merged_plots.png",
    units="in", width=15, height = 7, res = 800)
bray + jaccard + plot_layout(guides = "collect")
dev.off()

png("Results/beta_diversity/bray_pairwise_adonis_table.png",
    units="in", width=7, height = 3, res = 800)
bray_pairwise_adonis_summary_p  
dev.off()

png("Results/beta_diversity/jaccard_pairwise_adonis_table.png",
    units="in", width=7, height = 3, res = 800)
jaccard_dist_pairwise_adonis_summary_p  
dev.off()

#####################
# Alpha Diversity
#####################

alpha_diversity_df <- estimate_richness(phylo.rarefied)

# adding metadata to the alpha diversity df
alpha_diversity_df <- alpha_diversity_df %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(., 
            metadata %>% rownames_to_column(var = "SampleID"), 
            by = "SampleID") %>%
  column_to_rownames(var = "SampleID")

# write_tsv(alpha_diversity_df %>% rownames_to_column("SampleID"),
#           "Results/alpha_diversity/all_alpha_div_all_indices_calculated.tsv")


###########################
# GLM with Gamma log link
###########################

###################
# MR or MS vs HC
###################

chao1_three_groups_glm_coef <- summary(glm(Chao1 ~ Group + Age + Sex + Carbohydrate + Protein + Fat + Fiber, data = alpha_diversity_df, 
                                           family = Gamma(link = "log")))$coefficients %>% 
  as.data.frame %>%
  rownames_to_column("Predictor") %>%
  mutate(Response = rep("Chao1", nrow(.)), .before = 1)

shannon_three_groups_glm_coef <- summary(glm(Shannon ~ Group + Age + Sex + Carbohydrate + Protein + Fat + Fiber, data = alpha_diversity_df, 
                                             family = Gamma(link = "log")))$coefficients %>% 
  as.data.frame %>%
  rownames_to_column("Predictor") %>%
  mutate(Response = rep("Shannon Index", nrow(.)), .before = 1)

simpson_three_groups_glm_coef <- summary(glm(Simpson ~ Group + Age + Sex + Carbohydrate + Protein + Fat + Fiber, data = alpha_diversity_df, 
                                             family = Gamma(link = "log")))$coefficients %>% 
  as.data.frame %>%
  rownames_to_column("Predictor") %>%
  mutate(Response = rep("Simpson's Index", nrow(.)), .before = 1)


# combine all results into one df
alpha_all_three_groups_glm <- dplyr::bind_rows(chao1_three_groups_glm_coef,
                                               shannon_three_groups_glm_coef,
                                               simpson_three_groups_glm_coef) %>%
  mutate(P.adj = p.adjust(`Pr(>|t|)`, method = "fdr")) %>% # adjusting p values
  mutate(Significance = case_when(P.adj >= 0.1 ~ "NS", # adding significance stars
                                  P.adj >= 0.05 & P.adj < 0.1 ~ ".",
                                  P.adj >= 0.01 & P.adj < 0.05 ~ "*",
                                  P.adj >= 0.001 & P.adj < 0.01 ~ "**",
                                  P.adj < 0.001 ~ "***"))

cbind(alpha_all_three_groups_glm$Predictor, exp(alpha_all_three_groups_glm$Estimate))

openxlsx::write.xlsx(alpha_all_three_groups_glm, "Results/alpha_diversity/Alpha_Diversity_MR_or_MS_vs_Healthy_GLM/alpha_all_three_groups_glm.xlsx")


#############
# MR vs MS
#############

alpha_diversity_df_epi <- alpha_diversity_df %>%
  filter(Group %in% c("MS", "MR"))

chao1_epi_groups_glm_coef <- summary(glm(Chao1 ~ Group + Age + Sex + Carbohydrate + Protein + Fat + Fiber, data = alpha_diversity_df_epi, 
                                           family = Gamma(link = "log")))$coefficients %>% 
  as.data.frame %>%
  rownames_to_column("Predictor") %>%
  mutate(Response = rep("Chao1", nrow(.)), .before = 1)

shannon_epi_groups_glm_coef <- summary(glm(Shannon ~ Group + Age + Sex + Carbohydrate + Protein + Fat + Fiber, data = alpha_diversity_df_epi, 
                                             family = Gamma(link = "log")))$coefficients %>% 
  as.data.frame %>%
  rownames_to_column("Predictor") %>%
  mutate(Response = rep("Shannon Index", nrow(.)), .before = 1)

simpson_epi_groups_glm_coef <- summary(glm(Simpson ~ Group + Age + Sex + Carbohydrate + Protein + Fat + Fiber, data = alpha_diversity_df_epi, 
                                             family = Gamma(link = "log")))$coefficients %>% 
  as.data.frame %>%
  rownames_to_column("Predictor") %>%
  mutate(Response = rep("Simpson's Index", nrow(.)), .before = 1)


# combine all results into one df
alpha_all_epi_groups_glm <- dplyr::bind_rows(chao1_epi_groups_glm_coef,
                                             shannon_epi_groups_glm_coef,
                                             simpson_epi_groups_glm_coef) %>%
  mutate(P.adj = p.adjust(`Pr(>|t|)`, method = "fdr")) %>% # adjusting p values
  mutate(Significance = case_when(P.adj >= 0.1 ~ "NS", # adding significance stars
                                  P.adj >= 0.05 & P.adj < 0.1 ~ ".",
                                  P.adj >= 0.01 & P.adj < 0.05 ~ "*",
                                  P.adj >= 0.001 & P.adj < 0.01 ~ "**",
                                  P.adj < 0.001 ~ "***"))

openxlsx::write.xlsx(alpha_all_epi_groups_glm, "Results/alpha_diversity/Alpha_Diversity_MR_vs_MS_GLM/alpha_all_epi_groups_glm.xlsx")

cbind(alpha_all_three_groups_glm$Predictor, exp(alpha_all_three_groups_glm$Estimate))
