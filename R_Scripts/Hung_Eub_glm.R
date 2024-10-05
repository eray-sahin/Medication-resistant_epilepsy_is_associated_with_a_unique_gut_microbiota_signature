rm(list = ls())

library(Maaslin2)
library(tidyverse)
library(readxl)
library(openxlsx)
library(ggplot2)

# setwd(PATH)

# import genus ra table
genus.ra.raw <- read.table("asv.genus.ra.tsv", 
                           header = T, row.names = 1, sep = "\t")

# import metadata, and sort based on the order of samples on genus ra table
metadata <- readxl::read_xlsx("Epi_Diet_Final_Metadata.xlsx") %>% column_to_rownames(var = "SampleID")
metadata <- metadata[match(rownames(genus.ra.raw), rownames(metadata)), ]
# factor the groups in desired order
metadata$Group <- factor(metadata$Group, levels = c("HC", "MS", "MR"))
# number of samples in each group
table(metadata$Group)

# checking order of the samples
all(rownames(metadata) %in% rownames(genus.ra.raw)) #TRUE
all(rownames(metadata) == rownames(genus.ra.raw)) #TRUE

# adding group information to genus ra
genus.ra.raw.meta <- genus.ra.raw %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(., metadata %>% rownames_to_column(var = "SampleID") %>% select(SampleID, Group), by = "SampleID") %>%
  column_to_rownames(var = "SampleID")

# filtering

filtered <- genus.ra.raw.meta %>% 
  # across the selected columns, check if at least 3 elements in each group are above 0.01
  summarise(across(where(is.numeric), ~ sum(.x > 0.01) >= 0.5*length(Group)), .by = Group) %>% 
  # Select the columns where any of the groups satisfy the conditions and extract the names
  select(where(~ is.logical(.x) && any(.x))) %>% names()

genus.ra <- genus.ra.raw.meta %>%
  select(all_of(filtered)) %>%
  mutate(hungatella_eubacterium_ratio = Hungatella/X.Eubacterium._siraeum_group) %>% 
  mutate(erysipelatoclostridium_eubacterium_ratio = Erysipelatoclostridium/X.Eubacterium._siraeum_group)

###############################
# Hung - Eub ratio from Greta
###############################

ratio_table <- read.xlsx("Hung_Eub_ratios_arranged_eray.xlsx", colNames = T)

ratio_table_meta <- ratio_table %>%
  left_join(., metadata %>% rownames_to_column("SampleID"), join_by("SampleID"))

# Three groups, HC as reference

ratio_table_meta_0_1 <- ratio_table_meta %>%
  mutate(Hung_ratio_0_1 = case_when(Hung_ratio < 50 ~ 0,
                                    Hung_ratio >= 50 ~ 1)) %>%
  mutate(Eub_ratio_0_1 = case_when(Eub_ratio < 50 ~ 0,
                                   Eub_ratio >= 50 ~ 1))

# openxlsx::write.xlsx(ratio_table_meta_0_1, "Results/genus/hung_eub_dominance_ratio_table_meta_0_1.xlsx")

m_hung <- glm(Hung_ratio_0_1 ~ Group + Age + Sex + Carbohydrate + Protein + Fat + Fiber, 
              data = ratio_table_meta_0_1,
              family = "binomial")

hung_dominance_all_three <- summary(m_hung)$coefficients %>% 
  as.data.frame %>%
  rownames_to_column("Predictor") %>%
  mutate(Response = rep("Hung/Eub Dominance", nrow(.)), .before = 1) %>%
  mutate(Significance = case_when(`Pr(>|z|)` >= 0.1 ~ "NS", # adding significance stars
                                   `Pr(>|z|)` >= 0.05 & `Pr(>|z|)` < 0.1 ~ ".",
                                   `Pr(>|z|)` >= 0.01 & `Pr(>|z|)` < 0.05 ~ "*",
                                   `Pr(>|z|)` >= 0.001 & `Pr(>|z|)` < 0.01 ~ "**",
                                   `Pr(>|z|)` < 0.001 ~ "***"))



# Two epilepsy group, MS as reference

ratio_table_meta_0_1_epi <- ratio_table_meta_0_1 %>%
  filter(Group %in% c("MS", "MR"))

m_hung_epi <- glm(Hung_ratio_0_1 ~ Group + Age + Sex + Carbohydrate + Protein + Fat + Fiber, 
                  data = ratio_table_meta_0_1_epi,
                  family = binomial)

hung_dominance_epi <- summary(m_hung_epi)$coefficients %>% 
  as.data.frame %>%
  rownames_to_column("Predictor") %>%
  mutate(Response = rep("Hung/Eub Dominance", nrow(.)), .before = 1) %>%
  mutate(Significance = case_when(`Pr(>|z|)` >= 0.1 ~ "NS", # adding significance stars
                                  `Pr(>|z|)` >= 0.05 & `Pr(>|z|)` < 0.1 ~ ".",
                                  `Pr(>|z|)` >= 0.01 & `Pr(>|z|)` < 0.05 ~ "*",
                                  `Pr(>|z|)` >= 0.001 & `Pr(>|z|)` < 0.01 ~ "**",
                                  `Pr(>|z|)` < 0.001 ~ "***"))

# openxlsx::write.xlsx(hung_dominance_all_three, "Results/genus/hung_eub_dominance_glm_all_three.xlsx")
# openxlsx::write.xlsx(hung_dominance_epi, "Results/genus/hung_eub_dominance_glm_epi.xlsx")

#################
# Visualizations
##################

########################################################################
# Creating a sample-wise distribution of Hung/Eub in isolated community
########################################################################

hung_eub_ratio_samples <- ratio_table_meta_0_1 %>%
 dplyr:: select(SampleID, Hung_ratio, Eub_ratio) %>%
  reshape2::melt() %>%
  left_join(., metadata %>% rownames_to_column("SampleID") %>% dplyr::select(SampleID, Group), join_by("SampleID")) %>%
  dplyr::rename("Ratio" = value)

p <- ggplot(hung_eub_ratio_samples, aes(x = SampleID, y = Ratio, fill = variable)) + 
  geom_bar(stat = "identity") + 
  facet_grid(. ~ Group, scale="free") + 
  theme(axis.line = element_line(colour = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 90),
        axis.text.y = element_text(size = 15),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 15, face = "italic"),
        legend.title = element_blank(),
        strip.text = element_text(size=15),
        legend.position='bottom') +
  labs(y = "Ratio (%)") +
  scale_fill_manual(values = c("#41AEAC", "#0F6A73"), 
                    breaks = c("Hung_ratio", "Eub_ratio"),
                    labels = c("Hungatella", "[Eubacterium] siraeum group"))

# using group colours for facet title boxes
g <- ggplot_gtable(ggplot_build(p))
strip_t <- which(grepl('strip-t', g$layout$name))
fills <- c("#F4B942", "#97D8C4", "#4059AD")
k <- 1
for (i in strip_t) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

png("Results/genus/hung_eub_abund_in_samples_plot.png", 
    units = "in", width = 15, height = 6, res = 800)
grid::grid.draw(g)
dev.off()


ratio_table_meta_0_1 %>%
  dplyr:: select(SampleID, Hung_ratio_0_1, Eub_ratio_0_1, Group)  %>%
  group_by(Group) %>%
  summarise_at(c("Hung_ratio_0_1", "Eub_ratio_0_1"), median, na.rm = TRUE)

library(gghalves)
library(ggpubr)

png("Results/genus/hung_abund_in_isol_comm_plot.png", 
    units = "in", width = 5, height = 6, res = 800)

ratio_table_meta_0_1 %>%
  dplyr:: select(SampleID, Hung_ratio, Eub_ratio) %>%
  left_join(., metadata %>% rownames_to_column("SampleID") %>% dplyr::select(SampleID, Group), join_by("SampleID")) %>%
  ggplot(., aes(x = Group, y = Hung_ratio, fill = Group, colour = Group)) + 
  geom_half_boxplot(side = "r", alpha = 0.8, outlier.shape = "triangle", colour = "black") +
  geom_half_point(side = "l", alpha = 0.8, size = 3, transformation = position_jitter(width = 0.1, height = 0, seed = 1)) +
  # theme_light() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        axis.title.y = element_text(size = 15),
        axis.title.x =element_blank(),
        axis.text.x = element_text(size = 15, vjust = 0.5, hjust=0.5),
        axis.text.y = element_text(size = 15),
        plot.title = element_text(size = 20),
        strip.text = element_text(size=15),
        legend.position='none') +
  theme(strip.text = element_text(colour = 'black')) +
  scale_fill_manual(values = c("#F4B942", "#97D8C4", "#4059AD")) + 
  scale_colour_manual(values = c("#F4B942", "#97D8C4", "#4059AD")) + 
  labs(y = "Hungatella Ratio (%) in Hung+Eub community") +
  geom_signif(comparisons = list(c("HC", "MS")), 
              y_position = 105,
              map_signif_level=function(p) {print("NS")}, 
              size = 1, textsize = 5, color = "black", ) +
  geom_signif(comparisons = list(c("HC", "MR")), 
              y_position = 116,
              map_signif_level=function(p) {print("**")}, 
              size = 1, textsize = 5, color = "black", ) +
  geom_signif(comparisons = list(c("MS", "MR")), 
              y_position = 110,
              map_signif_level=function(p) {print("*")}, 
              size = 1, textsize = 5, color = "black") +
  theme(plot.margin = margin(70, 10, 10, 10)) +
  coord_cartesian(clip = "off", ylim = c(0, 100))

dev.off()
