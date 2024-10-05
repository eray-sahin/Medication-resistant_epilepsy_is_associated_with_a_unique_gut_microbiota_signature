rm(list = ls())

library(tidyverse)
library(readxl)
library(readr)
library(ggplot2)
library(gghalves)
library(ggsignif)


# setwd(PATH)

###########
# Import
###########

# pathway relative abundances
pathway_ra <- read_tsv("picrust_kegg_pathway_ra.tsv", col_names = T) %>%
  column_to_rownames("...1")

# pathway descriptions
kegg.descr <- read_xlsx("picrust2_kegg_pathway_descriptions.xlsx", col_names = T, na = "")

# selected pathways
pathways_of_interest <- c("Fructose and mannose metabolism", 
                          "Fatty acid biosynthesis",
                          "D-Glutamine and D-glutamate metabolism")
# assigning kegg IDs as names
names(pathways_of_interest) <- kegg.descr$pathway[match(pathways_of_interest, kegg.descr$description)]

# metadata, and sort based on the order of samples on kegg ra table
metadata <- readxl::read_xlsx("Epi_Diet_Final_Metadata.xlsx") %>% column_to_rownames(var = "SampleID")
# factor the groups in desired order
metadata$Group <- factor(metadata$Group, levels = c("HC", "MS", "MR"))

##################
# DF Preparations
##################

# filtering the pathways of interest, adding group info to pathway ra, melting table, formatting
poi_ra_meta_melted <- pathway_ra %>%
  dplyr::select(all_of(names(pathways_of_interest))) %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(., metadata %>% rownames_to_column(var = "SampleID") %>% select(SampleID, Group), by = "SampleID") %>%
  reshape2::melt(., id = c("SampleID", "Group")) %>%
  dplyr::rename("Pathway_ID" = "variable") %>%
  dplyr::rename("Relative_Abundance" = "value") %>%
  left_join(., kegg.descr, join_by("Pathway_ID" == "pathway"))

###########
# Plotting
###########

# preparing annotation df for ggsignif as described on https://github.com/const-ae/ggsignif

# import maaslin2 result table and filter
poi_significance <- readxl::read_xlsx("PiCRUST_Maaslin2_POI_Sig.xlsx", col_names = T)

annotation_df <- data.frame(
  description = c(rep("Fructose and mannose metabolism", 2),
            "Fatty acid biosynthesis",
            "D-Glutamine and D-glutamate metabolism"),
  start = c("HC", "MS", "HC", "MS"),
  end = c("MR", "MR", "MS", "MR"),
  y = c(1.75, 1.5, 2, 2.6),
  Comparison = c("DR_vs_HC", 
            "MR_vs_MS",
            "MS_vs_HC",
            "MR_vs_DMS")) %>% 
  left_join(., poi_significance %>% dplyr::select(Comparison, description, pval, qval), 
            by = c("Comparison", "description")) %>%
  mutate(label = paste0("pval=", round(pval, 3), "\np.adj=", round(qval, 3)))

plot <- ggplot(poi_ra_meta_melted,
               aes(x = Group, y = Relative_Abundance)) + 
  geom_half_boxplot(side = "r", alpha = 0.8, outlier.shape = "triangle", colour = "black",
                    aes(colour = Group, fill = Group)) +
  geom_half_point(side = "l", alpha = 0.8, size = 3, transformation = position_jitter(width = 0.1, height = 0, seed = 1),
                  aes(colour = Group)) +
  ggsignif::geom_signif(
    data = annotation_df,
    aes(xmin = start, xmax = end, annotations = label, y_position = y),
    textsize = 4, vjust = -0.2,
    manual = TRUE
  ) +
  ylim(NA, 3) + 
  facet_grid(~ description, scale="free", labeller = label_wrap_gen(width = 16)) + 
  theme(axis.line = element_line(colour = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15, angle = 90),
        axis.text.y = element_text(size = 15),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 15, face = "italic"),
        legend.title = element_blank(),
        strip.text = element_text(size=15),
        legend.position='bottom') +
  labs(y = "Relative Abundance (%)") +
  scale_fill_manual(values = c("#F4B942", "#97D8C4", "#4059AD")) + 
  scale_colour_manual(values = c("#F4B942", "#97D8C4", "#4059AD"))

# using group colours for facet title boxes
g <- ggplot_gtable(ggplot_build(plot))
strip_t <- which(grepl('strip-t', g$layout$name))
fills <- c("white", "white", "white")
k <- 1
for (i in strip_t) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

png("Results/picrust2_kegg/picrust_poi_plot.png", 
    units = "in", width = 12, height = 8, res = 800)
grid::grid.draw(g)
dev.off()
