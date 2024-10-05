rm(list = ls())

library(phyloseq)
library(tidyverse)
library(ggplot2)
library(readr)
library(Maaslin2)
library(readxl)
library(egg)

# setwd(PATH)

# import alpha diversity indices 
alpha_diversity_df <-read_tsv("all_alpha_div_all_indices_calculated.tsv", col_names = T)
# factor and order of groups
alpha_diversity_df$Group <- factor(alpha_diversity_df$Group, levels = c("HC", "MS", "MR"))

# import GLM output for all three groups, significant only
sign_alpha <- readxl::read_xlsx("Alpha_Diversity_MR_or_MS_vs_Healthy_GLM/alpha_all_three_groups_glm.xlsx", col_names = T)

######################
# Correlation Plots
######################

# Age vs Chao1
age_Chao1 <- ggplot(alpha_diversity_df, 
                    aes(x = Age, 
                        y = Chao1, 
                        colour = Group)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_smooth(method = "glm", colour = "black") +
  theme_light() +
  theme(axis.line = element_line(colour = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        strip.text = element_text(size=15),
        legend.position='right',
        legend.title = element_text(size=15),
        legend.text = element_text(size=15)) +
  theme(strip.background =element_rect(fill="white", colour = "black")) +
  theme(strip.text = element_text(colour = 'black')) +
  # Adding text box for maaslin2 results
  annotate('text', 
           x = Inf,
           y = -Inf,
           hjust = 1,
           vjust = -0.25,
           label = paste0("Coefficient: ", 
                          round(sign_alpha$Estimate[sign_alpha$Response == "Chao1" & sign_alpha$Predictor == "Age"], 5),
                          "\np.adj: ",
                          round(sign_alpha$P.adj[sign_alpha$Response == "Chao1" & sign_alpha$Predictor == "Age"], 5)), 
           fontface = 'bold',
           size = 5) +
  scale_fill_manual(values = c("#F4B942", "#97D8C4", "#4059AD")) + 
  scale_colour_manual(values = c("#F4B942", "#97D8C4", "#4059AD"))



# Age vs Shannon
age_Shannon <- ggplot(alpha_diversity_df, 
                    aes(x = Age, 
                        y = Shannon, 
                        colour = Group)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_smooth(method = "glm", colour = "black") +
  theme_light() +
  theme(axis.line = element_line(colour = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        strip.text = element_text(size=15),
        legend.position='right',
        legend.title = element_text(size=15),
        legend.text = element_text(size=15)) +
  theme(strip.background =element_rect(fill="white", colour = "black")) +
  theme(strip.text = element_text(colour = 'black')) +
  # Adding text box for maaslin2 results
  annotate('text', 
           x = Inf,
           y = -Inf,
           hjust = 1,
           vjust = -0.25,
           label = paste0("Coefficient: ", 
                          round(sign_alpha$Estimate[sign_alpha$Response == "Shannon Index" & sign_alpha$Predictor == "Age"], 5),
                          "\np.adj: ",
                          round(sign_alpha$P.adj[sign_alpha$Response == "Shannon Index" & sign_alpha$Predictor == "Age"], 5)), 
           fontface = 'bold',
           size = 5) +
  scale_fill_manual(values = c("#F4B942", "#97D8C4", "#4059AD")) + 
  scale_colour_manual(values = c("#F4B942", "#97D8C4", "#4059AD"))


# Age vs Simpson
age_Simpson <- ggplot(alpha_diversity_df, 
                      aes(x = Age, 
                          y = Simpson, 
                          colour = Group)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_smooth(method = "glm", colour = "black") +
  theme_light() +
  theme(axis.line = element_line(colour = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        strip.text = element_text(size=15),
        legend.position='right',
        legend.title = element_text(size=15),
        legend.text = element_text(size=15)) +
  theme(strip.background =element_rect(fill="white", colour = "black")) +
  theme(strip.text = element_text(colour = 'black')) +
  # Adding text box for maaslin2 results
  annotate('text', 
           x = Inf,
           y = -Inf,
           hjust = 1,
           vjust = -0.25,
           label = paste0("Coefficient: ", 
                          round(sign_alpha$Estimate[sign_alpha$Response == "Simpson's Index" & sign_alpha$Predictor == "Age"], 5),
                          "\np.adj: ",
                          round(sign_alpha$P.adj[sign_alpha$Response == "Simpson's Index" & sign_alpha$Predictor == "Age"], 5)), 
           fontface = 'bold',
           size = 5) +
  scale_fill_manual(values = c("#F4B942", "#97D8C4", "#4059AD")) + 
  scale_colour_manual(values = c("#F4B942", "#97D8C4", "#4059AD"))

# Combine

library(patchwork)
png("Plots/Age_and_three_alpha_indices_corr_plots_glm_merged.png",
    units="in", width=15, height = 7, res = 800)

age_Chao1 + age_Shannon + age_Simpson + plot_layout(guides = "collect")

dev.off()
