rm(list = ls())

library(tidyverse)
library(readxl)
library(readr)
library(Maaslin2)

# setwd(PATH)

###########
# Import
###########

# pathway abundance table
kegg.abundance <- read_xlsx("picrust2_kegg_path_abundance.xlsx", col_names = T) %>% 
  column_to_rownames("pathway")

# pathway descriptions
kegg.descr <- read_xlsx("picrust2_kegg_pathway_descriptions.xlsx", col_names = T, na = "")

###########
# Pathway RA
###########

# calculate relative abundance of pathways
kegg.ra <- as.data.frame(vegan::decostand(kegg.abundance, MARGIN = 2, method = "total")*100) %>%
  as.data.frame()

write.table(t(kegg.ra),
            "picrust_kegg_pathway_ra.tsv",
            quote = F, sep = "\t", row.names = T, col.names = NA)


# import metadata, and sort based on the order of samples on kegg ra table
metadata <- readxl::read_xlsx("Epi_Diet_Final_Metadata_Jan24.xlsx") %>% column_to_rownames(var = "SampleID")
metadata <- metadata[match(colnames(kegg.ra), rownames(metadata)), ]
# factor the groups in desired order
metadata$Group <- factor(metadata$Group, levels = c("HC", "MS", "MR"))
# number of samples in each group
table(metadata$Group)

# checking order of the samples
all(rownames(metadata) %in% colnames(kegg.ra)) #TRUE
all(rownames(metadata) == colnames(kegg.ra)) #TRUE

# adding group information to pathway ra
kegg.ra.meta <- as.data.frame(t(kegg.ra)) %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(., metadata %>% rownames_to_column(var = "SampleID") %>% select(SampleID, Group), by = "SampleID") %>%
  column_to_rownames(var = "SampleID")

# filtering

filtered <- kegg.ra.meta %>% 
  # across the selected columns, check if at least 3 elements in each group are above 0.01
  summarise(across(where(is.numeric), ~ sum(.x > 0.001) >= 0.5*length(Group)), .by = Group) %>% 
  # Select the columns where any of the groups satisfy the conditions and extract the names
  select(where(~ is.logical(.x) && any(.x))) %>% names()

kegg.ra.filtered <- kegg.ra.meta %>%
  select(all_of(filtered))

kegg.ra.filtered_to_export <- kegg.ra.filtered %>%
  rownames_to_column("SampleID")

openxlsx::write.xlsx(kegg.ra.filtered_to_export,
                     "picrust_filtered_kegg_pathway_ra.xlsx")

# Maaslin2 analysis using all three groups by assigning HC as the reference group
kegg_three_groups = Maaslin2(input_data = kegg.ra.filtered,
                              input_metadata = metadata %>%
                                select(Group, Age, Sex, Carbohydrate, Protein, Fat, Fiber),
                              min_prevalence = -Inf,
                              normalization  = "NONE",
                              transform = "NONE",
                              max_significance = 0.05,
                              output         = "Results/picrust2_kegg/KEGG_MR_or_MS_vs_Healthy_Maaslin2", 
                              fixed_effects  = c("Group", "Age", "Sex", "Carbohydrate", "Protein", "Fat", "Fiber"),
                              standardize = TRUE,
                              reference      = c("Group,HC"),
                              analysis_method = "CPLM")


# Maaslin2 analysis for epilepsy groups by assigning DS as the reference group
kegg_maaslin_drug = Maaslin2(input_data     = kegg.ra.filtered %>%
                               rownames_to_column("SampleID") %>%
                                filter(SampleID %in% rownames(kegg.ra.meta[kegg.ra.meta$Group %in% c("MS", "MR"),])) %>%
                                column_to_rownames("SampleID"),
                              input_metadata = metadata %>%
                                filter(Group %in% c("MS", "MR")) %>%
                                select(Group, Age, Sex, Carbohydrate, Protein, Fat, Fiber),
                              min_prevalence = -Inf,
                              normalization  = "NONE",
                              transform = "NONE",
                              max_significance = 0.05,
                              output         = "Results/picrust2_kegg/KEGG_MR_vs_MS_Maaslin2", 
                              fixed_effects  = c("Group", "Age", "Sex", "Carbohydrate", "Protein", "Fat", "Fiber"),
                              standardize = TRUE,
                              reference      = c("group,MS"),
                              analysis_method = "CPLM")


# adding descriptions to the Maaslin2 output tables
three_groups_all_results <- read_tsv("Results/picrust2_kegg/KEGG_MR_or_MS_vs_Healthy_Maaslin2/all_results.tsv", col_names = T)

three_groups_all_results_descr <- three_groups_all_results %>%
  left_join(., kegg.descr, join_by("feature" == "pathway")) %>%
  relocate(description, .after = feature)

openxlsx::write.xlsx(three_groups_all_results_descr,
                     "Results/picrust2_kegg/KEGG_MR_or_MS_vs_Healthy_Maaslin2/all_results_with_path_descriptions.xlsx")

# no sign results, skipping that table for three groups


# Epi group tables
epi_groups_all_results <- read_tsv("Results/picrust2_kegg/KEGG_MR_vs_MS_Maaslin2/all_results.tsv", col_names = T)

epi_groups_all_results_descr <- epi_groups_all_results %>%
  left_join(., kegg.descr, join_by("feature" == "pathway")) %>%
  relocate(description, .after = feature)

openxlsx::write.xlsx(epi_groups_all_results_descr,
                     "Results/picrust2_kegg/KEGG_MR_vs_MS_Maaslin2/all_results_with_path_descriptions.xlsx")

# no sign results, skipping that table for epi groups
