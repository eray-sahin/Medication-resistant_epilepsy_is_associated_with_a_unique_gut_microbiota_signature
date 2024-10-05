rm(list = ls())

library(Maaslin2)
library(tidyverse)

# setwd(PATH)

# import phylum ra table
phylum.ra.raw <- read.table("asv.phylum.ra.tsv", 
                         header = T, row.names = 1, sep = "\t")

# import metadata, and sort based on the order of samples on phylum ra table
metadata <- readxl::read_xlsx("Epi_Diet_Final_Metadata.xlsx") %>% column_to_rownames(var = "SampleID")
metadata <- metadata[match(rownames(phylum.ra.raw), rownames(metadata)), ]
# factor the groups in desired order
metadata$Group <- factor(metadata$Group, levels = c("HC", "MS", "MR"))
# number of samples in each group
table(metadata$Group)

# checking order of the samples
all(rownames(metadata) %in% rownames(phylum.ra.raw)) #TRUE
all(rownames(metadata) == rownames(phylum.ra.raw)) #TRUE

# adding group information to phylum ra
phylum.ra.raw.meta <- phylum.ra.raw %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(., metadata %>% rownames_to_column(var = "SampleID") %>% select(SampleID, Group), by = "SampleID") %>%
  column_to_rownames(var = "SampleID")

# filtering phyla based on abundance distribution in groups

filtered <- phylum.ra.raw.meta %>% 
  # across the selected columns, check if at least 3 elements in each group are above 0.1
  summarise(across(where(is.numeric), ~ sum(.x > 0.1) >= 0.5*length(Group)), .by = Group) %>% 
  # Select the columns where any of the groups satisfy the conditions and extract the names
  select(where(~ is.logical(.x) && any(.x))) %>% names()

phylum.ra <- phylum.ra.raw.meta %>%
  select(all_of(filtered)) %>%
  mutate(bf.ratio = Firmicutes/Bacteroidota) # adding bf ratio


openxlsx::write.xlsx(phylum.ra, "Results/phylum/filtered_phylum_rel.abundance_with_bf.ratio.xlsx")


# Maaslin2 analysis using all three groups by assigning HC as the reference group
phylum_three_groups = Maaslin2(input_data = phylum.ra, 
                     input_metadata = metadata %>%
                       select(Group, Age, Sex, Carbohydrate, Protein, Fat, Fiber), 
                     min_prevalence = -Inf,
                     normalization  = "NONE",
                     transform = "NONE",
                     max_significance = 0.05,
                     output         = "Results/phylum/Phylum_MR_or_MS_vs_Healthy_Maaslin2", 
                     fixed_effects  = c("Group", "Age", "Sex", "Carbohydrate", "Protein", "Fat", "Fiber"),
                     reference      = c("Group,HC"),
                     standardize = TRUE,
                     analysis_method = "CPLM")


# Maaslin2 analysis for epilepsy groups by assigning DS as the reference group

phylum_maaslin_drug = Maaslin2(input_data     = phylum.ra %>%
                                 rownames_to_column("SampleID") %>%
                                 filter(SampleID %in% rownames(phylum.ra.raw.meta[phylum.ra.raw.meta$Group %in% c("MS", "MR"),])) %>%
                                 column_to_rownames("SampleID"),
                               input_metadata = metadata %>%
                                 filter(Group %in% c("MS", "MR")) %>%
                                 select(Group, Age, Sex, Carbohydrate, Protein, Fat, Fiber),
                               min_prevalence = -Inf,
                               normalization  = "NONE",
                               transform = "NONE",
                               max_significance = 0.05,
                               output         = "Results/phylum/Phylum_MR_vs_MS_Maaslin2",
                               fixed_effects  = c("Group", "Age", "Sex", "Carbohydrate", "Protein", "Fat", "Fiber"),
                               standardize = TRUE,
                               reference      = c("group,MS"),
                               analysis_method = "CPLM")


