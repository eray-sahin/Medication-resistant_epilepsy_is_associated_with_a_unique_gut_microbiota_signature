rm(list = ls())

library(Maaslin2)
library(tidyverse)

# setwd(PATH)

# import family ra table
family.ra.raw <- read.table("asv.family.ra.tsv", 
                            header = T, row.names = 1, sep = "\t")

# import metadata, and sort based on the order of samples on family ra table
metadata <- readxl::read_xlsx("Epi_Diet_Final_Metadata.xlsx") %>% column_to_rownames(var = "SampleID")
metadata <- metadata[match(rownames(family.ra.raw), rownames(metadata)), ]
# factor the groups in desired order
metadata$Group <- factor(metadata$Group, levels = c("HC", "MS", "MR"))
# number of samples in each group
table(metadata$Group)

## checking order of the samples
all(rownames(metadata) %in% rownames(family.ra.raw)) #TRUE
all(rownames(metadata) == rownames(family.ra.raw)) #TRUE

# adding group information to family ra
family.ra.raw.meta <- family.ra.raw %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(., metadata %>% rownames_to_column(var = "SampleID") %>% select(SampleID, Group), by = "SampleID") %>%
  column_to_rownames(var = "SampleID")

# filtering

filtered <- family.ra.raw.meta %>% 
  # across the selected columns, check if at least 3 elements in each group are above 0.01
  summarise(across(where(is.numeric), ~ sum(.x > 0.01) >= 0.5*length(Group)), .by = Group) %>% 
  # Select the columns where any of the groups satisfy the conditions and extract the names
  select(where(~ is.logical(.x) && any(.x))) %>% names()

family.ra <- family.ra.raw.meta %>%
  select(all_of(filtered))

openxlsx::write.xlsx(family.ra, "Results/family/filtered_family_rel.abundance.xlsx")


# Maaslin2 analysis using all three groups by assigning HC as the reference group
family_three_groups = Maaslin2(input_data = family.ra, 
                               input_metadata = metadata %>%
                                 select(Group, Age, Sex, Carbohydrate, Protein, Fat, Fiber), 
                               min_prevalence = -Inf,
                               normalization  = "NONE",
                               transform = "NONE",
                               max_significance = 0.05,
                               output         = "Results/family/Family_MR_or_MS_vs_Healthy_Maaslin2", 
                               fixed_effects  = c("Group", "Age", "Sex", "Carbohydrate", "Protein", "Fat", "Fiber"),
                               reference      = c("Group,HC"),
                               standardize = TRUE,
                               analysis_method = "CPLM")



# Maaslin2 analysis for epilepsy groups by assigning DS as the reference group
family_maaslin_drug = Maaslin2(input_data     = family.ra %>%
                                 rownames_to_column("SampleID") %>%
                                 filter(SampleID %in% rownames(family.ra.raw.meta[family.ra.raw.meta$Group %in% c("MS", "MR"),])) %>%
                                 column_to_rownames("SampleID"),
                               input_metadata = metadata %>%
                                 filter(Group %in% c("MS", "MR")) %>%
                                 select(Group, Age, Sex, Carbohydrate, Protein, Fat, Fiber),
                               min_prevalence = -Inf,
                               normalization  = "NONE",
                               transform = "NONE",
                               max_significance = 0.05,
                               output         = "Results/family/Family_MR_vs_MS_Maaslin2",
                               fixed_effects  = c("Group", "Age", "Sex", "Carbohydrate", "Protein", "Fat", "Fiber"),
                               standardize = TRUE,
                               reference      = c("group,MS"),
                               analysis_method = "CPLM")
