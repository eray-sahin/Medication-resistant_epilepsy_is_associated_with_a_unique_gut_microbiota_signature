rm(list = ls())

library(Maaslin2)
library(tidyverse)

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


genus.ra.no.zero <- genus.ra.raw.meta %>%
  select(all_of(filtered))
genus.ra.no.zero[genus.ra.no.zero == 0] <- min(genus.ra.no.zero[genus.ra.no.zero != 0])/2

genus.ra.no.zero <- genus.ra.no.zero %>%
  mutate(hungatella_eubacterium_no.zero_ratio = Hungatella/X.Eubacterium._siraeum_group) %>% 
  mutate(erysipelatoclostridium_no.zero_eubacterium_ratio = Erysipelatoclostridium/X.Eubacterium._siraeum_group)

genus.ra.and.ratio <- genus.ra.raw.meta %>%
  select(all_of(filtered)) %>%
  rownames_to_column("sample.id") %>%
  left_join(.,
            genus.ra.no.zero %>%
              rownames_to_column("sample.id") %>%
              dplyr::select(sample.id, hungatella_eubacterium_no.zero_ratio, erysipelatoclostridium_no.zero_eubacterium_ratio),
            join_by("sample.id")) %>%
  column_to_rownames("sample.id") 

openxlsx::write.xlsx(genus.ra %>% rownames_to_column("Sample_ID"), 
                     "Results/genus/Genus_rel_abundance_with_ratio.xlsx")

openxlsx::write.xlsx(genus.ra.and.ratio %>% rownames_to_column("Sample_ID"), 
                     "Results/genus/Genus_rel_abundance_with_no.zero_ratio.xlsx")

# Maaslin2 analysis using all three groups by assigning HC as the reference group
genus_three_groups = Maaslin2(input_data = genus.ra,
                              input_metadata = metadata %>%
                                select(Group, Age, Sex, Carbohydrate, Protein, Fat, Fiber),
                              min_prevalence = -Inf,
                              normalization  = "NONE",
                              transform = "NONE",
                              max_significance = 0.05,
                              output         = "Results/genus/Genus_MR_or_MS_vs_Healthy_Maaslin2", 
                              fixed_effects  = c("Group", "Age", "Sex", "Carbohydrate", "Protein", "Fat", "Fiber"),
                              standardize = TRUE,
                              reference      = c("Group,HC"),
                              analysis_method = "CPLM")

genus_three_groups = Maaslin2(input_data = genus.ra.and.ratio,
                              input_metadata = metadata %>%
                                select(Group, Age, Sex, Carbohydrate, Protein, Fat, Fiber),
                              min_prevalence = -Inf,
                              normalization  = "NONE",
                              transform = "NONE",
                              max_significance = 0.05,
                              output         = "Results/genus/Genus_and_No.Zero_Ratio_MR_or_MS_vs_Healthy_Maaslin2_no.method", 
                              fixed_effects  = c("Group", "Age", "Sex", "Carbohydrate", "Protein", "Fat", "Fiber"),
                              standardize = TRUE,
                              reference      = c("Group,HC"))


# Maaslin2 analysis for epilepsy groups by assigning DS as the reference group
genus_maaslin_drug = Maaslin2(input_data     = genus.ra %>%
                                rownames_to_column("SampleID") %>%
                                filter(SampleID %in% rownames(genus.ra.raw.meta[genus.ra.raw.meta$Group %in% c("MS", "MR"),])) %>%
                                column_to_rownames("SampleID"),
                              input_metadata = metadata %>%
                                filter(Group %in% c("MS", "MR")) %>%
                                select(Group, Age, Sex, Carbohydrate, Protein, Fat, Fiber),
                              min_prevalence = -Inf,
                              normalization  = "NONE",
                              transform = "NONE",
                              max_significance = 0.05,
                              output         = "Results/genus/Genus_MR_vs_MS_Maaslin2", 
                              fixed_effects  = c("Group", "Age", "Sex", "Carbohydrate", "Protein", "Fat", "Fiber"),
                              standardize = TRUE,
                              reference      = c("group,MS"),
                              analysis_method = "CPLM")
