rm(list = ls())

library(tidyverse)
library(readxl)
library(readr)

# setwd(PATH)

###########
# Import
###########

# feature count table
asv.counts <- read_xlsx("epi_diet_table_feature_table_arranged.xlsx", col_names = T) %>% 
  column_to_rownames("ASV")

# ASV taxonomy table
asv.tax <- read_xlsx("epi_diet_taxonomy_arranged.xlsx", col_names = T, na = "")
# Replacing NA with Unidentified
asv.tax[is.na(asv.tax)] <- "Unidentified"

###########
# ASV RA
###########

# calculate relative abundance of ASVs, then aggregate into taxa levels, and remove Eukaryota, Archaea, Mitochondria
asv.ra <- as.data.frame(vegan::decostand(asv.counts, MARGIN = 2, method = "total")*100) %>%
  rownames_to_column("ASV") %>%
  left_join(., 
            asv.tax,
            by = "ASV") %>%
  filter(Kingdom == "Bacteria") %>%
  filter(Genus != "Mitochondria") %>%
  column_to_rownames(var = "ASV") %>%
  as.data.frame()

write.table(t(asv.ra),
            "asv.ra.tsv",
            quote = F, sep = "\t", row.names = T, col.names = NA)

###########
# Phylum
###########

# aggregate into phylum level and remove unidentified 
asv.phylum.ra <- asv.ra %>%
  rownames_to_column("ASV") %>%
  group_by(Phylum) %>%
  summarise_if(is.numeric, sum) %>%
  filter(Phylum != "Unidentified") %>%
  column_to_rownames(var = "Phylum") %>%
  as.data.frame()

write.table(t(asv.phylum.ra),
            "asv.phylum.ra.tsv",
            quote = F, sep = "\t", row.names = T, col.names = NA)


###########
# Family
###########

# aggregate into family level and remove unidentified and uncultured
asv.family.ra <- asv.ra %>%
  rownames_to_column("ASV") %>%
  group_by(Family) %>%
  summarise_if(is.numeric, sum) %>%
  filter(Family != "Unidentified") %>%
  filter(Family != "uncultured") %>%
  column_to_rownames(var = "Family") %>%
  as.data.frame()

write.table(t(asv.family.ra),
            "asv.family.ra.tsv",
            quote = F, sep = "\t", row.names = T, col.names = NA)

###########
# Genus
###########

# aggregate into genus level and remove unidentified 
asv.genus.ra <- asv.ra %>%
  rownames_to_column("ASV") %>%
  group_by(Genus) %>%
  summarise_if(is.numeric, sum) %>%
  filter(Genus != "Unidentified") %>%
  filter(Genus != "uncultured") %>%
  column_to_rownames(var = "Genus") %>%
  as.data.frame()

write.table(t(asv.genus.ra),
            "asv.genus.ra.tsv",
            quote = F, sep = "\t", row.names = T, col.names = NA)


