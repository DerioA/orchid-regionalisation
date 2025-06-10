# ==============================================================================
# 2.- Data wrangling
# ==============================================================================
## 2.1.- Data collection: 
# We used georeferenced occurrences derived from herbarium data from the:
# - Global Biodiversity Information Facility (GBIF.ORG, 2024),
# - speciesLink (speciesLink network, 2025),
# - RAINBIO (Dauby et al., 2016), and
# - Atlas of Living Australia (ALA, 2025).
# ==============================================================================
# PACKAGES
library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
# ==============================================================================
# 1.- GBIF database
# ==============================================================================

orchid_gbif <- fread("raw-data/gbif/0052485-241126133413365.csv") # 2,172,761 records

# Exclude records at family and genus level, as the analyses are species-based.
table(orchid_gbif$taxonRank)
orchid_gbif <- orchid_gbif[taxonRank %in% c("SPECIES", "SUBSPECIES", "VARIETY", "FORM")] # 1,910,260 records
# 262,501 records removed
table(orchid_gbif$taxonRank)

# Keep only Orchidaceae
table(orchid_gbif$family)
orchid_gbif <- orchid_gbif[taxonRank %in% c("SPECIES", "SUBSPECIES", "VARIETY", "FORM")]
table(orchid_gbif$taxonRank)

# Keep only records with coordinates
orchid_gbif <- orchid_gbif %>% 
  filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) # 907,965 records

# Eliminate inconsistencies
names(orchid_gbif)
orchid_gbif <- orchid_gbif %>% filter(!is.na(genus) & genus != "") # 907,957 records
orchid_gbif <- orchid_gbif %>% filter(!str_detect(species, "×\\s?")) # 907,735 records
orchid_gbif <- orchid_gbif %>% filter(!str_detect(verbatimScientificName, "×\\s?")) # 907,032 records
orchid_gbif <- orchid_gbif %>% filter(!str_detect(scientificName, "×\\s?")) # 903,659 records

# ==============================================================================
# 2.- RAINBIO database
# ==============================================================================

orchid_rainbio <- fread("raw-data/rainbio_database-2/RAINBIO.csv") # 157,885 records

# Keep only Orchidaceae records
table(orchid_rainbio$family)
orchid_rainbio <- orchid_rainbio[family %in% c("Orchidaceae")] # 18,602 records

# Keep only records with coordinates
orchid_rainbio <- orchid_rainbio %>% 
  filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) # 18,602 records

# ==============================================================================
# 3.- speciesLink database
# ==============================================================================

orchid_links <- fread("raw-data/speciesLink/speciesLink-20250429025336-0006433.txt") # 306,822 records

# Keep only records with coordinates
names(orchid_links)
orchid_links <- orchid_links %>% 
  filter(!is.na(longitude) & !is.na(latitude)) # 306,822 records

# Keep only Orchidaceae
table(orchid_links$family)
orchid_links <- orchid_links[family %in% c("Orchidaceae")] # 306,179 records

# Keep only preserved specimens
table(orchid_links$basisofrecord)
orchid_links <- orchid_links[basisofrecord %in% c("PreservedSpecimen")] # 290,470 records

# Keep only species-level records
names(orchid_links)
orchid_links <- orchid_links %>% filter(!is.na(species) & species != "") # 244,774 records

# Eliminate hybrids
orchid_links <- orchid_links %>% 
  filter(!str_detect(species, regex("\\(h[ií]brido\\)", ignore_case = TRUE))) # 244,757 records
orchid_links <- orchid_links %>% 
  filter(!grepl("× ", species, fixed = TRUE))

# Additional filters for inconsistent names
orchid_links <- orchid_links %>%  filter(species != "×altaflorestense")
orchid_links <- orchid_links %>%  filter(species != "×apolloi")
orchid_links <- orchid_links %>%  filter(species != "×canaense")
orchid_links <- orchid_links %>%  filter(species != "×freitasii")
orchid_links <- orchid_links %>%  filter(species != "×pabstii")
orchid_links <- orchid_links %>%  filter(species != "×perimii")
orchid_links <- orchid_links %>%  filter(species != "×roseo-album")
orchid_links <- orchid_links %>%  filter(species != "×sgarbii")
orchid_links <- orchid_links %>%  filter(species != "×tapiriceps")
orchid_links <- orchid_links %>%  filter(species != "×valdisonianum")
orchid_links <- orchid_links %>%  filter(species != "3 tipos")
orchid_links <- orchid_links %>%  filter(species != "s")
orchid_links <- orchid_links %>% 
  filter(!str_detect(species, "\\d{1,3}°\\d{1,2}'\\d{1,2}\""))
orchid_links <- orchid_links %>%  filter(species != "sp.") # 237,469 records

# ==============================================================================
# 4.- ALA database
# ==============================================================================

orchid_ala <- fread("raw-data/records-ALA-2025-04-29/records-2025-04-29.csv", 
                    drop = "images") # 573,330 records

# Keep only records with coordinates
names(orchid_ala)
orchid_ala <- orchid_ala %>% 
  filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) # 535,264 records

# Keep only Orchidaceae
table(orchid_ala$family)

# Keep only species-level records
table(orchid_ala$taxonRank)
orchid_ala <- orchid_ala[taxonRank %in% c("species", "subspecies", "variety")] # 519,681 records

# Keep only preserved specimens
table(orchid_ala$basisOfRecord)
orchid_ala <- orchid_ala[basisOfRecord %in% c("PRESERVED_SPECIMEN")] # 155,175 records
