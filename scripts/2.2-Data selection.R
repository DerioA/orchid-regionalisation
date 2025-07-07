# ==============================================================================
# 2.- Data wrangling
# ==============================================================================
# ==============================================================================
## **2.2.- Data selection**
# ==============================================================================
# ==============================================================================
### *Geographic inaccuracies:*
# We identified and discarded occurrence records meeting one or more of the following criteria:
# - Equal latitude and longitude
# - Zero latitude and longitude
# - Coordinates within 5-km radius of country centroids
# - Coordinates within 10-km radius of capital cities
# - Coordinates within 2-km radius of biodiversity institutions
# - Marine coordinates (using CoordinateCleaner package; Zizka et al., 2019)
# ==============================================================================
### *Taxonomic standard:*
#After identifying and discarding occurrence records with spatial errors,
#we updated the taxonomy and revised the native distribution of the species.
#==============================================================================
################################################################################

## Load packages
library(CoordinateCleaner) #Automated Cleaning of Occurrence Records from Biological Collections
options("sp_evolution_status" = 2)
library(sp) #Spatial Data
library(ggplot2) #Data Visualisations 
library(dplyr) #Data Manipulation
library(stringr) #Wrappers for common string operations
library(tidyr) #Data Manipulation
library(data.table) #Fast aggregation of large data
library(rWCVP) #It includes functions to generate maps and species lists, as well as match names to the WCVP
library(rWCVPdata)
library(purrr) #Functional Programming Tools
library(furrr) #Apply Mapping Functions in Parallel using Futures
library(sf) #Simple Features for R
library(readr) #Read Rectangular Text Data

################################################################################
# We analyse each dataset separately and then join them together remove spatial errors
## Process GBIF database
flags_gbif <- clean_coordinates(
  x = orchid_gbif,#This database comes from "2.1.- Data collection"
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  countries = "countryCode",
  species = "species",
  tests = c("capitals", "centroids", "equal", "gbif", "institutions", "seas", "zeros"))

summary(flags_gbif)
plot(flags_gbif, lon = "decimalLongitude", lat = "decimalLatitude")

# Exclude problematic records
clear_gbif <- flags_gbif[flags_gbif$.summary,]
# The flagged records
Notclear_gbif <- flags_gbif[!flags_gbif$.summary,] # 72,667 excluded records

# Visualise the clean data on a map
wm <- borders("world", colour = "gray50", fill = "gray50")
ggplot() +
  coord_fixed() +
  wm +
  geom_point(data = clear_gbif,aes(x = decimalLongitude, y = decimalLatitude),
    colour = "darkgreen",
    size = 0.5) +
  theme_bw()

################################################################################
## Process RAINBIO database
# Remove spatial error - rainbio database
names(orchid_rainbio)
flags_rainbio <- clean_coordinates(
  x = orchid_rainbio,#This database comes from "2.1.- Data collection"
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  countries = "countryCode",
  species = "tax_sp_level",
  tests = c("capitals", "centroids", "equal", "gbif", "institutions", "seas", "zeros"))

summary(flags_rainbio)
plot(flags_rainbio, lon = "decimalLongitude", lat = "decimalLatitude")

# Exclude problematic records
clear_rainbio <- flags_rainbio[flags_rainbio$.summary,]
# The flagged records
Notclear_rainbio <- flags_rainbio[!flags_rainbio$.summary,] # 530 excluded records

# Visualise the clean data on a map
ggplot() +
  coord_fixed() +
  wm +
  geom_point(data = clear_rainbio,aes(x = decimalLongitude, y = decimalLatitude),
    colour = "darkgreen",
    size = 0.5) +
  theme_bw()

################################################################################
## Process SpeciesLink database
# Remove spatial error - specieslinks database
# Create a column with the specific epithet
orchid_links <- orchid_links %>% #This database comes from "2.1.- Data collection"
  mutate(
    scientificName.new = str_glue("{genus} {species}") %>%  # Combine with format
      as.character(),                                     # Convert to text
    scientificName.new = str_replace_all(scientificName.new, "  ", " "))  # Remove double spaces

names(orchid_links)
flags_links <- clean_coordinates(
  x = orchid_links,
  lon = "longitude",
  lat = "latitude",
  species = "scientificName.new",
  tests = c("capitals", "centroids", "equal", "gbif", "institutions", "seas", "zeros"))

# Clear invalid coordinates
inv <- c(10526, 10535, 10553, 10555, 10565, 10569, 10577, 10588, 
         10595, 10607, 10615, 10626, 10658, 10687, 10691, 11653, 
         38564, 149425, 149449, 149485, 149517, 149518, 149570, 
         149587, 149593, 149599, 152584, 152585, 152592, 152609, 
         152627, 152630)
orchid_links <- orchid_links[-inv,] # 237,437 records
# re-run lines 107-112

summary(flags_links)
plot(flags_links, lon = "longitude", lat = "latitude")

# Exclude problematic records
clear_links <- flags_links[flags_links$.summary,]
# The flagged records
Notclear_links <- flags_links[!flags_links$.summary,] # 108,729 excluded records

# Visualise the clean data on a map
ggplot() +
  coord_fixed() +
  wm +
  geom_point(data = clear_links,aes(x = longitude, y = latitude),
    colour = "darkgreen",
    size = 0.5) +
  theme_bw()

################################################################################
## Process ALA database
# Remove spatial error - ALA database
names(orchid_ala)
flags_ala <- clean_coordinates(
  x = orchid_ala,#This database comes from "2.1.- Data collection"
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  countries = "countryCode",
  species = "scientificName",
  tests = c("capitals", "centroids", "equal", "gbif", "institutions", "seas", "zeros"))

summary(flags_ala)
plot(flags_ala, lon = "decimalLongitude", lat = "decimalLatitude")

# Exclude problematic records
clear_ala <- flags_ala[flags_ala$.summary,]
# The flagged records
Notclear_ala <- flags_ala[!flags_ala$.summary,] # 10,510 excluded records

# Visualise the clean data on a map
ggplot() +
  coord_fixed() +
  wm +
  geom_point(data = clear_ala,aes(x = decimalLongitude, y = decimalLatitude),
    colour = "darkgreen",
    size = 0.5) +
  theme_bw()

# After checking for spatial errors, we will merge the databases
# with the help of the bdc package, and we will keep the columns
# that will later help us in the subsequent analyses.

# First, save the databases
write.csv(clear_gbif, file = "processed-data/dataselection_spatial-errors/orchid-clear-gbif.csv")
write.csv(clear_rainbio, file = "processed-data/dataselection_spatial-errors/orchid-clear-rainbio.csv")
write.csv(clear_links, file = "processed-data/dataselection_spatial-errors/orchid-clear-links.csv")
write.csv(clear_ala, file = "processed-data/dataselection_spatial-errors/orchid-clear-ala.csv")

#######################join the databases#######################
################################################################
# The bdc package needs a .csv file with the metadata to join the databases.
# Let's create it and add the columns we want to get from all the databases. 
# The metadata file has to be saved in the same folder where the databases to be joined are hosted.
rm(list = ls()) # Since we are saving our databases that we will use from now on, we can delete everything we generated in R. 
metadata <- read.csv("processed-data/dataselection_spatial-errors/metadata.csv")
metadata <- metadata %>% mutate_all(na_if,"")
occ_orchid <-
  bdc_standardize_datasets(metadata = metadata,
                           format = "csv",
                           overwrite = FALSE,
                           save_database = TRUE)

write.csv(occ_orchid, file = "processed-data/dataselection_spatial-errors/final-occ-orchid.csv")

# At this point, it is important to make decisions.
# We work with 'orchid_clear', although authors can manually check
# the records with spatial problems, as well as add or remove
# the tests they want to analyse within "clean_coordinates function". 

#######################duplicates#######################
################################################################
# Duplicated records were removed to avoid redundancy
# Check the number of species before and after removing duplicates
# to be sure that we have removed duplicates but not species.
names(occ_orchid)
occ_orchid %>%
  dplyr::select(scientificName) %>%
  unique() %>% # selecting columns with species names
  as.data.frame(.) -> splist_orchid # 24,132 species

# We will first standardize the names of collectors, collector number and year
# collector
occ_orchid$recordedBy[is.na(occ_orchid$recordedBy)] <- "without-collector"
occ_orchid$recordedBy[occ_orchid$recordedBy == ""] <- "without-collector"
table(occ_orchid$recordedBy)

# collector number
occ_orchid$recordNumber[is.na(occ_orchid$recordNumber)] <- "without-number"
occ_orchid$recordNumber[occ_orchid$recordNumber == ""] <- "without-number"
table(occ_orchid$recordNumber)

# year
occ_orchid$year[is.na(occ_orchid$year)] <- "without-year"
occ_orchid$year[occ_orchid$year == ""] <- "without-year"
table(occ_orchid$year)

# Duplicates
names(occ_orchid)
occ_orchid_dups <- occ_orchid %>%
  distinct(scientificName, recordNumber, recordedBy, year, .keep_all = TRUE) # 750,050

# Number of species before removing duplicates
occ_orchid_dups %>%
  dplyr::select(scientificName) %>%
  unique() %>% # selecting columns with species names
  as.data.frame(.) -> splist_orchid_2 # 24,132 species

write.csv(occ_orchid_dups, file = "processed-data/dataselection_spatial-errors/occ_orchid_dups.csv")

######################Taxonomic standardization#######################
####################################################################
# We standardized botanical nomenclature according to the World Checklist of Vascular Plants (WCVP; Govaerts et al., 2021),
# using the rWCVP package (Brown et al. 2023).

rm(list=ls())
orchids_dups<- fread("processed-data/dataselection_spatial-errors/occ_orchid_dups.csv",
                     na.strings = "") # 750,050 records
orchids_dups <- select(orchids_dups, -V1)

# Remove affinities (aff.) and closeness of (cf.) from original nouns before updating
orchids_dups <- orchids_dups %>% 
  filter(!str_detect(scientificName, " aff\\.")) %>% 
  filter(!str_detect(scientificName, " cf\\."))

# update
names(orchids_dups)
orchids.wcvp <- wcvp_match_names(
  orchids_dups,
  wcvp_names = NULL,
  name_col = "scientificName",
  id_col = "database_id",
  author_col = "authority",
  fuzzy = FALSE)

# Subsequent review
table(orchids.wcvp$wcvp_status)
# Finally, we only use accepted (705,989 records) and synonym (37,099) records.
# We do not use species marked as "illegitimate" (14,067), 'artificial hybrid' (14),
# 'invalid' (782), 'misapplied' (4,490), 'orthographic' (39) and 'unplaced' (986).
occ_accepted <- orchids.wcvp %>% filter(wcvp_status %in% c("Accepted"))
table(occ_accepted$match_type)

# We will update the records marked as synonyms in WCVP using U.Taxonstand to get their accepted names.
occ_synonym <- orchids.wcvp %>% filter(wcvp_status %in% c("Synonym"))

# Because we had problems with the U.taxonstand package,
# we decided to use the online shiny application to update
# the synonym records: https://ecoinfor.shinyapps.io/UTaxonstandOnline/
synonyms_updated <- read.csv(file = "processed-data/dataselection_taxonomic-standard/synonyms_updated.csv")
names(synonyms_updated)
synonyms_updated <- synonyms_updated %>%
  rename(scientificName = Submitted_Name)

# Update the names in occ_synonym with the new names updated with U.taxonstand (base synonyms_updated).
# Step 1: Join the databases keeping all records of occ_synonym
occ_updated <- occ_synonym %>%
  left_join(synonyms_updated %>% 
              select(scientificName, Accepted_SPNAME),
            by = "scientificName")

# Step 2: Update wcvp_name column with accepted names where available
occ_updated <- occ_updated %>%
  mutate(wcvp_name = ifelse(!is.na(Accepted_SPNAME), 
                            Accepted_SPNAME, 
                            wcvp_name))
# Step 3: Deleted Accepted_SPNAME column
occ_updated <- occ_updated %>%
  select(-Accepted_SPNAME)

# Verify results
head(occ_updated)
table(occ_updated$wcvp_name == occ_synonym$wcvp_name, useNA = "always")

# Finally, delete some records that after the update are species belonging to other families.
names <- c("Didymopanax selloi", "Didymopanax morototoni", 
           "Didymopanax vinosus", "Ixeris polycephala", 
           "Erucastrum gallicum", "Stigmatodon pseudoliganthus", 
           "Miconia stenopetala", "Eugenia monticola")
# Deleted records
occ_updated <- occ_updated %>%
  filter(!wcvp_name %in% names)

# Join the base of accepted names with the base of synonymous names that have been updated.
occ_final <- bind_rows(occ_accepted, occ_updated)
dim(occ_final)
table(occ_final$wcvp_status)
occ_final <- occ_final %>%
  mutate(wcvp_status = ifelse(wcvp_status == "Synonym", "Accepted", wcvp_status))
table(occ_final$wcvp_status)

# There are still hybrids in the wcvp_name column, they must be removed.
# Finally, delete some records that after the update are species belonging to other families.
names <- c("× Cyanthera glossodioides","× Dactylocamptis drudei",
           "× Laeliocattleya calimaniana", "× Laeliocattleya menezesiana",
           "× Phelodia tutelata","× Catyclia intermedia","× Brassocattleya arauji")
# Deleted records
occ_final <- occ_final %>%
  filter(!wcvp_name %in% names)
# More hybrids
occ_final <- occ_final %>%
  filter(!str_detect(wcvp_name, "\\s×\\s"))  # × between spaces

# IF there are "NAs" remove
occ_final <- occ_final[!is.na(occ_final$wcvp_name),]

# Separate "wcvp_name" column to access the genera
occ_final <- occ_final %>% separate(wcvp_name, c('wcvp_genus', 'wcvp_epithet'), remove = FALSE)

# Review genus 
# genus
occ_final <- occ_final[!(occ_final$wcvp_genus == "Platylepis"),]
occ_final <- occ_final[!(occ_final$wcvp_genus == "Phaius"),]

# Number of species in clear database
occ_final %>%
  dplyr::select(wcvp_name) %>%
  unique() %>% # selecting columns with species names
  as.data.frame(.) -> splist_orchid_final # 19,297 species
# number of genera
occ_final %>%
  dplyr::select(wcvp_genus) %>%
  unique() %>% # selecting columns with species names
  as.data.frame(.) -> genuslist_orchid_final # 690

# We are ready! Save our final database with revised coordinates, no duplicates and updated taxonomy.
write.csv(occ_final, file = "processed-data/dataselection_taxonomic-standard/occ_final.csv")

#####################################################################
#####################################################################
############# Cleaning coordinates to species' native range ###########
##### At genus level ########
rm(list=ls())
orchids <- fread("processed-data/dataselection_taxonomic-standard/occ_final.csv")
# 740,837 records

# Set up parallel processing - use all available cores minus one for stability (8 cores in my computer)
plan(multisession, workers = availableCores() - 1)

# Create output directory structure
output_dir <- "processed-data/dataselection_taxonomic-standard/native.range"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Define main processing function
process_genus <- function(genus_name, data) {
  # Filter occurrences for the current genus
  genus_data <- data %>% 
    filter(wcvp_genus == genus_name) %>% 
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
             crs = st_crs(4326), remove = FALSE)
  # Skip if no records found for this genus
  if(nrow(genus_data) == 0) {
    message(sprintf("No records found for genus %s", genus_name))
    return(NULL)
  }
  # Attempt to get native range from WCVP with error handling
  genus_range <- tryCatch({
    wcvp_distribution(genus_name, taxon_rank = "genus",
                      introduced = FALSE, extinct = FALSE, 
                      location_doubtful = FALSE)
  }, error = function(e) {
    message(sprintf("Error processing genus %s: %s", genus_name, e$message))
    return(NULL)
  })
  # Discard completely if genus not in WCVP or error occurred
  if(is.null(genus_range) || nrow(genus_range) == 0) {
    message(sprintf("Genus %s not found in WCVP - discarding %d records", 
                    genus_name, nrow(genus_data)))
    return(NULL)
  }
  # Create 1km buffer around native range (approx. 0.009 degrees)
  buffered_dist <- genus_range %>%
    st_union() %>%
    st_buffer(0.009)
  # Identify records within native range (with buffer)
  genus_data$native_buffer <- st_intersects(genus_data, buffered_dist, sparse = FALSE)[,1]
  # Filter to keep only native occurrences
  filtered_data <- genus_data %>% 
    filter(native_buffer)
  # Only save if we have valid records
  if(nrow(filtered_data) > 0) {
    output_file <- file.path(output_dir, paste0(genus_name, ".csv"))
    write_csv(filtered_data, file = output_file)
    return(filtered_data)
  } else {
    message(sprintf("Genus %s: no records within native range", genus_name))
    return(NULL)
  }}

# Get unique genus list (excluding NA values)
unique_genera <- unique(orchids$wcvp_genus)
unique_genera <- unique_genera[!is.na(unique_genera)]
message(sprintf("Processing %d orchid genera...", length(unique_genera)))

# Process in parallel with progress tracking
results <- future_map(
  unique_genera,
  ~{
    message(sprintf("Processing genus: %s", .x))
    safely(process_genus)(.x, orchids) },
  .progress = TRUE,
  .options = furrr_options(seed = TRUE))

# Generate comprehensive processing report
processed_log <- tibble(
  genus = unique_genera,
  success = map_lgl(results, ~is.null(.x$error)),
  error_msg = map_chr(results, ~{
    if(is.null(.x$error)) "" else .x$error$message
  }),
  original_records = map_int(results, ~{
    if(is.null(.x$result)) {
      0L
    } else {
      genus_name <- .x$result$wcvp_genus[1]
      nrow(orchids[orchids$wcvp_genus == genus_name, ])
    }
  }),
  filtered_records = map_int(results, ~{
    if(is.null(.x$result)) 0L else nrow(.x$result)
  }),
  status = case_when(
    !success ~ "Error in processing",
    filtered_records == 0 ~ "No valid records",
    TRUE ~ "Success"))

# Save detailed log
write_csv(processed_log, file.path(output_dir, "processing_summary.csv"))

# Calculate summary statistics
summary_stats <- processed_log %>% 
  summarise(
    total_genera = n(),
    processed_successfully = sum(status == "Success"),
    genera_with_errors = sum(status == "Error in processing"),
    genera_with_no_valid_records = sum(status == "No valid records"),
    total_original_records = sum(original_records),
    total_filtered_records = sum(filtered_records),
    percentage_retained = round(total_filtered_records/total_original_records*100, 2))

# Print summary to console
message("\nProcessing summary:")
message(sprintf("- Total genera processed: %d", summary_stats$total_genera))
message(sprintf("- Successfully processed: %d", summary_stats$processed_successfully))
message(sprintf("- Genera with errors: %d", summary_stats$genera_with_errors))
message(sprintf("- Genera with no valid records: %d", summary_stats$genera_with_no_valid_records))
message(sprintf("- Original records: %d", summary_stats$total_original_records))
message(sprintf("- Filtered records retained: %d (%.2f%%)", 
                summary_stats$total_filtered_records, summary_stats$percentage_retained))

# Clean up parallel workers
plan(sequential)

#######################join the databases#######################
################################################################
# Using data.table for efficiency with large datasets
genus_files <- list.files(
  path = "processed-data/dataselection_taxonomic-standard/native.range/",
  pattern = "\\.csv$",
  full.names = TRUE)

# Merge all files
orchids_native <- rbindlist(
  lapply(genus_files, function(f) {
    dt <- fread(f)
    dt[, source_genus := gsub("\\.csv$", "", basename(f))]
    dt
  }),
  fill = TRUE) # In case there is any minimal difference between files

################################################################################
# Manually check some of the genera to ensure correct processing
################## Aa genus ##################
##############################################
Aa_genus <- orchids_native[orchids_native$wcvp_genus == "Aa", ]
# genus, lat, long
occs_Aa <- 
  Aa_genus %>% 
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
           crs = st_crs(4326), remove = FALSE)
# Getting the native range of WCVP data
Aa_range <- wcvp_distribution("Aa", taxon_rank = "genus",
                              introduced = FALSE, extinct = FALSE, 
                              location_doubtful = FALSE)
(Aa <- wcvp_distribution_map(Aa_range, crop_map = TRUE) + 
    theme(legend.position = "none"))
# Add our records in the native distribution area
Aa + geom_sf(data = occs_Aa, fill = "#6e6ad9", col = "black", shape = 21) 

################## Epidendrum ##################
##############################################
Epidendrum_genus <- orchids_native[orchids_native$wcvp_genus == "Epidendrum", ]
# genus, lat, long
occs_Epidendrum <- 
  Epidendrum_genus %>% 
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
           crs = st_crs(4326), remove = FALSE)
# Getting the native range of WCVP data
Epidendrum_range <- wcvp_distribution("Epidendrum", taxon_rank = "genus",
                                      introduced = FALSE, extinct = FALSE, 
                                      location_doubtful = FALSE)
(Epidendrum <- wcvp_distribution_map(Epidendrum_range, crop_map = TRUE) + 
    theme(legend.position = "none"))
# Add our records in the native distribution area
Epidendrum + geom_sf(data = occs_Epidendrum, fill = "#6e6ad9", col = "black", shape = 21)

# Eliminate records that do not correspond to orchidaceae.
list <- c("Endlicheria", "Ouratea", "Scybalium",
          "Quiina", "Sauvagesia",
          "Selaginella", "Vriesea")
# Deleted records
orchids_native <- orchids[!wcvp_genus %in% list]

# Final database with native records only
write_csv(orchids_native, file = "processed-data/dataselection_taxonomic-standard/orchids_native_range.csv")
# Processing summary:
# - Total genera processed: 690
# - Successfully processed: 679
# - Genera with errors: 0
# - Genera with no valid records: 11
# - Original records: 740667
# - Filtered records retained: 732367 (98.88%)

### Taxonomic level:
# We carried out biogeographic regionalization at the species level.
# We will use the column 'wcvp_name' in our final database for analysis.
# Then, we transformed this occurrence into incidence-based data matrix
# where rows represent grid cells and columns represent species using
# four grain sizes (100 × 100, 200 × 200 and 400 × 400 km, and 800 x 800 km).