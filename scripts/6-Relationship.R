# ==============================================================================
# Relationship Between Bioregions - Species Sharing Analysis
# ==============================================================================
# Clear environment and load required packages
rm(list = ls())  # Clear workspace

# Load libraries:
library(sf)          # For spatial data handling
library(dplyr)       # For data manipulation
library(circlize)    # For circular visualisations
library(paletteer)   # For colour palettes
library(ggplot2)     # For additional plotting
library(ggpubr)      # For publication-ready graphics
library(indicspecies)# For indicator species
library(tidyr)
library(Matrix)
library(tibble)

# ==============================================================================
# 1. Load spatial and community data for 200km hexagonal grid
# 2. Define six biogeographical regions based on phylogenetic beta diversity
# 3. Calculate unique and shared orchid species between regions
# 4. Visualise species sharing patterns using chord diagrams
# 5. Use spectral colour scheme matching regional biogeographical affinities
# ==============================================================================
# 1. Load data
shape200 <- st_read("processed-data/community_matrix/pam_shape/grid_200km.gpkg")
comm200 <- readRDS("processed-data/community_matrix/pam/pam_200km.rds")
load("processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_200.RData")
rm(beta_sne_mean200, beta_sor_mean200)

# 2. Define cluster names
regions <- c("Holarctic","Indo-Malaysian","Australian",
             "Chile-Patagonian","Neotropical","Afrotropical")


# 3. Cluster analysis
hc200 <- stats::hclust(beta_sim_mean200, method = "average")
clusters200 <- cutree(hc200, k = 6)
# Assigning region names to groups
shape200$cluster <- factor(clusters200, labels = regions)
# Verify that comm200 has cells in the same order.
comm200 <- as.matrix(comm200)
# Ensure correspondence between cells and groups
group_vector <- shape200$cluster[match(rownames(comm200), shape200$idcell)]

# 4. Calculate unique and shared species
species_by_region <- lapply(unique(clusters200), function(cluster) {
  colnames(comm200)[colSums(comm200[clusters200 == cluster, , drop = FALSE]) > 0] %>%
    unique()  # Remove duplicates
  })
names(species_by_region) <- regions

# Unique species calculation
unique_species <- lapply(regions, function(region) {
  setdiff(species_by_region[[region]], 
          unlist(species_by_region[names(species_by_region) != region]))
})
names(unique_species) <- regions

# Shared species matrix
shared_matrix <- matrix(0, nrow = 6, ncol = 6, dimnames = list(regions, regions))

for (i in 1:6) {
  for (j in 1:6) {
    sp_i <- species_by_region[[i]]
    sp_j <- species_by_region[[j]]
    shared_matrix[i,j] <- length(intersect(sp_i, sp_j))
  }}

# 5. Visualisation settings
region_colours <- c(
  "Holarctic"="#ED820AFF",
  "Indo-Malaysian" = "#FCDE85FF",
  "Australian" = "#A71B4BFF",
  "Chile-Patagonian"="#BAEEAEFF",  
  "Neotropical"="#00B1B5FF",
  "Afrotropical"="#584B9FFF")  # Corrected spelling

# Generate chord diagram
png("orchid-regionalisation/results/figures/shared_sp.png", width = 2500, height = 2500, res = 400)
par(family = "sans")

circos.par(
  gap.after = c(rep(4, 5), 7),
  start.degree = 5,
  track.margin = c(0, 0))  # Remove label space
colors200<-c("#ED820AFF","#FCDE85FF","#A71B4BFF",
             "#BAEEAEFF","#00B1B5FF","#584B9FFF")
chordDiagram(shared_matrix,
             grid.col = colors200,
             transparency = 0,
             annotationTrack = NULL,  # Remove axes/labels
             directional = 0,
             link.sort = TRUE,
             link.decreasing = TRUE,
             link.lwd = 1,
             link.border = "gray40")

dev.off()

# Reset graphical parameters
circos.clear()
while (!is.null(dev.list())) {
  dev.off()
}

# ==============================================================================
#Bioregion Statistics
# ==============================================================================
# 1. Calculate species statistics ----
# Total species per region
total_species_counts <- sapply(species_by_region, length)
# Unique species counts (from your previous calculation)
unique_species_counts <- sapply(unique_species, length)

region_stats <- data.frame(
  Region = regions,
  Total_Species = total_species_counts,
  Unique_Species = unique_species_counts,
  stringsAsFactors = FALSE)

# 2. Calculate spatial statistics
# 1. Cell counting
cell_counts <- shape200 %>%
  st_drop_geometry() %>%
  count(cluster) %>%
  rename(Region = cluster, Hexagonal_Cells = n)

# 2. Calculating areas
region_areas <- shape200 %>%
  mutate(area_celda = as.numeric(sf::st_area(.)) / 1e6) %>% 
  group_by(cluster) %>%  
  summarise(Area_km2 = sum(area_celda)) %>% 
  st_drop_geometry() %>%  
  rename(Region = cluster)
            
# 3. Aggregate PD values ----
load("processed-data/community_matrix/phylogenetic_metrics/PD_site_means_200.RData")
            
# Verify and merge data (corrected version)
analysis_data <- shape200 %>%
              st_drop_geometry() %>%
              select(idcell, cluster) %>%
              mutate(idcell = as.character(idcell)) %>%
              left_join(
                PD_summary200,
                by = "idcell") %>%
              rename(Region = cluster)
            
# 4. Calculate region-level statistics ----
final_table <- region_stats %>% 
              left_join(cell_counts, by = "Region") %>%
              left_join(region_areas, by = "Region") %>%
              # Add PD metrics from analysis_data
              left_join(
                analysis_data %>%
                  group_by(Region) %>%
                  summarise(
                    Mean_Species_Richness = round(mean(Species_Richness), 1),
                    Mean_PD = round(mean(PD_mean), 2),
                    SD_PD = round(sd(PD_mean), 2),
                    PD_Range = paste(
                      round(min(PD_mean), 1),
                      "-",
                      round(max(PD_mean), 1))),
                by = "Region") %>%
              # Calculate derived metrics
              mutate(
                Species_Uniqueness = round(Unique_Species / Total_Species * 100, 1),
                PD_CV = round(SD_PD / Mean_PD * 100, 1) ) %>%
              select(
                Region,
                Hexagonal_Cells,
                Area_km2,
                Total_Species,
                Unique_Species,
                Species_Uniqueness,
                Mean_Species_Richness,
                Mean_PD,
                SD_PD) %>%
              arrange(desc(Total_Species))
            
# 5. Verification step ----
# Check that total species matches sum of species richness
species_consistency_check <- final_table %>%
              mutate(
                Calculated_Total = analysis_data %>%
                  group_by(Region) %>%
                  summarise(Sum_Richness = sum(Species_Richness)) %>%
                  pull(Sum_Richness),
                Check = Total_Species == Calculated_Total)
            
if(!all(species_consistency_check$Check)) {
warning("Discrepancy between total species counts and sum of species richness values")
print(species_consistency_check)
}
            
# 6. Format and export ----
colnames(final_table) <- c(
              "Region",
              "Hex cells",
              "Area (km²)",
              "Total species",
              "Unique species",
              "Uniqueness (%)",
              "Mean sichness",
              "Mean PD",
              "± SD")
write.csv(final_table,"results/table/bioregion_stats.csv", row.names = FALSE)

#############################################################################
#############################################################################
#Analysis of indicator species
#############################################################################
# 1. Preparación de datos -------------------------------------------------
# Usar tus grupos existentes (k=6)
rownames(comm200) <- shape200$idcell  # Ensure matrix and spatial data match
hc200 <- stats::hclust(beta_sim_mean200, method = "average")
clusters200 <- cutree(hc200, k = 6)

# 2. Pre-filtering of species:
# Why consider pre-filtering?
# Species with too few occurrences may generate spurious associations.
# - Statistical analysis requires sufficient power (usually n ≥ 5)
# - Reference literature:
# - Dufrêne & Legendre (1997): Recommend excluding rare species
# - De Cáceres et al. (2010): Suggest minimum thresholds for occurrence
# Pre-filtered for dgCMatrix
class(comm200)
dim(comm200)
occurrences <- Matrix::colSums(comm200 != 0)
# Conteo de ocurrencias
species_to_keep <- names(which(occurrences >= 5))
comm200_filtrado <- comm200[, species_to_keep]
# 3. Análisis de especies indicadoras
set.seed(123)  # Para reproducibilidad
# Versión robusta de IndVal
indval_result <- multipatt(
  x = as.matrix(comm200_filtrado),  # Convert the dense matrix
  cluster = clusters200,
  func = "IndVal.g",  # IIndVal for unequal groups
  duleg = TRUE,       # Consider only species with maximum in one group only
  control = how(nperm = 999))

# 3. Process ALL identified species
all_species_df <- indval_result$sign %>%
  as.data.frame() %>%
  rownames_to_column("Species") %>%
  mutate(
    Group = colnames(indval_result$comb)[index],
    # Calculate the difference between the maximum IndVal value and the second maximum
    IndVal_diff = apply(indval_result$str, 1, function(x) {
      sorted <- sort(x, decreasing = TRUE)
      if(length(sorted) > 1) sorted[1] - sorted[2] else NA
    })) %>%
  select(Species, Group, Stat = stat, P_value = p.value, IndVal_diff)

# 4. Calculate group distribution ----------------------------
# Attendance/absence matrix by group
group_presence <- t(apply(as.matrix(comm200_filtrado), 2, function(x) {
  tapply(x > 0, clusters200, sum)
}))

# Convert to percentage
group_percentage <- t(apply(group_presence, 1, function(x) {
  round(x / sum(x) * 100, 1)
}))

# Rename columns
colnames(group_percentage) <- paste0("Perc_Group_", colnames(group_percentage))
colnames(group_presence) <- paste0("Count_Group_", colnames(group_presence))

# 5. Combining all information
REVIEW:
  # - Consistency between assigned Group and highest Perc_Group_X
  # - Species with Max_Percentage >80% as best candidates
  # - Species with Max_Percentage <50% as potential false positives
  left_join(
    as.data.frame(group_presence) %>% 
      tibble::rownames_to_column("Species"),
    by = "Species") %>%
  left_join(
    as.data.frame(group_percentage) %>% 
      tibble::rownames_to_column("Species"),
    by = "Species") %>%
  mutate(
    Total_Occurrences = rowSums(group_presence[all_species_df$Species, ]),
    Max_Percentage = apply(group_percentage[all_species_df$Species, ], 1, max)) %>%
  arrange(Group, -Stat)
write.csv(full_results, "results/table/all_indicator_species_full.csv", row.names = FALSE)

# subsets
# Define thresholds
p_threshold <- 0.05
stat_threshold <- 0.3
concentration_threshold <- 80
indval_diff_threshold <- 0.1

# Selecting the best indicator species
best_indicators <- full_results %>%
  filter(
    P_value < p_threshold,
    Stat > stat_threshold,
    Max_Percentage > concentration_threshold,
    IndVal_diff > indval_diff_threshold) %>%
  arrange(Group, -Stat) %>%
  select(Species, Group, Stat, P_value, Max_Percentage, IndVal_diff)

# Exporting results
write.csv(best_indicators, "results/table/best_indicator_species.csv", row.names = FALSE)

# Optional: Top 5 per group
top5_per_group <- best_indicators %>%
  group_by(Group) %>%
  slice_max(order_by = Stat, n = 5) %>%
  ungroup()

write.csv(top5_per_group, "results/table/top5_indicators_per_group.csv", row.names = FALSE)# 7. Crear resumen por grupo -------------------------------------
# REVIEW:
# - Groups with mean_Percentage <60% contain poorly constrained species
# - Groups with n_species <10 may be biogeographically insignificant
group_summary <- full_results %>%
  group_by(Group) %>%
  summarise(
    n_species = n(),
    mean_Stat = mean(Stat, na.rm = TRUE),
    mean_Percentage = mean(Max_Percentage, na.rm = TRUE),
    .groups = "drop")
write.csv(group_summary, "results/table/indicator_species_group_summary.csv", row.names = FALSE)

##################################################################
# Cell Affiliation Analysis Across Bioregions
##################################################################
rm(list = ls())  # Clear workspace

library(sf)            # Spatial data handling
library(dplyr)         # Data manipulation
library(circlize)      # Circular visualisations
library(paletteer)     # Colour palettes
library(ggplot2)       # Data visualisation
library(ggpubr)        # Publication-ready graphics
require(devtools)
library(Herodotools)   # Biogeographical analyses
library(rnaturalearth) # Base map data
library(ggthemes)      # Map themes
library(scales)        # Colour scaling
library(dendextend)
library(tibble)

# ==============================================================================
# Methodology Summary:
# 1. Load spatial and community data
# 2. Perform cluster analysis and calculate cell affiliation
# 3. Merge affiliation data with spatial polygons
# 4. Prepare base map and colour scheme
# 5. Visualise affiliation patterns with gradient colours
# ==============================================================================

# 1. Load data ----------------------------------------------------------------
shape200 <- st_read("processed-data/community_matrix/pam_shape/grid_200km.gpkg")
comm200 <- readRDS("processed-data/community_matrix/pam/pam_200km.rds")
load("processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_200.RData")
rm(beta_sne_mean200, beta_sor_mean200)

# 2. Cluster analysis and affiliation -----------------------------------------
hc200 <- stats::hclust(beta_sim_mean200, method = "average")
clusters200 <- cutree(hc200, k = 6)
groups_factor <- factor(clusters200)

# Calculate cell-region affiliation
afi <- Herodotools::calc_affiliation_evoreg(beta_sim_mean200, groups_factor)
afi <- as.data.frame(afi) %>%
  mutate(idcell = rownames(.))

# 3. Merge affiliation data with spatial polygons -----------------------------
shape200_affiliation <- shape200 %>%
  mutate(idcell = as.character(idcell)) %>%
  left_join(afi, by = "idcell")

# 4. Prepare base map ---------------------------------------------------------
continents <- ne_countries(
  continent = c("Africa", "Asia", "North America", "Europe", "Oceania", "South America"),
  returnclass = "sf",
  scale = "medium")

# Set Behrmann projection
behrmann <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
map <- st_transform(continents, behrmann)
shape200_affiliation <- st_transform(shape200_affiliation, behrmann)

# 5. Define colour scheme -----------------------------------------------------
# to make sure that the colours are assigned to exactly the right group,
#we generate the same dendrogram as in step 5.
hc200 <- stats::hclust(beta_sim_mean200, method = "average")
clusters200 <- cutree(hc200, k = 6)
colors200 <- as.character(paletteer_c("grDevices::Spectral", 6))

# Colouring dendrogram
dend200 <- as.dendrogram(hc200)
dend200 <- color_branches(dend200, k = 6, col = colors200)

# Extract order, labels and true colours
dend_order <- order.dendrogram(dend200)
leaves <- labels(dend200)
leaf_colors <- get_leaves_branches_col(dend200)

# Create table with idcell, group and colour assigned according to dendrogram
group_df <- tibble(
  idcell = leaves,
  group = clusters200[leaves],
  color = leaf_colors
) %>% distinct(group, .keep_all = TRUE)

# Asignar colores a shape200_affiliation por left_join
shape200_affiliation <- shape200_affiliation %>%
  mutate(group = clusters200[as.character(idcell)]) %>%
  left_join(group_df %>% select(group, color), by = "group")

# 8. Identify low-affiliation cells -------------------------------------------
q1 <- quantile(shape200_affiliation$afilliation, probs = 0.25, na.rm = TRUE)
shape200_affiliation$low_affiliation <- shape200_affiliation$afilliation <= q1

# 9. Create and save affiliation map ------------------------------------------
affiliation_map <- ggplot() +
  geom_sf(data = map, fill = "grey60", colour = "grey60", linewidth = 0.2) +
  geom_sf(data = shape200_affiliation, aes(fill = color), colour = NA, size = 0.1) +
  geom_sf(data = filter(shape200_affiliation, low_affiliation), fill = "black", colour = NA, size = 0.1) +
  scale_fill_identity() +
  theme_map() +
  theme(
    legend.position = "none",
    text = element_text(size = 16),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()) +
  coord_sf(expand = FALSE)

affiliation_map

ggsave("results/figures/affiliation_200km.png",
       affiliation_map,dpi = 400,width = 10,height = 6,bg = "white")
