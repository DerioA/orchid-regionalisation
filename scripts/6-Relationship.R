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

# ==============================================================================
# 1. Load spatial and community data for 200km hexagonal grid
# 2. Define six biogeographical regions based on phylogenetic beta diversity
# 3. Calculate unique and shared orchid species between regions
# 4. Visualise species sharing patterns using chord diagrams
# 5. Use spectral colour scheme matching regional biogeographical affinities
# ==============================================================================
# 1. Load data
shape200 <- st_read("processed-data/community_matrix/pam_shape/shape_reduce200.shp")
load("processed-data/community_matrix/pam/pam200_reduce.RData")
load("processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_200.RData")
rm(beta_sne_mean200, beta_sor_mean200)

# 2. Define cluster names
regions <- c("Chile-Patagonian","Australian","Neotropical", 
             "Afrotropical","Indo-Malaysian","Holartic")

# 3. Cluster analysis
hc200 <- stats::hclust(beta_sim_mean200, method = "average")
clusters200 <- cutree(hc200, k = 6)
shape200$cluster <- factor(clusters200, labels = regions)
comm200 <- as.matrix(comm200)

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
  "Chile-Patagonian" = "#ED820AFF",
  "Australian" = "#A71B4BFF",
  "Neotropical" = "#FCDE85FF",          
  "Afrotropical" = "#BAEEAEFF",   
  "Indo-Malaysian" = "#584B9FFF",   
  "Holarctic" = "#00B1B5FF")  # Corrected spelling

# Generate chord diagram
png("results/figures/shared_sp.png", width = 2500, height = 2500, res = 400)
par(family = "sans")

circos.par(
  gap.after = c(rep(4, 5), 7),
  start.degree = 5,
  track.margin = c(0, 0))  # Remove label space
colors200<-c("#ED820AFF","#A71B4BFF","#FCDE85FF","#BAEEAEFF",
             "#584B9FFF","#00B1B5FF" )
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
# Total species per region (from your original community matrix)
total_species_counts <- sapply(species_by_region, length)

# Unique species counts (from your previous calculation)
unique_species_counts <- sapply(unique_species, length)

region_stats <- data.frame(
  Region = regions,
  Total_Species = total_species_counts,
  Unique_Species = unique_species_counts,
  stringsAsFactors = FALSE
)

# 2. Calculate spatial statistics ----
# (This part remains correct as is)
cell_counts <- shape200 %>%
  st_drop_geometry() %>%
  count(cluster) %>%
  rename(Region = cluster, Hexagonal_Cells = n)

region_areas <- shape200 %>%
  group_by(cluster) %>%  # Group by bioregion
  summarise(
    Area_km2 = sum(as.numeric(st_area(geometry)))/1e6  # Calculate area in km²
    ) %>% 
      st_drop_geometry() %>%  # Remove geometry column
      rename(Region = cluster)  # Rename cluster to Region
            
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
# Make sure the sites are in the same order.
group_vector <- shape200$cluster[match(rownames(comm200), shape200$idcell)]
# Run indicator species analysis
indval_result <- multipatt(comm200, cluster = group_vector, func = "r.g", control = how(nperm = 999))
# Convertir a data.frame y agregar nombre de especie
indval_table <- as.data.frame(indval_result$sign)
indval_table$Species <- rownames(indval_table)
# Identificar columnas que son grupos
group_cols <- grep("^s\\.", names(indval_table))

# Filtrar y obtener top 5 especies indicadoras por grupo (con p-value)
top_species_per_region <- indval_table %>%
  filter(p.value <= 0.05) %>%
  mutate(Group = apply(.[, group_cols], 1, function(x) {
    groups <- which(x == 1)
    if (length(groups) == 1) return(groups) else return(NA)
  })) %>%
  filter(!is.na(Group)) %>%
  group_by(Group) %>%
  slice_max(order_by = stat, n = 5, with_ties = FALSE) %>%
  mutate(Species_Info = paste0(
    Species,
    " (stat=", signif(stat, 3),
    ", p=", signif(p.value, 3), ")"
    )) %>%
  summarise(Indicator_Species = paste(Species_Info, collapse = ", "))
# Asignar nombres a los grupos
region_names <- c(
  "Chile-Patagonian","Australian","Neotropical", 
  "Afrotropical","Indo-Malaysian","Holartic")

# Agregar nombre de la región
top_species_per_region$Region <- region_names[as.numeric(top_species_per_region$Group)]

# Reordenar columnas
top_species_per_region <- top_species_per_region %>%
  select(Region, Top_Indicator_Species)
write.csv(top_species_per_region,"results/table/bio_species.csv", row.names = FALSE)
###################################################################
# Cell Affiliation Analysis Across Bioregions
##################################################################
rm(list = ls())  # Clear workspace

library(sf)            # Spatial data handling
library(dplyr)         # Data manipulation
library(circlize)      # Circular visualisations
library(paletteer)     # Colour palettes
library(ggplot2)       # Data visualisation
library(ggpubr)        # Publication-ready graphics
library(Herodotools)   # Biogeographical analyses
library(rnaturalearth) # Base map data
library(ggthemes)      # Map themes
library(scales)        # Colour scaling

# ==============================================================================
# Methodology Summary:
# 1. Load spatial and community data
# 2. Perform cluster analysis and calculate cell affiliation
# 3. Merge affiliation data with spatial polygons
# 4. Prepare base map and colour scheme
# 5. Visualise affiliation patterns with gradient colours
# ==============================================================================

# 1. Load data ----------------------------------------------------------------
shape200 <- st_read("processed-data/community_matrix/pam_shape/shape_reduce200.shp")
load("processed-data/community_matrix/pam/pam200_reduce.RData")
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
  scale = "medium"
)

# Set Behrmann projection
behrmann <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
map <- st_transform(continents, behrmann)
shape200_affiliation <- st_transform(shape200_affiliation, behrmann)

# 5. Define colour scheme -----------------------------------------------------
group_colours <- as.character(paletteer_c("grDevices::Spectral", 6))

# 6. Create gradient colour function ------------------------------------------
get_gradient_colour <- function(aff, group, colours) {
  if (is.na(aff) || is.na(group)) return(NA)
  
  group_num <- as.numeric(group)
  if (group_num < 1 || group_num > length(colours)) return(NA)
  
  base_rgb <- col2rgb(colours[group_num]) / 255
  final_rgb <- aff * base_rgb + (1 - aff) * 1  # Blend with white
  rgb(final_rgb[1], final_rgb[2], final_rgb[3], alpha = aff)
}

# 7. Assign gradient colours to cells -----------------------------------------
shape200_affiliation$fill_colour <- mapply(
  get_gradient_colour,
  aff = shape200_affiliation$afilliation,
  group = shape200_affiliation$group,
  MoreArgs = list(colours = group_colours)
)

# 8. Identify low-affiliation cells -------------------------------------------
q1 <- quantile(shape200_affiliation$afilliation, probs = 0.25, na.rm = TRUE)
shape200_affiliation$low_affiliation <- shape200_affiliation$afilliation <= q1

# 9. Create and save affiliation map ------------------------------------------
affiliation_map <- ggplot() +
  geom_sf(data = map, fill = "grey60", colour = "grey60", linewidth = 0.2) +
  #geom_sf(data = shape200_affiliation, aes(fill = fill_colour), colour = NA, size = 0.1) +
  scale_fill_identity() +
  geom_sf(
    data = filter(shape200_affiliation, low_affiliation),
    fill = NA, colour = "black", size = 0.01
  ) +
  theme_bw()+
  theme(
    panel.background = element_rect(fill = "white",colour = NA),
    plot.background = element_rect(fill = "white",colour = NA),
    text = element_text(size = 20, family = "sans"),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )
affiliation_map

affiliation_map <- ggplot() +
  geom_sf(data = map, fill = "grey60", colour = "grey60", linewidth = 0.2) +
  geom_sf(
    data = shape200_affiliation,
    fill = "black",
    colour = NA,
    size = 0.01) +
  geom_sf(
    data = shape200_affiliation,
    aes(fill = fill_colour,
        alpha = ifelse(low_affiliation, 0.1, 1)),  # Transparencia condicional
    colour = NA,
    size = 0.1) +
  scale_fill_identity() +
  scale_alpha_identity() +  # Añade escala para alpha
  theme_bw() +
  theme(
    legend.position = "none",
    text = element_text(size = 16),
    axis.title.x = element_blank(),        
    axis.text.x = element_blank(),     
    axis.ticks.x = element_blank(),    
    axis.line.x = element_blank(),    
    axis.title.y = element_blank(),    
    axis.text.y = element_blank(),    
    axis.ticks.y = element_blank(),    
    axis.line.y = element_blank()) +
  coord_sf(expand = FALSE)

affiliation_map

ggsave("results/figures/affiliation_200km.png",affiliation_map,dpi = 400,width = 10,height = 6)
