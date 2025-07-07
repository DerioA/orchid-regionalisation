####### Explore the data #######
### 1.- MAPPING BETA PHYLOGENETIC DIVERSITY
### 2.- SAMPLE COMPLETENESS ANALYSIS

#------------------------------------------
### 1.- MAPPING BETA PHYLOGENETIC DIVERSITY
#------------------------------------------
## Here we use the 'beta_sim_mean' file generated in the 'community matrix' section
rm(list = ls()) # Clear environmental

# Load all required packages at the beginning
library(sf)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(RColorBrewer)
library(extrafont)
library(ggthemes)
library(data.table)
library(tidyr)
library(KnowBR)
library(gridExtra)
library(grid)

### PHYLOGENETIC BETA DIVERSITY MAPPING ----

# 1. Load the data
shape200 <- st_read(dsn = "processed-data/community_matrix/pam_shape/grid_200km.gpkg", quiet = TRUE)
load("processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_200.RData") # Phylogenetic beta diversity components
comm200 <- readRDS("processed-data/community_matrix/pam/pam_200km.rds") # Species presence-absence matrix

# 2. Data processing
# Converting distance matrix and calculating averages per site
class(beta_sim_mean200)
beta_sim_matrix_200 <- as.matrix(beta_sim_mean200)

# Verify the converted matrix
dim(beta_sim_matrix_200)  # Show square dimensions
isSymmetric(beta_sim_matrix_200)  # Should be TRUE
diag(beta_sim_matrix_200)  # Will display the diagonal (normally 0)

# 3. Calculate row averages
beta_sim_site200 <- rowMeans(beta_sim_matrix_200, na.rm = TRUE)

# 4. Create data.frame with the results
beta_df200 <- data.frame(
  site_id = rownames(beta_sim_matrix_200),
  beta_sim_mean = beta_sim_site200,
  stringsAsFactors = FALSE)

# 4. Joining with shapefile
beta_df200 <- rename(beta_df200, idcell = site_id)
beta_df200$idcell <- as.character(beta_df200$idcell)
shape200$idcell <- as.character(shape200$idcell)
shape200 <- left_join(shape200, beta_df200, by = "idcell")

# 5. Behrmann projection
behrmann <- "+proj=cea +lon_0=0 +lat_ts=30 +datum=WGS84 +units=m +no_defs"
shape200 <- st_transform(shape200, crs = behrmann)

# 6. Base map
continents <- ne_countries(
  continent = c("Africa", "Asia", "North America", "Europe", "Oceania", "South America"),
  returnclass = "sf",
  scale = "medium") %>% st_transform(behrmann)

# 7. Visualisation ----
color_palette <- rev(brewer.pal(11, "Spectral"))

pb <- ggplot() +
  geom_sf(data = continents, fill = "grey60", colour = "grey60") +
  geom_sf(data = shape200, aes(fill = beta_sim_mean), color = NA, lwd = 0) + 
  scale_fill_gradientn(
    colours = color_palette,
    name = expression("Mean " * beta * "sim"),
    na.value = "transparent",
    breaks = c(0.40,0.50,0.60,0.70,0.80),
    labels = function(x) ifelse(x %in% c(0, 1), as.character(x), sprintf("%.2f", x)),
    limits = c(0.40, 0.81), 
    guide = guide_colorbar(
      title.position = "left",
      title.hjust = 0.5,
      barwidth = unit(7, "cm"),
      barheight = unit(0.7, "cm"),
      frame.colour = "black",
      ticks.colour = "black",
      direction = "horizontal")) +
  coord_sf(datum = NA) +
  theme_bw() +
  theme(
    text = element_text(family = "sans", size = 14),
    legend.title = element_text(family = "sans", size = 14),
    legend.text = element_text(family = "sans", size = 12),
    legend.position = "bottom",
    legend.justification = "center",
    legend.box = "horizontal",
    legend.margin = margin(t = 0, r = 0, b = 5, l = 0),
    plot.title = element_text(family = "sans"),
    plot.subtitle = element_text(family = "sans"))

pb

# Export
ggsave("results/Figures/phylogenetic.beta.diversity_200km.png", pb,
       width = 7.7, height = 5, units = "in", dpi = 400)

#------------------------------------------
### 2.- SAMPLE COMPLETENESS ANALYSIS
#------------------------------------------
rm(list = ls())

# 100KM analysis ----
shape100 <- st_read(dsn = "processed-data/community_matrix/pam_shape/grid_100km.gpkg")
shape100 <- st_transform(shape100, crs = 4326)
shape100 <- as(shape100, "Spatial")

orchids <- fread("processed-data/dataselection_taxonomic-standard/orchids_native_range.csv",
                 na.strings = "") # 732,359 records

# Prepare species list
species.list <- orchids %>% 
  dplyr::select(wcvp_name, decimalLongitude, decimalLatitude) %>%
  as.data.frame() %>%
  mutate(Counts = 1)

data(adworld)
# Calculate completeness
completeness100km <- KnowBPolygon(
  species.list, format = "A", shape = shape100,
  shapenames = "idcell", admAreas = FALSE,
  Area = "World", curve = "Rational", estimator = 1,
  cutoff = 1, cutoffCompleteness = 0,
  cutoffSlope = 1, extent = TRUE, inc = 0.005, exclude = NULL,
  colexc = NULL,
  save = "RData", 
  file1 = "processed-data/data_exploration/Completeness/Species per site100",
  file2 = "processed-data/data_exploration/Completeness/Estimators100", 
  file3 = "processed-data/data_exploration/Completeness/Standard error of the estimators100",
  Maps = FALSE, jpg = FALSE)

################################################################################
# 200KM analysis ----
shape200 <- st_read(dsn = "processed-data/community_matrix/pam_shape/grid_200km.gpkg")
shape200 <- st_transform(shape200, crs = 4326)
shape200 <- as(shape200, "Spatial")

completeness200km <- KnowBPolygon(
  species.list, format = "A", shape = shape200,
  shapenames = "idcell", admAreas = FALSE,
  Area = "World", curve = "Rational", estimator = 1,
  cutoff = 1, cutoffCompleteness = 0,
  cutoffSlope = 1, extent = TRUE, inc = 0.005, exclude = NULL,
  colexc = NULL,
  save = "RData", 
  file1 = "processed-data/data_exploration/Completeness/Species per site200",
  file2 = "processed-data/data_exploration/Completeness/Estimators200", 
  file3 = "processed-data/data_exploration/Completeness/Standard error of the estimators200",
  Maps = FALSE, jpg = FALSE)

################################################################################
# 400KM analysis ----
shape400 <- st_read(dsn = "processed-data/community_matrix/pam_shape/grid_400km.gpkg")
shape400 <- st_transform(shape400, crs = 4326)
shape400 <- as(shape400, "Spatial")

completeness400km <- KnowBPolygon(
  species.list, format = "A", shape = shape400,
  shapenames = "idcell", admAreas = FALSE,
  Area = "World", curve = "Rational", estimator = 1,
  cutoff = 1, cutoffCompleteness = 0,
  cutoffSlope = 1, extent = TRUE, inc = 0.005, exclude = NULL,
  colexc = NULL,
  save = "RData", 
  file1 = "processed-data/data_exploration/Completeness/species per site400",
  file2 = "processed-data/data_exploration/Completeness/Estimators400", 
  file3 = "processed-data/data_exploration/Completeness/Standard error of the estimators400",
  Maps = FALSE, jpg = FALSE)

################################################################################
# 800KM analysis ----
shape800 <- st_read(dsn = "orchid-regionalisation/processed-data/community_matrix/pam_shape/grid_800km.gpkg")
shape800 <- st_transform(shape800, crs = 4326)
shape800 <- as(shape800, "Spatial")

completeness800km <- KnowBPolygon(
  species.list, format = "A", shape = shape800,
  shapenames = "idcell", admAreas = FALSE,
  Area = "World", curve = "Rational", estimator = 1,
  cutoff = 1, cutoffCompleteness = 0,
  cutoffSlope = 1, extent = TRUE, inc = 0.005, exclude = NULL,
  colexc = NULL,
  save = "RData", 
  file1 = "processed-data/data_exploration/Completeness/species per site800",
  file2 = "processed-data/data_exploration/Completeness/Estimators800", 
  file3 = "processed-data/data_exploration/Completeness/Standard error of the estimators800",
  Maps = FALSE, jpg = FALSE)

#####################################################################
### Create completeness maps
#####################################################################

# Behrmann projection setup
behrmann <- "+proj=cea +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
continents <- ne_countries(continent = c("Africa", "Asia", "North America",
                                         "Europe", "Oceania", "South America"),
                           returnclass = "sf", scale = "medium")

behrmann <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

map <- st_transform(continents, behrmann)
# Color palette for quality categories
colors <- c("#3288BD", "#FEE08B", "#D7191C")

#####################################################################
### 100KM completeness map ----
#####################################################################
load("processed-data/data_exploration/Completeness/Estimators100.RData")
shape100 <- st_read(dsn = "processed-data/community_matrix/pam_shape/grid_100km.gpkg")
colnames(estimators)[1] <- "idcell"

completeness100 <- shape100 %>% 
  left_join(estimators, by = "idcell")

without_estimation100 <- completeness100 %>% filter(is.na(Completeness)) # 405 cells
completeness100 <- completeness100 %>% filter(!is.na(Completeness))

# Classify survey quality
completeness100 <- completeness100 %>%
  mutate(Quality = case_when(
    Slope < 0.02 & Completeness > 90 & Ratio > 15 ~ "High",
    Slope > 0.3 & Completeness < 50 & Ratio < 3 ~ "Poor",
    TRUE ~ "Fair"
  )) %>%
  mutate(Quality = factor(Quality, levels = c("High", "Fair", "Poor")))

table(completeness100$Quality)

completenes_100km <- ggplot() +
  geom_sf(data = map, fill = "grey60", colour = "grey60") +
  geom_sf(data = without_estimation100, fill = "gray10", color = "gray10", linewidth = 0.05) +
  geom_sf(data = completeness100, aes(fill = Quality), linewidth = 0.05, color = NA) +
  scale_fill_manual(values = colors, name = "") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = c(0.08, 0.7),
    legend.background = element_rect(fill = alpha("white", 0), color = "black", linewidth = NA),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.5, "cm") ) +
    coord_sf(expand = FALSE)
completenes_100km

#####################################################################
### 200KM completeness map ----
#####################################################################

load("processed-data/data_exploration/Completeness/Estimators200.RData")
shape200 <- st_read(dsn = "processed-data/community_matrix/pam_shape/grid_200km.gpkg")
colnames(estimators)[1] <- "idcell"

completeness200 <- shape200 %>% 
  left_join(estimators, by = "idcell")

without_estimation200 <- completeness200 %>% filter(is.na(Completeness))
completeness200 <- completeness200 %>% filter(!is.na(Completeness))

completeness200 <- completeness200 %>%
  mutate(Quality = case_when(
    Slope < 0.02 & Completeness > 90 & Ratio > 15 ~ "High",
    Slope > 0.3 & Completeness < 50 & Ratio < 3 ~ "Poor",
    TRUE ~ "Fair"
  )) %>%
  mutate(Quality = factor(Quality, levels = c("High", "Fair", "Poor")))

table(completeness200$Quality)

completenes_200km <- ggplot() +
  geom_sf(data = map, fill = "grey60", colour = "grey60") +
  geom_sf(data = without_estimation200, fill = "gray10", color = "gray10", linewidth = 0.05) +
  geom_sf(data = completeness200, aes(fill = Quality), linewidth = 0.05, color = NA) +
  scale_fill_manual(values = colors, name = "", guide = "none") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()) +
  coord_sf(expand = FALSE)
completenes_200km

#####################################################################
### 400KM completeness map ----
#####################################################################

load("processed-data/data_exploration/Completeness/Estimators400.RData")
shape400 <- st_read(dsn = "processed-data/community_matrix/pam_shape/grid_400km.gpkg")
colnames(estimators)[1] <- "idcell"

completeness400 <- shape400 %>% 
  left_join(estimators, by = "idcell")

without_estimation400 <- completeness400 %>% filter(is.na(Completeness))
completeness400 <- completeness400 %>% filter(!is.na(Completeness))

completeness400 <- completeness400 %>%
  mutate(Quality = case_when(
    Slope < 0.02 & Completeness > 90 & Ratio > 15 ~ "High",
    Slope > 0.3 & Completeness < 50 & Ratio < 3 ~ "Poor",
    TRUE ~ "Fair"
  )) %>%
  mutate(Quality = factor(Quality, levels = c("High", "Fair", "Poor")))

table(completeness400$Quality)

completenes_400km <- ggplot() +
  geom_sf(data = map, fill = "grey60", colour = "grey60") +
  geom_sf(data = without_estimation400, fill = "gray10", color = "gray10", linewidth = 0.05) +
  geom_sf(data = completeness400, aes(fill = Quality), linewidth = 0.05, color = NA) +
  scale_fill_manual(values = colors, name = "", guide = "none") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()) +
  coord_sf(expand = FALSE)
completenes_400km

#####################################################################
### 800KM completeness map ----
#####################################################################

load("processed-data/data_exploration/Completeness/Estimators800.RData")
shape800 <- st_read(dsn = "processed-data/community_matrix/pam_shape/grid_800km.gpkg")
colnames(estimators)[1] <- "idcell"

completeness800 <- shape800 %>% 
  left_join(estimators, by = "idcell")

without_estimation800 <- completeness800 %>% filter(is.na(Completeness))
completeness800 <- completeness800 %>% filter(!is.na(Completeness))

completeness800 <- completeness800 %>%
  mutate(Quality = case_when(
    Slope < 0.02 & Completeness > 90 & Ratio > 15 ~ "High",
    Slope > 0.3 & Completeness < 50 & Ratio < 3 ~ "Poor",
    TRUE ~ "Fair"
  )) %>%
  mutate(Quality = factor(Quality, levels = c("High", "Fair", "Poor")))
         
table(completeness800$Quality)
         
completenes_800km <- ggplot() +
  geom_sf(data = map, fill = "grey60", colour = "grey60") +
    geom_sf(data = without_estimation800, fill = "gray10", color = "gray10", linewidth = 0.05) +
    geom_sf(data = completeness800, aes(fill = Quality), linewidth = 0.05, color = NA) +
    scale_fill_manual(values = colors, name = "", guide = "none") +
    theme_bw() +
    theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()) +
    coord_sf(expand = FALSE)
completenes_800km      

#####################################################################
### Combine all completeness maps ----
#####################################################################
library(patchwork)
combined_maps <- (completenes_100km | completenes_200km) /
  (completenes_400km | completenes_800km)
combined_maps
ggsave("results/Figures/completeness.png", plot = combined_maps, width = 10, height = 8,units = "in",dpi = 400)
