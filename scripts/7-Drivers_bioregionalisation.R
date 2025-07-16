###############################################################################
#       Climate and Topographic Variables Processing at 200 km Resolution     #
#           Includes CHELSA layers + Quaternary climate velocity              #
###############################################################################
rm(list = ls())

# Load libraries -------------------------------------------------------------
library(terra)       # Raster data processing
library(sf)          # Vector spatial data
library(dplyr)       # Data manipulation
library(ggplot2)     # Plotting
library(ggthemes)    # Additional ggplot themes
library(tidyr)       # Data tidying
library(tibble)      # Enhanced data frames
library(readr)       # Data input/output
library(reshape2)    # Data reshaping
library(viridis)     # Colour scales
library(purrr)       # Functional programming (used for map2)
library(patchwork)   # Combine ggplots

# 1. Load CHELSA raster files ------------------------------------------------
chelsa_path <- "raw-data/CHELSA/"
temp_mean <- rast(paste0(chelsa_path, "CHELSA_bio1_1981-2010_V.2.1.tif"))
prec_annual <- rast(paste0(chelsa_path, "CHELSA_bio12_1981-2010_V.2.1.tif"))
temp_seasonality <- rast(paste0(chelsa_path, "CHELSA_bio4_1981-2010_V.2.1.tif"))
prec_seasonality <- rast(paste0(chelsa_path, "CHELSA_bio15_1981-2010_V.2.1.tif"))
elevation <- rast(paste0(chelsa_path, "elevation_1KMmn_GMTEDmn.tif"))
LGM_temp <- rast(paste0(chelsa_path, "CHELSA_PMIP_CCSM4_BIO_01.tif"))
LGM_prec <- rast(paste0(chelsa_path, "CHELSA_PMIP_CCSM4_BIO_12.tif"))

# 2. Define projection (Behrmann equal-area) ---------------------------------
behrmann_proj4 <- "+proj=cea +lon_0=0 +lat_ts=30 +datum=WGS84 +units=m +no_defs"

# 3. Function to project and resample to 200 km resolution in Behrmann projection
project_resample <- function(r, proj, res_km = 200) {
  projected <- project(r, proj, method = "bilinear")
  ref <- rast(ext(projected), resolution = res_km * 1000, crs = crs(projected))
  resampled <- resample(projected, ref, method = "bilinear")
  return(resampled)
}

# 4. Project and resample all layers ------------------------------------------
temp_mean_rescaled <- project_resample(temp_mean, behrmann_proj4)
prec_annual_rescaled <- project_resample(prec_annual, behrmann_proj4)
temp_seasonality_rescaled <- project_resample(temp_seasonality, behrmann_proj4)
prec_seasonality_rescaled <- project_resample(prec_seasonality, behrmann_proj4)
elevation_rescaled <- project_resample(elevation, behrmann_proj4)
LGM_temp_rescaled <- project_resample(LGM_temp, behrmann_proj4)
LGM_prec_rescaled <- project_resample(LGM_prec, behrmann_proj4)

# 5. Check CRS match ----------------------------------------------------------
all_same_crs <- all(
  crs(temp_mean_rescaled) == crs(prec_annual_rescaled),
  crs(temp_mean_rescaled) == crs(temp_seasonality_rescaled),
  crs(temp_mean_rescaled) == crs(prec_seasonality_rescaled),
  crs(temp_mean_rescaled) == crs(elevation_rescaled),
  crs(temp_mean_rescaled) == crs(LGM_temp_rescaled),
  crs(temp_mean_rescaled) == crs(LGM_prec_rescaled))
print(all_same_crs)  # TRUE if all match

# 6. Load polygon grid and project to Behrmann -------------------------------
shape200_vect <- vect("processed-data/community_matrix/pam_shape/grid_200km.gpkg")
shape200_vect <- terra::project(shape200_vect, behrmann_proj4)

# 7. Extract mean values per polygon -------------------------------------------
cv_temp_vals <- terra::extract(temp_mean_rescaled, shape200_vect, fun = mean, na.rm = TRUE)
cv_prec_vals <- terra::extract(prec_annual_rescaled, shape200_vect, fun = mean, na.rm = TRUE)
cv_temp_seas_vals <- terra::extract(temp_seasonality_rescaled, shape200_vect, fun = mean, na.rm = TRUE)
cv_prec_seas_vals <- terra::extract(prec_seasonality_rescaled, shape200_vect, fun = mean, na.rm = TRUE)
hetero_vals <- terra::extract(elevation_rescaled, shape200_vect, fun = mean, na.rm = TRUE)
LGM_temp_vals <- terra::extract(LGM_temp_rescaled, shape200_vect, fun = mean, na.rm = TRUE)
LGM_prec_vals <- terra::extract(LGM_prec_rescaled, shape200_vect, fun = mean, na.rm = TRUE)

shape200_vect$temp_mean <- cv_temp_vals[, 2]
shape200_vect$prec_annual <- cv_prec_vals[, 2]
shape200_vect$temp_seasonality <- cv_temp_seas_vals[, 2]
shape200_vect$precipitation_seasonality <- cv_prec_seas_vals[, 2]
shape200_vect$heterogeneity_topo <- hetero_vals[, 2]
shape200_vect$LGM_temp <- LGM_temp_vals[, 2]
shape200_vect$LGM_prec <- LGM_prec_vals[, 2]

# 9. Convert to sf ----------------------------------------
shape200_sf <- st_as_sf(shape200_vect)

# 10. Prepare data frame for plotting ------------------------------------------
climate_data <- shape200_sf %>%
  rename(
    "Mean temperature" = temp_mean,
    "Annual precipitation" = prec_annual,
    "Temperature seasonality" = temp_seasonality,
    "Precipitation seasonality" = precipitation_seasonality,
    "Elevation" = heterogeneity_topo,
    "LGM mean temperature" = LGM_temp,
    "LGM annual precipitation" = LGM_prec)

# 11. Function for thematic maps -----------------------------------------------
plot_map <- function(sf_data, variable) { 
  ggplot(sf_data) +
    geom_sf(aes(fill = !!sym(variable)), color = NA) +
    scale_fill_viridis_c(option = "C", direction = -1) +
    theme_bw(base_family = "sans",base_size = 16) +
    theme(plot.title = element_blank())
}

# 12. Create all maps ---------------------------------------------------------
map_vars <- c("Mean temperature", "Annual precipitation", "Temperature seasonality", 
              "Precipitation seasonality", "Elevation", 
              "LGM mean temperature", "LGM annual precipitation")

map_plots <- purrr::map(map_vars, ~ plot_map(climate_data, .x))

# 13. Combine all maps into one figure ---------------------------------------
combined_map <- patchwork::wrap_plots(map_plots, ncol = 2) + 
  patchwork::plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text())
combined_map
# 14. Save combined figure ----------------------------------------------------
ggsave("results/Figures/all_variables.png", 
       plot = combined_map,
       dpi = 400, 
       width = 15, 
       height = 10)

######################################################################################
## Effects of climate, topography and historic climate on floristic divisions ##
######################################################################################

# Required libraries
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
library(car)
library(relaimpo)
library(broom)

# Projection (Behrmann example)
shape200_sf <- st_transform(shape200_sf, crs = behrmann_proj4)

# Hierarchical clustering for phyloregions
load("processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_200.RData")
hc200 <- hclust(beta_sim_mean200, method = "average")
# 1. Get vector of groups with names
groups_vector <- cutree(hc200, k = 6)
groups_df <- data.frame(
  idcell = names(groups_vector),
  phyloregion = factor(groups_vector))

# 2. Convert idcell in shape200_sf to character if not a character
shape200_sf <- shape200_sf %>%
  mutate(idcell = as.character(idcell))

# 3. Join by idcell
shape200_sf <- left_join(shape200_sf, groups_df, by = "idcell")

# Phyloregions map (numeric codes)
ggplot(shape200_sf) +
  geom_sf(aes(fill = phyloregion), color = NA) +
  geom_sf_text(aes(label = as.factor(phyloregion)), size = 2, color = "black") +
  scale_fill_viridis_d(name = "Phyloregion (number)", option = "H") +
  theme_minimal() +
  labs(title = "Phyloregions Map (numeric codes)")

# --- 1. Load and prepare environmental data ---
vars <- c("temp_mean", "prec_annual", "temp_seasonality",
          "precipitation_seasonality", "heterogeneity_topo",
          "LGM_temp", "LGM_prec")

env_data <- shape200_sf %>% 
  dplyr::select(idcell, phyloregion, all_of(vars)) %>%
  filter(if_all(all_of(vars), ~ !is.na(.))) %>%
  mutate(across(all_of(vars), scale)) %>%
  mutate(idcell = as.character(idcell))  # Secure character type


# --- 2. Create cell pairs and calculate distances + log-transform + scale ---
pairs_df <- expand.grid(i = env_data$idcell, j = env_data$idcell, stringsAsFactors = FALSE) %>%
  filter(i < j) %>%
  left_join(env_data, by = c("i" = "idcell")) %>%
  left_join(env_data, by = c("j" = "idcell"), suffix = c(".i", ".j")) %>%
  mutate(
    dist_climate_mean_raw = sqrt((temp_mean.i - temp_mean.j)^2 + (prec_annual.i - prec_annual.j)^2),
    dist_climate_var_raw  = sqrt((temp_seasonality.i - temp_seasonality.j)^2 + (precipitation_seasonality.i - precipitation_seasonality.j)^2),
    dist_topo_raw = sqrt((heterogeneity_topo.i - heterogeneity_topo.j)^2),
    dist_LGM_raw          = sqrt((LGM_temp.i - LGM_temp.j)^2 + (LGM_prec.i - LGM_prec.j)^2),
    # Apply log1p and scale
    dist_climate_mean = scale(log1p(dist_climate_mean_raw)),
    dist_climate_var  = scale(log1p(dist_climate_var_raw)),
    dist_topo         = scale(log1p(dist_topo_raw)),
    dist_LGM          = scale(log1p(dist_LGM_raw)))
# --- 3. Add floristic dissimilarity ---
# Convert beta_sim_mean200 in long format without coercion to integer
beta_long <- as.data.frame(as.table(as.matrix(beta_sim_mean200)))
colnames(beta_long) <- c("i", "j", "beta_sim")

# Ensure i and j are characters (match idcell)
beta_long <- beta_long %>%
  mutate(i = as.character(i), j = as.character(j)) %>%
  filter(i < j)
# Pair with pairs
pairs_df <- inner_join(pairs_df, beta_long, by = c("i", "j"))

# --- 4. Add phyloregions ---
pairs_df <- pairs_df %>%
  left_join(env_data %>% dplyr::select(idcell, phyloregion) %>% rename(phyloregion_i = phyloregion), by = c("i" = "idcell")) %>%
  left_join(env_data %>% dplyr::select(idcell, phyloregion) %>% rename(phyloregion_j = phyloregion), by = c("j" = "idcell"))

# --- 5. Label phyloregions ---
# See how many cells each group has and manually add group names
table(shape200_sf$phyloregion)
# Create manual correlation table after visual inspection
correspondencia <- tibble::tibble(
  phyloregion = factor(1:6),
  region_name = c("Holartic","Indo-Malaysian","Australian",
                  "Chile-Patagonian","Neotropical","Afrotropical"))
# Assign names by join, not by factor directly
pairs_df <- pairs_df %>%
  left_join(correspondencia, by = c("phyloregion_i" = "phyloregion")) %>%
  rename(phyloregion_i_name = region_name) %>%
  left_join(correspondencia, by = c("phyloregion_j" = "phyloregion")) %>%
  rename(phyloregion_j_name = region_name)
# --- 6. Save
save(pairs_df, file = "processed-data/drivers_bioregionalisation/pairs_df.RData")

# --- 6. Define phyloregional contexts
# 1. Australian vs all other groups
regA1 <- "Australian"
regB1 <- c("Chile-Patagonian", "Neotropical", "Afrotropical", "Holartic", "Indo-Malaysian")
# 2. Chile-Patagonian vs Neotropical y Afrotropical
regA2 <- "Chile-Patagonian"
regB2 <- c("Neotropical", "Afrotropical")
# 3. Neotropical vs Afrotropical
regA3 <- "Neotropical"
regB3 <- "Afrotropical"
# 4. Holartic vs Indo-Malaysian
regA4 <- "Holartic"
regB4 <- "Indo-Malaysian"

#revie data
summary(pairs_df$beta_sim)
hist(pairs_df$beta_sim, breaks = 30, main = "Histogram of beta_sim", xlab = "beta_sim")
summary(pairs_df$beta_sim)
library(ggplot2)
ggplot(pairs_df, aes(x = beta_sim)) +
  geom_density(fill = "lightblue") +
  labs(title = "Density plot of beta_sim")
table(pairs_df$beta_sim == 0)
table(pairs_df$beta_sim == 1)
shapiro.test(sample(pairs_df$beta_sim, 5000))

#################################################################
# --- 7. Fit models by context ---
# Contexto 1: Australian vs Others
# ================================================================

# Intra A (Australian)
cat("\n--- Intra A (Australian) ---\n")
context1_intraA <- pairs_df %>% 
  filter(phyloregion_i_name == "Australian", 
         phyloregion_j_name == "Australian")
# 1. Correlation analysis
cat("\n[1] Matriz de correlación:\n")
cor_matrix_intraA1 <- cor(context1_intraA[, c("dist_climate_mean", "dist_climate_var", "dist_topo", "dist_LGM")])
print(round(cor_matrix_intraA1, 3))
# Correlation display
library(ggplot2)
library(reshape2)
cor_melted <- melt(cor_matrix_intraA1)
ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1)) +
  coord_fixed() +
  ggtitle("Context 1 Intra A")

# 2. Modelo OLS
model1_intraA <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo,
                     data = context1_intraA, family = gaussian)
# 3. Summary of the model
cat("\n[2] Summary of the model:\n")
print(summary(model1_intraA))
# 4. VIF analysis
cat("\n[3] VIF:\n")
print(vif(model1_intraA))
# 5. Model assumptions
cat("\n[4] Diagnostic charts:\n")
par(mfrow = c(2, 2))
plot(model1_intraA)
par(mfrow = c(1, 1))

###########################################################
# Intra B (Others)
cat("\n\n--- Intra B (Others) ---\n")
context1_intraB <- pairs_df %>% 
  filter(phyloregion_i_name == phyloregion_j_name,
         phyloregion_i_name %in% c("Chile-Patagonian", "Neotropical", "Afrotropical", "Holartic", "Indo-Malaysian"))

# 1. Correlation analysis
cat("\n[1] Correlation matrix:\n")
cor_matrix_intraB1 <- cor(context1_intraB[, c("dist_climate_mean", "dist_climate_var", "dist_topo", "dist_LGM")])
print(round(cor_matrix_intraB1, 3))

# Correlation display
cor_melted <- melt(cor_matrix_intraB1)
ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1)) +
  coord_fixed() +
  ggtitle("Context 1 Intra B")

# 2. OLS model
model1_intraB <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo,
                     data = context1_intraB, family = gaussian)
# 3. Summary of the model
cat("\n[2] Summary of the model:\n")
print(summary(model1_intraB))
# 4. VIF analysis
cat("\n[3] VIF:\n")
print(vif(model1_intraB))
# 5. Model assumptions
cat("\n[4] Diagnostic charts:\n")
par(mfrow = c(2, 2))
plot(model1_intraB)
par(mfrow = c(1, 1))

# Between A and B (Australian vs the remaining groups)
cat("\n\n--- Between A and B (Australian vs Others) ---\n")
context1_inter <- pairs_df %>% 
  filter((phyloregion_i_name == "Australian" & phyloregion_j_name %in% c("Chile-Patagonian", "Neotropical", "Afrotropical", "Holartic", "Indo-Malaysian")) |
           (phyloregion_j_name == "Australian" & phyloregion_i_name %in% c("Chile-Patagonian", "Neotropical", "Afrotropical", "Holartic", "Indo-Malaysian")))

# 1. Correlation analysis
cat("\n[1] Correlation matrix:\n")
cor_matrix_inter1 <- cor(context1_inter[, c("dist_climate_mean", "dist_climate_var", "dist_topo", "dist_LGM")])
print(round(cor_matrix_inter1, 3))
# Correlation display
cor_melted <- melt(cor_matrix_inter1)
ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1)) +
  coord_fixed() +
  ggtitle("Correlation Matrix - Context 1 Between A-B")

# 2. OLS model
model1_inter <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo,
                    data = context1_inter, family = gaussian)
# 3. Summary of the model
cat("\n[2] Summary of the model:\n")
print(summary(model1_inter))
# 4. VIF analisys
cat("\n[3] VIF:\n")
print(vif(model1_inter))
# 5. Model assumptions
cat("\n[4] Diagnostic charts:\n")
par(mfrow = c(2, 2))
plot(model1_inter)
par(mfrow = c(1, 1))

# Context 2: Chile-Patagonian vs. Neotropical + Afrotropical
# ================================================================

# Intra A (Chile-Patagonian)
cat("\n--- Intra A (Chile-Patagonian) ---\n")
context2_intraA <- pairs_df %>% 
  filter(phyloregion_i_name == "Chile-Patagonian", 
         phyloregion_j_name == "Chile-Patagonian")

# 1. Correlation analisys
cat("\n[1] Correlation matrix:\n")
cor_matrix_intraA2 <- cor(context2_intraA[, c("dist_climate_mean", "dist_climate_var", "dist_topo", "dist_LGM")])
print(round(cor_matrix_intraA2, 3))

# Correlation display
cor_melted <- melt(cor_matrix_intraA2)
ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1)) +
  coord_fixed() +
  ggtitle("Correlation Matrix - Context 2 Intra A")

# 2. OLS model
model2_intraA <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo,
                     data = context2_intraA, family = gaussian)
# 3. Summary of the model
cat("\n[2] Summary of the model:\n")
print(summary(model2_intraA))
# 4. Vif analisys
cat("\n[3] VIF:\n")
print(vif(model2_intraA))
# 5. Model assumptions
cat("\n[4] Diagnostic charts:\n")
par(mfrow = c(2, 2))
plot(model2_intraA)
par(mfrow = c(1, 1))

# Intra B (Neotropical + Afrotropical)
cat("\n\n--- Intra B (Neotropical + Afrotropical) ---\n")
context2_intraB <- pairs_df %>% 
  filter(phyloregion_i_name == phyloregion_j_name,
         phyloregion_i_name %in% c("Neotropical", "Afrotropical"))

# 1. Correlation analisys
cat("\n[1] Correlation matrix:\n")
cor_matrix_intraB2 <- cor(context2_intraB[, c("dist_climate_mean", "dist_climate_var", "dist_topo", "dist_LGM")])
print(round(cor_matrix_intraB2, 3))

# Correlation display
cor_melted <- melt(cor_matrix_intraB2)
ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1)) +
  coord_fixed() +
  ggtitle("Correlation Matrix - Context 2 Intra B")

# 2. OLS model
model2_intraB <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo,
                     data = context2_intraB, family = gaussian)
# 3. Summary of the model
cat("\n[2] Summary of the model:\n")
print(summary(model2_intraB))
# 4. VIF analisys
cat("\n[3] VIF:\n")
print(vif(model2_intraB))
# 5. Model assumptions
cat("\n[4] Diagnostic charts:\n")
par(mfrow = c(2, 2))
plot(model2_intraB)
par(mfrow = c(1, 1))

# Between A and B (Chile-Patagonian vs Neotropical + Afrotropical)
cat("\n\n--- Between A and B (Chile-Patagonian vs Neotropical + Afrotropical) ---\n")
context2_inter <- pairs_df %>% 
  filter((phyloregion_i_name == "Chile-Patagonian" & phyloregion_j_name %in% c("Neotropical", "Afrotropical")) |
           (phyloregion_j_name == "Chile-Patagonian" & phyloregion_i_name %in% c("Neotropical", "Afrotropical")))

# 1. Correlation analisys
cat("\n[1] Correlation matrix:\n")
cor_matrix_inter2 <- cor(context2_inter[, c("dist_climate_mean", "dist_climate_var", "dist_topo", "dist_LGM")])
print(round(cor_matrix_inter2, 3))

# Correlation display
cor_melted <- melt(cor_matrix_inter2)
ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1)) +
  coord_fixed() +
  ggtitle("Correlation Matrix - Context 2 Between A-B")

# 2. OLS model
model2_inter <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo,
                    data = context2_inter, family = gaussian)
# 3. Summary of the model
cat("\n[2] Summary of the model:\n")
print(summary(model2_inter))
# 4. VIF analisys
cat("\n[3] VIF:\n")
print(vif(model2_inter))
# 5. Model assumptions
cat("\n[4] Diagnostic charts:\n")
par(mfrow = c(2, 2))
plot(model2_inter)
par(mfrow = c(1, 1))

# Context 3: Neotropical vs Afrotropical
# ================================================================
cat("CONTEXT 3: Neotropical vs Afrotropical\n")

# Intra A (Neotropical)
cat("\n--- Intra A (Neotropical) ---\n")
context3_intraA <- pairs_df %>% 
  filter(phyloregion_i_name == "Neotropical", 
         phyloregion_j_name == "Neotropical")

# 1. Correlation analisys
cat("\n[1] Correlation matrix:\n")
cor_matrix_intraA3 <- cor(context3_intraA[, c("dist_climate_mean", "dist_climate_var", "dist_topo", "dist_LGM")])
print(round(cor_matrix_intraA3, 3))

# Correlation display
cor_melted <- melt(cor_matrix_intraA3)
ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1)) +
  coord_fixed() +
  ggtitle("Correlation Matrix - Context 3 Intra A")

# 2. OLS model
model3_intraA <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo,
                     data = context3_intraA, family = gaussian)
# 3. Summary of the model
cat("\n[2] Summary of the modelo:\n")
print(summary(model3_intraA))
# 4. VIF analisys
cat("\n[3] VIF:\n")
print(vif(model3_intraA))
# 5. Model assuptions
cat("\n[4] Diagnostic charts:\n")
par(mfrow = c(2, 2))
plot(model3_intraA)
par(mfrow = c(1, 1))

# Intra B (Afrotropical)
cat("\n\n--- Intra B (Afrotropical) ---\n")
context3_intraB <- pairs_df %>% 
  filter(phyloregion_i_name == "Afrotropical", 
         phyloregion_j_name == "Afrotropical")

# 1. Correlation analisys
cat("\n[1] Correlation matrix:\n")
cor_matrix_intraB3 <- cor(context3_intraB[, c("dist_climate_mean", "dist_climate_var", "dist_topo", "dist_LGM")])
print(round(cor_matrix_intraB3, 3))

# Correlation display
cor_melted <- melt(cor_matrix_intraB3)
ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1)) +
  coord_fixed() +
  ggtitle("Correlation Matrix - Context 3 Intra B")

# 2. OLS model
model3_intraB <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo,
                     data = context3_intraB, family = gaussian)
# 3. Summary of the model
cat("\n[2] Summary of the model:\n")
print(summary(model3_intraB))
# 4. VIF analisys
cat("\n[3] VIF:\n")
print(vif(model3_intraB))
# 5. Model assuptions
cat("\n[4] Diagnostic charts:\n")
par(mfrow = c(2, 2))
plot(model3_intraB)
par(mfrow = c(1, 1))

# Between A and B (Neotropical vs Afrotropical)
cat("\n\n--- Entre A y B (Neotropical vs Afrotropical) ---\n")
context3_inter <- pairs_df %>% 
  filter((phyloregion_i_name == "Neotropical" & phyloregion_j_name == "Afrotropical") |
           (phyloregion_j_name == "Neotropical" & phyloregion_i_name == "Afrotropical"))

# 1. Correlation analisys
cat("\n[1] Correlation matrix:\n")
cor_matrix_inter3 <- cor(context3_inter[, c("dist_climate_mean", "dist_climate_var", "dist_topo", "dist_LGM")])
print(round(cor_matrix_inter3, 3))

# Correlation display
cor_melted <- melt(cor_matrix_inter3)
ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1)) +
  coord_fixed() +
  ggtitle("Correlation Matrix - Context 3 Between A-B")

# 2. OLS model
model3_inter <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo,
                    data = context3_inter, family = gaussian)

# 3. Summary model
cat("\n[2] Summary of the model:\n")
print(summary(model3_inter))
# 4. VIF analisys
cat("\n[3] VIF:\n")
print(vif(model3_inter))
# 5. Model assuptions
cat("\n[4] Diagnostic charts:\n")
par(mfrow = c(2, 2))
plot(model3_inter)
par(mfrow = c(1, 1))

# Context 4: Holartic vs Indo-Malaysian
# ================================================================
cat("CONTEXT 4: Holartic vs Indo-Malaysian\n")

# Intra A (Holartic)
cat("\n--- Intra A (Holartic) ---\n")
context4_intraA <- pairs_df %>% 
  filter(phyloregion_i_name == "Holartic", 
         phyloregion_j_name == "Holartic")

# 1. Correlation analisys
cat("\n[1] Correlation matrix:\n")
cor_matrix_intraA4 <- cor(context4_intraA[, c("dist_climate_mean", "dist_climate_var", "dist_topo", "dist_LGM")])
print(round(cor_matrix_intraA4, 3))

# Correlation display
cor_melted <- melt(cor_matrix_intraA4)
ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1)) +
  coord_fixed() +
  ggtitle("Correlation Matrix - Context 4 Intra A")

# 2. OLS model
model4_intraA <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo,
                     data = context4_intraA, family = gaussian)
# 3. Summary of the model
cat("\n[2] Summary of the model:\n")
print(summary(model4_intraA))
# 4. VIF analisys
cat("\n[3] VIF:\n")
print(vif(model4_intraA))
# 5. Model assumptions
cat("\n[4] Diagnostic charts:\n")
par(mfrow = c(2, 2))
plot(model4_intraA)
par(mfrow = c(1, 1))

# Intra B (Indo-Malaysian)
cat("\n\n--- Intra B (Indo-Malaysian) ---\n")
context4_intraB <- pairs_df %>% 
  filter(phyloregion_i_name == "Indo-Malaysian", 
         phyloregion_j_name == "Indo-Malaysian")
# 1. Correlation analisys
cat("\n[1] Correlation matrix:\n")
cor_matrix_intraB4 <- cor(context4_intraB[, c("dist_climate_mean", "dist_climate_var", "dist_topo", "dist_LGM")])
print(round(cor_matrix_intraB4, 3))
# Correlation display
cor_melted <- melt(cor_matrix_intraB4)
ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1)) +
  coord_fixed() +
  ggtitle("Correlation Matrix - Context 4 Intra B")
# 2. OLS model
model4_intraB <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo,
                     data = context4_intraB, family = gaussian)
# 3. Summary of the model
cat("\n[2] Summary of the model:\n")
print(summary(model4_intraB))
# 4. VIF analisys
cat("\n[3] VIF:\n")
print(vif(model4_intraB))
# 5. Model assumptions
cat("\n[4] Diagnostic charts:\n")
par(mfrow = c(2, 2))
plot(model4_intraB)
par(mfrow = c(1, 1))

# Between A y B (Holartic vs Indo-Malaysian)
cat("\n\n--- Between A y B (Holartic vs Indo-Malaysian) ---\n")
context4_inter <- pairs_df %>% 
  filter((phyloregion_i_name == "Holartic" & phyloregion_j_name == "Indo-Malaysian") |
           (phyloregion_j_name == "Holartic" & phyloregion_i_name == "Indo-Malaysian"))
# 1. Correlation analisys
cat("\n[1] Correlation matrix:\n")
cor_matrix_inter4 <- cor(context4_inter[, c("dist_climate_mean", "dist_climate_var", "dist_topo", "dist_LGM")])
print(round(cor_matrix_inter4, 3))
# Correlation display
cor_melted <- melt(cor_matrix_inter4)
ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1)) +
  coord_fixed() +
  ggtitle("Correlation Matrix - Context 4 Between A-B")

# 2. OLS model
model4_inter <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo,
                    data = context4_inter, family = gaussian)
# 3. Summary of the model
cat("\n[2] Summary of the model:\n")
print(summary(model4_inter))
# 4. VIF analisys
cat("\n[3] VIF:\n")
print(vif(model4_inter))
# 5. Model assumptions
cat("\n[4] Diagnostic charts:\n")
par(mfrow = c(2, 2))
plot(model4_inter)
par(mfrow = c(1, 1))

################################
#######Summary of Models########
library(relaimpo)
library(broom)
library(dplyr)

# Create list of models
model_list <- list(
  "Australian" = model1_intraA,
  "The remaining groups" = model1_intraB,
  "Australian vs. The remaining groups" = model1_inter,
  "Chile-Patagonian" = model2_intraA,
  "Neotropical - Afrotropical" = model2_intraB,
  "Chile-Patagonian vs. Neotropical - Afrotropical" = model2_inter,
  "Neotropical" = model3_intraA,
  "Afrotropical" = model3_intraB,
  "Neotropical vs. Afrotropical " = model3_inter,
  "Holartic" = model4_intraA,
  "Indo-Malaysian" = model4_intraB,
  "Holartic vs. Indo-Malaysian" = model4_inter)

# Function to obtain the relative importancea
get_coef_df <- function(model, context_name) {
  if (!inherits(model, "glm")) {
    return(data.frame(Context = context_name,
                      Variable = "(Not glm)",
                      Estimate = NA, StdError = NA,
                      CI_low = NA, CI_high = NA,
                      p_value = NA, Signif = NA,
                      stringsAsFactors = FALSE))
  }
  
  tidy_df <- tryCatch({
    broom::tidy(model, conf.int = TRUE)
  }, error = function(e) {
    message("❌ Error en tidy para ", context_name, ": ", e$message)
    return(NULL)
  })
  
  if (is.null(tidy_df)) {
    return(data.frame(Context = context_name,
                      Variable = "(tidy error)",
                      Estimate = NA, StdError = NA,
                      CI_low = NA, CI_high = NA,
                      p_value = NA, Signif = NA,
                      stringsAsFactors = FALSE))
  }
  
  # Add NA columns if missing
  for (col in c("term", "estimate", "std.error", "conf.low", "conf.high", "p.value")) {
    if (!col %in% names(tidy_df)) {
      tidy_df[[col]] <- NA
    }
  }
  
  # Add context and rename terms
  tidy_df$Context <- context_name
  tidy_df$Variable <- dplyr::case_when(
    tidy_df$term == "dist_climate_mean" ~ "Average climate",
    tidy_df$term == "dist_climate_var"  ~ "Climate variation",
    tidy_df$term == "dist_topo"         ~ "Topographical heterogeneity",
    tidy_df$term == "(Intercept)"       ~ "Intercept",
    TRUE ~ tidy_df$term)
  
  # Significance labels
  tidy_df$Signif <- dplyr::case_when(
    is.na(tidy_df$p.value)             ~ " ",
    tidy_df$p.value < 0.001 ~ "***",
    tidy_df$p.value < 0.01  ~ "**",
    tidy_df$p.value < 0.05  ~ "*",
    tidy_df$p.value < 0.1   ~ ".",
    TRUE ~ " ")
  # Select columns
result <- tidy_df[, c("Context", "Variable", "estimate", "std.error",
                        "conf.low", "conf.high", "p.value", "Signif")]
  # Rename
names(result) <- c("Context", "Variable", "Estimate", "StdError",
                     "CI_low", "CI_high", "p_value", "Signif")
  return(result)
}

# Create table of coefficients
tabla_resultados <- bind_rows(
  lapply(names(model_list), function(x) get_coef_df(model_list[[x]], x)))

# Function for calculating relative importance
get_relimp_df <- function(model, context_name) {
  if (!inherits(model, "lm") && !inherits(model, "glm")) return(NULL)
  
  tryCatch({
    rel <- calc.relimp(model, type = "lmg", rela = TRUE)
    data.frame(
      Context = context_name,
      Variable = names(rel@lmg),
      R2_percent = 100 * rel@lmg,
      stringsAsFactors = FALSE)
  }, error = function(e) {
    message("❌ Error en relaimpo para ", context_name, ": ", e$message)
    return(NULL)
  })
}

# Create the complete relative importance table
relimp_results <- bind_rows(
  lapply(names(model_list), function(x) get_relimp_df(model_list[[x]], x)))

# Rename variables in relimp_results to match those in results_table
relimp_results <- relimp_results %>%
  mutate(Variable = dplyr::case_when(
    Variable == "dist_climate_mean" ~ "Average climate",
    Variable == "dist_climate_var"  ~ "Climate variation",
    Variable == "dist_topo"         ~ "Topographical heterogeneity",
    TRUE ~ Variable))

# Combine with relative importance
tabla_combinada <- tabla_resultados %>%
  left_join(relimp_results, by = c("Context", "Variable"))

# SAVE
write.csv(tabla_combinada, "results/table/results_OLS_models.csv", row.names = FALSE)

################################################################################
# view the results
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggthemes)
library(stringr)
library(tidyr)
library(cowplot)

# 1. Data preparation with 3 variables)
# 12 combinations: 4 contexts × 3 comparisons

tabla_combinada <- as_tibble(tabla_combinada)

context_labels <- tibble::tibble(
  Context = rep(c("Australian", "The remaining groups", "Between groups",
                  "Chile-Patagonian", "Neotropical + Afrotropical", "Between groups",
                  "Neotropical", "Afrotropical", "Between groups",
                  "Holarctic", "Indo-Malaysian", "Between groups"), each = 3),
  Biogeo_Group = rep(c("Australia", "Chile-Patagonia", "Neotropical", "Holarctic"), each = 9),
  Context_Type = rep(rep(c("IntraA", "IntraB", "Inter"), each = 3), times = 4))

# Add row index to original table
tabla_combinada <- tabla_combinada %>%
  filter(Variable %in% c("Average climate", "Climate variation", "Topographical heterogeneity")) %>%
  mutate(row_id = row_number())

# Join by position
plot_data <- tabla_combinada %>%
  dplyr::left_join(context_labels %>% mutate(row_id = row_number()), by = "row_id") %>%
  dplyr::select(-row_id) %>%
  dplyr::mutate(
    Biogeo_Group = factor(Biogeo_Group, levels = c("Australia", "Chile-Patagonia", "Neotropical", "Holarctic")),
    Context_Type = factor(Context_Type, levels = c("IntraA", "IntraB", "Inter")),
    Variable = factor(Variable, levels = c("Average climate", "Climate variation", "Topographical heterogeneity")),
    R2_real = R2_percent  ) %>%  # <<--- Here you store the actual values
  dplyr::group_by(Biogeo_Group, Context_Type) %>%
  dplyr::mutate(
    R2_percent = ifelse(is.na(R2_real), 0, R2_real),
    Total_R2 = sum(R2_real, na.rm = TRUE),
    R2_percent = ifelse(Total_R2 > 0, (R2_real / Total_R2) * 100, 0)) %>%
  dplyr::ungroup()

library(forcats)

plot_data <- plot_data %>%
  mutate(
    Context_Type = fct_recode(Context_Type,
                              "Within A" = "IntraA",
                              "Within B" = "IntraB",
                              "Between groups" = "Inter"))


#2. Function to create individual graphs
levels(plot_data$Biogeo_Group) <- c(
  "Australian vs. The remaining groups",
  "Chile-Patagonian vs. Neotropical + Afrotropical",
  "Neotropical vs. Afrotropical",
  "Holarctic vs. Indo-Malaysian")

final_plot <- ggplot(plot_data, aes(x = Context_Type, y = R2_percent, fill = Variable)) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(
    aes(label = ifelse(R2_real > 0.5, paste0(round(R2_real, 1), "%"), "")),
    position = position_stack(vjust = 0.80),
    size = 4.5,
    color = "grey75") +
  facet_wrap(~ Biogeo_Group, ncol = 2) +
  scale_fill_manual(
    values = c(
      "Average climate" = "#440154",
      "Climate variation" = "#3B528B",
      "Topographical heterogeneity" = "#21908C"),
    name = "") +
  labs(x = NULL, y = expression(R^2~" (%)")) +
  theme_bw(base_size = 16, base_family = "sans") +
  theme(legend.position = "bottom")

final_plot
# Save
ggsave("results/Figures/Drivers.png", final_plot, 
       width = 25, height =18 , units = "cm", dpi = 400, bg = "white")
