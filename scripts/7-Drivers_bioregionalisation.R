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
elevation <- rast(paste0(chelsa_path, "elevation_1KMmd_GMTEDmd-2.tif"))
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

# 6. Plot raw raster layers ---------------------------------------------------
png("results/Figures/all_layers_highres.png", width = 3000, height = 2000, res = 400)
par(mfrow = c(4, 2))
plot(temp_mean_rescaled, main = "Average annual temperature (°C)")
plot(prec_annual_rescaled, main = "Total annual rainfall (mm/year)")
plot(temp_seasonality_rescaled, main = "Temperature seasonality (°C)")
plot(prec_seasonality_rescaled, main = "Precipitation seasonality (mm/year)")
plot(elevation_rescaled, main = "Elevation m a.s.l.")
plot(LGM_temp_rescaled, main = "LGM annual temperature (°C)")
plot(LGM_prec_rescaled, main = "LGM annual precipitation (mm/year)")
dev.off()

# 7. Load polygon grid and project to Behrmann -------------------------------
shape200_vect <- vect("processed-data/community_matrix/pam_shape/shape_reduce200.shp")
shape200_vect <- terra::project(shape200_vect, behrmann_proj4)

# 8. Extract mean values per polygon -------------------------------------------
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

# 9. Convert to sf and scale variables ----------------------------------------
shape200_sf <- st_as_sf(shape200_vect)

# 10. Prepare data frame for plotting ------------------------------------------
climate_data <- shape200_sf %>%
  st_drop_geometry() %>%
  rename(
    "Mean Temperature" = temp_mean,
    "Annual Precipitation" = prec_annual,
    "Temperature Seasonality" = temp_seasonality,
    "Precipitation Seasonality" = precipitation_seasonality,
    "Topographic Heterogeneity" = heterogeneity_topo,
    "LGM Mean Temperature" = LGM_temp,
    "LGM Annual Precipitation" = LGM_prec)

# 11. Function for thematic maps -----------------------------------------------
plot_map <- function(sf_data, variable, title) {
  ggplot(sf_data) +
    geom_sf(aes(fill = !!sym(variable)), color = NA) +
    scale_fill_viridis_c(option = "C") +
    labs(title = title) +
    theme_minimal() +
    theme(plot.title = element_text(size = 11))
}

# 12. Create all maps ---------------------------------------------------------
map_vars <- c("temp_mean", "prec_annual", "temp_seasonality", 
              "precipitation_seasonality", "heterogeneity_topo", "LGM_temp", "LGM_prec")

titles <- c("Mean Temperature (°C)", "Annual Precipitation (mm)", 
            "Temperature Seasonality (CV)", "Precipitation Seasonality (CV)",
            "Topographic Heterogeneity (m)", 
            "LGM Mean Temperature (°C)", "LGM Annual Precipitation (mm)")

map_plots <- purrr::map2(map_vars, titles, ~ plot_map(shape200_sf, .x, .y))

# 13. Combine all maps into one figure ---------------------------------------
combined_map <- patchwork::wrap_plots(map_plots, ncol = 2) + 
  patchwork::plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = 'bold'))

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
hc200 <- hclust(as.dist(beta_sim_mean200), method = "average")
shape200_sf$phyloregion <- factor(cutree(hc200, k = 6))

# Phyloregions map (numeric codes)
ggplot(shape200_sf) +
  geom_sf(aes(fill = as.factor(phyloregion)), color = NA) +
  geom_sf_text(aes(label = as.factor(phyloregion)), size = 2, color = "black") +
  scale_fill_viridis_d(name = "Phyloregion (number)", option = "H") +
  theme_minimal() +
  labs(title = "Phyloregions Map (numeric codes)")

# --- 1. Load and prepare environmental data ---
vars <- c("temp_mean", "prec_annual", 
          "temp_seasonality", "prec_seasonality", 
          "heterogeneity_topo", 
          "LGM_temp", "LGM_prec")
env_data <- shape200_sf %>% 
  dplyr::select(idcell, phyloregion, all_of(vars)) %>%
  filter(if_all(all_of(vars), ~ !is.na(.))) %>%
  mutate(across(all_of(vars), scale))

# --- 2. Create cell pairs and calculate distances ---
pairs_df <- expand.grid(i = env_data$idcell, j = env_data$idcell) %>%
  filter(i < j) %>%
  left_join(env_data, by = c("i" = "idcell")) %>%
  left_join(env_data, by = c("j" = "idcell"), suffix = c(".i", ".j")) %>%
  mutate(
    dist_climate_mean = sqrt((temp_mean.i - temp_mean.j)^2 + (prec_annual.i - prec_annual.j)^2),
    dist_climate_var  = sqrt((temp_seasonality.i - temp_seasonality.j)^2 + (prec_seasonality.i - prec_seasonality.j)^2),
    dist_topo         = abs(heterogeneity_topo.i - heterogeneity_topo.j),
    dist_LGM          = sqrt((LGM_temp.i - LGM_temp.j)^2 + (LGM_prec.i - LGM_prec.j)^2)
  )

# --- 3. Add floristic dissimilarity ---
beta_long <- as.data.frame(as.table(as.matrix(beta_sim_mean200)))
colnames(beta_long) <- c("i", "j", "beta_sim")
beta_long <- beta_long %>% mutate(i = as.integer(as.character(i)), j = as.integer(as.character(j)))
pairs_df <- inner_join(pairs_df, beta_long, by = c("i", "j"))

# --- 4. Add phyloregions ---
pairs_df <- pairs_df %>%
  left_join(env_data %>% dplyr::select(idcell, phyloregion) %>% rename(phyloregion_i = phyloregion), by = c("i" = "idcell")) %>%
  left_join(env_data %>% dplyr::select(idcell, phyloregion) %>% rename(phyloregion_j = phyloregion), by = c("j" = "idcell"))

# --- 5. Label phyloregions ---
regiones_nombres <- c("Chile-Patagonian","Australian","Neotropical", 
                      "Afrotropical","Indo-Malaysian","Holartic")

pairs_df$phyloregion_i <- factor(pairs_df$phyloregion_i,
                                 levels = 1:6,
                                 labels = regiones_nombres)

pairs_df$phyloregion_j <- factor(pairs_df$phyloregion_j,
                                 levels = 1:6,
                                 labels = regiones_nombres)
save(pairs_df, file = "processed-data/drivers_bioregionalisation/pairs_df.RData")

# --- 6. Define phyloregional contexts ---
regA1 <- "Australian"; regB1 <- c("Chile-Patagonian","Neotropical","Afrotropical","Holartic","Indo-Malaysian")
regA2 <- "Australian"; regB2 <- c("Chile-Patagonian","Neotropical","Afrotropical")
regA3 <- "Australian"; regB3 <- c("Holartic","Indo-Malaysian")
regA4 <- "Chile-Patagonian"; regB4 <- c("Neotropical","Afrotropical")
regA5 <- c("Chile-Patagonian","Neotropical","Afrotropical"); regB5 <- c("Holartic","Indo-Malaysian")
regA6 <- "Holartic"; regB6 <- "Indo-Malaysian"

# --- 7. Fit models by context ---

# Context 1
context1_intraA <- pairs_df %>% filter(phyloregion_i == regA1 & phyloregion_j == regA1)
model1_intraA <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo + dist_LGM,
                     data = context1_intraA, family = gaussian)
summary(model1_intraA)
vif(model1_intraA)
par(mfrow = c(2, 2))
plot(model1_intraA)
par(mfrow = c(1, 1))

context1_intraB <- pairs_df %>% filter(phyloregion_i %in% regB1 & phyloregion_j %in% regB1)
model1_intraB <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo + dist_LGM,
                     data = context1_intraB, family = gaussian)
summary(model1_intraB)
vif(model1_intraB)
par(mfrow = c(2, 2))
plot(model1_intraB)
par(mfrow = c(1, 1))

context1_inter <- pairs_df %>% filter((phyloregion_i == regA1 & phyloregion_j %in% regB1) |
                                        (phyloregion_j == regA1 & phyloregion_i %in% regB1))
model1_inter <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo + dist_LGM,
                    data = context1_inter, family = gaussian)
summary(model1_inter)
vif(model1_inter)
par(mfrow = c(2, 2))
plot(model1_inter)
par(mfrow = c(1, 1))

# Context 2
context2_intraA <- pairs_df %>% filter(phyloregion_i == regA2 & phyloregion_j == regA2)
model2_intraA <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo + dist_LGM,
                     data = context2_intraA, family = gaussian)
summary(model2_intraA)
vif(model2_intraA)
par(mfrow = c(2, 2))
plot(model2_intraA)
par(mfrow = c(1, 1))

context2_intraB <- pairs_df %>% filter(phyloregion_i %in% regB2 & phyloregion_j %in% regB2)
model2_intraB <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo + dist_LGM,
                     data = context2_intraB, family = gaussian)
summary(model2_intraB)
vif(model2_intraB)
par(mfrow = c(2, 2))
plot(model2_intraB)
par(mfrow = c(1, 1))

context2_inter <- pairs_df %>% filter((phyloregion_i == regA2 & phyloregion_j %in% regB2) |
                                        (phyloregion_j == regA2 & phyloregion_i %in% regB2))
model2_inter <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo + dist_LGM,
                    data = context2_inter, family = gaussian)
summary(model2_inter)
vif(model2_inter)
par(mfrow = c(2, 2))
plot(model2_inter)
par(mfrow = c(1, 1))

# Context 3
context3_intraA <- pairs_df %>% filter(phyloregion_i == regA3 & phyloregion_j == regA3)
model3_intraA <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo + dist_LGM,
                     data = context3_intraA, family = gaussian)
summary(model3_intraA)
vif(model3_intraA)
par(mfrow = c(2, 2))
plot(model3_intraA)
par(mfrow = c(1, 1))

context3_intraB <- pairs_df %>% filter(phyloregion_i %in% regB3 & phyloregion_j %in% regB3)
model3_intraB <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo + dist_LGM,
                     data = context3_intraB, family = gaussian)
summary(model3_intraB)
vif(model3_intraB)
par(mfrow = c(2, 2))
plot(model3_intraB)
par(mfrow = c(1, 1))

context3_inter <- pairs_df %>% filter((phyloregion_i == regA3 & phyloregion_j %in% regB3) |
                                        (phyloregion_j == regA3 & phyloregion_i %in% regB3))
model3_inter <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo + dist_LGM,
                    data = context3_inter, family = gaussian)
summary(model3_inter)
vif(model3_inter)
par(mfrow = c(2, 2))
plot(model3_inter)
par(mfrow = c(1, 1))

# Context 4
context4_intraA <- pairs_df %>% filter(phyloregion_i == regA4 & phyloregion_j == regA4)
model4_intraA <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo + dist_LGM,
                     data = context4_intraA, family = gaussian)
summary(model4_intraA)
vif(model4_intraA)
par(mfrow = c(2, 2))
plot(model4_intraA)
par(mfrow = c(1, 1))

context4_intraB <- pairs_df %>% filter(phyloregion_i %in% regB4 & phyloregion_j %in% regB4)
model4_intraB <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo + dist_LGM,
                     data = context4_intraB, family = gaussian)
summary(model4_intraB)
vif(model4_intraB)
par(mfrow = c(2, 2))
plot(model4_intraB)
par(mfrow = c(1, 1))

context4_inter <- pairs_df %>% filter((phyloregion_i == regA4 & phyloregion_j %in% regB4) |
                                        (phyloregion_j == regA4 & phyloregion_i %in% regB4))
model4_inter <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo + dist_LGM,
                    data = context4_inter, family = gaussian)
summary(model4_inter)
vif(model4_inter)
par(mfrow = c(2, 2))
plot(model4_inter)
par(mfrow = c(1, 1))

# Context 5: Gondwana vs Laurasia
context5_intraA <- pairs_df %>% filter(phyloregion_i %in% regA5 & phyloregion_j %in% regA5)
model5_intraA <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo + dist_LGM,
                     family = gaussian, data = context5_intraA)
summary(model5_intraA)
vif(model5_intraA)
par(mfrow = c(2,2))
plot(model5_intraA)
par(mfrow = c(1,1))

context5_intraB <- pairs_df %>% filter(phyloregion_i %in% regB5 & phyloregion_j %in% regB5)
model5_intraB <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo + dist_LGM,
                     family = gaussian, data = context5_intraB)
summary(model5_intraB)
vif(model5_intraB)
par(mfrow = c(2,2))
plot(model5_intraB)
par(mfrow = c(1,1))

context5_inter <- pairs_df %>% filter((phyloregion_i %in% regA5 & phyloregion_j %in% regB5) |
                                        (phyloregion_j %in% regA5 & phyloregion_i %in% regB5))
model5_inter <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo + dist_LGM,
                    family = gaussian, data = context5_inter)
summary(model5_inter)
vif(model5_inter)
par(mfrow = c(2,2))
plot(model5_inter)
par(mfrow = c(1,1))

# Context 6: Intra-Laurasia
context6_intraA <- pairs_df %>% filter(phyloregion_i == regA6 & phyloregion_j == regA6)
model6_intraA <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo + dist_LGM,
                     family = gaussian, data = context6_intraA)
summary(model6_intraA)
vif(model6_intraA)
par(mfrow = c(2,2))
plot(model6_intraA)
par(mfrow = c(1,1))

context6_intraB <- pairs_df %>% filter(phyloregion_i == regB6 & phyloregion_j == regB6)
model6_intraB <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo + dist_LGM,
                     family = gaussian, data = context6_intraB)
summary(model6_intraB)
vif(model6_intraB)
par(mfrow = c(2,2))
plot(model6_intraB)
par(mfrow = c(1,1))

context6_inter <- pairs_df %>% filter((phyloregion_i == regA6 & phyloregion_j == regB6) |
                                        (phyloregion_j == regA6 & phyloregion_i == regB6))
model6_inter <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo + dist_LGM,
                    family = gaussian, data = context6_inter)
summary(model6_inter)
vif(model6_inter)
par(mfrow = c(2,2))
plot(model6_inter)
par(mfrow = c(1,1))

# --- Global model: all regions ---
context_global <- pairs_df
model_global <- glm(beta_sim ~ dist_climate_mean + dist_climate_var + dist_topo + dist_LGM,
                    data = context_global, family = gaussian)
summary(model_global)
vif(model_global)
par(mfrow = c(2, 2))
plot(model_global)
par(mfrow = c(1, 1))

# --- VIF analysis ---
model_list <- list(
  model1_intraA, model1_intraB, model1_inter,
  model2_intraA, model2_intraB, model2_inter,
  model3_intraA, model3_intraB, model3_inter,
  model4_intraA, model4_intraB, model4_inter,
  model5_intraA, model5_intraB, model5_inter,
  model6_intraA, model6_intraB, model6_inter,
  model_global
)

model_names <- c(
  "1_intraA", "1_intraB", "1_inter",
  "2_intraA", "2_intraB", "2_inter",
  "3_intraA", "3_intraB", "3_inter",
  "4_intraA", "4_intraB", "4_inter",
  "5_intraA", "5_intraB", "5_inter",
  "6_intraA", "6_intraB", "6_inter",
  "global"
)

vif_table <- lapply(model_list, function(mod) {
  vif_values <- vif(mod)
  data.frame(
    Variable = names(vif_values),
    VIF = as.numeric(vif_values)
  )
})

vif_table <- Map(function(df, name) {
  df$model <- name
  df
}, vif_table, model_names)

vif_df <- bind_rows(vif_table)
vif_df_sorted <- vif_df %>% arrange(model, Variable)
print(vif_df_sorted)

# --- Relative importance analysis ---
model_list <- list(
  "1_intraA" = model1_intraA,
  "1_intraB" = model1_intraB,
  "1_inter"  = model1_inter,
  "2_intraA" = model2_intraA,
  "2_intraB" = model2_intraB,
  "2_inter"  = model2_inter,
  "3_intraA" = model3_intraA,
  "3_intraB" = model3_intraB,
  "3_inter"  = model3_inter,
  "4_intraA" = model4_intraA,
  "4_intraB" = model4_intraB,
  "4_inter"  = model4_inter,
  "5_intraA" = model5_intraA,
  "5_intraB" = model5_intraB,
  "5_inter"  = model5_inter,
  "6_intraA" = model6_intraA,
  "6_intraB" = model6_intraB,
  "6_inter"  = model6_inter,
  "global"   = model_global
)

get_relimp_df <- function(model, context_name) {
  relimp <- try(calc.relimp(model, type = "lmg", rela = TRUE), silent = TRUE)
  if(inherits(relimp, "try-error")) {
    data.frame(Context = context_name,
               Variable = c("Average climate", "Climate variation", "Altitude variation", "Climate of LGM"),
               R2_percent = NA)
  } else {
    df <- as.data.frame(relimp$lmg)
    colnames(df) <- "R2_fraction"
    df$Variable <- rownames(df)
    df$Context <- context_name
    df$R2_percent <- df$R2_fraction * 100
    dplyr::select(df, Context, Variable, R2_percent)
  }
}

model_stats_detailed <- bind_rows(
  lapply(names(model_list), function(x) get_relimp_df(model_list[[x]], x))
)

# --- Results table ---
get_coef_df <- function(model, context_name) {
  tidy_df <- try(broom::tidy(model), silent = TRUE)
  
  if(inherits(tidy_df, "try-error")) {
    data.frame(Context = context_name,
               Variable = c("Average climate", "Climate variation", "Altitude variation", "Climate of LGM"),
               Estimate = NA,
               StdError = NA,
               t_value = NA,
               p_value = NA,
               Signif = NA)
  } else {
    tidy_df %>%
      dplyr::mutate(Context = context_name,
                    Variable = term,
                    Estimate = estimate,
                    StdError = std.error,
                    t_value = statistic,
                    p_value = p.value,
                    Signif = dplyr::case_when(
                      p_value < 0.001 ~ "***",
                      p_value < 0.01  ~ "**",
                      p_value < 0.05  ~ "*",
                      p_value < 0.1   ~ ".",
                      TRUE            ~ " "
                    )) %>%
      dplyr::select(Context, Variable, Estimate, StdError, t_value, p_value, Signif)
  }
}

tabla_resultados <- bind_rows(
  lapply(names(model_list), function(x) get_coef_df(model_list[[x]], x)))

tabla_resultados <- tabla_resultados %>%
  mutate(Context = gsub("^Context", "", Context),
         Variable = case_when(
           Variable == "dist_climate_mean" ~ "Average climate",
           Variable == "dist_climate_var"  ~ "Climate variation",
           Variable == "dist_topo"         ~ "Altitude variation",
           Variable == "dist_LGM"          ~ "Climate of LGM",
           Variable == "(Intercept)"       ~ "(Intercept)",
           TRUE                           ~ Variable))

tabla_combinada <- tabla_resultados %>%
  left_join(model_stats_detailed, by = c("Context", "Variable"))
write.csv(tabla_combinada, "results/table/tabla_resultados_modelos_OLS.csv", row.names = FALSE)

# --- Final plot ---
model_stats_detailed <- model_stats_detailed %>%
  mutate(Group = ifelse(Context == "global", "7",
                        sub("_.*", "", Context)),
         Variable = dplyr::recode(Variable,
                                  "dist_climate_mean" = "Average climate",
                                  "dist_climate_var" = "Climate variation",
                                  "dist_topo" = "Altitude variation",
                                  "dist_LGM" = "Climate of LGM"),
         Group = case_when(
           Context == "global" ~ "1",
           grepl("^1_", Context) ~ "2",
           grepl("^2_", Context) ~ "3",
           grepl("^3_", Context) ~ "4",
           grepl("^4_", Context) ~ "5",
           grepl("^5_", Context) ~ "6",
           grepl("^6_", Context) ~ "7",
           TRUE ~ "NA"
         ))

model_stats_detailed$Group <- factor(model_stats_detailed$Group,
                                     levels = c("1", "2", "3", "4", "5", "6", "7"))

drivers <- ggplot(model_stats_detailed, aes(x = Context, y = R2_percent, fill = Variable)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "R² (%)", fill = "") +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  facet_wrap(~ Group, scales = "free_x") +
  scale_fill_viridis_d(direction = -1)

ggsave("results/figures/Drivers.png", drivers, dpi = 400, width = 10, height = 6)
