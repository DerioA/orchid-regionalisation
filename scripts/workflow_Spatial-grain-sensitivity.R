################################################################################
# Spatial Grain Sensitivity Analysis
# Quantifying bioregional fragmentation and fusion relative
# to the 200 km reference regionalisation
################################################################################

# Load required packages for spatial data handling, wrangling, and plotting
library(sf)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rnaturalearth)
library(scico)

################################################################################
# 1. Data Importation
################################################################################
# Read spatial datasets containing bioregional definitions at multiple scales
g100 <- st_read("results/SIG/bioregion100km.gpkg", quiet = TRUE)
g200 <- st_read("results/SIG/bioregion200km.gpkg", quiet = TRUE)
g400 <- st_read("results/SIG/bioregion400km.gpkg", quiet = TRUE)
g800 <- st_read("results/SIG/bioregion800km.gpkg", quiet = TRUE)

################################################################################
# 2. Bioregional Scale-Sensitivity Function
################################################################################
analyze_fragmentation <- function(alternative, 
                                  reference = g200, 
                                  cluster_col = "group") {
message("Calculating spatial intersections...")
  
# Step 2.1: Dissolve spatial grid cells into continuous bioregional polygons
# for both the baseline reference layer and the alternative grain layer.
ref_regions <- reference %>%
  group_by(ref_group = .data[[cluster_col]]) %>%
  summarise()
  
alt_regions <- alternative %>%
  group_by(alt_group = .data[[cluster_col]]) %>%
  summarise()
  
# Step 2.2: Compute geometric intersections to trace which alternative clusters 
# nest within or overlap with the 200 km baseline reference regions.
overlap <- st_intersection(ref_regions, alt_regions)
overlap$area_m2 <- as.numeric(st_area(overlap))
  
# Step 2.3: Filter out negligible geometric slivers or topological artefacts 
# caused by slight boundary mismatches along the spatial coastline gridding.
overlap <- overlap %>% filter(area_m2 > 100000)
  
# Step 2.4: Quantify the number of distinct alternative clusters present 
# within each singular baseline 200 km reference macro-region.
fragmentation_stats <- overlap %>%
  st_drop_geometry() %>%
  group_by(ref_group) %>%
  summarise(n_fragments = n_distinct(alt_group),
            fragment_ids = paste(alt_group, collapse = ", ")) %>%
  ungroup()
  
# Step 2.5: Map the resulting metrics back onto the fixed 200 km baseline geometry.
# This operation ensures that the geographic grid structure remains completely constant,
# allowing for direct visual comparison of scale vulnerability across all panels.
mapped_cells <- reference %>%
  select(group) %>%
  left_join(fragmentation_stats, by = c("group" = "ref_group")) %>%
  mutate(bioregional_status = case_when(
        n_fragments == 1 ~ "Cohesive",
        n_fragments == 2 ~ "Minor",
        n_fragments == 3 ~ "Moderate",
        n_fragments == 4 ~ "Advanced",
        n_fragments >= 5 ~ "Extreme",
        TRUE             ~ "Cohesive"),
        bioregional_status = factor(
        bioregional_status, 
        levels = c(
          "Cohesive", 
          "Minor", 
          "Moderate", 
          "Advanced", 
          "Extreme")))
  return(list(cells = mapped_cells, stats = fragmentation_stats))
}

################################################################################
# 3. Coordinate Projections & Spatial Transformations
################################################################################
# Fetch base continental landmasses for background contextual mapping
continents <- ne_countries(continent = c("Africa", "Asia", "North America",
                                         "Europe", "Oceania", "South America"),
                           returnclass = "sf", scale = "medium")
# Define the equal-area Behrmann cylindrical projection to safeguard 
# accurate geographic proportions during geometric calculation and layout design
behrmann <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
# Transform all vector objects to the synchronised Behrmann coordinate reference system (CRS)
map  <- st_transform(continents, behrmann)
g100 <- st_transform(g100, behrmann)
g200 <- st_transform(g200, behrmann)
g400 <- st_transform(g400, behrmann)
g800 <- st_transform(g800, behrmann)

################################################################################
# 4. Execute Scale Analyses
################################################################################
# Run the cross-scale sensitivity assessment against the 200 km reference system
analysis100 <- analyze_fragmentation(g100, reference = g200)
analysis400 <- analyze_fragmentation(g400, reference = g200)
analysis800 <- analyze_fragmentation(g800, reference = g200)

################################################################################
# 5. Text-Based Reporting for Manuscript Compilation
################################################################################
cat("\n--- BIOME FRAGMENTATION REPORT ---\n")
for(res in c("100", "400", "800")) {
  df <- get(paste0("analysis", res))$stats
  cat(paste0("\nAt ", res, " km grain relative to 200 km reference:\n"))
  for(i in 1:nrow(df)) {
    cat(paste0("  - Reference Region '", df$ref_group[i], "' contains ", df$n_fragments[i], " cluster(s).\n"))
  }
}

################################################################################
# 6. Cartographic Production & Composition
################################################################################
# Generate dynamic hex colours using the scico package for the 200 km baseline configuration
colors200 <- as.character(scico(6, palette = "batlow"))

# Define a standardised, uniform cmocean thermal gradient for fragmentation tracking
frag_palette <- c(
  "Cohesive"               = "#F9C378", 
  "Minor"      = "#f7aa41", 
  "Moderate"   = "#ea5e36", 
  "Advanced" = "#b82255", 
  "Extreme" = "#470b6a")

# Build a neat, publication-ready minimal theme with distinct outer panel frames
biom_theme <- theme_minimal() +
  theme(
    text = element_text(family = "Gill Sans", size = 16),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank(), 
    legend.position = "none" )

# Panel A: Baseline 200 km Bioregional Schema (Pencilled via scico batlow)
p200 <- ggplot(g200) +
  geom_sf(data = map, fill = "gray90", colour = "gray90") +
  geom_sf(aes(fill = factor(group)), colour = NA) +
  scale_fill_manual(values = colors200) +
  guides(fill = "none") + 
  biom_theme

# Panel B: 100 km Grain Projection Mapping
p100 <- ggplot(analysis100$cells) +
  geom_sf(data = map, fill = "gray90", colour = "gray90") +
  geom_sf(aes(fill = bioregional_status), colour = NA) +
  scale_fill_manual(
    values = frag_palette, 
    drop = FALSE,
    name = "Fragmentation status:") +
  guides(fill = "none") + 
  biom_theme

# Panel C: 400 km Grain Projection Mapping (Serves as the master for legend collection)
p400 <- ggplot(analysis400$cells) +
  geom_sf(data = map, fill = "gray90", colour = "gray90") +
  geom_sf(aes(fill = bioregional_status), colour = NA) +
  scale_fill_manual(
    values = frag_palette, 
    drop = FALSE, 
    name = "Fragmentation status:") +
  biom_theme # Kept active to let patchwork scrape and build the global guide

# Panel D: 800 km Grain Projection Mapping
p800 <- ggplot(analysis800$cells) +
  geom_sf(data = map, fill = "gray90", colour = "gray90") +
  geom_sf(aes(fill = bioregional_status), colour = NA) +
  scale_fill_manual(
    values = frag_palette, 
    drop = FALSE, 
    name = "Fragmentation status:") +
  guides(fill = "none") + 
  biom_theme

################################################################################
# 7. Final Layout Composition & Plot Annotations
################################################################################
# Stitch panels into a 2x2 multi-panel layout, isolate the single legend bar,
# and programmatically append elegant typographic labels (A, B, C, D)
final_plot <- (p200 | p100) / (p400 | p800) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') & 
  theme(
    text = element_text(family = "Gill Sans", size = 15),
    plot.tag = element_text(size = 15), 
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(),
    legend.key.width = unit(1.5, "cm"))

# Render output in the local plotting window
final_plot

################################################################################
# 8. File Exportation
################################################################################

# 8.1. Save as high-resolution raster PNG (ideal for quick previews and presentations)
ggsave("results/Figures/Figure1_BI.png",
  final_plot,
  width = 10,    
  height = 6,   
  dpi = 350)

# 8.2. Save as production-ready vector PDF using the robust Cairo engine
# 8.2. Save as production-ready vector PDF using the robust Cairo engine
# Ajustamos las dimensiones a un lienzo más panorámico y compacto (8 x 4.5)
ggsave(
  "results/Figures/Figure1_BI.pdf",
  final_plot + plot_layout(guides = "collect") & 
    theme(
      plot.margin = margin(1, 1, 1, 1, "mm"),       # Eliminamos márgenes al absoluto cero
      panel.spacing = unit(-2, "mm")),                 # Fuerza a los paneles a colapsar espacio entre sí
  width = 9,    
  height = 5,                                      # Al bajar el alto, el PDF ya no tiene espacio para dejar aire arriba/abajo
  device = cairo_pdf)

# 8.3. Export calculated statistics for appendices or external tabulations
write.csv(analysis100$stats, "results/table/frag_stats_100km.csv", row.names = FALSE)
write.csv(analysis400$stats, "results/table/frag_stats_400km.csv", row.names = FALSE)
write.csv(analysis800$stats, "results/table/frag_stats_800km.csv", row.names = FALSE)

################################################################################
# METHODOLOGICAL PIPELINE: SAMPLING COMPLETENESS SENSITIVITY BENCHMARK
# Project: Global Orchid Bioregionalisation Sensitivity Analysis
# Component: Step 5/6 - Data-Quality and Survey-Completeness Stress Testing
# Purpose: Quantifying spatial boundary shifts, macro-regional fragmentation, 
#          and cluster fusion under strict completeness stratification.
# Framework: Vector-based topological intersection analysis using 'sf'
################################################################################

# Ensure core spatial, statistical, and visualization packages are loaded
library(sf)          # Simple Features for R: handles topological relations and geometries
library(dplyr)       # A Grammar for Data Manipulation: handles attribute aggregation
library(ggplot2)     # Declarative Graphics: controls cartographic rendering
library(patchwork)   # The Composer of Plots: manages multi-panel figure compilation
library(rnaturalearth) # World Map Data: provides clean physical landmass vectors

################################################################################
# STEP 1: GEOSPATIAL DATA IMPORTATION & LAYER INTEGRATION
################################################################################
# Import the baseline reference grid compiled from the unrestricted comprehensive dataset
# Geometries represent a fixed 200 km fishnet mesh where 'group' denotes the 6 optimal realms
g200_full <- st_read("results/SIG/bioregion200km.gpkg", quiet = TRUE)
# Import alternative bioregional partitions generated using data restricted to high/fair grid networks
g200_well <- st_read("results/SIG/bioregion200km_WELL-SAMPLED_CELLS.gpkg", quiet = TRUE)
# Import alternative bioregional partitions generated using data restricted to undersampled networks
g200_poor <- st_read("results/SIG/bioregion200km_POORLY-SAMPLED_CELLS.gpkg", quiet = TRUE)

################################################################################
# STEP 2: CARTOGRAPHIC PROJECTION & EQUAL-AREA TRANSFORMATION
################################################################################
# Load background physical landmasses to provide geographical context for global maps
continents <- ne_countries(
  continent = c("Africa", "Asia", "North America", "Europe", "Oceania", "South America"),
  returnclass = "sf", 
  scale = "medium")
# Define the Behrmann cylindrical equal-area projection string
# This standard CRS preserves true areal proportions, which is critical for spatial overlap metrics
behrmann <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# Transform all spatial vectors to the target projection to ensure geometric compatibility
map       <- st_transform(continents, behrmann)
g200_full <- st_transform(g200_full, behrmann)
g200_well <- st_transform(g200_well, behrmann)
g200_poor <- st_transform(g200_poor, behrmann)

################################################################################
# STEP 3: REPLICABLE TOPOLOGICAL INTERSECTION & FRAGMENTATION FUNCTION
################################################################################
# This reusable function executes a vector-based intersection to track cross-scale 
# and cross-stratum partition shifts, standardising metrics back onto the fixed reference grid.
analyze_completeness_effect <- function(stratified_layer, 
                                        reference = g200_full, 
                                        cluster_col = "group") {
# Step 3.1: Dissolve discrete grid cells into large continuous regional polygons.
# This step simplifies individual pixel structures into consolidated bioregional zones.
ref_regions <- reference %>%
  group_by(ref_group = .data[[cluster_col]]) %>%
  summarise()

strat_regions <- stratified_layer %>%
  group_by(strat_group = .data[[cluster_col]]) %>%
  summarise()
  
# Step 3.2: Perform a vector spatial intersection between layers.
# This geometrically breaks overlapping regions down into intersection polygons.
overlap <- st_intersection(ref_regions, strat_regions)
  
# Step 3.3: Calculate absolute surface area for every intersecting piece.
overlap$area_m2 <- as.numeric(st_area(overlap))
  
# Step 3.4: Apply a strict geometric filter to remove boundary slivers.
# Small overlapping shapes along coastline borders are statistical noise caused by snapping 
# variations rather than true macro-regional division. We drop fragments smaller than 100,000 m2.
overlap <- overlap %>% filter(area_m2 > 100000)
  
# Step 3.5: Group by baseline region and count distinct intersecting stratified clusters.
# This assesses how many independent regional entities are spawned inside the original boundaries.
fragmentation_stats <- overlap %>%
  st_drop_geometry() %>%
  group_by(ref_group) %>%
  summarise(n_fragments = n_distinct(strat_group),
            fragment_ids = paste(strat_group, collapse = ", ")) %>%
  ungroup()
  
# Step 3.6: Index categorical Cohesion Status metrics back onto the fixed baseline 200 km cells.
# Maps maintain a constant spatial resolution whilst grid attributes dynamically track stability.
mapped_cells <- reference %>%
  select(group) %>%
  left_join(fragmentation_stats, by = c("group" = "ref_group")) %>%
  mutate(
      cohesion_status = case_when(
        n_fragments == 1 ~ "Cohesive",
        n_fragments == 2 ~ "Minor",
        n_fragments == 3 ~ "Moderate",
        n_fragments == 4 ~ "Advanced",
        n_fragments >= 5 ~ "Extreme",
        TRUE             ~ "Cohesive"),
      cohesion_status = factor(
        cohesion_status, 
        levels = c(
          "Cohesive", "Minor", "Moderate", 
          "Advanced", "Extreme")))
  return(list(cells = mapped_cells, stats = fragmentation_stats))
}

################################################################################
# STEP 4: EXECUTE WORKFLOW STRATIFICATION
################################################################################
# Benchmark high/fair survey quality configurations against the comprehensive model
analysis_well <- analyze_completeness_effect(g200_well, reference = g200_full)
# Benchmark poorly sampled configurations against the comprehensive model
analysis_poor <- analyze_completeness_effect(g200_poor, reference = g200_full)

################################################################################
# STEP 5: CARTOGRAPHIC THEMING & COLOR SPECIFICATION
################################################################################

# Standardised colour palette - Clean short names
frag_palette <- c(
  "Cohesive" = "#F9C378", 
  "Minor"    = "#f7aa41", 
  "Moderate" = "#ea5e36", 
  "Advanced" = "#b82255", 
  "Extreme"  = "#470b6a")

# Build a streamlined theme template targeting academic publication standards
# CRITICAL: We remove 'legend.position = "none"' from here so it doesn't conflict
biom_theme <- theme_minimal() +
  theme(
    text = element_text(family = "Gill Sans", size = 16),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, hjust = 0.5))

################################################################################
################################################################################
# STEP 6: PREPARE UNIFIED DATASET FOR FACETING (CRS-Safe Row Injection)
################################################################################
cells_well <- analysis_well$cells %>% mutate(stratum = "A")
cells_poor <- analysis_poor$cells %>% mutate(stratum = "B")

combined_cells <- rbind(cells_well, cells_poor)

combined_cells$stratum <- factor(
  combined_cells$stratum, 
  levels = c("A", "B"))

# Force strict evaluation of factor levels
combined_cells$cohesion_status <- factor(
  combined_cells$cohesion_status,
  levels = c("Cohesive", "Minor", "Moderate", "Advanced", "Extreme"))

# CRS-SAFE HACK: We duplicate the row but we DO NOT alter the geometry mathematically.
# Instead, we just assign its value to "Advanced". To prevent it from rendering 
# on the map, we set its coordinates to an impossible empty state (st_crs preservation).
dummy_row <- combined_cells[1, ]
dummy_row$stratum         <- "A"
dummy_row$cohesion_status <- "Advanced"

# This hides the geometry completely without destroying its CRS metadata
st_geometry(dummy_row) <- st_geometry(combined_cells[1, ])[0] 

# Now rbind will execute smoothly without CRS conflicts
combined_cells <- rbind(combined_cells, dummy_row)


################################################################################
# STEP 7: CARTOGRAPHIC PRODUCTION (External tagging and forced legend graphics)
################################################################################
final_completeness_plot <- ggplot(combined_cells) +
  # Render the global landmass background
  geom_sf(data = map, fill = "gray90", colour = "gray90") +
  # Render the 200 km reference grids colored by their cohesion status
  geom_sf(aes(fill = cohesion_status), colour = NA) +
  # Split into two vertical panels based on the stratum (A and B letters)
  facet_wrap(~stratum, ncol = 1) + 
  # Inject the custom color palette forcing all categories to stay
  scale_fill_manual(
    values = frag_palette, 
    drop = FALSE, 
    name = "Fragmentation status:",
    # TOTAL VICTORY OVER GGPLOT: override.aes forces the legend grid to draw
    # every single color block manually, bypassing any empty layer constraints.
    guide = guide_legend(
      ncol = 1,
      override.aes = list(
        fill = c("#F9C378", "#f7aa41", "#ea5e36", "#b82255", "#470b6a")))) +
  # Apply publication themes and configure the single right-hand legend
  biom_theme +
  theme(
    text = element_text(family = "Gill Sans", size = 16),
    # EXTERNAL TAGGING: Turn on facet headers to serve as external labels (A and B)
    strip.text = element_text(
      size = 16,                  # Letras grandes por fuera del mapa
      hjust = 0,                  # Alineado perfecto a la izquierda
      margin = margin(b = 5, t = 10)), 
    strip.background = element_blank(), # Remueve cualquier caja gris de fondo
    # Layout configuration for a single clean right-hand legend vertical stack
    legend.position = "right", 
    legend.box = "vertical",
    legend.title = element_text(vjust = 1),
    legend.text = element_text(size = 12),
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.5, "cm"))

# Render the layout onto the active plotting device
final_completeness_plot

################################################################################
# STEP 8: EXPORTATION OF PRODUCTION-READY REPRODUCIBLE ASSETS
################################################################################
# 8.1. Save as a crisp raster PNG (Optimized for quick sharing and presentations)
ggsave("results/Figures/Figure3_BI.png",
  final_completeness_plot,
  width = 10,    
  height = 5,   
  dpi = 350,
  bg = "white")

# 8.2. Save as a production-ready vector PDF using the robust Cairo engine.
# Notice the negative panel spacing and zero-margin configuration: this prevents
# Cairo from injecting excessive canvas padding between panels on vector plots.
ggsave("results/Figures/Figure3_BI.pdf",
  final_completeness_plot + plot_layout(guides = "collect") & 
    theme(
      plot.margin = margin(1, 1, 1, 1, "mm"),      
      panel.spacing = unit(-2, "mm")),                 
  height = 5,                                   
  device = cairo_pdf)

# 8.3. Export backing statistics as CSV logs to enable absolute transparency
write.csv(analysis_well$stats, "results/table/completeness_stats_well.csv", row.names = FALSE)
write.csv(analysis_poor$stats, "results/table/completeness_stats_poor.csv", row.names = FALSE)
################################################################################
################################################################################
################################################################################

################################################################################
# BENCHMARK: STRUCTURAL VARIATION OF BIOREGIONS ACROSS CLUSTERING ALGORITHMS
################################################################################

rm(list = ls())  # Limpiar entorno para ejecución limpia

library(sf)
library(phyloregion)
library(cluster)
library(dplyr)

# ------------------------------------------------------------------------------
# STEP 1: LOAD 200KM DATA ASSETS
# ------------------------------------------------------------------------------
message("Step 1: Loading 200 km spatial and phylogenetic matrices...")
comm200  <- readRDS("processed-data/community_matrix/pam/pam_200km.rds")
load(file = "processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_200.RData")
shape200 <- st_read(dsn = "processed-data/community_matrix/pam_shape/grid_200km.gpkg", quiet = TRUE)
rownames(comm200) <- shape200$idcell

# ------------------------------------------------------------------------------
# STEP 2: MULTI-ALGORITHM EVALUATION LOOP (Cophenetic, K-optimal, and EV)
# ------------------------------------------------------------------------------
message("Step 2: Evaluating classification shifts across all 8 algorithms...")

methods <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")

# Create a master dataframe to host the structural comparison results
algorithm_benchmark <- data.frame(
  Method = methods,
  Cophenetic_Correlation = NA,
  Optimal_K = NA,
  Explained_Variance = NA
)

for (i in seq_along(methods)) {
  m <- methods[i]
  message(paste("Processing method:", m))
  
  # A. Calculate cluster hierarchy
  hc_local <- stats::hclust(beta_sim_mean200, method = m)
  
  # B. Calculate cophenetic correlation (Data fidelity)
  algorithm_benchmark$Cophenetic_Correlation[i] <- round(cor(cophenetic(hc_local), beta_sim_mean200), 3)
  
  # C. Run optimal evaluation up to k = 30 to extract structural behavior
  opt_eval <- tryCatch({
    optimal_phyloregion(beta_sim_mean200, method = m, k = 30)
  }, error = function(e) { NULL })
  
  if (!is.null(opt_eval) && !is.null(opt_eval$optimal) && length(opt_eval$optimal) > 0) {
    
    # [CRITICAL FIX]: Force selection of the FIRST optimal K value in case 
    # the algorithm suggests multiple tie-breaking inflection points.
    k_selected <- head(opt_eval$optimal, 1)
    algorithm_benchmark$Optimal_K[i] <- k_selected
    
    # Safely extract Explained Variance for that exact single K-value
    ev_value <- opt_eval$df$ev[opt_eval$df$k == k_selected]
    
    if (length(ev_value) > 0) {
      algorithm_benchmark$Explained_Variance[i] <- round(ev_value[1], 4)
    } else {
      algorithm_benchmark$Explained_Variance[i] <- NA
    }
    
  } else {
    algorithm_benchmark$Optimal_K[i]          <- NA
    algorithm_benchmark$Explained_Variance[i]  <- NA
  }
}

# Order the final table by data fidelity (Cophenetic Correlation)
algorithm_benchmark <- algorithm_benchmark %>% arrange(desc(Cophenetic_Correlation))

# ------------------------------------------------------------------------------
# STEP 3: PRINT AND EXPORT BENCHMARK TABLE
# ------------------------------------------------------------------------------
print(algorithm_benchmark)
str(algorithm_benchmark)
algorithm_benchmark$Optimal_K <-
  sapply(algorithm_benchmark$Optimal_K, function(x) x[1])

str(algorithm_benchmark)
write.csv(algorithm_benchmark, "results/table/Algorithm_Variation_IB.csv", row.names = FALSE)
message("Success! Multi-algorithm comparison matrix successfully compiled.")


# ------------------------------------------------------------------------------
# QUICK FIX: RE-LOAD SPATIAL BACKGROUND ASSETS AND FIX TYPO
# ------------------------------------------------------------------------------
library(rnaturalearth)

# Re-defining the background and projection that were wiped from the environment
continents <- ne_countries(
  continent = c("Africa", "Asia", "North America", "Europe", "Oceania", "South America"),
  returnclass = "sf", scale = "medium"
)
behrmann <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
map_bg   <- st_transform(continents, behrmann)

# ------------------------------------------------------------------------------
# STEP 4: COMPILING MULTI-PANEL MAP CONFIGURATION (Filtering out NA Methods)
# ------------------------------------------------------------------------------
message("Step 4: Initialising multi-panel map compilation (Valid methods only)...")

cluster_colors <- as.character(scico(6, palette = "batlow"))
plot_list <- list()

for (i in seq_along(methods)) {
  m <- methods[i]
  
  # [CRITICAL FILTER]: Skip centroid and median as they did not resolve stable clusters (NA)
  if (m %in% c("centroid", "median")) {
    message(paste("Skipping method due to unresolved groups (NA):", m))
    next
  }
  
  message(paste("Generating map layout for method:", m))
  
  hc_local <- stats::hclust(beta_sim_mean200, method = m)
  clusters_local <- cutree(hc_local, k = 6)
  
  dend_local <- as.dendrogram(hc_local)
  dend_local <- color_branches(dend_local, k = 6, values = cluster_colors)
  leaves_local <- labels(dend_local)
  leaf_colors_local <- get_leaves_branches_col(dend_local)
  
  group_df <- tibble(
    idcell = leaves_local,
    group = clusters_local[leaves_local],
    color = leaf_colors_local) %>% 
    distinct(group, .keep_all = TRUE) %>%
    mutate(color = ifelse(is.na(color), cluster_colors[as.numeric(group)], color))
  
  bioregion_local <- shape200 %>% 
    mutate(group = clusters_local[as.character(idcell)]) %>% 
    left_join(group_df %>% select(group, color), by = "group") %>% 
    st_transform(behrmann)
  
  # Build clean individual map canvas
  p_map <- ggplot() +
    geom_sf(data = map_bg, fill = "gray90", colour = "gray90") +
    geom_sf(data = bioregion_local, aes(fill = color), colour = NA) +
    scale_fill_identity() +
    labs(subtitle = paste("Method:", m)) +
    theme_void() +
    theme(
      plot.subtitle = element_text(family = "Gill Sans", size = 16, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(2, 2, 2, 2, "mm"))
  
  plot_list[[m]] <- p_map
}

# ------------------------------------------------------------------------------
# STITCH THE 6 VALID MAPS INTO A PERFECT 3x2 MATRIX
# ------------------------------------------------------------------------------
message("Stitching 6 valid maps together into a 3x2 canvas via patchwork...")
library(patchwork)

# ncol = 2 over 6 plots automatically creates an elegant 3-row, 2-column grid
combined_maps_figure <- wrap_plots(plot_list, ncol = 2) + 
  plot_annotation(
    tag_levels = 'A',
    title = "",
    caption = "") & 
  theme(
    plot.title = element_text(family = "Gill Sans", size = 16, hjust = 0.5),
    plot.caption = element_text(family = "Gill Sans", size = 10, hjust = 0.5),
    plot.tag = element_text(family = "Gill Sans", size = 14),
    plot.background = element_rect(fill = "white", color = NA))
combined_maps_figure
# ------------------------------------------------------------------------------
# EXPORT MASTER FIGURE (Optimized Dimensions for 3x2 Grid)
# ------------------------------------------------------------------------------
message("Exporting final 3x2 multi-panel image to results folder...")
ggsave(
  filename = "results/Figures/Figure5_Valid_Algorithms_IB.png", 
  plot = combined_maps_figure, 
  width = 11, 
  height = 11, # Dimensiones cuadradas perfectas para 3 filas x 2 columnas de mapas
  dpi = 350, 
  bg = "white")

ggsave(
  filename = "results/Figures/Figure5_Valid_Algorithms_IB.pdf", 
  plot = combined_maps_figure, 
  width = 11, 
  height = 11, 
  device = cairo_pdf, 
  bg = "white")

message("Success! 3x2 multi-panel figure compiled and exported beautifully.")
# ------------------------------------------------------------------------------
# STITCH ALL 8 MAPS INTO A SINGLE MULTI-PANEL MATRIX
# ------------------------------------------------------------------------------
message("Stitching all 8 maps together into a unified canvas via patchwork...")
library(patchwork)

combined_maps_figure <- wrap_plots(plot_list, ncol = 2) + 
  plot_annotation(
    tag_levels = 'A',
    title = "",
    caption = "") & 
  theme(
    plot.title = element_text(family = "Gill Sans", size = 16, hjust = 0.5),
    plot.caption = element_text(family = "Gill Sans", face = "italic", size = 10, hjust = 0.5),
    plot.tag = element_text(family = "Gill Sans", size = 14),
    plot.background = element_rect(fill = "white", color = NA))
combined_maps_figure
# Export Master Figure
ggsave("results/Figures/Figure5_All_Algorithms_Maps.png", plot = combined_maps_figure, width = 11, height = 14, dpi = 350, bg = "white")
ggsave("results/Figures/Figure5_All_Algorithms_Maps.pdf", plot = combined_maps_figure, width = 11, height = 14, device = cairo_pdf, bg = "white")

message("Success! Complete multi-panel algorithm figure exported clean.")
# ------------------------------------------------------------------------------
# STITCH ALL 8 MAPS INTO A SINGLE MULTI-PANEL MATRIX
# ------------------------------------------------------------------------------
message("Stitching all 8 maps together into a unified canvas via patchwork...")

# wrap_plots takes the list and automatically arranges it into a grid
# ncol = 2 creates a balanced 4x2 vertical matrix of maps
combined_maps_figure <- wrap_plots(plot_list, ncol = 2) + 
  plot_annotation(
    tag_levels = 'A',
    title = "Comparison of Global Bioregional Boundary Behaviour Across 8 Clustering Algorithms",
    caption = "Framework demonstration at the 200 km spatial grain with a constant k = 6 partition filter."
  ) & 
  theme(
    plot.title = element_text(family = "Gill Sans", face = "bold", size = 16, hjust = 0.5),
    plot.caption = element_text(family = "Gill Sans", face = "italic", size = 10, hjust = 0.5),
    plot.tag = element_text(family = "Gill Sans", face = "bold", size = 14),
    plot.background = element_rect(fill = "white", color = NA)
  )

# ------------------------------------------------------------------------------
# EXPORT MASTER FIGURE (Solid White Canvas)
# ------------------------------------------------------------------------------
message("Exporting final multi-panel image to results folder...")

ggsave(
  filename = "results/Figures/Figure5_All_Algorithms_Maps.png", 
  plot = combined_maps_figure, 
  width = 11, 
  height = 14, # Proporciones optimizadas para una cuadrícula vertical de 4 filas x 2 columnas
  dpi = 350, 
  bg = "white"
)

ggsave(
  filename = "results/Figures/Figure5_All_Algorithms_Maps.pdf", 
  plot = combined_maps_figure, 
  width = 11, 
  height = 14, 
  device = cairo_pdf, 
  bg = "white"
)

message("Success! Complete multi-panel algorithm figure exported clean with zero text distortions.")
