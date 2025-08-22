# ==============================================================================
# Relationship Between Bioregions - Species Sharing Analysis
# ==============================================================================
# Clear environment and load required packages
rm(list = ls())  # Clear workspace

# Load libraries:
library(sf)          # For spatial data handling
library(dplyr)       # For data manipulation
library(paletteer)   # For colour palettes
library(ggplot2)     # For additional plotting
library(ggpubr)      # For publication-ready graphics
library(tidyr)
library(Matrix)
library(tibble)
library(bioregion)

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
# 1. Convertir disimilitud a formato bioregion
# ======================================================
convert_dist_to_bioregion <- function(dist_obj, metric_name = "Simpson", nb_species = NA) {
  if (!inherits(dist_obj, "dist")) stop("Objeto debe ser de clase 'dist'.")
  labels <- attr(dist_obj, "Labels")
  mat <- as.matrix(dist_obj)
  df <- data.frame(
    Site1 = rep(labels, times = length(labels)),
    Site2 = rep(labels, each = length(labels)),
    value = as.vector(mat),
    stringsAsFactors = FALSE
  )
  df <- df[df$Site1 < df$Site2, ]
  colnames(df)[3] <- metric_name
  class(df) <- c("bioregion.pairwise.metric", "data.frame")
  attr(df, "type") <- "dissimilarity"
  attr(df, "nb_sites") <- attr(dist_obj, "Size")
  attr(df, "nb_species") <- nb_species
  return(df)
}

dissim_br <- convert_dist_to_bioregion(beta_sim_mean200, metric_name = "Simpson")

# ======================================================
# 2. Ajustar regionalizacion para obtener métricas
# ======================================================
tree4 <- hclu_hierarclust(dissim_br, n_clust = 2:100)
eval_tree4 <- bioregionalization_metrics(tree4, dissimilarity = dissim_br, eval_metric = "pc_distance")
realms.regions <- find_optimal_n(eval_tree4, metrics_to_use = "pc_distance", criterion = "cutoff", metric_cutoffs = c(.8,.85))

optimal_realms <- realms.regions$optimal_nb_clusters$pc_distance[1]  # 6
optimal_regions <- realms.regions$optimal_nb_clusters$pc_distance[2] # 18

clusters_realm200 <- tree4$clusters[[paste0("K_", optimal_realms)]]
clusters_region200 <- tree4$clusters[[paste0("K_", optimal_regions)]]

# Asignar nombres de celdas
names(clusters_realm200) <- rownames(tree4$clusters)
names(clusters_region200) <- rownames(tree4$clusters)

# Assigning realms names to groups
# Crear un data.frame de asignación
clusters_df <- data.frame(
  idcell = names(clusters_realm200),
  cluster_realm = as.numeric(clusters_realm200))

# Unir con shape200 usando idcell
shape200 <- shape200 %>%
  left_join(clusters_df, by = "idcell")
#map
ggplot(shape200) +
  geom_sf(fill = NA, color = "grey40", size = 0.1) +  # solo los polígonos sin relleno
  geom_sf_text(aes(label = cluster_realm), size = 3, color = "black") + # números de reino
  theme_bw() +
  labs(title = "Mapa de celdas con números de reinos") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank())

# Definir los nombres en el mismo orden que los clusters (1 a 6)
realms <- c("Holartic","Australian","Chile-Patagonian",
            "Neotropical","Afrotropical","Indo-Malaysian")

# Crear un factor con etiquetas
shape200$realm_name <- factor(shape200$cluster_realm,
                              levels = 1:6,
                              labels = realms)

# Verificar
table(shape200$realm_name)

#Metrics
# Hierarchical bioregionalization
set.seed(1)
realms_orquis <- hclu_hierarclust(dissimilarity = dissim_br,
                                          index = names(dissim_br)[3],
                                          method = "average", n_clust = 6,
                                          optimal_tree_method = "iterative_consensus_tree")
realms_orquis$cluster_info
# 1. Convertir dgCMatrix a matrix numérica
comm200_mat <- as.matrix(comm200)  # sigue siendo numérica
# 2. Crear dataframe solo si quieres una columna explícita 'Site'
comm_df <- data.frame(
  Site = rownames(comm200_mat),
  comm200_mat,
  check.names = FALSE)
# 3. Convertir de nuevo a matrix numérica
# Seleccionamos solo las columnas de abundancia/presencia (sin 'Site')
comm200_numeric <- as.matrix(comm_df[, -1])
# Asignamos los nombres de fila según la columna 'Site'
rownames(comm200_numeric) <- comm_df$Site
# 4. Verificar
str(comm200_numeric)

# Ahora sí podemos correr las métricas
# Revisar coincidencia de nombres de filas
all(names(realms_orquis) %in% rownames(comm200_numeric))
# Revisar si hay nombres faltantes
setdiff(names(realms_orquis), rownames(comm200_numeric))
#Obtener metricas
summary <- bioregion_metrics(
  bioregionalization = realms_orquis,
  comat = comm200_numeric)
summary
write.csv(summary,"results/table/bioregion_metrics.csv", row.names = FALSE)
contrib_orquis <- site_species_metrics(realms_orquis, comm200_numeric,
                  indices = c("rho", "affinity", "fidelity", "indicator_value"))
contrib_orquis
write.csv(contrib_orquis,"results/table/contribution_species.csv", row.names = FALSE)

##################################################################
# Cell Affiliation Analysis Across Bioregions
##################################################################

library(sf)            # Spatial data handling
library(dplyr)         # Data manipulation
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
groups_factor <- factor(clusters_realm200) #Obtenido anteriormente para obtener las metricas

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
best.hclust <- tree4$algorithm$final.tree
best.hclust$labels <- names(clusters_region200)
dend200 <- as.dendrogram(best.hclust)

# Paleta fija para 6 realms
colors200 <- as.character(paletteer_c("grDevices::Spectral", 6))
# Colouring dendrogram
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
