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
# 2. Define six biogeographical realms based on phylogenetic beta diversity
# 3. Calculate unique orchid species between regions
# ==============================================================================
# 1. Load data
shape200 <- st_read("processed-data/community_matrix/pam_shape/grid_200km.gpkg")
comm200 <- readRDS("processed-data/community_matrix/pam/pam_200km.rds")
load("processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_200.RData")
rm(beta_sne_mean200, beta_sor_mean200)

# Define cluster names
# Convert dissimilarity to bioregion format
# ======================================================
#This is the same as in step 5 (building regionalisation) for creating realms.
convert_dist_to_bioregion <- function(dist_obj, metric_name = "Simpson", nb_species = NA) {
  if (!inherits(dist_obj, "dist")) stop("Objeto debe ser de clase 'dist'.")
  labels <- attr(dist_obj, "Labels")
  mat <- as.matrix(dist_obj)
  df <- data.frame(
    Site1 = rep(labels, times = length(labels)),
    Site2 = rep(labels, each = length(labels)),
    value = as.vector(mat),
    stringsAsFactors = FALSE)
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
# Adjust regionalisation to obtain metrics
# ======================================================
tree4 <- hclu_hierarclust(dissim_br, n_clust = 2:100)
eval_tree4 <- bioregionalization_metrics(tree4, dissimilarity = dissim_br, eval_metric = "pc_distance")
realms.regions <- find_optimal_n(eval_tree4, metrics_to_use = "pc_distance", criterion = "cutoff", metric_cutoffs = c(.8,.85))

optimal_realms <- realms.regions$optimal_nb_clusters$pc_distance[1]  # 6
optimal_regions <- realms.regions$optimal_nb_clusters$pc_distance[2] # 18

clusters_realm200 <- tree4$clusters[[paste0("K_", optimal_realms)]]
clusters_region200 <- tree4$clusters[[paste0("K_", optimal_regions)]]

# Assign cell names
names(clusters_realm200) <- rownames(tree4$clusters)
names(clusters_region200) <- rownames(tree4$clusters)

# Assigning realms names to groups
# Create an allocation data frame
clusters_df <- data.frame(
  idcell = names(clusters_realm200),
  cluster_realm = as.numeric(clusters_realm200))

# Join with shape200 using idcell
shape200 <- shape200 %>%
  left_join(clusters_df, by = "idcell")
#map
ggplot(shape200) +
  geom_sf(fill = NA, color = "grey40", size = 0.1) + 
  geom_sf_text(aes(label = cluster_realm), size = 3, color = "black") +
  theme_bw() +
  labs(title = "Mapa de celdas con números de reinos") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank())

# Define the names in the same order as the clusters (1 to 6)
realms <- c("Holartic","Australian","Andean-Patagonian",
            "Neotropical","Afrotropical","Indo-Malaysian")

# Create a factor with labels
shape200$realm_name <- factor(shape200$cluster_realm,
                              levels = 1:6,
                              labels = realms)

# Verify
table(shape200$realm_name)

#Metrics
# Hierarchical bioregionalisation
set.seed(1)
realms_orquis <- hclu_hierarclust(dissimilarity = dissim_br,
                                          index = names(dissim_br)[3],
                                          method = "average", n_clust = 6,
                                          optimal_tree_method = "iterative_consensus_tree")
realms_orquis$cluster_info
# Convert dgCMatrix to numerical matrix
comm200_mat <- as.matrix(comm200)
# Create a dataframe with an explicit “Site” column.
comm_df <- data.frame(
  Site = rownames(comm200_mat),
  comm200_mat,
  check.names = FALSE)
# Convert back to numerical matrix
comm200_numeric <- as.matrix(comm_df[, -1])
# We assign row names based on the “Site” column.
rownames(comm200_numeric) <- comm_df$Site
# 4. Verify
str(comm200_numeric)
# Check row name matches
all(names(realms_orquis) %in% rownames(comm200_numeric))
# Check for missing names
setdiff(names(realms_orquis), rownames(comm200_numeric))
# Obtain metrics
summary <- bioregion_metrics(
  bioregionalization = realms_orquis,
  comat = comm200_numeric)
str(comm200_numeric)
summary
write.csv(summary,"results/table/bioregion_metrics.csv", row.names = FALSE)
contrib_orquis <- site_species_metrics(realms_orquis, comm200_numeric,
                  indices = c("affinity", "fidelity", "indicator_value"))
contrib_orquis
str(contrib_orquis)
write.csv(contrib_orquis,"results/table/contribution_species.csv", row.names = FALSE)

#Count how many indicator species there are in each bioregion.
# Filter only indicator species (indval > 0 and not NA)
#
indicadoras <- contrib_orquis %>%
  filter(!is.na(indval), indval > 0)
# Contar número de especies únicas por bioregion
conteo_indicadoras <- indicadoras %>%
  group_by(Bioregion) %>%
  summarise(n_species_indicadoras = n_distinct(Species)) %>%
  arrange(Bioregion)

conteo_indicadoras
##################################################################
##################################################################
#In the review, we have been asked to perform this analysis for genders,
#families, and even tribes. We have decided to do so for genders;
#here is the change.
#Extract the genera(everything before the first space)
species <- colnames(comm200_numeric)
# Extraer el género (todo antes del primer espacio)
genera <- sub(" .*", "", species)
#Collapse all species of the same genus (sum or convert to presence/absence)
# Convertir a data.frame para manipular
comm_df <- as.data.frame(comm200_numeric)
# Split columns by gender and apply OR function (1 if any species is present)
comm_genus <- sapply(split(colnames(comm_df), genera), function(cols) {
  as.numeric(rowSums(comm_df[, cols, drop = FALSE]) > 0)
})
dim(comm_genus)
head(colnames(comm_genus))
rownames(comm_genus) <- rownames(comm200_numeric)
comm_genus <- as.matrix(comm_genus)
# Convert back to numerical matrix
comm_genus <- as.matrix(comm_genus[, -1])
str(comm_genus)
##Run again species indicators
contrib_orquis_genus <- site_species_metrics(realms_orquis, comm_genus,
                                       indices = c("affinity", "fidelity", "indicator_value"))
contrib_orquis_genus
str(contrib_orquis_genus)
write.csv(contrib_orquis_genus,"results/table/contribution_genera.csv", row.names = FALSE)

#Count how many indicator species there are in each bioregion.
# Filter only indicator species (indval > 0 and not NA)
indicadoras_genus <- contrib_orquis_genus %>%
  filter(!is.na(indval), indval > 0)
# Count number of unique species per bioregion
conteo_indicadoras_genus <- indicadoras_genus %>%
  group_by(Bioregion) %>%
  summarise(n_species_indicadoras = n_distinct(Species)) %>%
  arrange(Bioregion)

conteo_indicadoras_genus
#Bioregion n_species_indicadoras
#1                           161
#2                            67
#3                             6
#4                           337
#5                           100
#6                           231

##################################################################
##################################################################


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
rownames(comm200) <- shape200$idcell  # Ensure matrix and spatial data match
hc200 <- stats::hclust(beta_sim_mean200, method = "average")
clusters200 <- cutree(hc200, k = 6)
# Paleta fija para 6 realms
colors200 <- as.character(scico(6, palette = "batlow"))
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

ggplot(shape200_affiliation) +
  geom_sf(fill = NA, color = "grey40", size = 0.1) + 
  geom_sf_text(aes(label = realm_name), size = 3, color = "black") +
  theme_bw() +
  labs(title = "Mapa de celdas con números de reinos") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank())

# Definimos el diccionario (Holarctic con 'c' para New Phytologist)
cluster_to_realm <- c(
  "1" = "Holarctic",
  "2" = "Indo-Malayan",
  "3" = "Australian",
  "4" = "Andean-Patagonian",
  "5" = "Neotropical",
  "6" = "Afrotropical")

# Realizamos la asignación usando la columna 'group'
shape200_affiliation <- shape200_affiliation %>%
  mutate(realm_name = cluster_to_realm[as.character(group)])

# Verificamos que se haya asignado correctamente
table(shape200_affiliation$realm_name)
# 9. Create and save affiliation map ------------------------------------------
# Creamos un vector de colores nombrado para asegurar consistencia
# Usamos los nombres correctos (British English)
realm_colors <- c(
  "Holarctic"         = "#185461", 
  "Indo-Malayan"      = "#577646", 
  "Australian"        = "#001959", 
  "Andean-Patagonian" = "#B28D2E", 
  "Neotropical"       = "#FAA588", 
  "Afrotropical"      = "#F9CCF9") 

library(ggplot2)
library(dplyr)
library(ggthemes) # Para theme_map

affiliation_map <- ggplot() +
  geom_sf(data = map, fill = "grey85", colour = "white", linewidth = 0.1) +
  geom_sf(data = filter(shape200_affiliation, !low_affiliation), 
          aes(fill = realm_name), colour = NA) +
  geom_sf(data = filter(shape200_affiliation, low_affiliation), 
          fill = "gray25", colour = NA) +
  scale_fill_manual(values = realm_colors, guide = "none") +
  theme_map() +
  theme(
    text = element_text(family = "Gill Sans")) +
  coord_sf(expand = FALSE)

affiliation_map

ggsave("results/figures/Figure2.png",
       affiliation_map,dpi = 400,width = 10,height = 6,bg = "white")
