#############################################################
### 1.- Biogeographical Regionalisation using Hierarchical
#.      Clustering for 100km, 200km, 400km, and 800km resolution with men phylogenetic beta diversity
### 2.- Sensitivity analyses
### 3.- Validation
#############################################################
rm(list = ls())  # Clear environment
## ============================================================
## Step 1: Load Required Packages
## ============================================================
library(sf)              # Handle spatial vector data
library(phyloregion)     # Biogeographic regionalisation
library(ape)             # Tools for phylogenetic analysis and tree manipulation
library(cluster)         # Clustering algorithms and validation methods
library(recluster)       # Community clustering using resampling techniques
library(rnaturalearth)   # Download and use natural earth base maps
library(ggplot2)         # Create plots
library(cowplot)         # Combine multiple ggplot2 plots into one
library(ggthemes)        # Additional themes and style options for ggplot2
library(factoextra)      # Visualise and interpret clustering and ordination
library(dplyr)           # Efficient data manipulation
library(dendextend)      # Enhanced dendrogram customisation and comparison
library(paletteer)       # Unified access to many colour palettes
library(gridExtra)       # Arrange multiple grid-based plots
library(vegan)           # Ecological analysis, including diversity and dissimilarity indices
library(patchwork)       # Flexible system to compose ggplot2 plots
library(sabre)           # Spatial Association Between Regionalizations analysis
library(viridis)         # Color scales
library(vegan)
## ============================================================
## Step 2: Load Input Files
## Load pre-processed phylogenetic beta diversity and site data
## ============================================================
comm100 <- readRDS("processed-data/community_matrix/pam/pam_100km.rds")
load(file = "processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_100.RData")
shape100 <- st_read(dsn = "processed-data/community_matrix/pam_shape/grid_100km.gpkg")
rm(beta_sne_mean100, beta_sor_mean100)  # These metrics will not be used

## ============================================================
## Step 3: Evaluate Clustering Methods
## Use cophenetic correlation to identify best clustering algorithm
## ============================================================
methods <- c("ward.D", "ward.D2", "single", "complete", "average", 
             "mcquitty", "median", "centroid")
cophenetic_results <- data.frame(Method = methods, Cophenetic_Correlation = NA)
for (i in seq_along(methods)) {
  hc <- hclust(beta_sim_mean100, method = methods[i])
  cophenetic_results$Cophenetic_Correlation[i] <- cor(cophenetic(hc), beta_sim_mean100)
}

cophenetic_results <- cophenetic_results[order(-cophenetic_results$Cophenetic_Correlation), ]
cophenetic_results

# Visualise correlation
corr <- ggplot(cophenetic_results, aes(x = reorder(Method, -Cophenetic_Correlation), y = Cophenetic_Correlation)) +
  geom_bar(stat = "identity", fill = "grey30") +
  geom_text(aes(label = round(Cophenetic_Correlation, 2)), vjust = -0.5, size = 5) +
  labs(x = "Clustering method", y = "Cophenetic correlation") +
  theme_bw() +
  ylim(0, 1) +
  theme(
    axis.text.x = element_text(size = 16, angle = 20, hjust = 1, family = "sans"),
    axis.text.y = element_text(size = 16, family = "sans"),
    axis.title = element_text(size = 16, family = "sans"),
    legend.text = element_text(size = 16, family = "sans"),
    legend.title = element_text(size = 16, family = "sans"))
corr

## ============================================================
## Step 4: Determine Optimal Number of Clusters (k)
## Using explained variance to guide number of regions
## ============================================================
optimal_result <- optimal_phyloregion(beta_sim_mean100, method = "average", k = 20)
print(optimal_result$optimal)
quality_values <- optimal_result$df#ev = 0.46255309

optimal <- ggplot(quality_values, aes(x = k, y = ev)) +
  geom_line(color = "grey30", size = 1) +
  geom_point(color = "grey30", size = 3) +
  geom_vline(xintercept = 5, linetype = "dashed", colour = "black") +
  geom_point(aes(x = 5, y = 0.46255309), colour = "red", size = 5) +
  labs(x = "Number of clusters", y = "Explained variance") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 16, hjust = 1, family = "sans"),
    axis.text.y = element_text(size = 16, family = "sans"),
    axis.title = element_text(size = 16, family = "sans"),
    legend.text = element_text(size = 16, family = "sans"),
    legend.title = element_text(size = 16, family = "sans"))
optimal

# Combine visual comparisons
e <- grid.arrange(corr, optimal, ncol = 2)

## ============================================================
## Step 5: Cluster Sites and Visualise Dendrogram
## ============================================================
rownames(comm100) <- shape100$idcell  # Ensure matrix and spatial data match
hc100 <- stats::hclust(beta_sim_mean100, method = "average")
clusters <- cutree(hc100, k = 5)

# Colour palette
colors100 <- as.character(paletteer_c("grDevices::Spectral", 5))

# Constructing dendrograms and extracting colours
dend100 <- as.dendrogram(hc100)
dend100 <- color_branches(dend100, k = 5, col = colors100)

# Get leaf order and true colours from the dendrogram
dend_order <- order.dendrogram(dend100)
leaves <- labels(dend100)
leaf_colors <- get_leaves_branches_col(dend100)

# Associate correct cell ID, group and colour
group_df <- tibble(
  idcell = leaves,
  group = clusters[leaves],
  color = leaf_colors) %>%
  distinct(group, .keep_all = TRUE)

## Assign colours to spatial data
bioregion100 <- shape100 %>% 
  mutate(group = clusters[as.character(idcell)]) %>% 
  left_join(group_df %>% select(group, color), by = "group")

# Plot the dendrogram
dend_100 <- fviz_dend(dend100, k = 5, show_labels = FALSE, k_colors = colors100,
                      rect = FALSE, horiz = FALSE, main = "") + 
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, family = "sans"),   # Y-axis labels
    axis.title.y = element_text(size = 16, family = "sans"),  # Y-axis title
    axis.text.x = element_text(size = 16, family = "sans"),   # X-axis labels
    axis.title.x = element_text(size = 16, family = "sans")) + # X-axis title
  ylab("Mean βsim") + 
  xlab("Regions")
dend_100

# Plot the map
map <- ggplot() +
  geom_sf(data = bioregion100, aes(fill = color), colour = NA, size = 0.1) +
  scale_fill_identity() +
  theme_bw() 
map

# Save the shape
st_write(bioregion100, "results/SIG/bioregion100km.gpkg")

## ============================================================
## Step 6: Map with Behrmann Projection
## ============================================================
continents <- ne_countries(continent = c("Africa", "Asia", "North America",
                                         "Europe", "Oceania", "South America"),
                           returnclass = "sf", scale = "medium")

behrmann <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

map <- st_transform(continents, behrmann)
bioregion100 <- st_transform(bioregion100, behrmann)

map.regions100 <- ggplot() +
  geom_sf(data = map, fill = "grey60", colour = "grey60") +
  geom_sf(data = bioregion100, aes(fill = color), colour = NA, size = 0.1) +
  scale_fill_identity() +
  theme_map() +
  theme(
    legend.position = "none",
    text = element_text(size = 16),
    # Remove X-axis
    axis.title.x = element_blank(),        # X-axis title
    axis.text.x = element_blank(),         # X-axis labels
    axis.ticks.x = element_blank(),        # X-axis ticks
    axis.line.x = element_blank(),         # X-axis line
    # Remove Y-axis
    axis.title.y = element_blank(),        # Y-axis title
    axis.text.y = element_blank(),         # Y-axis labels
    axis.ticks.y = element_blank(),        # Y-axis ticks
    axis.line.y = element_blank()) +       # Y-axis line
  coord_sf(expand = FALSE)                 # Prevent automatic axis expansion
map.regions100

## ============================================================
## Step 7: NMDS Ordination of Phylogenetic Beta Diversity
## ============================================================
set.seed(123)
nmds_100 <- metaMDS(beta_sim_mean100, k = 2, trymax = 100,
                    engine = "monoMDS", autotransform = FALSE)

nmds_points100 <- as.data.frame(nmds_100$points)
colnames(nmds_points100) <- c("NMDS1", "NMDS2")

# Assignment of groups
nmds_points100$group <- factor(
  clusters[row.names(nmds_points100)],
  levels = 1:5)

group_color_map <- group_df %>% 
  select(group, dend_color = color) %>%
  mutate(group = as.character(group))

nmds100 <- ggplot(nmds_points100, aes(x = NMDS1, y = NMDS2, colour = factor(group))) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(values = setNames(group_color_map$dend_color, group_color_map$group)) +
  labs(x = "NMDS 1", y = "NMDS 2") +
  annotate("text", 
           x = min(nmds_points100$NMDS1), 
           y = min(nmds_points100$NMDS2), 
           label = paste("Stress =", round(nmds_100$stress, 3)),
           hjust = 0, size = 5, fontface = "italic") +
  theme_bw(base_size = 16, base_family = "sans") +
  theme(legend.position = "none",
        axis.title = element_text(size = 16, family = "sans"),
        axis.text = element_text(size = 16, family = "sans"))
nmds100

## ============================================================
## Step 8: Combine Outputs into Single Figure
## ============================================================
a <- map.regions100 / (nmds100 + dend_100 + plot_layout(widths = c(1, 1))) +
  plot_layout(heights = c(5, 3))
a
###################################################################
### ============================================================
##                             200 KM
## ============================================================
## Step 2: Load Input Files
## ============================================================
comm200 <- readRDS("processed-data/community_matrix/pam/pam_200km.rds")
load(file = "processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_200.RData")
shape200 <- st_read(dsn = "processed-data/community_matrix/pam_shape/grid_200km.gpkg")
rm(beta_sne_mean200, beta_sor_mean200)

## ============================================================
## Step 3: Evaluate Clustering Methods
## ============================================================
methods <- c("ward.D", "ward.D2", "single", "complete", "average", 
             "mcquitty", "median", "centroid")
cophenetic_results200 <- data.frame(Method = methods, Cophenetic_Correlation = NA)
for (i in seq_along(methods)) {
  hc <- hclust(beta_sim_mean200, method = methods[i])
  cophenetic_results200$Cophenetic_Correlation[i] <- cor(cophenetic(hc), beta_sim_mean200)
}

cophenetic_results200 <- cophenetic_results200[order(-cophenetic_results200$Cophenetic_Correlation), ]
cophenetic_results200

# Visualise correlation
corr200 <- ggplot(cophenetic_results200, aes(x = reorder(Method, -Cophenetic_Correlation), y = Cophenetic_Correlation)) +
  geom_bar(stat = "identity", fill = "grey30") +
  geom_text(aes(label = round(Cophenetic_Correlation, 2)), vjust = -0.5, size = 5) +
  labs(x = "Clustering method", y = "Cophenetic correlation") +
  theme_bw() +
  ylim(0, 1) +
  theme(
    axis.text.x = element_text(size = 16, angle = 20, hjust = 1, family = "sans"),
    axis.text.y = element_text(size = 16, family = "sans"),
    axis.title = element_text(size = 16, family = "sans"),
    legend.text = element_text(size = 16, family = "sans"),
    legend.title = element_text(size = 16, family = "sans"))
corr200

## ============================================================
## Step 4: Determine Optimal Number of Clusters (k)
## ============================================================
optimal_result200 <- optimal_phyloregion(beta_sim_mean200, method = "average", k = 20)
optimal_result200$optimal#EV=0.51
quality_values200 <- optimal_result200$df

optimal200 <- ggplot(quality_values200, aes(x = k, y = ev)) +
  geom_line(color = "grey30", size = 1) +
  geom_point(color = "grey30", size = 3) +
  geom_vline(xintercept = 6, linetype = "dashed", colour = "black") +
  geom_point(aes(x = 6, y = 0.51), colour = "red", size = 5) +
  labs(x = "Number of clusters", y = "Explained variance") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 16, hjust = 1, family = "sans"),
    axis.text.y = element_text(size = 16, family = "sans"),
    axis.title = element_text(size = 16, family = "sans"),
    legend.text = element_text(size = 16, family = "sans"),
    legend.title = element_text(size = 16, family = "sans"))
optimal200

# Combine visual comparisons
f <- grid.arrange(corr200, optimal200, ncol = 2)

## ============================================================
## Step 5: Cluster Sites and Visualise Dendrogram
## ============================================================
rownames(comm200) <- shape200$idcell  # Ensure matrix and spatial data match
hc200 <- stats::hclust(beta_sim_mean200, method = "average")
clusters200 <- cutree(hc200, k = 6)

# Colour palette
colors200 <- as.character(paletteer_c("grDevices::Spectral", 6))

# Constructing dendrograms and extracting colours
dend200 <- as.dendrogram(hc200)
dend200 <- color_branches(dend200, k = 6, col = colors200)

# Get leaf order and true colours from the dendrogram
dend_order200 <- order.dendrogram(dend200)
leaves200 <- labels(dend200)
leaf_colors200 <- get_leaves_branches_col(dend200)

# Associate correct cell ID, group and colour
group_df200 <- tibble(
  idcell = leaves200,
  group = clusters200[leaves200],
  color = leaf_colors200) %>% 
  distinct(group, .keep_all = TRUE)  # One colour per group

## Assign colours to spatial data
bioregion200 <- shape200 %>% 
  mutate(group = clusters200[as.character(idcell)]) %>% 
  left_join(group_df200 %>%
  select(group, color), by = "group")

# Plot the dendrogram
dend_200 <- fviz_dend(dend200, k = 6, show_labels = FALSE, k_colors = colors200,
                      rect = FALSE, horiz = FALSE, main = "") + 
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, family = "sans"),   # Y-axis labels
    axis.title.y = element_text(size = 16, family = "sans"),  # Y-axis title
    axis.text.x = element_text(size = 16, family = "sans"),   # X-axis labels
    axis.title.x = element_text(size = 16, family = "sans")) + # X-axis title
  ylab("Mean βsim") + 
  xlab("Regions")
dend_200

# Plot the map
map200 <- ggplot() +
  geom_sf(data = bioregion200, aes(fill = color), colour = NA, size = 0.1) +
  scale_fill_identity() +
  theme_bw() 
map200

# Save the shape
st_write(bioregion200, "results/SIG/bioregion200km.gpkg")

## ============================================================
## Step 6: Map with Behrmann Projection
## ============================================================
bioregion200 <- st_transform(bioregion200, behrmann)

# Without names
map.regions200 <- ggplot() +
  geom_sf(data = map, fill = "grey60", colour = "grey60") +
  geom_sf(data = bioregion200, aes(fill = color), colour = NA, size = 0.1) +
  scale_fill_identity() + 
  theme_map() +
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
map.regions200

# With names
groups <- c("Australian","Holarctic","Indo-Malaysian","Chile-Patagonian",
            "Neotropical","Afrotropical")

map.regions200_names <- ggplot() +
  geom_sf(data = map, fill = "grey60", colour = "grey60") +
  geom_sf(data = bioregion200, aes(fill = color), colour = NA, size = 0.1) +
  scale_fill_identity(
    name = NULL,
    breaks = colors200,
    labels = groups,
    guide = guide_legend(
      nrow = 1, # Horizontal legend
      title.position = "top",
      label.position = "bottom",
      keywidth = unit(1.2, "cm"),
      keyheight = unit(0.5, "cm"),
      direction = "horizontal")) +
  theme_map() +
  theme(
    text = element_text(family = "sans", size = 30),
    legend.position = "bottom",
    legend.justification = "center",
    legend.box.just = "center",
    legend.margin = margin(t = 10, b = 5),
    legend.spacing.x = unit(1, "cm")) +
  theme(
    text = element_text(size = 30),
    axis.title.x = element_blank(),        
    axis.text.x = element_blank(),     
    axis.ticks.x = element_blank(),    
    axis.line.x = element_blank(),    
    axis.title.y = element_blank(),    
    axis.text.y = element_blank(),    
    axis.ticks.y = element_blank(),    
    axis.line.y = element_blank()) +
  coord_sf(expand = FALSE)
map.regions200_names

## ============================================================
## Step 7: NMDS Ordination of Phylogenetic Beta Diversity
## ============================================================
set.seed(123)
nmds_200 <- metaMDS(beta_sim_mean200, k = 2, trymax = 100,
                    engine = "monoMDS", autotransform = FALSE)

nmds_points200 <- as.data.frame(nmds_200$points)
colnames(nmds_points200) <- c("NMDS1", "NMDS2")

# Assignment of groups
nmds_points200$group <- factor(
  clusters200[row.names(nmds_points200)],
  levels = 1:6)

group_color_map200 <- group_df200 %>% 
  select(group, dend_color = color) %>%
  mutate(group = as.character(group))

nmds200 <- ggplot(nmds_points200, aes(x = NMDS1, y = NMDS2, colour = factor(group))) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(values = setNames(group_color_map200$dend_color, group_color_map200$group)) +
  labs(x = "NMDS 1", y = "NMDS 2") +
  annotate("text", 
           x = min(nmds_points200$NMDS1), 
           y = min(nmds_points200$NMDS2), 
           label = paste("Stress =", round(nmds_200$stress, 3)),
           hjust = 0, size = 5, fontface = "italic") +
  theme_bw(base_size = 16, base_family = "sans") +
  theme(legend.position = "none",
        axis.title = element_text(size = 16, family = "sans"),
        axis.text = element_text(size = 16, family = "sans"))
nmds200

## ============================================================
## Step 8: Combine Outputs into Single Figure
## ============================================================

b1 <- map.regions200_names / (nmds200 + dend_200 + plot_layout(widths = c(1, 1))) +
  plot_layout(heights = c(5, 3))
b1
ggsave("results/Figures/phyloregions200km_all_names.png", b1, dpi = 400, width = 15, height = 12)

b <- map.regions200 / (nmds200 + dend_200 + plot_layout(widths = c(1, 1))) +
  plot_layout(heights = c(5, 3))
b
ggsave("results/Figures/phyloregions200km_all.png", b, dpi = 400, width = 15, height = 12)

###################################################################
### ============================================================
##                             400 KM
## ============================================================
comm400 <- readRDS("processed-data/community_matrix/pam/pam_400km.rds")
load(file = "processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_400.RData")
shape400 <- st_read(dsn = "processed-data/community_matrix/pam_shape/grid_400km.gpkg")
rm(beta_sne_mean400, beta_sor_mean400)

## ============================================================
## Step 3: Evaluate Clustering Methods
## ============================================================
methods <- c("ward.D", "ward.D2", "single", "complete", "average", 
             "mcquitty", "median", "centroid")
cophenetic_results400 <- data.frame(Method = methods, Cophenetic_Correlation = NA)
for (i in seq_along(methods)) {
  hc <- hclust(beta_sim_mean400, method = methods[i])
  cophenetic_results400$Cophenetic_Correlation[i] <- cor(cophenetic(hc), beta_sim_mean400)
}

cophenetic_results400 <- cophenetic_results400[order(-cophenetic_results400$Cophenetic_Correlation), ]
cophenetic_results400

# Visualise correlation
corr400 <- ggplot(cophenetic_results400, aes(x = reorder(Method, -Cophenetic_Correlation), y = Cophenetic_Correlation)) +
  geom_bar(stat = "identity", fill = "grey30") +
  geom_text(aes(label = round(Cophenetic_Correlation, 2)), vjust = -0.5, size = 5) +
  labs(x = "Clustering method", y = "Cophenetic correlation") +
  theme_bw() +
  ylim(0, 1) +
  theme(
    axis.text.x = element_text(size = 16, angle = 20, hjust = 1, family = "sans"),
    axis.text.y = element_text(size = 16, family = "sans"),
    axis.title = element_text(size = 16, family = "sans"),
    legend.text = element_text(size = 16, family = "sans"),
    legend.title = element_text(size = 16, family = "sans"))
corr400

## ============================================================
## Step 4: Determine Optimal Number of Clusters (k)
## ============================================================
optimal_result400 <- optimal_phyloregion(beta_sim_mean400, method = "average", k = 20)
optimal_result400$optimal
quality_values400 <- optimal_result400$df

optimal400 <- ggplot(quality_values400, aes(x = k, y = ev)) +
  geom_line(color = "grey30", size = 1) +
  geom_point(color = "grey30", size = 3) +
  geom_vline(xintercept = 20, linetype = "dashed", colour = "black") +
  geom_point(aes(x = 20, y = 0.6757274), colour = "red", size = 5) +
  labs(x = "Number of clusters", y = "Explained variance") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 16, hjust = 1, family = "sans"),
    axis.text.y = element_text(size = 16, family = "sans"),
    axis.title = element_text(size = 16, family = "sans"),
    legend.text = element_text(size = 16, family = "sans"),
    legend.title = element_text(size = 16, family = "sans"))
optimal400

# Combine visual comparisons
g <- grid.arrange(corr400, optimal400, ncol = 2)

## ============================================================
## Step 5: Cluster Sites and Visualise Dendrogram
## ============================================================
rownames(comm400) <- shape400$idcell  # Ensure matrix and spatial data match
hc400 <- stats::hclust(beta_sim_mean400, method = "average")
clusters400 <- cutree(hc400, k = 20)
# Colour palette
colors400 <- as.character(paletteer_c("grDevices::Spectral", 20))
# Constructing dendrograms and extracting colours
dend400 <- as.dendrogram(hc400)
dend400 <- color_branches(dend400, k = 20, col = colors400)
# Get leaf order and true colours from the dendrogram
dend_order400 <- order.dendrogram(dend400)
leaves400 <- labels(dend400)
leaf_colors400 <- get_leaves_branches_col(dend400)
# Associate correct cell ID, group and colour
group_df400 <- tibble(
  idcell = leaves400,
  group = clusters400[leaves400],
  color = leaf_colors400) %>% 
  distinct(group, .keep_all = TRUE)  # One colour per group
## Assign colours to spatial data
bioregion400 <- shape400 %>% 
  mutate(group = clusters400[as.character(idcell)]) %>% 
  left_join(group_df400 %>% select(group, color), by = "group")
# Plot the dendrogram
dend_400 <- fviz_dend(dend400, k = 20, show_labels = FALSE, k_colors = colors400,
                      rect = FALSE, horiz = FALSE, main = "") + theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, family = "sans"),   # Y-axis labels
    axis.title.y = element_text(size = 16, family = "sans"),  # Y-axis title
    axis.text.x = element_text(size = 16, family = "sans"),   # X-axis labels
    axis.title.x = element_text(size = 16, family = "sans")) + # X-axis title
  ylab("Mean βsim") + 
  xlab("Regions")
dend_400

# Plot the map
map400 <- ggplot() +
  geom_sf(data = bioregion400, aes(fill = color), colour = NA, size = 0.1) +
  scale_fill_identity() +
  theme_bw() 
map400

# Save the shape
st_write(bioregion400, "results/SIG/bioregion400km.gpkg")

## ============================================================
## Step 6: Map with Behrmann Projection
## ============================================================
bioregion400 <- st_transform(bioregion400, behrmann)

map.regions400 <- ggplot() +
  geom_sf(data = map, fill = "grey60", colour = "grey60") +
  geom_sf(data = bioregion400, aes(fill = color), colour = NA, size = 0.1) +
  scale_fill_identity() +
  theme_map() +
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
map.regions400

## ============================================================
## Step 7: NMDS Ordination of Phylogenetic Beta Diversity
## ============================================================
set.seed(123)
nmds_400 <- metaMDS(beta_sim_mean400, k = 2, trymax = 100,
                    engine = "monoMDS", autotransform = FALSE)

nmds_points400 <- as.data.frame(nmds_400$points)
colnames(nmds_points400) <- c("NMDS1", "NMDS2")

# Assignment of groups
nmds_points400$group <- factor(
  clusters400[row.names(nmds_points400)],
  levels = 1:20)

group_color_map400 <- group_df400 %>% 
  select(group, dend_color = color) %>%
  mutate(group = as.character(group))

nmds400 <- ggplot(nmds_points400, aes(x = NMDS1, y = NMDS2, colour = factor(group))) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(values = setNames(group_color_map400$dend_color, group_color_map400$group)) +
  labs(x = "NMDS 1", y = "NMDS 2") +
  annotate("text", 
           x = min(nmds_points400$NMDS1), 
           y = min(nmds_points400$NMDS2), 
           label = paste("Stress =", round(nmds_400$stress, 3)),
           hjust = 0, size = 5, fontface = "italic") +
  theme_bw(base_size = 16, base_family = "sans") +
  theme(legend.position = "none",
        axis.title = element_text(size = 16, family = "sans"),
        axis.text = element_text(size = 16, family = "sans"))
nmds400

## ============================================================
## Step 8: Combine Outputs into Single Figure
## ============================================================
c <- map.regions400 / (nmds400 + dend_400 + plot_layout(widths = c(1, 1))) +
  plot_layout(heights = c(5, 3))
c

###################################################################
### ============================================================
##                              800 KM
## ============================================================
comm800 <- readRDS("processed-data/community_matrix/pam/pam_800km.rds")
load(file = "processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_800.RData")
shape800 <- st_read(dsn = "processed-data/community_matrix/pam_shape/grid_800km.gpkg")
rm(beta_sne_mean800, beta_sor_mean800)

## ============================================================
## Step 3: Evaluate Clustering Methods
## ============================================================
methods <- c("ward.D", "ward.D2", "single", "complete", "average", 
             "mcquitty", "median", "centroid")
cophenetic_results800 <- data.frame(Method = methods, Cophenetic_Correlation = NA)
for (i in seq_along(methods)) {
  hc <- hclust(beta_sim_mean800, method = methods[i])
  cophenetic_results800$Cophenetic_Correlation[i] <- cor(cophenetic(hc), beta_sim_mean800)
}

cophenetic_results800 <- cophenetic_results800[order(-cophenetic_results800$Cophenetic_Correlation), ]
cophenetic_results800

# Visualise correlation
corr800 <- ggplot(cophenetic_results800, aes(x = reorder(Method, -Cophenetic_Correlation), y = Cophenetic_Correlation)) +
  geom_bar(stat = "identity", fill = "grey30") +
  geom_text(aes(label = round(Cophenetic_Correlation, 2)), vjust = -0.5, size = 5) +
  labs(x = "Clustering method", y = "Cophenetic correlation") +
  theme_bw() +
  ylim(0, 1) +
  theme(
    axis.text.x = element_text(size = 16, angle = 20, hjust = 1, family = "sans"),
    axis.text.y = element_text(size = 16, family = "sans"),
    axis.title = element_text(size = 16, family = "sans"),
    legend.text = element_text(size = 16, family = "sans"),
    legend.title = element_text(size = 16, family = "sans"))
corr800

## ============================================================
## Step 4: Determine Optimal Number of Clusters (k)
## ============================================================
optimal_result800 <- optimal_phyloregion(beta_sim_mean800, method = "average", k = 20)
optimal_result800$optimal
quality_values800 <- optimal_result800$df

optimal800 <- ggplot(quality_values800, aes(x = k, y = ev)) +
  geom_line(color = "grey30", size = 1) +
  geom_point(color = "grey30", size = 3) +
  geom_vline(xintercept = 16, linetype = "dashed", colour = "black") +
  geom_point(aes(x = 16, y = 0.74903708), colour = "red", size = 5) +
  labs(x = "Number of clusters", y = "Explained variance") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 16, hjust = 1, family = "sans"),
    axis.text.y = element_text(size = 16, family = "sans"),
    axis.title = element_text(size = 16, family = "sans"),
    legend.text = element_text(size = 16, family = "sans"),
    legend.title = element_text(size = 16, family = "sans"))
optimal800

# Combine visual comparisons
h <- grid.arrange(corr800, optimal800, ncol = 2)

## ============================================================
## Step 5: Cluster Sites and Visualise Dendrogram
## ============================================================
rownames(comm800) <- shape800$idcell  # Ensure matrix and spatial data match
hc800 <- stats::hclust(beta_sim_mean800, method = "average")
clusters800 <- cutree(hc800, k = 16)

# Colour palette
colors800 <- as.character(paletteer_c("grDevices::Spectral", 16))

# Constructing dendrograms and extracting colours
dend800 <- as.dendrogram(hc800)
dend800 <- color_branches(dend800, k = 16, col = colors800)

# Get leaf order and true colours from the dendrogram
dend_order800 <- order.dendrogram(dend800)
leaves800 <- labels(dend800)
leaf_colors800 <- get_leaves_branches_col(dend800)

# Associate correct cell ID, group and colour
group_df800 <- tibble(
  idcell = leaves800,
  group = clusters800[leaves800],
  color = leaf_colors800) %>% 
  distinct(group, .keep_all = TRUE)  # One colour per group

## Assign colours to spatial data
bioregion800 <- shape800 %>% 
  mutate(group = clusters800[as.character(idcell)]) %>% 
  left_join(group_df800 %>% select(group, color), by = "group")

# Plot the dendrogram
dend_800 <- fviz_dend(dend800, k = 16, show_labels = FALSE, k_colors = colors800,
                      rect = FALSE, horiz = FALSE, main = "") + theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, family = "sans"),   # Y-axis labels
    axis.title.y = element_text(size = 16, family = "sans"),  # Y-axis title
    axis.text.x = element_text(size = 16, family = "sans"),   # X-axis labels
    axis.title.x = element_text(size = 16, family = "sans")) + # X-axis title
  ylab("Mean βsim") + 
  xlab("Regions")
dend_800

# Plot the map
map800 <- ggplot() +
  geom_sf(data = bioregion800, aes(fill = color), colour = NA, size = 0.1) +
  scale_fill_identity() +
  theme_bw() 
map800

# Save the shape
st_write(bioregion800, "results/SIG/bioregion800km.gpkg")

## ============================================================
## Step 6: Map with Behrmann Projection
## ============================================================
bioregion800 <- st_transform(bioregion800, behrmann)

map.regions800 <- ggplot() +
  geom_sf(data = map, fill = "grey60", colour = "grey60") +
  geom_sf(data = bioregion800, aes(fill = color), colour = NA, size = 0.1) +
  scale_fill_identity() +
  theme_map() +
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
map.regions800

## ============================================================
## Step 7: NMDS Ordination of Phylogenetic Beta Diversity
## ============================================================
set.seed(123)
nmds_800 <- metaMDS(beta_sim_mean800, k = 2, trymax = 100,
                    engine = "monoMDS", autotransform = FALSE)

nmds_points800 <- as.data.frame(nmds_800$points)
colnames(nmds_points800) <- c("NMDS1", "NMDS2")

# Assignment of groups
nmds_points800$group <- factor(
  clusters800[row.names(nmds_points800)],
  levels = 1:16)

group_color_map800 <- group_df800 %>% 
  select(group, dend_color = color) %>%
  mutate(group = as.character(group))

nmds800 <- ggplot(nmds_points800, aes(x = NMDS1, y = NMDS2, colour = factor(group))) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(values = setNames(group_color_map800$dend_color, group_color_map800$group)) +
  labs(x = "NMDS 1", y = "NMDS 2") +
  annotate("text", 
           x = min(nmds_points800$NMDS1), 
           y = min(nmds_points800$NMDS2), 
           label = paste("Stress =", round(nmds_800$stress, 3)),
           hjust = 0, size = 5, fontface = "italic") +
  theme_bw(base_size = 16, base_family = "sans") +
  theme(legend.position = "none",
        axis.title = element_text(size = 16, family = "sans"),
        axis.text = element_text(size = 16, family = "sans"))
nmds800

## ============================================================
## Step 8: Combine Outputs into Single Figure
## ============================================================
d <- map.regions800 / (nmds800 + dend_800 + plot_layout(widths = c(1, 1))) +
  plot_layout(heights = c(5, 3))
d

#################################################################
## Join all figures
## ============================================================
# Maps
all <- wrap_plots(
  wrap_elements(a),
  wrap_elements(b),
  wrap_elements(c),
  wrap_elements(d),
  nrow = 2, ncol = 2
) &
  theme(plot.margin = margin(0, 0, 0, 0))
all
ggsave("results/Figures/all-maps.png", all, dpi = 400, width = 18, height = 15)

# Number of clusters and clustering algorithm
all_algorithm <- grid.arrange(e, f, g, h,
                              layout_matrix = rbind(c(1, 2),
                                                    c(3, 4)),
                              heights = c(5, 5))
ggsave("results/Figures/all-algorithm.png", all_algorithm, dpi = 400, width = 20, height = 13)

###################################################################
###################################################################
### ============================================================
##                   Sensitivity analyses   
##             (100 km, 200 km, 400 km, 800 km) 
## =============================================================
# Silhouette analysis to assess clustering performance
##############
# 100 km analysis
# Verify input data structure
class(beta_sim_mean100)  # Should be "dist"
class(hc100)             # Should be "hclust"
length(clusters)         # Must match observation count

# Calculate silhouette for k = 5 clusters
sil100 <- silhouette(clusters, beta_sim_mean100)
# Statistical summary
summary(sil100)  # Average silhouette per cluster and overall

# Prepare silhouette plot
sil100[,"cluster"] <- factor(sil100[,"cluster"], levels = 1:5)

sil_100 <- fviz_silhouette(sil100) +
  aes(fill = cluster, color = cluster) +
  scale_fill_manual(values = colors100, name = "Cluster") +
  scale_color_manual(values = colors100, guide = "none") +
  geom_hline(
    yintercept = mean(sil100[, "sil_width"]),
    linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(
    title = "a) 100 km (k = 5)",
    x = "",
    y = "Silhouette width") +
  theme_bw(base_size = 16, base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0, margin = margin(b = 15)),
    legend.position = "right",
    legend.key.size = unit(0.8, "cm"),
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank())
sil_100
# Cluster summary:
#cluster size ave.sil.width
#1      305          0.48
#2     1199          0.26
#3      631          0.26
#4      664          0.19
#5 2    204          0.38

##############
# 200 km analysis
class(beta_sim_mean200)  # Verify distance structure
class(hc200)             # Verify clustering object
length(clusters200)      # Check observation count

# Calculate silhouette for k = 6 clusters
sil200 <- silhouette(clusters200, beta_sim_mean200)
# Statistical summary
summary(sil200)  # Average silhouette per cluster and overall

# Prepare silhouette plot
sil200[,"cluster"] <- factor(sil200[,"cluster"], levels = 1:6)

sil_200 <- fviz_silhouette(sil200) +
  aes(fill = cluster, color = cluster) +
  scale_fill_manual(values = colors200, name = "Cluster") +
  scale_color_manual(values = colors200, guide = "none") +
  geom_hline(
    yintercept = mean(sil200[, "sil_width"]),
    linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(
    title = "b) 200 km (k = 6)",
    x = "",
    y = "Silhouette width") +
  theme_bw(base_size = 16, base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0, margin = margin(b = 15)),
    legend.position = "right",
    legend.key.size = unit(0.8, "cm"),
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank())
sil_200
# Cluster summary:
# cluster size ave.sil.width
#       1   955          0.41
#       2   335          0.19
#       3   109          0.53
#       4   10           0.90
#       5   473          0.28
#       6   318          0.30

##############
# 400 km analysis
class(beta_sim_mean400)  # Verify distance structure
class(hc400)             # Verify clustering object
length(clusters400)      # Check observation count

# Calculate silhouette for k = 14 clusters
sil400 <- silhouette(clusters400, beta_sim_mean400)
# Statistical summary
summary(sil400)  # Average silhouette per cluster and overall

# Prepare silhouette plot
sil400[,"cluster"] <- factor(sil400[,"cluster"], levels = 1:20)

sil_400 <- fviz_silhouette(sil400) +
  aes(fill = cluster, color = cluster) +
  scale_fill_manual(values = colors400, name = "Cluster") +
  scale_color_manual(values = colors400, guide = "none") +
  geom_hline(
    yintercept = mean(sil400[, "sil_width"]),
    linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(
    title = "c) 400 km (k = 20)",
    x = "",
    y = "Silhouette width") +
  theme_bw(base_size = 16, base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0, margin = margin(b = 15)),
    legend.position = "right",
    legend.key.size = unit(0.8, "cm"),
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank())
sil_400
# Cluster summary:
#cluster size ave.sil.width
#1        1  102          0.22
#2        2   13          0.47
#3        3   47          0.44
#4        4   15          0.48
#5        5   68         -0.08
#6        6   35          0.05
#7        7  145          0.04
#8        8    1          0.00
#9        9    4          0.07
#10      10    1          0.00
#11      11    1          0.00
#12      12   14          0.54
#13      13    1          0.00
#14      14   23          0.34
#15      15  309          0.19
#16      16    1          0.00
#17      17    1          0.00
#18      18   57          0.49
#19      19    7          0.94
#20      20    5          0.80

##############
# 800 km analysis
class(beta_sim_mean800)  # Verify distance structure
class(hc800)             # Verify clustering object
length(clusters800)      # Check observation count

# Calculate silhouette for k = 16 clusters
sil800 <- silhouette(clusters800, beta_sim_mean800)
# Statistical summary
summary(sil800)  # Average silhouette per cluster and overall

# Prepare silhouette plot
sil800[,"cluster"] <- factor(sil800[,"cluster"], levels = 1:16)

sil_800 <- fviz_silhouette(sil800) +
  aes(fill = cluster, color = cluster) +
  scale_fill_manual(values = colors800, name = "Cluster") +
  scale_color_manual(values = colors800, guide = "none") +
  geom_hline(
    yintercept = mean(sil800[, "sil_width"]),
    linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(
    title = "d) 800 km (k = 16)",
    x = "",
    y = "Silhouette width") +
  theme_bw(base_size = 16, base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0, margin = margin(b = 15)),
    legend.position = "right",
    legend.key.size = unit(0.8, "cm"),
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank())

sil_800
# Cluster summary:
#cluster size ave.sil.width
#1        1    5          0.82
#2        2   16          0.29
#3        3   38          0.16
#4        4    3          0.98
#5        5    5          0.82
#6        6   37          0.22
#7        7    6          0.57
#8        8   21          0.14
#9        9    6          0.69
#10      10   12          0.19
#11      11   11          0.23
#12      12   14          0.19
#13      13    5          0.52
#14      14    1          0.00
#15      15  120          0.34
#16      16   16          0.53
##############
# Composite plot
all_sil <- grid.arrange(sil_100, sil_200, sil_400, sil_800,
                        layout_matrix = rbind(c(1, 2),
                                              c(3, 4)),
                        heights = c(5, 5))

ggsave("results/Figures/all-silhouette.png", all_sil, 
       dpi = 400, width = 17, height = 15)

##############
# Comparative analysis
# Calculate average silhouette widths
sil100_avg <- mean(sil100[, "sil_width"])
sil200_avg <- mean(sil200[, "sil_width"])
sil400_avg <- mean(sil400[, "sil_width"])
sil800_avg <- mean(sil800[, "sil_width"])

# Generate comparative table
grain_comparison <- data.frame(
  "Tamaño_grano" = c(100, 200, 400, 800),
  "Silueta_promedio" = c(sil100_avg, sil200_avg, sil400_avg, sil800_avg))

# Display results
print(grain_comparison)
#Tamaño_grano Silueta_promedio
#1          100        0.3180056
#2          200        0.3398030
#3          400        0.2022927
#4          800        0.3166438

#####################################################################
# VALIDATION ANALYSIS
# Regionalization using only well-sampled cells (200km resolution)
#####################################################################

# ============================================================
# SECTION 1: DATA LOADING AND PREPARATION
# ============================================================
# Load completeness estimators data
load("processed-data/data_exploration/Completeness/Estimators200.RData")

# Read shapefile and rename column
shape200 <- st_read(dsn = "processed-data/community_matrix/pam_shape/grid_200km.gpkg")
colnames(estimators)[1] <- "idcell"

# Load phylogenetic beta diversity components
load("processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_200.RData")

# ============================================================
# SECTION 2: DATA QUALITY ASSESSMENT
# ============================================================
# Join estimators with spatial data
completeness200 <- shape200 %>% 
  left_join(estimators, by = "idcell")

# Filter NA values in completeness
names(completeness200)
without_estimation200 <- completeness200 %>% filter(is.na(Completeness))
completeness200 <- completeness200 %>% filter(!is.na(Completeness))

# Classify cells by sampling quality
completeness200 <- completeness200 %>%
  mutate(Quality = NA_character_) %>%
  mutate(Quality = case_when(
    Slope < 0.02 & Completeness > 90 & Ratio > 15 ~ "High",
    Slope > 0.3 & Completeness < 50 & Ratio < 3 ~ "Poor",
    TRUE ~ Quality
  )) %>%
  mutate(Quality = ifelse(is.na(Quality), "Fair", Quality))

# Convert Quality to factor
completeness200 <- completeness200 %>%
  mutate(Quality = factor(Quality, levels = c("High", "Fair", "Poor")))
table(completeness200$Quality)

# ============================================================
# WELL-SAMPLED CELLS ANALYSIS (HIGH + FAIR QUALITY)
# ============================================================
# Create subset of High and Fair quality cells
completeness_subset <- completeness200 %>%
  filter(Quality %in% c("High", "Fair")) %>%
  select(idcell) # 1,473 grid cells (66.95 %)

# Subset shapefile
shape200_subset <- shape200 %>%
  filter(idcell %in% completeness_subset$idcell)

# ============================================================
# SECTION 4: PHYLOGENETIC REGIONALIZATION - WELL-SAMPLED CELLS
# ============================================================
# Load phylogenetic beta diversity data
load(file = "processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_200.RData")
rm(beta_sne_mean200, beta_sor_mean200)
str(beta_sim_mean200)

# Prepare subset of beta diversity matrix
subset_ids <- as.character(shape200_subset$idcell)
labels_dist <- attr(beta_sim_mean200, "Labels")
keep_ids <- labels_dist[labels_dist %in% subset_ids]
beta_mat <- as.matrix(beta_sim_mean200)
beta_subset <- beta_mat[keep_ids, keep_ids]
beta_subset_dist <- as.dist(beta_subset)

# ============================================================
# SECTION 5: CLUSTERING METHOD EVALUATION
# ============================================================
methods <- c("ward.D", "ward.D2", "single", "complete", "average", 
             "mcquitty", "median", "centroid")
cophenetic_results200 <- data.frame(Method = methods, Cophenetic_Correlation = NA)

for (i in seq_along(methods)) {
  hc <- hclust(beta_subset_dist, method = methods[i])
  cophenetic_results200$Cophenetic_Correlation[i] <- cor(cophenetic(hc), beta_subset_dist)
}

cophenetic_results200 <- cophenetic_results200[order(-cophenetic_results200$Cophenetic_Correlation), ]
cophenetic_results200

# Visualization of cophenetic correlation
corr200_subset <- ggplot(cophenetic_results200, aes(x = reorder(Method, -Cophenetic_Correlation), y = Cophenetic_Correlation)) +
  geom_bar(stat = "identity", fill = "grey30") +
  geom_text(aes(label = round(Cophenetic_Correlation, 2)), vjust = -0.5, size = 5) +
  labs(x = "Clustering method", y = "Cophenetic correlation") +
  theme_bw() +
  ylim(0, 1) +
  theme(
    axis.text.x = element_text(size = 16, angle = 20, hjust = 1, family = "sans"),
    axis.text.y = element_text(size = 16, family = "sans"),
    axis.title = element_text(size = 16, family = "sans"),
    legend.text = element_text(size = 16, family = "sans"),
    legend.title = element_text(size = 16, family = "sans"))
corr200_subset
# ============================================================
# SECTION 6: OPTIMAL NUMBER OF CLUSTERS DETERMINATION
# ============================================================
optimal_result200 <- optimal_phyloregion(beta_subset_dist, method = "average", k = 20)
optimal_result200$optimal
quality_values200 <- optimal_result200$df

optimal200_subset <- ggplot(quality_values200, aes(x = k, y = ev)) +
  geom_line(color = "grey30", size = 1) +
  geom_point(color = "grey30", size = 3) +
  geom_vline(xintercept = 6, linetype = "dashed", colour = "black") +
  geom_point(aes(x = 6, y = 0.51), colour = "red", size = 5) +
  labs(x = "Number of clusters", y = "Explained variance") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 16, hjust = 1, family = "sans"),
    axis.text.y = element_text(size = 16, family = "sans"),
    axis.title = element_text(size = 16, family = "sans"),
    legend.text = element_text(size = 16, family = "sans"),
    legend.title = element_text(size = 16, family = "sans"))
optimal200_subset
# Combine visualizations
s <- grid.arrange(corr200_subset, optimal200_subset, ncol = 2)
# ============================================================
# SECTION 7: CLUSTER ANALYSIS AND VISUALIZATION
# ============================================================
rownames(comm200) <- shape200$idcell
hc200_subset <- stats::hclust(beta_subset_dist, method = "average")
clusters200_subset <- cutree(hc200_subset, k = 6)

# Color palette
colors200_subset <- as.character(paletteer_c("grDevices::Spectral", 6))

# Dendrogram construction
dend200_subset <- as.dendrogram(hc200_subset)
dend200_subset <- color_branches(dend200_subset, k = 6, col = colors200_subset)

# Extract dendrogram information
dend_order200_subset <- order.dendrogram(dend200_subset)
leaves200_subset <- labels(dend200_subset)
leaf_colors200_subset <- get_leaves_branches_col(dend200_subset)

# Create group-color mapping
group_df200_subset <- tibble(
  idcell = leaves200_subset,
  group = clusters200_subset[leaves200_subset],
  color = leaf_colors200_subset) %>% distinct(group, .keep_all = TRUE)

# Assign colors to spatial data
bioregion200_subset <- shape200_subset %>% 
  mutate(group = clusters200_subset[as.character(idcell)]) %>% 
  left_join(group_df200_subset %>% select(group, color), by = "group")

# ============================================================
# SECTION 8: DENDROGRAM AND MAP VISUALIZATION
# ============================================================
# Plot dendrogram
dend_200_subset <- fviz_dend(dend200_subset, k = 6, show_labels = FALSE, k_colors = colors200_subset,
                             rect = FALSE, horiz = FALSE, main = "") + theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, family = "sans"),
    axis.title.y = element_text(size = 16, family = "sans"),
    axis.text.x = element_text(size = 16, family = "sans"),
    axis.title.x = element_text(size = 16, family = "sans")) +
  ylab("Mean βsim") + 
  xlab("regions")
dend_200_subset
# Plot map
map200_subset <- ggplot() +
  geom_sf(data = bioregion200_subset, aes(fill = color), colour = NA, size = 0.1) +
  scale_fill_identity() +
  theme_map()
map200_subset

st_write(bioregion200_subset, "results/SIG/bioregion200km_WELL-SAMPLED_CELLS.gpkg")
# ============================================================
# SECTION 9: BEHRMANN PROJECTION MAP
# ============================================================
bioregion200_subset <- st_transform(bioregion200_subset, behrmann)

map.regions200_subset <- ggplot() +
  geom_sf(data = map, fill = "grey60", colour = "grey60") +
  geom_sf(data = bioregion200_subset, aes(fill = color), colour = NA, size = 0.1) +
  scale_fill_identity() + 
  theme_map() +
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
map.regions200_subset
# ============================================================
# SECTION 10: NMDS ORDINATION
# ============================================================
set.seed(123)
nmds_200_subset <- metaMDS(beta_subset_dist, k = 2, trymax = 100,
                           engine = "monoMDS", autotransform = FALSE)

nmds_points200_subset <- as.data.frame(nmds_200_subset$points)
colnames(nmds_points200_subset) <- c("NMDS1", "NMDS2")

nmds_points200_subset$group <- factor(
  clusters200_subset[row.names(nmds_points200_subset)],
  levels = 1:6)

group_color_map200_subset <- group_df200_subset %>% 
  select(group, dend_color = color) %>%
  mutate(group = as.character(group))

nmds200_subset <- ggplot(nmds_points200_subset, aes(x = NMDS1, y = NMDS2, colour = factor(group))) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(
    values = setNames(group_color_map200_subset$dend_color, group_color_map200_subset$group)) +
  labs(x = "NMDS 1", y = "NMDS 2") +
  annotate("text", 
           x = min(nmds_points200_subset$NMDS1), 
           y = min(nmds_points200_subset$NMDS2), 
           label = paste("Stress =", round(nmds_200_subset$stress, 3)),
           hjust = 0, size = 5, fontface = "italic") +
  theme_bw(base_size = 16, base_family = "sans") +
  theme(legend.position = "none",
        axis.title = element_text(size = 16, family = "sans"),
        axis.text = element_text(size = 16, family = "sans"))
nmds200_subset
# ============================================================
# SECTION 11: COMBINE OUTPUTS
# ============================================================
x <- map.regions200_subset / (nmds200_subset + dend_200_subset + plot_layout(widths = c(1, 1))) +
  plot_layout(heights = c(5, 3))
x
####################################################################
# POORLY-SAMPLED CELLS ANALYSIS (POOR QUALITY)
#####################################################################

# ============================================================
# POORLY-SAMPLED CELLS SUBSET
# ============================================================
completeness_subset2 <- completeness200 %>%
  filter(Quality %in% c("Poor")) %>%
  select(idcell) # 643 grid cells (30.7%)

shape200_subset2 <- shape200 %>%
  filter(idcell %in% completeness_subset2$idcell)

# ============================================================
# SECTION 13: PHYLOGENETIC REGIONALIZATION - POORLY-SAMPLED CELLS
# ============================================================
load(file = "processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_200.RData")
rm(beta_sne_mean200, beta_sor_mean200)

subset_ids <- as.character(shape200_subset2$idcell)
labels_dist <- attr(beta_sim_mean200, "Labels")
keep_ids <- labels_dist[labels_dist %in% subset_ids]
beta_mat <- as.matrix(beta_sim_mean200)
beta_subset <- beta_mat[keep_ids, keep_ids]
beta_subset_dist2 <- as.dist(beta_subset)

# ============================================================
# SECTION 14: CLUSTERING METHOD EVALUATION
# ============================================================
methods <- c("ward.D", "ward.D2", "single", "complete", "average", 
             "mcquitty", "median", "centroid")
cophenetic_results200 <- data.frame(Method = methods, Cophenetic_Correlation = NA)

for (i in seq_along(methods)) {
  hc <- hclust(beta_subset_dist2, method = methods[i])
  cophenetic_results200$Cophenetic_Correlation[i] <- cor(cophenetic(hc), beta_subset_dist2)
}

cophenetic_results200 <- cophenetic_results200[order(-cophenetic_results200$Cophenetic_Correlation), ]

corr200_subset2 <- ggplot(cophenetic_results200, aes(x = reorder(Method, -Cophenetic_Correlation), y = Cophenetic_Correlation)) +
  geom_bar(stat = "identity", fill = "grey30") +
  geom_text(aes(label = round(Cophenetic_Correlation, 2)), vjust = -0.5, size = 5) +
  labs(x = "Clustering method", y = "Cophenetic correlation") +
  theme_bw() +
  ylim(0, 1) +
  theme(
    axis.text.x = element_text(size = 16, angle = 20, hjust = 1, family = "sans"),
    axis.text.y = element_text(size = 16, family = "sans"),
    axis.title = element_text(size = 16, family = "sans"),
    legend.text = element_text(size = 16, family = "sans"),
    legend.title = element_text(size = 16, family = "sans"))
corr200_subset2
# ============================================================
# SECTION 15: OPTIMAL NUMBER OF CLUSTERS
# ============================================================
optimal_result200 <- optimal_phyloregion(beta_subset_dist2, method = "average", k = 20)
optimal_result200$optimal
quality_values200 <- optimal_result200$df

optimal200_subset2 <- ggplot(quality_values200, aes(x = k, y = ev)) +
  geom_line(color = "grey30", size = 1) +
  geom_point(color = "grey30", size = 3) +
  geom_vline(xintercept = 20, linetype = "dashed", colour = "black") +
  geom_point(aes(x = 20, y = 0.50121577), colour = "red", size = 5) +
  labs(x = "Number of clusters", y = "Explained variance") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 16, hjust = 1, family = "sans"),
    axis.text.y = element_text(size = 16, family = "sans"),
    axis.title = element_text(size = 16, family = "sans"),
    legend.text = element_text(size = 16, family = "sans"),
    legend.title = element_text(size = 16, family = "sans"))
optimal200_subset2
ww <- grid.arrange(corr200_subset2, optimal200_subset2, ncol = 2)

# ============================================================
# SECTION 16: CLUSTER ANALYSIS AND VISUALIZATION
# ============================================================
rownames(comm200) <- shape200$idcell
hc200_subset2 <- stats::hclust(beta_subset_dist2, method = "average")
clusters200_subset2 <- cutree(hc200_subset2, k = 20)

colors200_subset2 <- as.character(paletteer_c("grDevices::Spectral", 20))

dend200_subset2 <- as.dendrogram(hc200_subset2)
dend200_subset2 <- color_branches(dend200_subset2, k = 20, col = colors200_subset2)

dend_order200_subset2 <- order.dendrogram(dend200_subset2)
leaves200_subset2 <- labels(dend200_subset2)
leaf_colors200_subset2 <- get_leaves_branches_col(dend200_subset2)

group_df200_subset2 <- tibble(
  idcell = leaves200_subset2,
  group = clusters200_subset2[leaves200_subset2],
  color = leaf_colors200_subset2) %>% distinct(group, .keep_all = TRUE)

bioregion200_subset2 <- shape200_subset2 %>% 
  mutate(group = clusters200_subset2[as.character(idcell)]) %>% 
  left_join(group_df200_subset2 %>% select(group, color), by = "group")

# ============================================================
# SECTION 17: DENDROGRAM AND MAP VISUALIZATION
# ============================================================
dend_200_subset2 <- fviz_dend(dend200_subset2, k = 20, show_labels = FALSE, k_colors = colors200_subset2,
                              rect = FALSE, horiz = FALSE, main = "") + theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, family = "sans"),
    axis.title.y = element_text(size = 16, family = "sans"),
    axis.text.x = element_text(size = 16, family = "sans"),
    axis.title.x = element_text(size = 16, family = "sans")) +
  ylab("Mean βsim") + 
  xlab("regions")
dend_200_subset2

map200_subset2 <- ggplot() +
  geom_sf(data = bioregion200_subset2, aes(fill = color), colour = NA, size = 0.1) +
  scale_fill_identity() +
  theme_bw()
map200_subset2

st_write(bioregion200_subset2, "results/SIG/bioregion200km_POORLY-SAMPLED_CELLS.gpkg")

# ============================================================
# SECTION 18: BEHRMANN PROJECTION MAP
# ============================================================
bioregion200_subset2 <- st_transform(bioregion200_subset2, behrmann)

map.regions200_subset2 <- ggplot() +
  geom_sf(data = map, fill = "grey60", colour = "grey60") +
  geom_sf(data = bioregion200_subset2, aes(fill = color), colour = NA, size = 0.1) +
  scale_fill_identity() + 
  theme_map() +
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
map.regions200_subset2
# ============================================================
# SECTION 19: NMDS ORDINATION
# ============================================================
set.seed(123)
nmds_200_subset2 <- metaMDS(beta_subset_dist2, k = 2, trymax = 100,
                            engine = "monoMDS", autotransform = FALSE)

nmds_points200_subset2 <- as.data.frame(nmds_200_subset2$points)
colnames(nmds_points200_subset2) <- c("NMDS1", "NMDS2")

nmds_points200_subset2$group <- factor(
  clusters200_subset2[row.names(nmds_points200_subset2)],
  levels = 1:20)

group_color_map200_subset2 <- group_df200_subset2 %>% 
  select(group, dend_color = color) %>%
  mutate(group = as.character(group))

nmds200_subset2 <- ggplot(nmds_points200_subset2, aes(x = NMDS1, y = NMDS2, colour = factor(group))) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(
    values = setNames(group_color_map200_subset2$dend_color, group_color_map200_subset2$group)) +
  labs(x = "NMDS 1", y = "NMDS 2") +
  annotate("text", 
           x = min(nmds_points200_subset2$NMDS1), 
           y = min(nmds_points200_subset2$NMDS2), 
           label = paste("Stress =", round(nmds_200_subset2$stress, 3)),
           hjust = 0, size = 5, fontface = "italic") +
  theme_bw(base_size = 16, base_family = "sans") +
  theme(legend.position = "none",
        axis.title = element_text(size = 16, family = "sans"),
        axis.text = element_text(size = 16, family = "sans"))
nmds200_subset2
# ============================================================
# SECTION 20: COMBINE OUTPUTS
# ============================================================
yy <- map.regions200_subset2 / (nmds200_subset2 + dend_200_subset2 + plot_layout(widths = c(1, 1))) +
  plot_layout(heights = c(5, 3))
yy
x
#####################################################################
# SORENSEN INDEX ANALYSIS
#####################################################################

# ============================================================
# SECTION 21: SORENSEN INDEX REGIONALIZATION
# ============================================================
load(file = "processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_200.RData")
shape200 <- st_read(dsn = "processed-data/community_matrix/pam_shape/grid_200km.gpkg")

# ============================================================
# SECTION 22: CLUSTERING METHOD EVALUATION
# ============================================================
methods <- c("ward.D", "ward.D2", "single", "complete", "average", 
             "mcquitty", "median", "centroid")
cophenetic_results200 <- data.frame(Method = methods, Cophenetic_Correlation = NA)

for (i in seq_along(methods)) {
  hc <- hclust(beta_sor_mean200, method = methods[i])
  cophenetic_results200$Cophenetic_Correlation[i] <- cor(cophenetic(hc), beta_sor_mean200)
}

cophenetic_results200_sor <- cophenetic_results200[order(-cophenetic_results200$Cophenetic_Correlation), ]

corr200_sor <- ggplot(cophenetic_results200_sor, aes(x = reorder(Method, -Cophenetic_Correlation), y = Cophenetic_Correlation)) +
  geom_bar(stat = "identity", fill = "grey30") +
  geom_text(aes(label = round(Cophenetic_Correlation, 2)), vjust = -0.5, size = 5) +
  labs(x = "Clustering method", y = "Cophenetic correlation") +
  theme_bw() +
  ylim(0, 1) +
  theme(
    axis.text.x = element_text(size = 16, angle = 20, hjust = 1, family = "sans"),
    axis.text.y = element_text(size = 16, family = "sans"),
    axis.title = element_text(size = 16, family = "sans"),
    legend.text = element_text(size = 16, family = "sans"),
    legend.title = element_text(size = 16, family = "sans"))
corr200_sor
# ============================================================
# SECTION 23: OPTIMAL NUMBER OF CLUSTERS
# ============================================================
optimal_result200_sor <- optimal_phyloregion(beta_sor_mean200, method = "average", k = 20)
optimal_result200_sor$optimal
quality_values200_sor <- optimal_result200_sor$df

optimal200_sor <- ggplot(quality_values200_sor, aes(x = k, y = ev)) +
  geom_line(color = "grey30", size = 1) +
  geom_point(color = "grey30", size = 3) +
  geom_vline(xintercept = 9, linetype = "dashed", colour = "black") +
  geom_point(aes(x = 9, y = 0.34), colour = "red", size = 5) +
  labs(x = "Number of clusters", y = "Explained variance") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 16, hjust = 1, family = "sans"),
    axis.text.y = element_text(size = 16, family = "sans"),
    axis.title = element_text(size = 16, family = "sans"),
    legend.text = element_text(size = 16, family = "sans"),
    legend.title = element_text(size = 16, family = "sans"))
optimal200_sor
i <- grid.arrange(corr200_sor, optimal200_sor, ncol = 2)
# ============================================================
# SECTION 24: CLUSTER ANALYSIS AND VISUALIZATION
# ============================================================
rownames(comm200) <- shape200$idcell
hc200_sor <- stats::hclust(beta_sor_mean200, method = "average")
clusters200_sor <- cutree(hc200_sor, k = 9)

colors200_sor <- as.character(paletteer_c("grDevices::Spectral", 9))

dend200_sor <- as.dendrogram(hc200_sor)
dend200_sor <- color_branches(dend200_sor, k = 9, col = colors200_sor)

dend_order200_sor <- order.dendrogram(dend200_sor)
leaves200_sor <- labels(dend200_sor)
leaf_colors200_sor <- get_leaves_branches_col(dend200_sor)

group_df200_sor <- tibble(
  idcell = leaves200_sor,
  group = clusters200_sor[leaves200_sor],
  color = leaf_colors200_sor) %>% distinct(group, .keep_all = TRUE)

bioregion200_sor <- shape200 %>% 
  mutate(group = clusters200_sor[as.character(idcell)]) %>% 
  left_join(group_df200_sor %>% select(group, color), by = "group")

# ============================================================
# SECTION 25: DENDROGRAM AND MAP VISUALIZATION
# ============================================================
dend_200_sor <- fviz_dend(dend200_sor, k = 9, show_labels = FALSE, k_colors = colors200_sor,
                          rect = FALSE, horiz = FALSE, main = "") + theme_bw() +
  theme(
    axis.text.y = element_text(size = 16, family = "sans"),
    axis.title.y = element_text(size = 16, family = "sans"),
    axis.text.x = element_text(size = 16, family = "sans"),
    axis.title.x = element_text(size = 16, family = "sans")) +
  ylab("Mean βsim") + 
  xlab("regions")
dend_200_sor

map200_sor <- ggplot() +
  geom_sf(data = bioregion200_sor, aes(fill = color), colour = NA, size = 0.1) +
  scale_fill_identity() +
  theme_bw()
map200_sor

st_write(bioregion200_sor, "results/SIG/bioregion200km_sorensen.gpkg")

# ============================================================
# SECTION 26: BEHRMANN PROJECTION MAP
# ============================================================
bioregion200_sor <- st_transform(bioregion200_sor, behrmann)

map.regions200_sor <- ggplot() +
  geom_sf(data = map, fill = "grey60", colour = "grey60") +
  geom_sf(data = bioregion200_sor, aes(fill = color), colour = NA, size = 0.1) +
  scale_fill_identity() + 
  theme_map() +
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
map.regions200_sor
# ============================================================
# SECTION 27: NMDS ORDINATION
# ============================================================
set.seed(123)
nmds_200_sor <- metaMDS(beta_sor_mean200, k = 2, trymax = 50,
                        engine = "monoMDS", autotransform = FALSE)

nmds_points200_sor <- as.data.frame(nmds_200_sor$points)
colnames(nmds_points200_sor) <- c("NMDS1", "NMDS2")

nmds_points200_sor$group <- factor(
  clusters200_sor[row.names(nmds_points200_sor)],
  levels = 1:9)

group_color_map200_sor <- group_df200_sor %>% 
  select(group, dend_color = color) %>%
  mutate(group = as.character(group))

nmds200_sor <- ggplot(nmds_points200_sor, aes(x = NMDS1, y = NMDS2, colour = factor(group))) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(
    values = setNames(group_color_map200_sor$dend_color, group_color_map200_sor$group)) +
  labs(x = "NMDS 1", y = "NMDS 2") +
  annotate("text", 
           x = min(nmds_points200_sor$NMDS1), 
           y = min(nmds_points200_sor$NMDS2), 
           label = paste("Stress =", round(nmds_200_sor$stress, 3)),
           hjust = 0, size = 5, fontface = "italic") +
  theme_bw(base_size = 16, base_family = "sans") +
  theme(legend.position = "none",
        axis.title = element_text(size = 16, family = "sans"),
        axis.text = element_text(size = 16, family = "sans"))
nmds200_sor
# ============================================================
# SECTION 28: COMBINE OUTPUTS
# ============================================================
w <- map.regions200_sor / (nmds200_sor + dend_200_sor + plot_layout(widths = c(1, 1))) +
  plot_layout(heights = c(5, 3))
w
# ============================================================
# SECTION 29: COMBINE ALL MAPS
# ============================================================
final_plot <- wrap_elements(x) / 
  wrap_elements(yy) / 
  wrap_elements(w) +
  plot_layout(ncol = 1)
final_plot
ggsave("results/Figures/phyloregions200km_validation.png", final_plot, dpi = 400, width = 10, height = 18)

#####################################################################
# COMPARATIVE ANALYSIS
#####################################################################
# ============================================================
# SECTION 30: SIMPSON VS SORENSEN COMPARISON
# ============================================================
# Compare only the bioregionalizations between Simpson and Sorensen
# V-measure calculation between Simpson and Sorensen
v.measure<-vmeasure_calc(bioregion200, bioregion200_sor, group,group, B = 1, precision = NULL)
v.measure
#The SABRE results:
#V-measure: 0.83 
#Homogeneity: 0.86 
#Completeness: 0.81 
p1 <- ggplot(v.measure$map1) + geom_sf(aes(fill = rih)) +
  scale_fill_viridis_c(option = "B", direction = -1, name = "rih") +  # option C
  ggtitle("Bioregion mean βsim ") +
  theme_bw(base_size = 16,base_family = "sans")

p2 <- ggplot(v.measure$map2) + geom_sf(aes(fill = rih)) +
  scale_fill_viridis_c(option = "B", direction = -1, name = "rih") +  # option B
  ggtitle("Bioregion mean βsor") +
  theme_bw(base_size = 16,base_family = "sans")
final1<-p1 / p2
final1
# Goodness-of-fit between Simpson and Sorensen
# Mapcurves calculation: It calculates the Mapcurves's goodness-of-fit (GOF)
map.curv<-mapcurves_calc(bioregion200, bioregion200_sor, group,group, precision = NULL)
map.curv
#The MapCurves results:
#The goodness of fit: 0.91 
#Reference map: x 

# ============================================================
# SECTION 31: SIMPSON VS WELL-SAMPLED SUBSET COMPARISON
# ============================================================
# V-measure calculation between Simpson and subset (well and moderately sampled cells)
v.measure2<-vmeasure_calc(bioregion200, bioregion200_subset, group,group, B = 1, precision = NULL)
v.measure2
#The SABRE results:
#V-measure: 0.93 
#Homogeneity: 0.92 
#Completeness: 0.93 
p3 <- ggplot(v.measure2$map1) + geom_sf(aes(fill = rih)) +
  scale_fill_viridis_c(option = "B", direction = -1, name = "rih") +  # option C
  ggtitle("Bioregion mean βsim all data ") +
  theme_bw(base_size = 16,base_family = "sans")

p4 <- ggplot(v.measure2$map2) + geom_sf(aes(fill = rih)) +
  scale_fill_viridis_c(option = "B", direction = -1, name = "rih") +  # option B
  ggtitle("Bioregion mean βsim subset") +
  theme_bw(base_size = 16,base_family = "sans")
p3 / p4
final2<-p3 / p4
final2
# Goodness-of-fit between Simpson and subset
# Mapcurves calculation: It calculates the Mapcurves's goodness-of-fit (GOF)
map.curv2<-mapcurves_calc(bioregion200, bioregion200_subset, group,group, precision = NULL)
map.curv2
#The MapCurves results:
#The goodness of fit: 0.98
#Reference map: x 
# ============================================================
# SECTION 32: SIMPSON VS POOR-SAMPLED SUBSET COMPARISON
# ============================================================
# V-measure calculation between Simpson and subset2 (poorly sampled cells)
v.measure3<-vmeasure_calc(bioregion200, bioregion200_subset2, group,group, B = 1, precision = NULL)
v.measure3
#The SABRE results:
#V-measure: 0.82
#Homogeneity: 0.94 
#Completeness: 0.72 
p5 <- ggplot(v.measure3$map1) + geom_sf(aes(fill = rih)) +
  scale_fill_viridis_c(option = "B", direction = -1, name = "rih") +  # option C
  ggtitle("Bioregion mean βsim all data ") +
  theme_bw(base_size = 16,base_family = "sans")

p6 <- ggplot(v.measure3$map2) + geom_sf(aes(fill = rih)) +
  scale_fill_viridis_c(option = "B", direction = -1, name = "rih") +  # option B
  ggtitle("Bioregion mean βsim subset2") +
  theme_bw(base_size = 16,base_family = "sans")
p5 / p6
# Goodness-of-fit between Simpson and subset2
# Mapcurves calculation: It calculates the Mapcurves's goodness-of-fit (GOF)
map.curv3<-mapcurves_calc(bioregion200, bioregion200_subset2, group,group, precision = NULL)
map.curv3
#The MapCurves results:
#The goodness of fit: 0.87 
#Reference map: x