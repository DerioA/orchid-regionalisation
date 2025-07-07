###########################################################################
#                             # 3. Community Matrix
# In this section, we create the species presence-absence matrix (PAM). If using different grain sizes (ideal approach), # this is where the process is implemented.
# Our methodology:
# - Base sampling unit: 10×10 km grid cells
# - Aggregated scales:
#   1. 100×100 km
#   2. 200×200 km
#   3. 400×400 km
#   4. 800×800 km
#
# Rationale for multi-scale approach:
# Enables comparative analysis of bioregionalisation patterns across spatial scales.
# Additional steps:
# - Generate phylogenetic tree hypothesis
# - Cross-reference phylogenetic data with community matrices at each grain size
# Regionalisation strategy decision:
# Options include phylogenetic regionalisation, taxonomic regionalisation, or both.
# Our study implements: PHYLOGENETIC REGIONALISATION
###########################################################################
#working with the final base ‘orchid.final.tax.stand’ stored in taxonomic standardisation 
# Reproducible workflow to generate multiscale PAMs and visualisation
rm(list = ls())

# Packages
library(sf)
library(terra)
library(dplyr)
library(Matrix)
library(rnaturalearth)
library(ggplot2)
library(patchwork)
library(data.table)
library(purrr)

# Berhmann projection
behrmann <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# Occurrence data
occs <- fread("processed-data/dataselection_taxonomic-standard/orchids_native_range.csv") %>% 
  dplyr::select(wcvp_name, decimalLongitude, decimalLatitude) %>%
  rename(species = wcvp_name)# 732359 records

# Preparing spatial data
occs_sf <- st_as_sf(occs, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
  st_transform(crs = behrmann)

# CRS Verification
cat("CRS of occurrences:", st_crs(occs_sf)$input, "\n")

# Only cells for land surface
land <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(crs = behrmann)

# CRS verification and outreach
cat("Land CRS:", st_crs(land)$input, "\n")
cat("Extent of land:", st_bbox(land), "\n")

# 10 km baseline grid generation for sampling
res_m <- 10 * 1000  # 10 km in metres (Behrmann uses metres)
grid_10 <- st_make_grid(
  land,
  cellsize = c(res_m, res_m),
  square = TRUE,
  crs = behrmann) %>% 
  st_as_sf() %>%
  mutate(idcell = paste0("cell_", seq_len(nrow(.)))) %>%
  st_intersection(land)

# Basic verification
cat("Total number of cells generated:", nrow(grid_10), "\n")

# Assignment of records to a baseline grid
cat("Assigning records to a 10 km baseline grid...\n")
occs_10 <- st_join(occs_sf, grid_10["idcell"], left = FALSE)

# Verification
cat("Records assigned to base cells:", nrow(occs_10), "de", nrow(occs_sf), "\n")
if(nrow(occs_10) == 0) stop("No records assigned to base cells")

# Scales and results
scales <- c(100, 200, 400, 800)
all_results <- list()

# F10 km filter (minimum area)
grid_10_filtered <- grid_10 %>% 
  filter(as.numeric(st_area(.)) >= 0.5 * (10*1000)^2)

cat("=== Summary of initial filtering ===\n")
cat("Original 10 km cellss:", nrow(grid_10), "\n")
cat("Cells 10 km after filtering (<50% area):", nrow(grid_10_filtered), "\n\n")
#68,899 deleted grid cells
for(sc in scales) {
  cat("\n=== Processing scale:", sc, "km ===\n")
  res_m <- sc * 1000
  grid_sc <- st_make_grid(land, cellsize = c(res_m, res_m), square = TRUE, crs = behrmann) %>% 
    st_as_sf() %>%
    mutate(idcell_large = paste0("cell_", sc, "km_", seq_len(nrow(.))))
  
  cat("Large cells generated:", nrow(grid_sc), "\n")
  
  grid_10_centroids <- grid_10_filtered %>% select(idcell) %>% st_centroid()
  grid_assignment <- st_join(grid_10_centroids, grid_sc["idcell_large"], join = st_within, left = FALSE) %>%
    st_drop_geometry() %>% distinct()
  
  occs_joined <- occs_10 %>%
    st_drop_geometry() %>%
    inner_join(grid_assignment, by = "idcell") %>%
    select(idcell_large, species) %>%
    distinct()
  
  if(nrow(occs_joined) > 0) {
    total_species <- length(unique(occs_joined$species))
    cat("Unique species detected:", total_species, "\n")
    
    occs_joined$idcell_large <- factor(occs_joined$idcell_large)
    occs_joined$species <- factor(occs_joined$species)
    
    pam_sc <- Matrix::sparseMatrix(
      i = as.integer(occs_joined$idcell_large),
      j = as.integer(occs_joined$species),
      x = 1,
      dims = c(
        length(levels(occs_joined$idcell_large)),
        length(levels(occs_joined$species))),
      dimnames = list(levels(occs_joined$idcell_large), levels(occs_joined$species)))
    
    richness_sc <- Matrix::rowSums(pam_sc)
    keep_sc <- which(richness_sc >= 5)
    pam_sc_filtered <- pam_sc[keep_sc, , drop = FALSE]
    grid_sc_filtered <- grid_sc %>% filter(idcell_large %in% rownames(pam_sc_filtered))
    
    cat("Cells with species richness >=5 species:", length(keep_sc), "\n")
    cat("Percentage of cells retained:", round(length(keep_sc)/nrow(grid_sc)*100, 1), "%\n")
    
    all_results[[paste0(sc, "km")]] <- list(
      pam = pam_sc_filtered,
      grid = grid_sc_filtered,
      n_cells_total = nrow(grid_sc),
      n_cells_filtered = length(keep_sc),
      n_species = total_species,
      species_list = unique(occs_joined$species))
    
    print(
      ggplot() +
        geom_sf(data = land, fill = "gray90", color = NA) +
        geom_sf(data = grid_sc, fill = NA, color = "gray80", size = 0.1) +
        geom_sf(data = grid_sc_filtered, fill = "forestgreen", alpha = 0.5, color = NA) +
        labs(title = paste("Grid", sc, "km"),
             subtitle = paste("Cells with ≥5 species:", length(keep_sc), "de", nrow(grid_sc)),
             caption = paste("Unique species:", total_species)) +
        theme_void())
    
  } else {
    warning("There are no records assigned for the scale", sc, " km")
    all_results[[paste0(sc, "km")]] <- list(
      pam = Matrix(0,0,0), 
      grid = grid_sc[0,],
      n_cells_total = nrow(grid_sc),
      n_cells_filtered = 0,
      n_species = 0)
  }
}
#=== Processing scale: 100 km ====
# Large cells generated: 51156
#Unique species detected: 19123
#Cells with richness >=5 species: 5003
#Percentage of cells conserved: 9.8%.

#=== Processing scale: 200 km ====
# Large cells generated: 12876
#Unique species detected: 19123
#Cells with richness >=5 species: 2200
#Percentage of cells conserved: 17.1 %

#=== Processing scale: 400 km ===
# Large cells generated: 3219
#Unique species detected: 19123
#Cells with richness >=5 species: 850
#Percentage of cells conserved: 26.4%.

#=== Processing scale: 800 km ====
# Large cells generated: 836
#Unique species detected: 19123
#Cells with richness >=5 species: 316
#Percentage of cells conserved: 37.8% 

# Final report
cat("\n=== Summary ===\n")
summary_df <- map_dfr(names(all_results), function(nm) {
  data.frame(
    Escala = nm,
    Celdas_totales = all_results[[nm]]$n_cells_total,
    Celdas_filtradas = all_results[[nm]]$n_cells_filtered,
    Especies = all_results[[nm]]$n_species,
    Porcentaje_conservado = round(all_results[[nm]]$n_cells_filtered/all_results[[nm]]$n_cells_total*100, 1))
})
summary_df

# Save results
cat("\nSaving results...\n")
dir.create("processed-data/community_matrix/pam", recursive = TRUE, showWarnings = FALSE)
dir.create("processed-data/community_matrix/pam_shape", recursive = TRUE, showWarnings = FALSE)

summary_data <- data.frame(scale = character(), n_cells = integer(), n_species = integer(), stringsAsFactors = FALSE)

for(scale_name in names(all_results)) {
  res <- all_results[[scale_name]]
  if(nrow(res$grid) > 0 && nrow(res$pam) > 0) {
    saveRDS(res$pam, file = paste0("processed-data/community_matrix/pam/pam_", scale_name, ".rds"))
    grid_to_save <- res$grid %>% rename(idcell = idcell_large) %>% select(idcell)
    st_write(grid_to_save, paste0("processed-data/community_matrix/pam_shape/grid_", scale_name, ".gpkg"), delete_dsn = TRUE)
    summary_data <- rbind(summary_data, data.frame(scale = scale_name, n_cells = res$n_cells_filtered, n_species = res$n_species))
  } else {
    cat("Escala", scale_name, "has no valid data to save.\n")
  }
}

write.csv(summary_data, "processed-data/community_matrix/pam/summary_pams_grids.csv", row.names = FALSE)
cat("\nSummary of PAMs generateds:\n")
print(summary_data)

# Final validation
cat("\nFinal validation...\n")
pam_400 <- readRDS("processed-data/community_matrix/pam/pam_400km.rds")
grid_400 <- st_read("processed-data/community_matrix/pam_shape/grid_400km.gpkg", quiet = TRUE)
cat("First 5 IDs in PAM 400 km:", head(rownames(pam_400), 5), "\n")
cat("First 5 IDs in grid 400 km:", head(grid_400$idcell, 5), "\n")
matching_cells <- sum(rownames(pam_400) %in% grid_400$idcell)
cat("Overlapping cells between PAM 400 km and grid:", matching_cells, "de", nrow(pam_400), "\n")
if(matching_cells != nrow(pam_400)) warning("¡There is inconsistency in the IDs between PAM and grid.!")

##########################################################################
#                   Phylogenetic tree hypothesis             #
##########################################################################
# Clear workspace
rm(list = ls())
# Load required packages
library(data.table)     # For efficient reading of large CSV files
library(tidyr)          # For data tidying
library(dplyr) 
library(devtools)
#devtools::install_github("jinyizju/V.PhyloMaker2")# For data manipulation
library(V.PhyloMaker2)  # For phylogenetic tree construction
library(parallel)       # For parallel computing

# Load final orchid database
orchids <- fread("processed-data/dataselection_taxonomic-standard/orchids_native_range.csv",
                 na.strings = "")  # 732,359 records

# Extract list of unique species from the database
species.list <- orchids %>%
  dplyr::select(wcvp_name, wcvp_genus, family) %>%
  distinct() %>%
  as.data.frame()

# Rename columns to match expected input format
species.list <- species.list %>%
  rename(
    species = wcvp_name,
    genus = wcvp_genus)

#load community matrix
comm100 <- readRDS("processed-data/community_matrix/pam/pam_100km.rds")

# Check if there are differences in species between species.list and comm100.
# It must be true because cells with less than 50% of the area were eliminated in 10 x 10 km closures.
# Total species
cat("Total number of species in species.list:", nrow(species.list), "\n")
cat("Total number of species in comm100:", ncol(comm100), "\n")
# Species in each set
species_in_list <- species.list$species
species_in_comm <- colnames(comm100)
# Check differences
missing_in_comm <- setdiff(species_in_list, species_in_comm)
missing_in_list <- setdiff(species_in_comm, species_in_list)
cat("Species in species.list but not in comm100:", length(missing_in_comm), "\n")
cat("Species in comm100 but not in species.list:", length(missing_in_list), "\n")

# Filter species.list to match comm100, comm200, comm400 and comm800
species.list <- species.list %>%
  filter(species %in% species_in_comm)
cat("Species list filter:", nrow(species.list), "species\n")

# Replace spaces with underscores in species names (for matching with tree tips)
species <- gsub(" ", "_", species.list$species)
species <- as.data.frame(species)

# Combine and reorder columns: species, genus, family
species.list <- cbind(species.list[, c("genus", "family")], species)
species.list <- species.list[, c("species", "genus", "family")]
# Remove temporary objects
rm(species, orchids)

### Construct 100 phylogenetic trees using V.PhyloMaker2 in parallel ###
# Set up a parallel cluster
cl <- makeCluster(detectCores() - 1)

# Export necessary objects to the cluster
clusterExport(cl, c("species.list", "GBOTB.extended.WP", "nodes.info.1.WP"))

# Define function to generate a single phylogenetic tree
generate_named_tree <- function(i) {
  tree <- V.PhyloMaker2::phylo.maker(
    species.list,
    tree = GBOTB.extended.WP,
    nodes = nodes.info.1.WP,
    scenarios = "S2",
    r = 1
  )$scenario.2
  return(tree)}

# Generate 100 phylogenetic trees in parallel
trees_list <- parLapply(cl, 1:100, generate_named_tree)
# Approximately 7-8 hours on an 8-core Mac with 16 MB RAM.

# Assign names to each tree
names(trees_list) <- paste0("tree", 1:100)

# Stop the cluster
stopCluster(cl)

# Save the complete list of trees
saveRDS(trees_list, "processed-data/community_matrix/100_random_trees.rds")

##########################################################################
# Calculate average phylogenetic diversity and phylogenetic beta diversity
# from 100 posterior phylogenetic trees
##########################################################################
# Clean environment
rm(list = ls())
# Load required libraries
library(phyloregion)    # For phylogenetic diversity and beta diversity calculations
library(phytools)       # For phylogenetic tree manipulation
library(ape)            # For phylogenetic tree manipulation
library(vegan)          # For ecological community analysis
library(data.table)     # For efficient data handling
library(parallel)       # For parallel processing
library(Matrix)         # For handling sparse matrices

#load community matrix
comm100 <- readRDS("processed-data/community_matrix/pam/pam_100km.rds")
comm200 <- readRDS("processed-data/community_matrix/pam/pam_200km.rds")
comm400 <- readRDS("processed-data/community_matrix/pam/pam_400km.rds")
comm800 <- readRDS("processed-data/community_matrix/pam/pam_800km.rds")
# Corregir nombres de especies en todas las PAMs
colnames(comm100) <- gsub(" ", "_", colnames(comm100))# Add underscore between species nobmres
colnames(comm200) <- gsub(" ", "_", colnames(comm200))
colnames(comm400) <- gsub(" ", "_", colnames(comm400))
colnames(comm800) <- gsub(" ", "_", colnames(comm800))

# Check change in comm100
head(colnames(comm100))

#load posterior trees
posterior_trees<-readRDS("processed-data/community_matrix/phylogenetic_metrics/100_random_trees.rds")

##########################################################################
# Mean phylogenetic beta diversity
##########################################################################
#100km

#1.- Calculate phylogenetic beta diversity for each tree
beta_values100 <- lapply(seq_along(posterior_trees), function(i) {
  cat("Processing tree", i, "de 100\n")
  # Verify species match 
  sp_match <- colnames(comm100) %in% posterior_trees[[i]]$tip.label
  if(!all(sp_match)) {
    missing <- sum(!sp_match)
    warning(paste("tree", i, ":", missing, "species not found"))}
  # Calculate beta diversity only with matched species
  phyloregion::phylobeta(comm100[, sp_match], posterior_trees[[i]], index = "sorensen")})
# 2. extraction and calculation of averages 
components <- c("phylo.beta.sor", "phylo.beta.sim", "phylo.beta.sne")
mean_results <- lapply(components, function(comp) {
  Reduce("+", lapply(beta_values100, `[[`, comp)) / length(beta_values100)})
names(mean_results) <- c("beta_sor", "beta_sim", "beta_sne")
# 3. Name assignment
beta_sor_mean100 <- mean_results$beta_sor
beta_sim_mean100 <- mean_results$beta_sim
beta_sne_mean100 <- mean_results$beta_sne
# 4. Verification of results
cat("\n=== SUMMARY OF RESULTS ===\n")
cat("Dimension of the matrices:", dim(beta_sor_mean100), "\n")
cat("Average values:\n")
cat(" - SOR:", mean(beta_sor_mean100, na.rm = TRUE), "\n")
cat(" - SIM:", mean(beta_sim_mean100, na.rm = TRUE), "\n")
cat(" - SNE:", mean(beta_sne_mean100, na.rm = TRUE), "\n")
#Average values:
#- SOR: 0.755243 
#- SIM: 0.6226606
#- SNE: 0.1325824

# 5. saving the results
save(beta_sor_mean100, beta_sim_mean100, beta_sne_mean100,
     file ="processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_100.RData")
##########################################################
##########################################################
#200km
#1.- Calculate phylogenetic beta diversity for each tree
beta_values200 <- lapply(seq_along(posterior_trees), function(i) {
  cat("Processing tree", i, "de 100\n")
  # Verify species match 
  sp_match <- colnames(comm200) %in% posterior_trees[[i]]$tip.label
  if(!all(sp_match)) {
    missing <- sum(!sp_match)
    warning(paste("tree", i, ":", missing, "species not found"))}
  # Calculate beta diversity only with matched species
  phyloregion::phylobeta(comm200[, sp_match], posterior_trees[[i]], index = "sorensen")})
# 2. extraction and calculation of averages 
components <- c("phylo.beta.sor", "phylo.beta.sim", "phylo.beta.sne")
mean_results <- lapply(components, function(comp) {
  Reduce("+", lapply(beta_values200, `[[`, comp)) / length(beta_values200)})
names(mean_results) <- c("beta_sor", "beta_sim", "beta_sne")
# 3. Name assignment
beta_sor_mean200 <- mean_results$beta_sor
beta_sim_mean200 <- mean_results$beta_sim
beta_sne_mean200 <- mean_results$beta_sne
# 4. Verification of results
cat("\n=== SUMMARY OF RESULTS ===\n")
cat("Dimension of the matrices:", dim(beta_sor_mean200), "\n")
cat("Average values:\n")
cat(" - SOR:", mean(beta_sor_mean200, na.rm = TRUE), "\n")
cat(" - SIM:", mean(beta_sim_mean200, na.rm = TRUE), "\n")
cat(" - SNE:", mean(beta_sne_mean200, na.rm = TRUE), "\n")
#Average values:
#- SOR: 0.7649638  
#- SIM: 0.6141359
#- SNE: 0.1508278 

# 5. saving the results
save(beta_sor_mean200, beta_sim_mean200, beta_sne_mean200,
     file ="processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_200.RData")
#########################################################
##########################################################
#400km
#1.- Calculate phylogenetic beta diversity for each tree
beta_values400 <- lapply(seq_along(posterior_trees), function(i) {
  cat("Processing tree", i, "de 100\n")
  # Verify species match 
  sp_match <- colnames(comm400) %in% posterior_trees[[i]]$tip.label
  if(!all(sp_match)) {
    missing <- sum(!sp_match)
    warning(paste("tree", i, ":", missing, "species not found"))}
  # Calculate beta diversity only with matched species
  phyloregion::phylobeta(comm400[, sp_match], posterior_trees[[i]], index = "sorensen")})
# 2. extraction and calculation of averages 
components <- c("phylo.beta.sor", "phylo.beta.sim", "phylo.beta.sne")
mean_results <- lapply(components, function(comp) {
  Reduce("+", lapply(beta_values400, `[[`, comp)) / length(beta_values400)})
names(mean_results) <- c("beta_sor", "beta_sim", "beta_sne")
# 3. Name assignment
beta_sor_mean400 <- mean_results$beta_sor
beta_sim_mean400 <- mean_results$beta_sim
beta_sne_mean400 <- mean_results$beta_sne
# 4. Verification of results
cat("\n=== SUMMARY OF RESULTS ===\n")
cat("Dimension of the matrices:", dim(beta_sor_mean400), "\n")
cat("Average values:\n")
cat(" - SOR:", mean(beta_sor_mean400, na.rm = TRUE), "\n")
cat(" - SIM:", mean(beta_sim_mean400, na.rm = TRUE), "\n")
cat(" - SNE:", mean(beta_sne_mean400, na.rm = TRUE), "\n")
# Average values:
#- SOR: 0.7774449 
#- SIM: 0.6044446
#- SNE: 0.1730003

# 5. saving the results
save(beta_sor_mean400, beta_sim_mean400, beta_sne_mean400,
     file ="processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_400.RData")
#########################################################
##########################################################
#800
#1.- Calculate phylogenetic beta diversity for each tree
beta_values800 <- lapply(seq_along(posterior_trees), function(i) {
  cat("Processing tree", i, "de 100\n")
  # Verify species match 
  sp_match <- colnames(comm800) %in% posterior_trees[[i]]$tip.label
  if(!all(sp_match)) {
    missing <- sum(!sp_match)
    warning(paste("tree", i, ":", missing, "species not found"))}
  # Calculate beta diversity only with matched species
  phyloregion::phylobeta(comm800[, sp_match], posterior_trees[[i]], index = "sorensen")})
# 2. extraction and calculation of averages 
components <- c("phylo.beta.sor", "phylo.beta.sim", "phylo.beta.sne")
mean_results <- lapply(components, function(comp) {
  Reduce("+", lapply(beta_values800, `[[`, comp)) / length(beta_values800)})
names(mean_results) <- c("beta_sor", "beta_sim", "beta_sne")
# 3. Name assignment
beta_sor_mean800 <- mean_results$beta_sor
beta_sim_mean800 <- mean_results$beta_sim
beta_sne_mean800 <- mean_results$beta_sne
# 4. Verification of results
cat("\n=== SUMMARY OF RESULTS ===\n")
cat("Dimension of the matrices:", dim(beta_sor_mean800), "\n")
cat("Average values:\n")
cat(" - SOR:", mean(beta_sor_mean800, na.rm = TRUE), "\n")
cat(" - SIM:", mean(beta_sim_mean800, na.rm = TRUE), "\n")
cat(" - SNE:", mean(beta_sne_mean800, na.rm = TRUE), "\n")
# Average values:
#- SOR: 0.7938941 
#- SIM: 0.6011444
#- SNE: 0.1927497 

# 5. saving the results
save(beta_sor_mean800, beta_sim_mean800, beta_sne_mean800,
     file ="processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_800.RData")

##########################################################################
# Mean phylogenetic diversity
##########################################################################
#100km
library(phyloregion)
library(Matrix)
### 1. Calculate PD ----
PD_results100 <- vapply(posterior_trees, function(tree) {
  pd <- phyloregion::PD(comm100, tree)
  # Ensure order matches comm100
  pd[match(rownames(comm100), names(pd))]
}, FUN.VALUE = numeric(nrow(comm100)))
# Assign names
rownames(PD_results100) <- rownames(comm100)
colnames(PD_results100) <- paste0("tree", seq_along(posterior_trees))
### 2. Statistical calculation ----
PD_summary100 <- data.frame(
  idcell = rownames(comm100),
  PD_mean = rowMeans(PD_results100),
  PD_sd = apply(PD_results100, 1, sd),
  PD_min = apply(PD_results100, 1, min),
  PD_max = apply(PD_results100, 1, max),
  Species_Richness = Matrix::rowSums(comm100 > 0),
  stringsAsFactors = FALSE)

# 3. saving the results
save(PD_summary100,
     file ="processed-data/community_matrix/phylogenetic_metrics/PD_site_means_100.RData")
# Full PD matrix (binary efficient format)
save(PD_results100, file ="processed-data/community_matrix/phylogenetic_metrics/PD_all_trees100.RData")
### 4. PD-richness relationship graph ----
plot(PD_summary100$Species_Richness, PD_summary100$PD_mean,
     xlab = "Species number", ylab = "Mean PD", pch = 19, col = "gray20")
cat("✅ Analysis successfully completed!\n")
cat("📊 Sites processed:", nrow(comm100), "\n")
cat("🌳 Processed trees:", length(posterior_trees), "\n")
# Sites processed: 5003 
# Processed trees: 100 

#########################################################
##########################################################
#200km
### 1. Calculate PD ----
PD_results200 <- vapply(posterior_trees, function(tree) {
  pd <- phyloregion::PD(comm200, tree)
  # Ensure order matches comm100
  pd[match(rownames(comm200), names(pd))]
}, FUN.VALUE = numeric(nrow(comm200)))
# Assign names
rownames(PD_results200) <- rownames(comm200)
colnames(PD_results200) <- paste0("tree", seq_along(posterior_trees))
### 2. Statistical calculation ----
PD_summary200 <- data.frame(
  idcell = rownames(comm200),
  PD_mean = rowMeans(PD_results200),
  PD_sd = apply(PD_results200, 1, sd),
  PD_min = apply(PD_results200, 1, min),
  PD_max = apply(PD_results200, 1, max),
  Species_Richness = Matrix::rowSums(comm200 > 0),
  stringsAsFactors = FALSE)

# 3. saving the results
save(PD_summary200,
     file ="processed-data/community_matrix/phylogenetic_metrics/PD_site_means_200.RData")
# Full PD matrix (binary efficient format)
save(PD_results200, file ="processed-data/community_matrix/phylogenetic_metrics/PD_all_trees200.RData")
### 4. PD-richness relationship graph ----
plot(PD_summary200$Species_Richness, PD_summary200$PD_mean,
     xlab = "Species number", ylab = "Mean PD", pch = 19, col = "gray20")
cat("✅ Analysis successfully completed!\n")
cat("📊 Sites processed:", nrow(comm200), "\n")
cat("🌳 Processed trees:", length(posterior_trees), "\n")
# Sites processed: 2200 
# Processed trees: 100 

#########################################################
#########################################################
#400km
### 1. Calculate PD ----
PD_results400 <- vapply(posterior_trees, function(tree) {
  pd <- phyloregion::PD(comm400, tree)
  # Ensure order matches comm100
  pd[match(rownames(comm400), names(pd))]
}, FUN.VALUE = numeric(nrow(comm400)))
# Assign names
rownames(PD_results400) <- rownames(comm400)
colnames(PD_results400) <- paste0("tree", seq_along(posterior_trees))
### 2. Statistical calculation ----
PD_summary400 <- data.frame(
  idcell = rownames(comm400),
  PD_mean = rowMeans(PD_results400),
  PD_sd = apply(PD_results400, 1, sd),
  PD_min = apply(PD_results400, 1, min),
  PD_max = apply(PD_results400, 1, max),
  Species_Richness = Matrix::rowSums(comm400 > 0),
  stringsAsFactors = FALSE)
# 3. saving the results
save(PD_summary400,
     file ="processed-data/community_matrix/phylogenetic_metrics/PD_site_means_400.RData")
# Full PD matrix (binary efficient format)
save(PD_results400, file ="processed-data/community_matrix/phylogenetic_metrics/PD_all_trees400.RData")
### 4. PD-richness relationship graph ----
plot(PD_summary400$Species_Richness, PD_summary400$PD_mean,
     xlab = "Species number", ylab = "Mean PD", pch = 19, col = "gray20")
cat("✅ Analysis successfully completed!n")
cat("📊 Sites processed:", nrow(comm400), "\n")
cat("🌳 Processed trees:", length(posterior_trees), "\n")
# Sites processed: 850 
# processed trees: 100

#########################################################
#########################################################
#800km
### 1. Calculate PD ----
PD_results800 <- vapply(posterior_trees, function(tree) {
  pd <- phyloregion::PD(comm800, tree)
  # Ensure order matches comm100
  pd[match(rownames(comm800), names(pd))]
}, FUN.VALUE = numeric(nrow(comm800)))
# Assign names
rownames(PD_results800) <- rownames(comm800)
colnames(PD_results800) <- paste0("tree", seq_along(posterior_trees))
### 2. Statistical calculation ----
PD_summary800 <- data.frame(
  idcell = rownames(comm800),
  PD_mean = rowMeans(PD_results800),
  PD_sd = apply(PD_results800, 1, sd),
  PD_min = apply(PD_results800, 1, min),
  PD_max = apply(PD_results800, 1, max),
  Species_Richness = Matrix::rowSums(comm800 > 0),
  stringsAsFactors = FALSE)
# 3. saving the results
save(PD_summary800,
     file ="processed-data/community_matrix/phylogenetic_metrics/PD_site_means_800.RData")
# Full PD matrix (binary efficient format)
save(PD_results800, file ="processed-data/community_matrix/phylogenetic_metrics/PD_all_trees800.RData")
### 4. PD-richness relationship graph ----
plot(PD_summary800$Species_Richness, PD_summary800$PD_mean,
     xlab = "Species number", ylab = "Mean PD", pch = 19, col = "gray20")
cat("✅ Analysis successfully completed!n")
cat("📊 Sites processed:", nrow(comm800), "\n")
cat("🌳 Processed trees:", length(posterior_trees), "\n")
# Sites processed: 316 
# Processed trees: 100

# ==============================================================================
#matriz Statistics
# ==============================================================================
# List of community matrices at different grain sizes
matrices <- list(
  comm100 = comm100,
  comm200 = comm200,
  comm400 = comm400,
  comm800 = comm800)

# Function to extract metrics for each matrix
extract_metrics <- function(mat) {
  n_cells <- nrow(mat)
  n_species <- ncol(mat)
  species_per_cell <- Matrix::rowSums(mat > 0)
  records_per_cell <- Matrix::rowSums(mat)
  occupancy_per_species <- Matrix::colSums(mat > 0)
  n_singletons <- sum(occupancy_per_species == 1)
  n_doubletons <- sum(occupancy_per_species == 2)
  data.frame(
    N_cells = n_cells,
    Gamma_diversity = n_species,
    Mean_alpha_diversity = mean(species_per_cell),
    SD_alpha_diversity = sd(species_per_cell),
    Min_alpha_diversity = min(species_per_cell),
    Max_alpha_diversity = max(species_per_cell),
    Mean_records_per_cell = mean(records_per_cell),
    SD_records_per_cell = sd(records_per_cell),
    N_singleton_species = n_singletons,
    N_doubleton_species = n_doubletons
  )
}

# Apply to all matrices
results <- lapply(matrices, extract_metrics)
# Combine results
results_df <- bind_rows(results, .id = "Grain_size_km")
results_df$Grain_size_km <- gsub("comm", "", results_df$Grain_size_km)
results_df$Grain_size_km <- paste0(results_df$Grain_size_km, " km")
# View final table
results_df

# Adding average phylogenetic diversity at each grain size
pd_stats <- data.frame(
  Grain_size_km = c("100 km", "200 km", "400 km", "800 km"),
  Mean_PD = c(
    mean(PD_summary100$PD_mean, na.rm = TRUE),
    mean(PD_summary200$PD_mean, na.rm = TRUE),
    mean(PD_summary400$PD_mean, na.rm = TRUE),
    mean(PD_summary800$PD_mean, na.rm = TRUE)
  ),
  SD_PD = c(
    sd(PD_summary100$PD_mean, na.rm = TRUE),
    sd(PD_summary200$PD_mean, na.rm = TRUE),
    sd(PD_summary400$PD_mean, na.rm = TRUE),
    sd(PD_summary800$PD_mean, na.rm = TRUE)
  )
)

results_df <- dplyr::left_join(results_df, pd_stats, by = "Grain_size_km")

# Export to CSV
write.csv(results_df, "results/table/_metrics_by_grain.csv", row.names = FALSE)
