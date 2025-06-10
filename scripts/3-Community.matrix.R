###########################################################################
#                             3.- Community matrix
#in the ‘community matrix’ section you create the species presence-absence
#matrix (PAM), if you intend to use different grain sizes (which would be ideal)
#this step is where you do it, you also generate the hypothesis of the phylogenetic
#tree and use that information to match it with the community matrix. Finally,
#it is decided whether to perform biogeographical phylogenetic regionalization,
#taxonomic regionalization or both.
###########################################################################
#working with the final base ‘orchid.final.tax.stand’ stored in taxonomic standardisation 

#1.- We use only one grain size (200 x 200 km) to show our results,
#although in supplementary material we add three more grain sizes:
#200 x 200 km, 400 x 400 km, and 800 x 800 km.

#2.- We concentrated on phylogenetic biogeographic regionalization.

# Clear workspace
rm(list = ls())
# Load required packages
library(data.table)   # For efficient data handling (fread)
library(dplyr)        # For data manipulation
library(sf)           # For spatial data handling
library(ggplot2)      # For plotting
library(rnaturalearth) # For map data
library(ape)          # For phylogenetic tree handling
library(phyloregion)  # For dense2sparse

posterior_trees<-readRDS("processed-data/community_matrix/phylogenetic_metrics/100_random_trees.rds")#THow these phylogenetic trees have been generated is described below.
orchids<- fread("processed-data/dataselection_taxonomic-standard/orchids_native_range.csv")
#732,359 records
#subset of species and coordinates to create the presence-absence matrix
names(orchids)
occ.orchid <- orchids %>% 
  dplyr::select(wcvp_name,decimalLongitude,decimalLatitude)

#rename variables
names(occ.orchid)
occ.orchid = rename(occ.orchid, c(species="wcvp_name"))
names(occ.orchid)

# Add hyphen to species to match them later with the phylogenetic tree.
species<-gsub(" ", "_", occ.orchid$species)
species<-as.data.frame(species)
occ.orchid <- cbind(occ.orchid[, c(2,3)],species)
#Rearrange columns
occ.orchid = occ.orchid[ , c(3,1,2)]
occ.orchid<-as.data.frame(occ.orchid)
names(occ.orchid)
rm(species,orchids)

#convert our data into a ‘Simple Features’ object, 
#and reproject the coordinates to an equal-area projection= Behrmann projection.
# Verificación exhaustiva
occ.orchid <- st_as_sf(occ.orchid, coords=c("decimalLongitude", "decimalLatitude"), crs='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
behrmann <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
occ.orchid <- st_transform(occ.orchid, behrmann) 

#Now, let's look at the species on our rnaturalearth map,
#and we will also reproject to match the occurrence data.
continents<- ne_countries(continent = c("Africa","Asia","North America",
                                        "Europe","Oceania","South America"),
                          returnclass = "sf", scale = "medium")
str(continents)
map <- st_transform(continents, behrmann)
ggplot(map) +
  geom_sf(fill="grey40", color="grey40") +
  theme_void()

#Matching phylogenetic and occurrence data
tree <- posterior_trees[[1]]#extract one of the 100 trees that were generated
matched.names <- intersect(tree$tip.label, occ.orchid$species)
occ.orchid <- occ.orchid %>% filter(species %in% matched.names)
length(unique(occ.orchid$species))#check

#####Creating a grid######
#In order to make our species-community matrix,
#we need to make a grid to intersect with our occurrence data.
#Although we present our results in 100 x 100 grids. Here we generate 3 more grain sizes, see the supplementary material for the results.
grid.size <- 100000 ##specify grid dimensions (in this case, in meters)
grid.size2 <- 200000 
grid.size4 <- 400000 
grid.size8 <- 800000 
#100km
grid_shape1<- st_make_grid(x=map,what = "polygons", cellsize = grid.size,
                           square=TRUE,crs=behrmann,flat_topped=TRUE)#square=false will generate hexagonal cells
grid <- st_sf(idcell = 1:length(grid_shape1), geom = grid_shape1,
              crs=behrmann) %>% st_cast("POLYGON")#reproject
#200km
grid_shape2<- st_make_grid(x=map,what = "polygons", cellsize = grid.size2,
                           square=TRUE,crs=behrmann,flat_topped=TRUE)#square=false will generate hexagonal cells
grid2 <- st_sf(idcell = 1:length(grid_shape2), geom = grid_shape2,
               crs=behrmann) %>% st_cast("POLYGON")#reproject
#400km
grid_shape4<- st_make_grid(x=map,what = "polygons", cellsize = grid.size4,
                           square=TRUE,crs=behrmann,flat_topped=TRUE)#square=false will generate hexagonal cells
grid4 <- st_sf(idcell = 1:length(grid_shape4), geom = grid_shape4,
               crs=behrmann) %>% st_cast("POLYGON")#reproject
#800km
grid_shape8<- st_make_grid(x=map,what = "polygons", cellsize = grid.size8,
                           square=TRUE,crs=behrmann,flat_topped=TRUE)#square=false will generate hexagonal cells
grid8 <- st_sf(idcell = 1:length(grid_shape8), geom = grid_shape8,
               crs=behrmann) %>% st_cast("POLYGON")#reproject
#save grids
#100km
  st_write(grid,"processed-data/community_matrix/pam_shape/grid100km_species.shp")
plot(grid)
#200km
st_write(grid2,"processed-data/community_matrix/pam_shape/grid200km_species.shp")
plot(grid2)
#400km
st_write(grid4,"processed-data/community_matrix/pam_shape/grid400km_species.shp")
plot(grid4)
#800km
st_write(grid8,"processed-data/community_matrix/pam_shape/grid800km_species.shp")
plot(grid8)

grid.intersect100 <- st_intersection(map,grid)
grid.intersect200 <- st_intersection(map,grid2)
grid.intersect400 <- st_intersection(map,grid4)
grid.intersect800 <- st_intersection(map,grid8)
#Joining data to make the species-community grid.intersect <- st_intersection(map,grid) matrix
#Now that we have a grid, the next step is to join it with our species occurrence data.
#We’ll then convert this to a data frame and retain just the taxon names
#and the corresponding grid cell ID values. This will be the basis 
#for making our species-community matrix.
#100km
intersected.grid100 <- st_intersection(grid,occ.orchid)
intersected.gridDF <- as.data.frame(intersected.grid100)[,c(1,2)]
head(intersected.gridDF)

#200km
intersected.grid200 <- st_intersection(grid2,occ.orchid)
intersected.gridDF2 <- as.data.frame(intersected.grid200)[,c(1,2)]
head(intersected.gridDF2)

#400k
intersected.grid400 <- st_intersection(grid4,occ.orchid)
intersected.gridDF4 <- as.data.frame(intersected.grid400)[,c(1,2)]
head(intersected.gridDF4)

#800k
intersected.grid800 <- st_intersection(grid8,occ.orchid)
intersected.gridDF8 <- as.data.frame(intersected.grid800)[,c(1,2)]
head(intersected.gridDF8)

#species-community matrix 
#100km
grid.species.matrix <- merge.data.frame(x = data.frame(idcell = grid$idcell), 
                                        y = intersected.gridDF, 
                                        all.x = TRUE, by.x = TRUE, 
                                        sort = TRUE) %>% table()
grid.species.matrix[grid.species.matrix > 0] <- 1 #converts from frequence to P/A
#200km
grid.species.matrix2 <- merge.data.frame(x = data.frame(idcell = grid2$idcell), 
                                         y = intersected.gridDF2, 
                                         all.x = TRUE, by.x = TRUE, 
                                         sort = TRUE) %>% table()
grid.species.matrix2[grid.species.matrix2 > 0] <- 1 #converts from frequence to P/A
#400km
grid.species.matrix4 <- merge.data.frame(x = data.frame(idcell = grid4$idcell), 
                                         y = intersected.gridDF4, 
                                         all.x = TRUE, by.x = TRUE, 
                                         sort = TRUE) %>% table()
grid.species.matrix4[grid.species.matrix4 > 0] <- 1 #converts from frequence to P/A
#800km
grid.species.matrix8 <- merge.data.frame(x = data.frame(idcell = grid8$idcell), 
                                         y = intersected.gridDF8, 
                                         all.x = TRUE, by.x = TRUE, 
                                         sort = TRUE) %>% table()
grid.species.matrix8[grid.species.matrix8 > 0] <- 1 #converts from frequence to P/A

#Minimum 5 species per cell
#100km
com_reduce <- 
  which(rowSums(grid.species.matrix) <5) # remove cells with <5 especies
matrix.clean1 <- grid.species.matrix[-com_reduce, ]#save matrix.clean1
dim(matrix.clean1)#5085 grid cells
grid_reduce <- grid[-com_reduce,]
dim(grid_reduce)

shape1 <- st_read(dsn ="processed-data/community_matrix/pam_shape/grid100km_species.shp")
shape100<-shape1[-com_reduce,]
dim(shape100)
plot(shape100)
st_write(shape100, dsn = "processed-data/community_matrix/pam_shape/",layer="shape_reduce100", driver="ESRI Shapefile")
#save
comm100 <- dense2sparse(matrix.clean1)# 5 species
save(comm100, file = "processed-data/community_matrix/pam/pam100_reduce.RData")
write.csv(grid_reduce,file = "processed-data/community_matrix/pam/grid.pam100_reduce.csv")

#200km
com_reduce2 <- 
  which(rowSums(grid.species.matrix2) <5) # remove cells with <5 especies
matrix.clean2 <- grid.species.matrix2[-com_reduce2, ]#save matrix.clean1
dim(matrix.clean2)#2,201 grid cells
grid_reduce2 <- grid2[-com_reduce2,]
dim(grid_reduce2)

shape2 <- st_read(dsn ="processed-data/community_matrix/pam_shape/grid200km_species.shp")
shape200<-shape2[-com_reduce2,]
dim(shape200)
plot(shape200)
st_write(shape200, dsn = "processed-data/community_matrix/pam_shape/",layer="shape_reduce200", driver="ESRI Shapefile")
#save
comm200 <- dense2sparse(matrix.clean2)# 5 species
save(comm200, file = "processed-data/community_matrix/pam/pam200_reduce.RData")
write.csv(grid_reduce2,file = "processed-data/community_matrix/pam/grid.PAM200_reduce.csv")

#400km
com_reduce4 <- 
  which(rowSums(grid.species.matrix4) <5) # remove cells with <5 especies
matrix.clean4 <- grid.species.matrix4[-com_reduce4, ]#save matrix.clean1
dim(matrix.clean4)#837 grid cells
grid_reduce4 <- grid4[-com_reduce4,]
dim(grid_reduce4)

shape4 <- st_read(dsn ="processed-data/community_matrix/pam_shape/grid400km_species.shp")
shape400<-shape4[-com_reduce4,]
dim(shape400)
plot(shape400)
st_write(shape400, dsn = "processed-data/community_matrix/pam_shape/",layer="shape_reduce400", driver="ESRI Shapefile")
#save
comm400 <- dense2sparse(matrix.clean4)# 5 species
save(comm400, file = "processed-data/community_matrix/pam/pam400_reduce.RData")
write.csv(grid_reduce4,file = "processed-data/community_matrix/pam/grid.pam400_reduce.csv")

#800km
com_reduce8 <- 
  which(rowSums(grid.species.matrix8) <5) # remove cells with <5 especies
matrix.clean8 <- grid.species.matrix8[-com_reduce8, ]#save matrix.clean1
dim(matrix.clean8)#310 grid cells
grid_reduce8 <- grid8[-com_reduce8,]
dim(grid_reduce8)

shape800 <- st_read(dsn ="processed-data/community_matrix/pam_shape/grid800km_species.shp")
shape800<-shape800[-com_reduce8,]
dim(shape800)
plot(shape800)
st_write(shape800, dsn = "processed-data/community_matrix/pam_shape/",layer="shape_reduce800", driver="ESRI Shapefile")
#save
comm800 <- dense2sparse(matrix.clean8)# 5 species
save(comm800, file = "processed-data/community_matrix/pam/pam800_reduce.RData")
write.csv(grid_reduce8,file = "processed-data/community_matrix/pam/grid.pam800_reduce.csv")

##########################################################################
#                   Phylogenetic tree hypothesis             #
##########################################################################
# Clear workspace
rm(list = ls())
# Load required packages
library(data.table)     # For efficient reading of large CSV files
library(tidyr)          # For data tidying
library(dplyr)          # For data manipulation
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
load("processed-data/community_matrix/pam/pam100_reduce.RData")
load("processed-data/community_matrix/pam/pam200_reduce.RData")
load("processed-data/community_matrix/pam/pam400_reduce.RData")
load("processed-data/community_matrix/pam/pam800_reduce.RData")
#load posterior trees
posterior_trees<-readRDS("processed-data/community_matrix/phylogenetic_metrics/100_random_trees.rds")

##########################################################################
# Mean phylogenetic beta diversity
##########################################################################
#100km

#1.- Calculate phylogenetic beta diversity for each tree
beta_values100 <- lapply(seq_along(posterior_trees), function(i) {
  cat("Procesando árbol", i, "de 100\n")
  # Verify species match 
  sp_match <- colnames(comm100) %in% posterior_trees[[i]]$tip.label
  if(!all(sp_match)) {
    missing <- sum(!sp_match)
    warning(paste("Árbol", i, ":", missing, "especies no encontradas"))}
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
cat("\n=== RESUMEN DE RESULTADOS ===\n")
cat("Dimensión de las matrices:", dim(beta_sor_mean100), "\n")
cat("Valores promedio:\n")
cat(" - SOR:", mean(beta_sor_mean100, na.rm = TRUE), "\n")
cat(" - SIM:", mean(beta_sim_mean100, na.rm = TRUE), "\n")
cat(" - SNE:", mean(beta_sne_mean100, na.rm = TRUE), "\n")
#Valores promedio:
#- SOR: 0.7554568 
#- SIM: 0.6258285 
#- SNE: 0.1296283 
# 5. saving the results
save(beta_sor_mean100, beta_sim_mean100, beta_sne_mean100,
     file ="processed-data/community_matrix/phylogenetic_metrics//mean_beta_components_100.RData")
##########################################################
##########################################################
#200km
#1.- Calculate phylogenetic beta diversity for each tree
beta_values200 <- lapply(seq_along(posterior_trees), function(i) {
  cat("Procesando árbol", i, "de 100\n")
  # Verify species match 
  sp_match <- colnames(comm200) %in% posterior_trees[[i]]$tip.label
  if(!all(sp_match)) {
    missing <- sum(!sp_match)
    warning(paste("Árbol", i, ":", missing, "especies no encontradas"))}
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
cat("\n=== RESUMEN DE RESULTADOS ===\n")
cat("Dimensión de las matrices:", dim(beta_sor_mean200), "\n")
cat("Valores promedio:\n")
cat(" - SOR:", mean(beta_sor_mean200, na.rm = TRUE), "\n")
cat(" - SIM:", mean(beta_sim_mean200, na.rm = TRUE), "\n")
cat(" - SNE:", mean(beta_sne_mean200, na.rm = TRUE), "\n")
#Valores promedio:
#- SOR: 0.7646349 
#- SIM: 0.616163
#- SNE: 0.1484719 
# 5. saving the results
save(beta_sor_mean200, beta_sim_mean200, beta_sne_mean200,
     file ="processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_200.RData")
#########################################################
##########################################################
#400km
#1.- Calculate phylogenetic beta diversity for each tree
beta_values400 <- lapply(seq_along(posterior_trees), function(i) {
  cat("Procesando árbol", i, "de 100\n")
  # Verify species match 
  sp_match <- colnames(comm400) %in% posterior_trees[[i]]$tip.label
  if(!all(sp_match)) {
    missing <- sum(!sp_match)
    warning(paste("Árbol", i, ":", missing, "especies no encontradas"))}
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
cat("\n=== RESUMEN DE RESULTADOS ===\n")
cat("Dimensión de las matrices:", dim(beta_sor_mean400), "\n")
cat("Valores promedio:\n")
cat(" - SOR:", mean(beta_sor_mean400, na.rm = TRUE), "\n")
cat(" - SIM:", mean(beta_sim_mean400, na.rm = TRUE), "\n")
cat(" - SNE:", mean(beta_sne_mean400, na.rm = TRUE), "\n")
#Valores promedio:
#- SOR: 0.7775248 
#- SIM: 0.6083764
#- SNE: 0.1691484
# 5. saving the results
save(beta_sor_mean400, beta_sim_mean400, beta_sne_mean400,
     file ="processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_400.RData")
#########################################################
##########################################################
#800
#1.- Calculate phylogenetic beta diversity for each tree
beta_values800 <- lapply(seq_along(posterior_trees), function(i) {
  cat("Procesando árbol", i, "de 100\n")
  # Verify species match 
  sp_match <- colnames(comm800) %in% posterior_trees[[i]]$tip.label
  if(!all(sp_match)) {
    missing <- sum(!sp_match)
    warning(paste("Árbol", i, ":", missing, "especies no encontradas"))}
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
cat("\n=== RESUMEN DE RESULTADOS ===\n")
cat("Dimensión de las matrices:", dim(beta_sor_mean800), "\n")
cat("Valores promedio:\n")
cat(" - SOR:", mean(beta_sor_mean800, na.rm = TRUE), "\n")
cat(" - SIM:", mean(beta_sim_mean800, na.rm = TRUE), "\n")
cat(" - SNE:", mean(beta_sne_mean800, na.rm = TRUE), "\n")
#Valores promedio:
#- SOR: 0.7917338 
#- SIM: 0.6036937
#- SNE: 0.1880401
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
cat("✅ Análisis completado exitosamente!\n")
cat("📊 Sitios procesados:", nrow(comm100), "\n")
cat("🌳 Árboles procesados:", length(posterior_trees), "\n")
#Sitios procesados: 5460 
#Árboles procesados: 100 

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
cat("✅ Análisis completado exitosamente!\n")
cat("📊 Sitios procesados:", nrow(comm200), "\n")
cat("🌳 Árboles procesados:", length(posterior_trees), "\n")

#Sitios procesados: 2417 
#Árboles procesados: 100 
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
cat("✅ Análisis completado exitosamente!\n")
cat("📊 Sitios procesados:", nrow(comm400), "\n")
cat("🌳 Árboles procesados:", length(posterior_trees), "\n")
#Sitios procesados: 935 
#Árboles procesados: 100

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
cat("✅ Análisis completado exitosamente!\n")
cat("📊 Sitios procesados:", nrow(comm800), "\n")
cat("🌳 Árboles procesados:", length(posterior_trees), "\n")
#Sitios procesados: 342 
#Árboles procesados: 100

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
  species_per_cell <- rowSums(mat > 0)
  records_per_cell <- rowSums(mat)
  occupancy_per_species <- colSums(mat > 0)
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
    N_doubleton_species = n_doubletons)
}
# Apply to all matrices
results <- lapply(matrices, extract_metrics)
# Combine results
results_df <- bind_rows(results, .id = "Grain_size_km")
results_df$Grain_size_km <- gsub("comm", "", results_df$Grain_size_km)
results_df$Grain_size_km <- paste0(results_df$Grain_size_km, " km")
# View final table
results_df
# Export to CSV
write.csv(results_df, "results/table/_metrics_by_grain.csv", row.names = FALSE)
