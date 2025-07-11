# Assessing the distribution of beta diversity (beta_sim) by scale
# See if the dissimilarities between cells actually increase with grain size,
#which could justify more regions.
#############################################################
rm(list = ls())  # Clear environment
library(tidyverse)
library(data.table) 
library(dplyr)
library(ggplot2)
#############################################################
#load community matrix
load(file = "processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_100.RData")
load(file = "processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_200.RData")
load(file = "processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_400.RData")
load(file = "processed-data/community_matrix/phylogenetic_metrics/mean_beta_components_800.RData")

# Convert `dist` objects to data frames
dist_to_df <- function(d, escala) {
  df <- as.data.frame(as.table(as.matrix(d)))
  df <- df %>%
    filter(Var1 != Var2) %>%  # Eliminar diagonales
    rename(cell1 = Var1, cell2 = Var2, beta_sim = Freq) %>%
    mutate(scale = escala)
  return(df)
}

# Aplicar para cada escala
df_100 <- dist_to_df(beta_sim_mean100, "100 km")
df_200 <- dist_to_df(beta_sim_mean200, "200 km")
df_400 <- dist_to_df(beta_sim_mean400, "400 km")
df_800 <- dist_to_df(beta_sim_mean800, "800 km")

# Unir todo
beta_all_scales <- bind_rows(df_100, df_200, df_400, df_800)

# Función para resumir valores de beta_sim
resumen_beta <- function(df, escala) {
  df %>%
    summarise(
      scale = escala,
      n_pairs = n(),
      mean = mean(beta_sim),
      median = median(beta_sim),
      sd = sd(beta_sim),
      max = max(beta_sim),
      min = min(beta_sim))
}

# Aplicar para cada escala
summary_100 <- resumen_beta(df_100, "100 km")
summary_200 <- resumen_beta(df_200, "200 km")
summary_400 <- resumen_beta(df_400, "400 km")
summary_800 <- resumen_beta(df_800, "800 km")

# Unir todo
beta_summary <- bind_rows(summary_100, summary_200, summary_400, summary_800)
# Mostrar tabla
beta_summary

#############################################################
# ver si los clusters están más separados entre sí a gran escala
# ----------------------------
# Función generalizada
# ----------------------------
calc_between_cluster_beta <- function(beta_dist, clusters_named) {
  beta_mat <- as.matrix(beta_dist)
  cell_ids <- rownames(beta_mat)
  # Asegurar que los clusters estén en el mismo orden que las filas/columnas
  cluster_vec <- clusters_named[cell_ids]
  # Lista de clusters únicos
  cluster_ids <- sort(unique(cluster_vec))
  # Crear tabla con pares de clústeres y calcular disimilitud promedio
  between_df <- expand.grid(cluster1 = cluster_ids, cluster2 = cluster_ids) %>%
    filter(cluster1 < cluster2) %>%
    rowwise() %>%
    mutate(
      beta_mean = mean(beta_mat[
        names(cluster_vec[cluster_vec == cluster1]),
        names(cluster_vec[cluster_vec == cluster2])
      ], na.rm = TRUE)) %>%
    ungroup()
  return(between_df)
}

# Aplica la función a cada escala
between_100 <- calc_between_cluster_beta(beta_sim_mean100, clusters)
between_200 <- calc_between_cluster_beta(beta_sim_mean200, clusters200)
between_400 <- calc_between_cluster_beta(beta_sim_mean400, clusters400)
between_800 <- calc_between_cluster_beta(beta_sim_mean800, clusters800)

# Agrega la escala como columna
between_100$scale <- "100 km"
between_200$scale <- "200 km"
between_400$scale <- "400 km"
between_800$scale <- "800 km"

# Junta todo
between_all <- bind_rows(between_100, between_200, between_400, between_800)

# Función para resumir resultados de between_df
resumen_between <- function(df, escala) {
  df %>%
    summarise(
      scale = escala,
      n_pairs = n(),
      mean = mean(beta_mean, na.rm = TRUE),
      median = median(beta_mean, na.rm = TRUE),
      sd = sd(beta_mean, na.rm = TRUE),
      max = max(beta_mean, na.rm = TRUE),
      min = min(beta_mean, na.rm = TRUE))
}

# Aplicar a cada escala
summary_100 <- resumen_between(between_100, "100 km")
summary_200 <- resumen_between(between_200, "200 km")
summary_400 <- resumen_between(between_400, "400 km")
summary_800 <- resumen_between(between_800, "800 km")

# Unir todo en una tabla
between_summary <- bind_rows(summary_100, summary_200, summary_400, summary_800)

# Ver tabla
print(between_summary, n = Inf)

#############################################################
# Calcular Silhouette Score por tamaño de grano
#Evaluar si los clusters en cada escala están bien definidos internamente
#(alta similitud dentro del cluster) y bien separados externamente
#(alta disimilitud con otros clusters). Un valor alto de Silhouette apoyaría
#que el aumento de regiones a grano grueso es legítimo, no un artefacto.

library(cluster)
library(tidyverse)

# Función para calcular silhouette promedio por grupo
calc_silhouette_by_group <- function(beta_dist, clusters_named) {
  cell_ids <- attr(beta_dist, "Labels")
  clust_vec <- clusters_named[cell_ids]
  
  sil <- silhouette(clust_vec, beta_dist)
  sil_df <- as.data.frame(sil)
  
  # Agrupar por cluster y calcular promedio
  sil_summary <- sil_df %>%
    group_by(cluster) %>%
    summarise(
      mean_silhouette = mean(sil_width, na.rm = TRUE),
      n_cells = n()
    ) %>%
    mutate(cluster = as.character(cluster))  # Asegurar tipo texto
  
  return(sil_summary)
}

# Calcular silhouette por grupo para cada escala
sil_by_group_100 <- calc_silhouette_by_group(beta_sim_mean100, clusters) %>%
  mutate(scale = "100 km")
sil_by_group_200 <- calc_silhouette_by_group(beta_sim_mean200, clusters200) %>%
  mutate(scale = "200 km")
sil_by_group_400 <- calc_silhouette_by_group(beta_sim_mean400, clusters400) %>%
  mutate(scale = "400 km")
sil_by_group_800 <- calc_silhouette_by_group(beta_sim_mean800, clusters800) %>%
  mutate(scale = "800 km")

# Combinar todos en un solo dataframe
silhouette_all_groups <- bind_rows(
  sil_by_group_100,
  sil_by_group_200,
  sil_by_group_400,
  sil_by_group_800
)

# Revisar
print(silhouette_all_groups, n = Inf)

#############################################################
# Visualizar dendrogramas por escala
library(tidyverse)

# Función para graficar dendrograma
plot_dendrogram <- function(beta_dist, k, scale_label) {
  hc <- hclust(beta_dist, method = "average")  # Método UPGMA (average linkage)
  
  # Graficar dendrograma
  plot(hc, labels = FALSE, hang = -1,
       main = paste("Dendrograma -", scale_label),
       xlab = "", ylab = "Altura (β_sim)")
  
  # Línea horizontal para el corte (opcional)
  rect.hclust(hc, k = k, border = "red")
}

# Visualizar para cada escala con número de clusters conocido
par(mfrow = c(2, 2))  # Panel de 2x2

# 100 km (5 regiones)
plot_dendrogram(beta_sim_mean100, k = 5, scale_label = "100 km")
# 200 km (6 regiones)
plot_dendrogram(beta_sim_mean200, k = 6, scale_label = "200 km")
# 400 km (10-20 regiones, elige un valor intermedio)
plot_dendrogram(beta_sim_mean400, k = 20, scale_label = "400 km")
# 800 km (16 regiones)
plot_dendrogram(beta_sim_mean800, k = 16, scale_label = "800 km")

# Restaurar parámetros de gráfico
par(mfrow = c(1, 1))

# Altura de corte para cada k
library(tibble)
# Función para extraer alturas de corte de hclust
get_cut_heights <- function(dist_obj, method = "average", max_k = 20) {
  hc <- hclust(dist_obj, method = method)
  sapply(2:max_k, function(k) {
    hc$height[length(hc$height) - (k - 2)]
  })
}

# Aplicar a cada escala
heights_100 <- get_cut_heights(beta_sim_mean100)
heights_200 <- get_cut_heights(beta_sim_mean200)
heights_400 <- get_cut_heights(beta_sim_mean400)
heights_800 <- get_cut_heights(beta_sim_mean800)

# Unir en un solo data frame para comparar
cut_heights_df <- tibble(
  k = 2:20,
  `100 km` = heights_100,
  `200 km` = heights_200,
  `400 km` = heights_400,
  `800 km` = heights_800)

print(cut_heights_df, n = Inf)

#Exportar la argumentacion junto con las tablas y/o figuras
library(officer)
library(flextable)
library(dplyr)

# === Crear documento Word ===
doc <- read_docx()

# === Función auxiliar para añadir subtítulos en negrita ===
body_add_bold_par <- function(doc, title_text) {
  doc %>% body_add_fpar(fpar(ftext(title_text, prop = fp_text(bold = TRUE))))
}

# === Título principal ===
doc <- doc %>%
  body_add_par("Effect of Grain Size on Bioregionalisation of Orchids", style = "heading 1") %>%
  body_add_par("This report presents a series of analyses designed to understand why an increase in spatial grain size leads to a greater number of biogeographical regions in the clustering of global orchid phylogenetic beta diversity. Contrary to the prevailing assumption that coarser grains yield fewer, more homogeneous regions, our results reveal the opposite trend.", style = "Normal")

# === 1. Distribución de β_sim ===
doc <- doc %>%
  body_add_bold_par("Distribution of mean β_sim across spatial scales") %>%
  body_add_par("To interpret why broader grain sizes yielded more biogeographical regions, we first examined whether spatial aggregation increased or decreased floristic dissimilarity between grid cells. We calculated the distribution of mean β_sim across the four spatial resolutions (100, 200, 400, and 800 km), using the same distance matrices employed for clustering. This analysis was designed to test whether spatial aggregation reduces compositional differences, which would support fewer, broader regions. As expected, the average β_sim values slightly decreased with grain size (from 0.62 at 100 km to 0.60 at 800 km), indicating some spatial homogenisation. However, the overall variance in dissimilarity increased, with larger standard deviations and broader ranges at coarser grains. These results suggest that while many grid cells became more similar as grain size increased, others remained highly differentiated, contributing to greater heterogeneity.", style = "Normal") %>%
  body_add_par("Summary statistics of β_sim values by grain size.", style = "Normal") %>%
  body_add_flextable(
    flextable(beta_summary) %>%
      autofit() %>%
      fontsize(size = 8) %>%
      theme_booktabs()
  )

# === 2. Disimilitud entre clústeres ===
doc <- doc %>%
  body_add_bold_par("Between-cluster phylogenetic dissimilarity across spatial scales") %>%
  body_add_par("Based on our hypothesis, we expected that coarser grains would yield more internally homogeneous clusters—reflected in higher silhouette values. However, the number of clusters increased rather than decreased with scale, raising questions about their coherence. To investigate this, we calculated silhouette scores both overall and by cluster, at each grain size. The results indicate slightly higher between-cluster dissimilarities at 200 and 800 km, though this pattern is not linear. Higher variance at coarser grains suggests more fragmented regionalisation at broad scales.", style = "Normal") %>%
  body_add_par("Average β_sim between clusters per grain size.", style = "Normal") %>%
  body_add_flextable(
    flextable(between_summary) %>%
      autofit() %>%
      fontsize(size = 8) %>%
      theme_booktabs()
  )

# === 3. Evaluación de silueta ===
doc <- doc %>%
  body_add_bold_par("Clustering quality based on silhouette widths") %>%
  body_add_par("To assess the coherence of cluster assignments, silhouette widths were calculated for each grain size. The average silhouette peaked at 200 km (0.340), followed by 100 km (0.318), 800 km (0.317), with lower values at 400 km (0.202). Cluster-level silhouettes at 800 km showed a wide range, from 0.03 to 0.98, with several large clusters (e.g. clusters 1, 5, and 9) having values > 0.5—indicating good internal consistency. These results suggest that many of the additional clusters at 800 km are internally coherent, rather than artefacts. The 200 km grain appears optimal in balancing within-cluster homogeneity and regional distinctiveness, but 800 km still retains biologically meaningful structure.", style = "Normal") %>%
  body_add_par("Average silhouette width per grain size.", style = "Normal") %>%
  body_add_flextable(
    flextable(sil_summary) %>%
      autofit() %>%
      fontsize(size = 8) %>%
      theme_booktabs()
  ) %>%
  body_add_par("Silhouette width per cluster per scale.", style = "Normal") %>%
  body_add_flextable(
    flextable(silhouette_all_groups) %>%
      autofit() %>%
      fontsize(size = 7) %>%
      theme_booktabs()
  )

# === 4. Altura de corte en dendrogramas ===
doc <- doc %>%
  body_add_bold_par("Hierarchical depth of clusters (dendrogram analysis)") %>%
  body_add_par("To determine whether the additional regions at coarser scales arose from low clustering thresholds (i.e. shallow dendrogram cuts), we evaluated the height values at which clusters formed in UPGMA dendrograms across grain sizes. At a given number of clusters (e.g. k = 6 or k = 10), the cut heights were consistently lower at coarser grains: for instance, the height at k = 10 was 0.596 for 100 km but only 0.493 for 800 km. This indicates that regions at broader grains were defined at shallower levels of the dendrogram, supporting the idea that finer-scale structure becomes subdivided or re-partitioned differently with spatial aggregation.", style = "Normal") %>%
  body_add_par("Dendrogram cut heights by cluster number (k) and grain size.", style = "Normal") %>%
  body_add_flextable(
    flextable(cut_heights_df) %>%
      autofit() %>%
      fontsize(size = 7) %>%
      theme_booktabs()
  )

# === 5. Comentarios finales ===
doc <- doc %>%
  body_add_bold_par("Final remarks") %>%
  body_add_par("Across all analyses, we found that larger grain sizes did not reduce the number of regions or increase floristic homogeneity, as initially hypothesised. Instead, broader grains revealed additional, internally coherent regions that often corresponded to subdivisions of finer-grain clusters. These results suggest that spatial aggregation does not dilute biogeographical signal, but may instead highlight hierarchical structuring of orchid diversity. Consequently, coarse grains offer complementary perspectives on floristic organisation rather than obscured or artefactual views.", style = "Normal")

# === Guardar documento ===
print(doc, target = "results/table/Grain_size_Report.docx")
