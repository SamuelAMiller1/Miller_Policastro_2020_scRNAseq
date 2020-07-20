
library("tidyverse")
library("data.table")
library("Seurat")
library("SeuratWrappers")
library("monocle3")

seu <- readRDS("results/r_objects/seurat_expanded.RDS")

samples <- unique(seu$orig.ident)
seu <- map(setNames(samples, samples), function(x) {
  x <- subset(seu, subset = orig.ident == x)
  return(x)
})

cds <- map(seu, function(x) {

  x[["UMAP"]] <- x[["umap"]]
  x[["umap"]] <- NULL

  x[["PCA"]] <- x[["pca"]]
  x[["pca"]] <- NULL

  cds <- as.cell_data_set(x)
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)

  return(cds)

})

plot_cells(
  cds$H508_LSD1_KD, color_cells_by = "custom_clusters", label_groups_by_cluster=FALSE,
  label_leaves=FALSE, label_branch_points=FALSE
)
