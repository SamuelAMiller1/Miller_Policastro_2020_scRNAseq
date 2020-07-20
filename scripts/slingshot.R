
library("tidyverse")
library("data.table")
library("Seurat")
library("slingshot")
library("tradeSeq")

##################################
## Slingshot Trajectory Analysis ##
##################################

## Prepare the Data
## ----------

## Load the integrated data.

seu <- readRDS("results/r_objects/seurat_expanded.RDS")

## susbet the data.

seu <- subset(seurat_obj, subset = custom_clusters %in% c(
  "early EEC", "EEC", "early goblet", "goblet"
))

## Slingshot
## ----------

## Run slingshot.

traj <- slingshot(Embeddings(seu, "umap"), clusterLabels = seu$integrated_snn_res.0.6)

saveRDS(traj, file.path("results", "r_objects", "slingshot_trajectory.RDS"))

## Color palettes for plotting.

cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

cell_colors_ident <- cell_pal(seu$orig.ident, scales::hue_pal())
cell_colors_clust <- cell_pal(seu$custom_clusters, scales::hue_pal())

## Plot the overall results.

if (!dir.exists(file.path("results", "trajectory"))) {
  dir.create(file.path("results", "trajectory"))
}

pdf(file.path("results", "trajectory", "slingshot_overall.pdf"), height = 10, width = 12)
par(mfrow = c(3, 1))
  plot(reducedDim(traj), col = cell_colors_clust, pch = 16, cex = 0.5)
  lines(traj, lwd = 2, type = 'lineages', col = 'black')
dev.off()

## Calculate pseudotime over the various lineages.

pseudotime <- map(trajectory, slingPseudotime)

## Plot pseudotime over the various lineages.

curve_names <- colnames(pseudotime)
pal <- viridis::viridis(100, end = 0.95)

pdf(file.path("results", "trajectory", "slingshot_pseudotime.pdf"), height = 10, width = 10)
par(mfrow = c(3, 2))
for (i in curve_names) {
  colors <- pal[cut(pseudotime[,i], breaks = 100)]
  plot(reducedDim(trajectory), col = colors, pch = 16, cex = 0.5, main = i)
  lines(trajectory, lwd = 2, col = 'black', type = 'lineages')
}
dev.off()

