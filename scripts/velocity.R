
## Load singularity container.
##
## singularity shell -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_seurat_velocytor_0.2.sif
##
## . /opt/conda/etc/profile.d/conda.sh
## conda activate seurat; R

library("Seurat")
library("tidyverse")
library("SeuratWrappers")
library("velocyto.R")

###########################
## RNA Velocity Analysis ##
###########################

setwd("..")

## Load Integrated Data
## ----------

seurat_integrated <- readRDS("integrated_spliced.RDS")

## Plot RNA Velocity
## ----------

ident.colors <- (scales::hue_pal())(n = length(x = levels(x = seurat_integrated)))
names(x = ident.colors) <- levels(x = seurat_integrated)

cell.colors <- ident.colors[Idents(object = seurat_integrated)]
names(x = cell.colors) <- colnames(x = seurat_integrated)

p <- show.velocity.on.embedding.cor(
	emb = Embeddings(object = seurat_integrated, reduction = "umap"),
	vel = Tool(object = seurat_integrated, slot = "RunVelocity"), n = 200, scale = "sqrt",
	cell.colors = ac(x = cell.colors, alpha = 0.5), cex = 0.8, arrow.scale = 3,
	show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
	do.par = FALSE, cell.border.alpha = 0.1
)
