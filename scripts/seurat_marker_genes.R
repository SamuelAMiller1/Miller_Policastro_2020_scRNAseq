
library("Seurat")

##################
## Marker Genes ##
##################

## For raw clusters.

seurat_integrated <- readRDS(file.path("results", "r_objects", "seurat_integrated.RDS"))
markers <- FindAllMarkers(seurat_integrated, assay = "RNA", min.pct = 0.25)
saveRDS(markers, file.path("results", "r_objects", "markers.RDS"))

## For custom clusters.

seurat_expanded <- readRDS(file.path("results", "r_objects", "seurat_expanded.RDS"))
Idents(seurat_expanded) <- "custom_clusters"

custom_markers <- FindAllMarkers(seurat_expanded, assay = "RNA", min.pct = 0.25)
saveRDS(custom_markers, file.path("results", "r_objects", "custom_markers.RDS"))
