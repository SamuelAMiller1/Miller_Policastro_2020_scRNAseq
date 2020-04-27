
library("Seurat")

##################
## Marker Genes ##
##################

seurat_integrated <- readRDS(file.path("results", "r_objects", "seurat_integrated.RDS"))
markers <- FindAllMarkers(seurat_integrated, assay = "RNA", min.pct = 0.25)
saveRDS(markers, file.path("results", "r_objects", "markers.RDS"))
