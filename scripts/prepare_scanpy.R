
library("tidyverse")
library("Seurat")
library("SeuratWrappers")
library("SeuratDisk")
library("velocyto.R")

##############################
## Export Seurat for scanpy ##
##############################

## Seurat Spliced
## ----------

## Load samples.

seurat_obj <- readRDS(file.path("results", "r_objects", "seurat_integrated_spliced.RDS"))

samples <- setNames(unique(seurat_obj$orig.ident), unique(seurat_obj$orig.ident))
seurat_obj <- map(samples, function(x) {
  x <- subset(seurat_obj, subset = orig.ident == x)
  return(x)
})

## Save as H5Seurat.

if (!dir.exists(file.path("results", "py_objects"))) {
  dir.create(file.path("results", "py_objects"))
}

iwalk(seurat_obj, function(x, y) {
  x[["RNA"]] <- x[["spliced"]]
  DefaultAssay(x) <- "RNA"
  SaveH5Seurat(x, file.path("results", "py_objects", str_c(y, "seurat.h5seurat", sep = "_")))
})

## Convert to H5ad.

walk(names(seurat_obj), function(x) {
  file_name <- file.path("results", "py_objects", str_c(x, "seurat.h5seurat", sep = "_"))
  Convert(file_name, dest = "h5ad")
})
