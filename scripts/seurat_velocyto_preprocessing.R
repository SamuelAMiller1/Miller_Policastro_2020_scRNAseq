
## Load singularity container.
##
## singularity shell -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_seurat_velocytor_0.2.sif
##
## . /opt/conda/etc/profile.d/conda.sh
## conda activate seurat; R

## Libraries

library("Seurat")
library("tidyverse")
library("clustree")
library("future")
library("unixtools")
library("cerebroApp")
library("SeuratWrappers")
library("velocyto.R")

##########################
## Seurat Analysis LSD1 ##
##########################

setwd("..")

options(future.globals.maxSize = 10000 * 1024 ^2)
plan("multiprocess", workers = 3)

if (!dir.exists("results")) dir.create("results")

## Loading Data
## ----------

## Samples to analyze.

samples <- list(
	COLON_1 = file.path("velocity", "COLON_1.loom"),
	H508_EV = file.path("velocity", "H508_EV.loom"),
	H508_LSD1_KD = file.path("velocity", "H508_LSD1_KD.loom"),
	HT29_EV = file.path("velocity", "HT29_EV.loom"),
	HT29_LSD1_KD = file.path("velocity", "HT29_LSD1_KD.loom")
)

## Read in data.

data_velocyto <- map(samples, ReadVelocity)

## Convert to seurat object.

seurat_obj <- map(data_velocyto, as.Seurat)
seurat_obj <- imap(seurat_obj, function(x, y) {
	x@meta.data$orig.ident <- y
	return(x)
})

## Cell Quality Control
## ----------

## Add mitochondrial percentage.

seurat_obj <- map(seurat_obj, function(x) {
	x[["percent_mt_spliced"]] <- PercentageFeatureSet(x, "^MT-")
	return(x)
})

mt_data <- map(seurat_obj, function(x) {
	meta_data <- x@meta.data
	meta_data <- as_tibble(meta_data, .name_repair = "unique")
	return(meta_data)
})

mt_data <- bind_rows(mt_data)

p <- ggplot(mt_data, aes(x = percent_mt_spliced, y = nFeature_spliced)) +
	geom_point(size = 0.1) +
	facet_wrap(. ~ orig.ident, ncol = 3, scales = "free")

if (!dir.exists(file.path("results", "preprocessing"))) {
	dir.create(file.path("results", "preprocessing"))
}

## Plot mitochondrial percentage versus feature counts.

pdf(file.path("results", "preprocessing", "mt_content_spliced.pdf"), height = 8, width = 10)
p
dev.off()

## Filter the data based on number of features and mitochondrial content.

seurat_obj <- imap(seurat_obj, function(x, y) {
	if (y == "COLON_1") {
		x <- subset(x, subset = percent_mt_spliced <= 25 & nFeature_spliced >= 1000)
	} else {
		x <- subset(x, subset = percent_mt_spliced <= 25 & nFeature_spliced >= 2500)
	}
})

## Integrating Data
## ----------

## SCTransform to normalize data.

seurat_obj <- map(seurat_obj, ~ SCTransform(., assay = "spliced"))

## Prepare for integration.

integration_features <- SelectIntegrationFeatures(seurat_obj, nfeatures = 3000)
seurat_obj <- PrepSCTIntegration(seurat_obj, anchor.features = integration_features)

## Integrate the reference dataset.

reference_datasets <- which(names(seurat_obj) %in% c("HT29_EV", "H508_EV"))
integration_anchors <- FindIntegrationAnchors(
	seurat_obj, normalization.method = "SCT", anchor.features = integration_features,
	reference = reference_datasets
)

seurat_integrated <- IntegrateData(integration_anchors, normalization.method = "SCT")

## Dimensionality Reduction and Clustering
## ----------

## PCA dimension reduction for clustering.

seurat_integrated <- RunPCA(seurat_integrated, npcs = 100)

## Elbow plot of PCA dimensions.

p <- ElbowPlot(seurat_integrated, ndims = 100)

if (!dir.exists(file.path("results", "clustering"))) {
        dir.create(file.path("results", "clustering"))
}

pdf(file.path("results", "clustering", "pca_elbow_plot_spliced.pdf"), height = 5, width = 5)
p
dev.off()

## Clustering the data.

seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:30)
seurat_integrated <- FindClusters(seurat_integrated, resolution = seq(0.2, 2.0, 0.2))

## Plotting a cluster tree.

p <- clustree(seurat_integrated, prefix = "integrated_snn_res.")

pdf(file.path("results", "clustering", "cluster_tree_spliced.pdf"), height = 10, width = 10)
p
dev.off()

## Switch identity to a presumptive good clustering resolution.

Idents(seurat_integrated) <- "integrated_snn_res.0.8"

## UMAP dimension reduction for visualization.

if (!dir.exists("tempdir")) dir.create("tempdir")
set.tempdir("tempdir")

seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:30)

## Plot Clusters.

p <- DimPlot(seurat_integrated, group.by = "ident", split.by = "orig.ident", ncol = 2)

pdf(file.path("results", "clustering", "clusters_spliced.pdf"), height = 6, width = 10)
p
dev.off()

## RNA Velocity
## ----------

seurat_integrated <- RunVelocity(
	seurat_integrated, ambiguous = "ambiguous", ncores = 1,
	deltaT = 1, kCells = 25, fit.quantile = 0.02
)

## Save Integrated Data
## ----------

saveRDS(seurat_integrated, "integrated_spliced.RDS")


