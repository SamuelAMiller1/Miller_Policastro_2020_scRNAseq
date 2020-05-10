## Load singularity container.
##
## singularity shell -eCB `pwd` -H `pwd` scrnaseq_software_seurat_velocytor_0.3.sif
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
library("loomR")
library("velocyto.R")
library("SeuratWrappers")

##########################
## Seurat Analysis LSD1 ##
##########################

setwd("..")

# Commented out because it uses too much memory.
# options(future.globals.maxSize = 10000 * 1024 ^2)
# plan("multiprocess", workers = 2)

if (!dir.exists("results")) dir.create("results")

## Loading Data
## ----------

## Samples to analyze.

velocyto_samples <- list(
	HT29_EV = file.path("aligned", "HT29_EV", "velocyto", "HT29_EV.loom"),
	HT29_LSD1_KD = file.path("aligned", "HT29_LSD1_KD", "velocyto", "HT29_LSD1_KD.loom"),
	H508_EV = file.path("aligned", "H508_EV", "velocyto", "H508_EV.loom"),
	H508_LSD1_KD = file.path("aligned", "H508_LSD1_KD", "velocyto", "H508_LSD1_KD.loom"),
	COLON_1 = file.path("aligned", "COLON_1", "velocyto", "COLON_1.loom")
)

## Read in data.

velocyto_counts <- map(velocyto_samples, ReadVelocity)

## Convert to seurat object.

seurat_obj <- map(velocyto_counts, as.Seurat)

seurat_obj <- imap(seurat_obj, function(x, y) {
	x@meta.data$orig.ident <- y
	return(x)
})

## Cell Quality Control
## ----------

## Add mitochondrial percentage.

seurat_obj <- map(seurat_obj, function(x) {
	x[["percent.mt"]] <- PercentageFeatureSet(x, "^MT-")
	return(x)
})

## Plot mitochondrial percentage versus feature counts.

mt_data <- map(seurat_obj, function(x) {
	meta_data <- x@meta.data
	meta_data <- as_tibble(meta_data, .name_repair = "unique")
	return(meta_data)
})

mt_data <- bind_rows(mt_data)

p <- ggplot(mt_data, aes(x = percent.mt, y = nFeature_spliced)) +
	geom_point(size = 0.1) +
	facet_wrap(. ~ orig.ident, ncol = 3, scales = "free")

if (!dir.exists(file.path("results", "preprocessing"))) {
	dir.create(file.path("results", "preprocessing"))
}

png(
	file.path("results", "preprocessing", "mt_content_spliced.png"),
	height = 10, width = 10, res = 300, units = "in"
)
p
dev.off()

## Filter the data based on number of features and mitochondrial content.

seurat_obj <- imap(seurat_obj, function(x, y) {
	if (y == "COLON_1") {
		x <- subset(x, subset = percent.mt <= 30 & nFeature_spliced >= 750)
	} else {
		x <- subset(x, subset = percent.mt <= 25 & nFeature_spliced >= 2000)
	}
})

## Add cell-cycle scores.

# Gives error.
#seurat_obj <- map(seurat_obj, function(x) {
#	x <- CellCycleScoring(
#		x, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes,
#		set.ident = TRUE
#	)
#	return(x)
#})

## Save the object.

if (!dir.exists(file.path("results", "r_objects"))) {
	dir.create(file.path("results", "r_objects"))
}

saveRDS(seurat_obj, file.path("results", "r_objects", "seurat_obj_spliced.RDS"))

## Integrating Data
## ----------

## SCTransform to normalize data.

seurat_obj <- map(seurat_obj, ~ SCTransform(., assay = "spliced"))

saveRDS(seurat_obj, file.path("results", "r_objects", "seurat_obj_spliced.RDS"))

## Prepare for integration.

# Commented out because it uses too much memory.
# options(future.globals.maxSize = 10000 * 1024 ^2)

integration_features <- SelectIntegrationFeatures(seurat_obj, nfeatures = 3000)
seurat_obj <- PrepSCTIntegration(seurat_obj, anchor.features = integration_features)

## Integrate the reference dataset.

reference_datasets <- which(names(seurat_obj) %in% c("HT29_EV", "H508_EV"))
integration_anchors <- FindIntegrationAnchors(
	seurat_obj, normalization.method = "SCT", anchor.features = integration_features,
	reference = reference_datasets
)

seurat_integrated <- IntegrateData(integration_anchors, normalization.method = "SCT")

saveRDS(seurat_integrated, file.path("results", "r_objects", "seurat_integrated_spliced.RDS"))

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
seurat_integrated <- FindClusters(
	seurat_integrated, resolution = seq(0.2, 1.2, 0.1),
	method = "igraph", algorithm = 4, weights = TRUE
)

saveRDS(seurat_integrated, file.path("results", "r_objects", "seurat_integrated_spliced.RDS"))

## Plotting a cluster tree.

p <- clustree(seurat_integrated, prefix = "integrated_snn_res.")

pdf(file.path("results", "clustering", "cluster_tree_spliced.pdf"), height = 16, width = 12)
p
dev.off()

## Switch identity to a presumptive good clustering resolution.

Idents(seurat_integrated) <- "integrated_snn_res.0.6"

## UMAP dimension reduction for visualization.

if (!dir.exists("tempdir")) dir.create("tempdir")
set.tempdir("tempdir")

seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:30)

saveRDS(seurat_integrated, file.path("results", "r_objects", "seurat_integrated_spliced.RDS"))

## Plot Clusters.

p <- DimPlot(seurat_integrated, group.by = "ident", split.by = "orig.ident", ncol = 2)

pdf(file.path("results", "clustering", "clusters_spliced.pdf"), height = 12, width = 10)
p
dev.off()

## RNA Velocity
## ----------

# Takes a long time to run so will move to its own script.
#seurat_velocity <- RunVelocity(
#	seurat_integrated, ambiguous = "ambiguous", ncores = 1,
#	deltaT = 1, kCells = 25, fit.quantile = 0.02
#)

## Export to Cerebro
## ----------

## Preparation steps for cerebro.

# Already addressed by seurat.
#seurat_cerebro <- addPercentMtRibo(
#	seurat_integrated, assay = "SCT", organism = "hg",
#	gene_nomenclature = "name"
#)

# Crashes the cerebro browser.
#seurat_cerebro <- getMostExpressedGenes(
#	seurat_cerebro, assay = "SCT", column_sample = "orig.ident",
#	column_cluster = "integrated_snn_res.0.8"
#)

# Crashes the cerebro browser.
#seurat_cerebro <- getMarkerGenes(
#	seurat_cerebro, assay = "RNA", column_sample = "orig.ident",
#	column_cluster = "integrated_snn_res.0.8", only_pos = FALSE,
#	organism = "hg"
#)

# No marker genes calculated so can't export.
#seurat_cerebro <- getEnrichedPathways(
#	seurat_integrated, column_sample = "orig.ident",
#	column_cluster = "integrated_snn_res.0.8"
#)

## Export cerebro object.

if (!dir.exists(file.path("results", "cerebro"))) {
	dir.create(file.path("results", "cerebro"))
}

exportFromSeurat(
	seurat_integrated, assay = "SCT", file = file.path("results", "cerebro", "cerebro_spliced.crb"),
	experiment_name = "LSD1_KD", organism = "hg", column_sample = "orig.ident",
	column_cluster = "integrated_snn_res.0.6", column_nUMI = "nCount_spliced",
	column_nGene = "nFeature_spliced"
)

## Export as Loom for PAGA
## ----------

if (!dir.exists(file.path("results", "loom"))) {
	dir.create(file.path("results", "loom"))
}

seurat_loom <- as.loom(
	seurat_integrated,
	filename = file.path("results", "loom", "seurat_integrated_spliced.loom")
)

seurat_loom$close_all()
