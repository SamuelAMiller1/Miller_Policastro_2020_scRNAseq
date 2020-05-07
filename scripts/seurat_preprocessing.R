## Load singularity container.
##
## singularity shell -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_seurat_velocytor_0.3.sif
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

samples_10x <- list(
	HT29_EV = file.path("aligned", "HT29_EV", "outs", "filtered_feature_bc_matrix"),
	HT29_LSD1_KD = file.path("aligned", "HT29_LSD1_KD", "outs", "filtered_feature_bc_matrix"),
	H508_EV = file.path("aligned", "H508_EV", "outs", "filtered_feature_bc_matrix"),
	H508_LSD1_KD = file.path("aligned", "H508_LSD1_KD", "outs", "filtered_feature_bc_matrix"),
	COLON_1 = file.path("aligned", "COLON_1", "outs", "filtered_feature_bc_matrix")
)

## Read in data.

counts_10X <- map(samples_10x, Read10X)

## Convert to seurat object.

seurat_obj <- imap(counts_10X, ~ CreateSeuratObject(
	counts = .x, project = .y, min.cells = 10, min.features = 250
))

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

p <- ggplot(mt_data, aes(x = percent.mt, y = nFeature_RNA)) +
	geom_point(size = 0.1) +
	facet_wrap(. ~ orig.ident, ncol = 3, scales = "free")

if (!dir.exists(file.path("results", "preprocessing"))) {
	dir.create(file.path("results", "preprocessing"))
}

pdf(file.path("results", "preprocessing", "mt_content.pdf"), height = 10, width = 10)
p
dev.off()

## Filter the data based on number of features and mitochondrial content.

seurat_obj <- imap(seurat_obj, function(x, y) {
	if (y == "COLON_1") {
		x <- subset(x, subset = percent.mt <= 30 & nFeature_RNA >= 1000)
	} else {
		x <- subset(x, subset = percent.mt <= 25 & nFeature_RNA >= 2500)
	}
})

## Add cell-cycle scores.

seurat_obj <- map(seurat_obj, function(x) {
	x <- CellCycleScoring(
		x, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes,
		set.ident = TRUE
	)
	return(x)
})

## Save the object.

if (!dir.exists(file.path("results", "r_objects"))) {
	dir.create(file.path("results", "r_objects"))
}

saveRDS(seurat_obj, file.path("results", "r_objects", "seurat_obj.RDS"))

## Integrating Data
## ----------

## SCTransform to normalize data.

seurat_obj <- map(seurat_obj, SCTransform)

saveRDS(seurat_obj, file.path("results", "r_objects", "seurat_obj.RDS"))

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

saveRDS(seurat_integrated, file.path("results", "r_objects", "seurat_integrated.RDS"))

## Dimensionality Reduction and Clustering
## ----------

## PCA dimension reduction for clustering.

seurat_integrated <- RunPCA(seurat_integrated, npcs = 100)

## Elbow plot of PCA dimensions.

p <- ElbowPlot(seurat_integrated, ndims = 100)

if (!dir.exists(file.path("results", "clustering"))) {
	dir.create(file.path("results", "clustering"))
}

pdf(file.path("results", "clustering", "pca_elbow_plot.pdf"), height = 5, width = 5)
p
dev.off()

## Clustering the data.

seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:30)
seurat_integrated <- FindClusters(
	seurat_integrated, resolution = seq(0.2, 1.2, 0.1),
	algorithm = 4, method = "igraph", weights = TRUE
)

## Plotting a cluster tree.

p <- clustree(seurat_integrated, prefix = "integrated_snn_res.")

pdf(file.path("results", "clustering", "cluster_tree.pdf"), height = 16, width = 12)
p
dev.off()

## Switch identity to a presumptive good clustering resolution.

Idents(seurat_integrated) <- "integrated_snn_res.0.6"

## UMAP dimension reduction for visualization.

if (!dir.exists("tempdir")) dir.create("tempdir")
set.tempdir("tempdir")

seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:30)

## Plot Clusters.

p <- DimPlot(seurat_integrated, group.by = "ident", split.by = "orig.ident", ncol = 2)

pdf(file.path("results", "clustering", "clusters.pdf"), height = 12, width = 10)
p
dev.off()

## Save Integrated Data
## ----------

saveRDS(seurat_integrated, file.path("results", "r_objects", "seurat_integrated.RDS"))

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
	seurat_integrated, assay = "SCT", file = file.path("results", "cerebro", "cerebro.crb"),
	experiment_name = "LSD1_KD", organism = "hg", column_sample = "orig.ident",
	column_cluster = "integrated_snn_res.0.7", column_nUMI = "nCount_RNA",
	column_nGene = "nFeature_RNA", column_cell_cycle_seurat = "Phase"
)

## Export as Loom for PAGA
## ----------

if (!dir.exists(file.path("results", "loom"))) {
	dir.create(file.path("results", "loom"))
}

seurat_loom <- as.loom(
	seurat_integrated,
	filename = file.path("results", "loom", "seurat_integrated.loom")
)

seurat_loom$close_all()
