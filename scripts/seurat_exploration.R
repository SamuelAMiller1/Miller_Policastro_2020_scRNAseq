
## load singularity container.
##
## singularity shell -eCB `pwd` -H `pwd` scrnaseq_software_seurat_velocytor_0.3.sif
##
## . /opt/conda/etc/profile.d/conda.sh
## conda activate seurat; R

library("Seurat")
library("tidyverse")
library("data.table")
library("future")
library("wesanderson")
library("readxl")
library("scProportionTest")

## Variables.

###################################
## Exploration of scRNA-seq Data ##
###################################

setwd("..")

options(future.globals.maxSize = 10000 * 1024 ^2)
plan("multiprocess", workers = 2)

## Load Integrated Data
## ----------

seurat_integrated <- readRDS(file.path("results", "r_objects", "seurat_integrated.RDS"))

##################
## Marker Genes ##
##################

## Run the scripts seurat_marker_genes.sh first.

## Load and prepare marker list.

markers <- readRDS(file.path("results", "r_objects", "markers.RDS"))
setDT(markers)

markers <- markers[p_val_adj < 0.05]
markers[, c("avg_log2FC", "cluster") := list(log2(exp(avg_logFC)), str_c("cluster_", cluster))]
markers <- markers[order(p_val_adj)]

## Split the list based on cluster and save results.

markers <- split(markers, markers$cluster)

if (!dir.exists(file.path("results", "marker_tables"))) {
        dir.create(file.path("results", "marker_tables"), recursive = TRUE)
}

iwalk(markers, function(x, y) {
        file_name <- file.path("results", "marker_tables", str_c("markers_", y, ".tsv"))
        fwrite(x, file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}


##############################
## Defining Custom Clusters ##
##############################

## Read in the custom clusters.

custom_clusts <- read_xlsx(file.path("resources", "clust_classifications.xlsx")) %>%
	filter_all(any_vars(!is.na(.))) %>%
	dplyr::rename("cluster" = "Cluster(s)") %>%
	separate_rows(cluster, sep = ",") %>%
	dplyr::filter(
		!(cluster == 14 & Type == "TA"),
		cluster != 17
	) %>%
	dplyr::select(Type, cluster) %>%
	arrange(as.numeric(cluster)) %>%
	pull("Type") %>%
	set_names(levels(seurat_integrated))

## Add metadata column with new clusters.

seurat_integrated[["custom_clusters"]] <- seurat_integrated[["integrated_snn_res.0.6"]][[1]] %>%
	as.character %>%
	custom_clusts[.]

seurat_integrated[["custom_clusters"]] <- seurat_integrated[[c("custom_clusters", "orig.ident")]] %>%
	mutate(custom_clusters = ifelse(
		orig.ident == "COLON_1" & custom_clusters == "early EEC",
		"EEC", custom_clusters
	)) %>%
	pull("custom_clusters")

## Remove potential doublets and dying cells.

seurat_integrated <- subset(seurat_integrated, subset = custom_clusters != "dying" & custom_clusters != "delete doublets")

## Set factor order for custom clusters.

seurat_integrated[["custom_clusters"]] <- factor(seurat_integrated[["custom_clusters"]][[1]], levels = c(
	"TA", "early secretory/stem", "early EEC", "EEC", "early goblet", "goblet",
	"PLC", "Enterocyte", "Cancer/misc 1", "Cancer/misc 2"
))

## DimPlot of new clusters.

if (!dir.exists(file.path("results", "custom_clusters"))) {
	dir.create(file.path("results", "custom_clusters"))
}

p <- DimPlot(seurat_integrated, group.by = "custom_clusters", split.by = "orig.ident", ncol = 3)

pdf(file.path("results", "custom_clusters", "custom_clusters_dimplot.pdf"), height = 10, width = 16)
p
dev.off()

## Set custom clusters as default clusters and save object.

Idents(seurat_integrated) <- "custom_clusters"

saveRDS(seurat_integrated, file.path("results", "r_objects", "seurat_expanded.RDS"))

################
## Cell Cycle ##
################

if (!dir.exists(file.path("results", "cell_cycle"))) {
	dir.create(file.path("results", "cell_cycle"))
}

cell_cycle_palette <- wes_palette("Zissou1", 3, type = "continuous") %>%
	as.character

## Dim plots of cell cycle phase.

p <- DimPlot(
	seurat_integrated, group.by = "Phase", split.by = "orig.ident",
	ncol = 2, cols = cell_cycle_palette
)

pdf(file.path("results", "cell_cycle", "cell_cycle_dimplot.pdf"), height = 10, width = 10)
p
dev.off()

## Stacked bar plots of cell cycle phase.

cell_cycle_phase <- as.data.table(seurat_integrated@meta.data, keep.rownames = "cell_id")[,
	.(cell_id, orig.ident, Phase, custom_clusters)
]

p <- ggplot(cell_cycle_phase, aes(x = orig.ident, fill = Phase)) +
	geom_bar(position = "fill") +
	scale_fill_manual(values = cell_cycle_palette) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))


pdf(file.path("results", "cell_cycle", "cell_cycle_by_sample.pdf"), height = 4, width = 6)
p
dev.off()

p <- ggplot(cell_cycle_phase, aes(x = custom_clusters, fill = Phase)) +
	geom_bar(position = "fill") +
	scale_fill_manual(values = cell_cycle_palette) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path("results", "cell_cycle", "cell_cycle_by_cluster.pdf"), height = 4, width = 8)
p
dev.off()

p <- ggplot(cell_cycle_phase, aes(x = custom_clusters, fill = Phase)) +
	geom_bar(position = "fill") +
	scale_fill_manual(values = cell_cycle_palette) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	facet_grid(orig.ident ~ .)

pdf(file.path("results", "cell_cycle", "cell_cycle_by_sample_and_cluster.pdf"), height = 8, width = 8)
p
dev.off()

## Cell Cycle Counts
## ----------

## Prepare data for permutation test.

sc_utils_obj <- sc_utils(seurat_integrated)

comparisons <- list(
	c("COLON_1", "HT29_EV"),
	c("COLON_1", "H508_EV"),
	c("HT29_EV", "H508_EV"),
	c("HT29_EV", "HT29_LSD1_KD"),
	c("H508_EV", "H508_LSD1_KD"),
	c("HT29_LSD1_KD", "H508_LSD1_KD")
)

## Permutation tests amd plotting.

walk(comparisons, function(x) {
	sc_utils_obj <- permutation_test(
		sc_utils_obj, cluster_identity = "Phase",
		sample_1 = x[1], sample_2 = x[2]
	)

	p <- permutation_plot(sc_utils_obj, log2FD_threshold = 1)

	file_name <- str_c(x[1], "_vs_", x[2], ".pdf")
	pdf(file.path("results", "cell_cycle", file_name), height = 3, width = 8)
	print(p); dev.off()
})

#########################
## Cluster Cell Counts ##
#########################

## Prepare data for cluster count analysis.

sc_utils_obj <- sc_utils(seurat_integrated)

comparisons <- list(
        c("COLON_1", "HT29_EV"),
        c("COLON_1", "H508_EV"),
	c("HT29_EV", "H508_EV"),
        c("HT29_EV", "HT29_LSD1_KD"),
        c("H508_EV", "H508_LSD1_KD"),
        c("HT29_LSD1_KD", "H508_LSD1_KD")
)

## Permutation tests and plotting.

walk(comparisons, function(x) {
        sc_utils_obj <- permutation_test(
                sc_utils_obj, cluster_identity = "custom_clusters",
                sample_1 = x[1], sample_2 = x[2]
        )

        p <- permutation_plot(sc_utils_obj, log2FD_threshold = 1)

        file_name <- str_c(x[1], "_vs_", x[2], ".pdf")
        pdf(file.path("results", "cluster_counts", file_name), height = 3, width = 8)
        print(p); dev.off()
})
