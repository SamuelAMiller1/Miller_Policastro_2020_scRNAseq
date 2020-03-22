
## load singularity container.
##
## singularity shell -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_seurat_3.1.4.sif
##
## . /opt/conda/etc/profile.d/conda.sh
## conda activate seurat; R

library("Seurat")
library("tidyverse")
library("data.table")
library("future")
library("wesanderson")

###################################
## Exploration of scRNA-seq Data ##
###################################

setwd("..")

options(future.globals.maxSize = 10000 * 1024 ^2)
plan("multiprocess", workers = 2)

## Load Integrated Data
## ----------

seurat_integrated <- readRDS("integrated.RDS")

## Cell Cycle
## ----------

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
	.(cell_id, orig.ident, Phase, integrated_snn_res.0.8)
]

cell_cycle_phase[, integrated_snn_res.0.8 := factor(
	integrated_snn_res.0.8, levels = seq(0, max(as.numeric(integrated_snn_res.0.8)))
)]

p <- ggplot(cell_cycle_phase, aes(x = orig.ident, fill = Phase)) +
	geom_bar(position = "fill") +
	scale_fill_manual(values = cell_cycle_palette) +
	theme_bw()

pdf(file.path("results", "cell_cycle", "cell_cycle_by_sample.pdf"), height = 4, width = 6)
p
dev.off()

p <- ggplot(cell_cycle_phase, aes(x = integrated_snn_res.0.8, fill = Phase)) +
	geom_bar(position = "fill") +
	scale_fill_manual(values = cell_cycle_palette) +
	theme_bw()

pdf(file.path("results", "cell_cycle", "cell_cycle_by_cluster.pdf"), height = 4, width = 8)
p
dev.off()

p <- ggplot(cell_cycle_phase, aes(x = integrated_snn_res.0.8, fill = Phase)) +
	geom_bar(position = "fill") +
	scale_fill_manual(values = cell_cycle_palette) +
	theme_bw() +
	facet_grid(orig.ident ~ .)

pdf(file.path("results", "cell_cycle", "cell_cycle_by_sample_and_cluster.pdf"), height = 8, width = 8)
p
dev.off()

## Cluster Cell Counts
## ----------

## Prepare data for cluster count analysis.

cluster_counts <- as.data.table(seurat_integrated@meta.data, keep.rownames = "cell_id")[
	orig.ident != "COLON_1",
	.(cell_id, orig.ident, integrated_snn_res.0.8,
	line = str_extract(orig.ident, "^HT?\\d+"),
	condition = str_extract(orig.ident, "(EV|LSD1_KD)$"))
]
cluster_counts[, orig.ident := NULL]
cluster_counts <- split(cluster_counts, cluster_counts$line)

## Get observed cluster fractional differences.

observed_counts <- map(cluster_counts, function(x) {
	x <- x[,
		.(count = .N), by = .(condition, integrated_snn_res.0.8)
	]
	x[, fraction := count / sum(count), by = condition]
	x <- dcast(x, integrated_snn_res.0.8 ~ condition, value.var = "fraction") 
	x[, obs_log2_frac_diff := log2(LSD1_KD) - log2(EV)]

	return(x)
})

## Get permutation cluster fractional differences.

permuted_counts <- map(cluster_counts, function(x) {
	perm_samples <- modelr::permute(x, n = 10, .id = "resample") %>%
		pull(perm) %>%
		map(function(resample) {
			resampled_data <- resample$data
			resampled_data$condition <- resampled_data$condition[resample$idx]
			
			resampled_data <- resampled_data[,
				.(count = .N), by = .(condition, integrated_snn_res.0.8)
			]
			resampled_data[, fraction := count / sum(count), by = condition]
			resampled_data <- dcast(
				resampled_data, integrated_snn_res.0.8 ~ condition,
				value.var = "fraction"
			)
			resampled_data[, sim_log2_frac_diff := log2(LSD1_KD) - log2(EV)]
			resampled_data[, c("EV", "LSD1_KD") := NULL]

			return(resampled_data)
		})
	perm_samples <- rbindlist(perm_samples, idcol = "resample")
	perm_samples <- dcast(
		perm_samples, integrated_snn_res.0.8 ~ resample,
		value.var = "sim_log2_frac_diff"
	)

	return(perm_samples)
})
