
## load singularity container.
##
## singularity shell -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_seurat_velocytor_0.3.sif
##
## . /opt/conda/etc/profile.d/conda.sh
## conda activate seurat; R

library("Seurat")
library("tidyverse")
library("data.table")
library("future")
library("wesanderson")

## Variables.

###################################
## Exploration of scRNA-seq Data ##
###################################

setwd("..")

options(future.globals.maxSize = 10000 * 1024 ^2)
plan("multiprocess", workers = 2)

## Load Integrated Data
## ----------

seurat_integrated <- readRDS("integrated.RDS")

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
	.(cell_id, orig.ident, Phase, integrated_snn_res.0.7)
]

cell_cycle_phase[, integrated_snn_res.0.7 := factor(
	integrated_snn_res.0.7, levels = seq(0, max(as.numeric(integrated_snn_res.0.7)))
)]

p <- ggplot(cell_cycle_phase, aes(x = orig.ident, fill = Phase)) +
	geom_bar(position = "fill") +
	scale_fill_manual(values = cell_cycle_palette) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))


pdf(file.path("results", "cell_cycle", "cell_cycle_by_sample.pdf"), height = 4, width = 6)
p
dev.off()

p <- ggplot(cell_cycle_phase, aes(x = integrated_snn_res.0.7, fill = Phase)) +
	geom_bar(position = "fill") +
	scale_fill_manual(values = cell_cycle_palette) +
	theme_bw()

pdf(file.path("results", "cell_cycle", "cell_cycle_by_cluster.pdf"), height = 4, width = 8)
p
dev.off()

p <- ggplot(cell_cycle_phase, aes(x = integrated_snn_res.0.7, fill = Phase)) +
	geom_bar(position = "fill") +
	scale_fill_manual(values = cell_cycle_palette) +
	theme_bw() +
	facet_grid(orig.ident ~ .)

pdf(file.path("results", "cell_cycle", "cell_cycle_by_sample_and_cluster.pdf"), height = 8, width = 8)
p
dev.off()

#########################
## Cluster Cell Counts ##
#########################

## Prepare data for cluster count analysis.

cluster_counts <- as.data.table(seurat_integrated@meta.data, keep.rownames = "cell_id")[
	!str_detect(orig.ident, "colon"),
	.(cell_id, orig.ident, integrated_snn_res.0.7,
	line = str_extract(orig.ident, "^HT?\\d+"),
	condition = str_extract(orig.ident, "(EV|LSD1_KD)$"))
]
cluster_counts[, orig.ident := NULL]
cluster_counts <- split(cluster_counts, cluster_counts$line)

## Get observed cluster fractional differences.

observed_counts <- map(cluster_counts, function(x) {
	x <- x[,
		.(count = .N), by = .(condition, integrated_snn_res.0.7)
	]
	x[, fraction := count / sum(count), by = condition]
	x <- dcast(x, integrated_snn_res.0.7 ~ condition, value.var = "fraction") 
	x[, obs_log2_frac_diff := log2(LSD1_KD) - log2(EV)]

	return(x)
})

## Merge back the observed cluster fractions with the cluster counts.

merged <- map2(cluster_counts, observed_counts, function(x, y) {
	obs <- y[, .(integrated_snn_res.0.7, obs_log2_frac_diff)]
	setkey(obs, integrated_snn_res.0.7)

	clusts <- copy(x)
	setkey(clusts, integrated_snn_res.0.7)

	clusts <- merge(clusts, obs)
	return(clusts)
})

## Get permutation cluster fractional differences.

permuted_counts <- map(merged, function(x) {
	perm_samples <- modelr::permute(x, n = 10000, .id = "resample") %>%
		pull(perm) %>%
		map(function(resample) {
			resampled_data <- resample$data
			resampled_data[, condition := resampled_data$condition[resample$idx]]
			
			resampled_data <- resampled_data[,
				.(count = .N),
				by = .(condition, integrated_snn_res.0.7, obs_log2_frac_diff)
			]
			resampled_data[, fraction := count / sum(count), by = condition]
			resampled_data <- dcast(
				resampled_data, integrated_snn_res.0.7 + obs_log2_frac_diff ~ condition,
				value.var = "fraction"
			)
			resampled_data[, sim_log2_frac_diff := log2(LSD1_KD) - log2(EV)]
			resampled_data[, conditional := case_when(
				obs_log2_frac_diff > 0 & sim_log2_frac_diff >= obs_log2_frac_diff ~ TRUE,
				obs_log2_frac_diff < 0 & sim_log2_frac_diff <= obs_log2_frac_diff ~ TRUE,
				TRUE ~ FALSE
			)]
			resampled_data[, c(
				"EV", "LSD1_KD", "obs_log2_frac_diff", "sim_log2_frac_diff"
			) := NULL]

			return(resampled_data)
		})
	perm_samples <- rbindlist(perm_samples, idcol = "resample")
	perm_samples <- dcast(
		perm_samples, integrated_snn_res.0.7 ~ resample,
		value.var = "conditional"
	)
	perm_samples[,
		pvalue := (rowSums(.SD) + 1) / ((ncol(perm_samples) - 1) + 1),
		.SDcols = !"integrated_snn_res.0.7"
	]
	perm_samples <- perm_samples[, .(integrated_snn_res.0.7, pvalue)]
	perm_samples[, FDR := p.adjust(pvalue, "fdr")]

	return(perm_samples)
})

## Add p-values back to data.

perm_results <- map2(observed_counts, permuted_counts, function(x, y) {
	x <- setkey(x, integrated_snn_res.0.7)
	y <- setkey(y, integrated_snn_res.0.7)

	merged <- merge(x, y)
	return(merged)
})

## Bootstrap fraction difference.

boot_counts <- map(cluster_counts, function(x) {
	x <- split(x, x$condition)
	x <- map(x, function(y) {
		resampled_data <- modelr::bootstrap(y, n = 10000, id = "resample") %>%
			pull(strap) %>%
			map(function(resample) {
				resampled_data <- resample$data
				resampled_data <- resampled_data[resample$idx, ]
				resampled_data <- resampled_data[,
					.(count = .N), by = integrated_snn_res.0.7
				]
			})
		resampled_data <- rbindlist(resampled_data, idcol = "resample")
		resampled_data[, fraction := count / sum(count), by = resample]
		resampled_data[, count := NULL]
	})
	x <- rbindlist(x, idcol = "condition")
	x <- dcast(
		x, resample + integrated_snn_res.0.7 ~ condition,
		value.var = "fraction", fill = 0
	)
	x[, sim_log2_frac_diff := log2(LSD1_KD + 1E-10) - log2(EV + 1E-10)]
	x[, c("EV", "LSD1_KD") := NULL]
	x <- dcast(x, integrated_snn_res.0.7 ~ resample, value.var = "sim_log2_frac_diff")

	x[,
		c("boot_mean", "boot_lower_ci", "boot_upper_ci") := list(
			rowMeans(.SD, na.rm = TRUE),
			apply(.SD, 1, function(x) {
				x <- quantile(x, probs = 0.05, na.rm = TRUE)
				x <- as.numeric(x)
				return(x)
			}),
			apply(.SD, 1, function(x) {
				x <- quantile(x, probs = 0.95, na.rm = TRUE)
				x <- as.numeric(x)
				return(x)
			})
		),
		.SDcols = !"integrated_snn_res.0.7"
	]
	x <- x[, .(integrated_snn_res.0.7, boot_mean, boot_lower_ci, boot_upper_ci)]
	return(x)
})

## Add boostrap back to data.

boot_results <- map2(perm_results, boot_counts, function(x, y) {
	x <- setkey(x, integrated_snn_res.0.7)
	y <- setkey(y, integrated_snn_res.0.7)

	merged <- merge(x, y)
	return(merged)
})

## Export permutation and bootstrap results.

if (!dir.exists(file.path("results", "cluster_counts"))) {
	dir.create(file.path("results", "cluster_counts"))
}

write.table(
	boot_results, file.path("results", "cluster_counts", "cluster_counts_table.tsv"),
	sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

## Make cluster counts figure.

plot_results <- rbindlist(boot_results, idcol = "line")
plot_results[, integrated_snn_res.0.7 := fct_reorder(integrated_snn_res.0.7, boot_mean, .desc = TRUE)]
plot_results[, significant := ifelse(FDR < 0.05, "FDR < 0.05", "n.s.")]

p <- ggplot(plot_results, aes(x = integrated_snn_res.0.7, y = boot_mean, color = significant)) +
	geom_pointrange(aes(ymin = boot_lower_ci, ymax = boot_upper_ci)) +
	geom_hline(yintercept = 0, lty = 2) +
	scale_color_manual(values = c(cell_cycle_palette[3], "grey")) +
	theme_bw() +
	coord_flip() +
	ylim(c(-3.5,2)) +
	facet_wrap(line ~ ., ncol = 1, scales = "free")

pdf(file.path("results", "cluster_counts", "cluster_counts_pointrange.pdf"), height = 8, width = 6)
p
dev.off()

##################
## Marker Genes ##
##################

## Run the scripts seurat_marker_genes.sh first.

## Load and perpare marker list.

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
})
