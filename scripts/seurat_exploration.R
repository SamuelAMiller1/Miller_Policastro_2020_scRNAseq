
## load singularity container.
##
## singularity shell -eCB `pwd` -H `pwd` scrnaseq_software_seurat_velocytor_0.3.sif
##
## . /opt/conda/etc/profile.d/conda.sh && conda activate seurat && R

library("Seurat")
library("tidyverse")
library("data.table")
library("future")
library("wesanderson")
library("readxl")
library("scProportionTest")
library("clusterProfiler")
library("ReactomePA")
library("org.Hs.eg.db")

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

options(future.globals.maxSize = 10000 * 1024 ^2)
plan("multiprocess", workers = 6)

## Find Markers
## ----------

Idents(seurat_integrated) <- "integrated_snn_res.0.6"

markers <- FindAllMarkers(
	seurat_integrated, assay = "SCT", slot = "data", min.pct = 0.25,
	return.thresh = 0.05, logfc.threshold = log(1.5)
)

saveRDS(markers, file.path("results", "r_objects", "markers.RDS"))

## Load and prepare marker list.

setDT(markers)
markers <- markers[p_val_adj < 0.05]
markers[, c("avg_log2FC", "cluster") := list(log2(exp(avg_logFC)), str_c("cluster_", cluster))]
markers <- markers[order(cluster, -avg_log2FC)]

## Split the list based on cluster and save results.

markers <- split(markers, markers$cluster)

if (!dir.exists(file.path("results", "marker_tables"))) {
        dir.create(file.path("results", "marker_tables"), recursive = TRUE)
}

iwalk(markers, function(x, y) {
        file_name <- file.path("results", "marker_tables", str_c("markers_", y, ".tsv"))
        fwrite(x, file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
})


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

## Set custom clusters as default clusters and save object.

Idents(seurat_integrated) <- "custom_clusters"

saveRDS(seurat_integrated, file.path("results", "r_objects", "seurat_expanded.RDS"))


##########################
## Custom Cluster Plots ##
##########################

## DimPlot of new clusters.

if (!dir.exists(file.path("results", "custom_clusters"))) {
	dir.create(file.path("results", "custom_clusters"))
}

p <- DimPlot(seurat_integrated, group.by = "custom_clusters", split.by = "orig.ident", ncol = 3)

pdf(file.path("results", "custom_clusters", "custom_clusters_dimplot.pdf"), height = 10, width = 16)
p
dev.off()

## DimPlots of COLON_1, HT29_EV, and H508_EV colored by cluster.

select_cells <- which(seurat_integrated[["orig.ident"]][[1]] %in% c("COLON_1", "HT29_EV", "H508_EV"))
select_cells <- rownames(seurat_integrated[[]])[select_cells]

p <- DimPlot(seurat_integrated, cells = select_cells, group.by = "custom_clusters")

pdf(file.path("results", "custom_clusters", "COLON1_HT29EV_H508EV_dimplot.pdf"), height = 4, width = 6)
p
dev.off()

## DimPlots of COLON_1, HT29_EV, and H508_EV colored by cell cycle phase.

cell_cycle_palette <- wes_palette("Zissou1", 3, type = "continuous") %>%
        as.character

p <- DimPlot(seurat_integrated, cells = select_cells, group.by = "Phase", cols = cell_cycle_palette)

pdf(file.path("results", "custom_clusters", "COLON1_HT29EV_H508EV_cellcycle.pdf"), height = 4, width = 6)
p
dev.off()

## TFF3 expression plot.

meta_data <- as.data.table(seurat_integrated[[]], keep.rownames = "cell_id")

exp_data <- as.data.table(
	Assays(seurat_integrated, "SCT")@scale.data,
	keep.rownames = "gene"
)[
	gene == "TFF3"
]
exp_data <- melt(
	exp_data, id.vars = "gene",
	variable.name = "cell_id",
	value.name = "TFF3_SCT_scaled_UMI"
)
exp_data[, gene := NULL]

TFF3_merged <- merge(meta_data, exp_data, by = "cell_id", all.x = TRUE)
TFF3_merged <- TFF3_merged[orig.ident %in% c("HT29_EV", "H508_EV")]
TFF3_merged <- TFF3_merged[, TFF3_mean := mean(TFF3_SCT_scaled_UMI), by = custom_clusters]
TFF3_merged <- TFF3_merged[, custom_clusters := fct_reorder(custom_clusters, -TFF3_mean)][]

p <- ggplot(TFF3_merged, aes(x = custom_clusters, y = TFF3_SCT_scaled_UMI, fill = custom_clusters)) +
	geom_boxplot(outlier.size = 0.25) +
	scale_fill_manual(
		values = rev(as.character(wes_palette("Zissou1", 9, type = "continuous"))),
		guide = FALSE
	) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path("results", "custom_clusters", "TFF3_HT29EV_H508EV_violinplot.pdf"), height = 4.5, width = 4.5)
p; dev.off()

#######################
## Marker Gene Plots ##
#######################

## Pick genes to use for markers for each cluster.

marker_genes <- c(
	"TA" = "HMGN2", "early_secratory" = "MEX3A", "goblet" = "MUC2",
	"early_goblet" = "REG4", "early_EEC" = "PAX4", "PLC" = "SPIB",
	"cancer_misc_1" = "DUSP5", "cancer_misc_2" = "FOXQ1",
	"enterocyte" = "LGALS2", "EEC" = "CHGA"
)

## Make the feature plot.

DefaultAssay(seurat_integrated) <- "SCT"

select_cells <- which(seurat_integrated[["orig.ident"]][[1]] %in% c("COLON_1", "HT29_EV", "H508_EV"))
select_cells <- rownames(seurat_integrated[[]])[select_cells]

p <- FeaturePlot(
	seurat_integrated, features = marker_genes,
	cells = select_cells, ncol = 3
)

pdf(file.path("results", "custom_clusters", "cluster_markers.pdf"), height = 16, width = 16)
p; dev.off()

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

############################
## Custom Marker Analysis ##
############################

options(future.globals.maxSize = 10000 * 1024 ^2)
plan("multiprocess", workers = 4)

## Marker Analysis
## ----------

Idents(seurat_integrated) <- "custom_clusters"

custom_markers <- FindAllMarkers(
	seurat_integrated, assay = "SCT",
	slot = "data", min.pct = 0.25
)

saveRDS(custom_markers, file.path("results", "r_objects", "custom_markers.RDS"))

## Preparing Data
## ----------

setDT(custom_markers)
custom_markers[, avg_log2FC := log2(exp(avg_logFC))]
custom_markers <- custom_markers[p_val_adj < 0.05 & abs(avg_log2FC) >= log2(1.5)]
custom_markers <- custom_markers[order(cluster, p_val_adj)]

## Get cells to plot.

select_cells <- which(seurat_integrated[["orig.ident"]][[1]] %in% c("COLON_1", "HT29_EV", "H508_EV"))
select_cells <- rownames(seurat_integrated[[]])[select_cells]

## Export table.

fwrite(
	as.data.table(mutate_if(custom_markers, is.factor, as.character)),
	file.path("results", "custom_clusters", "marker_table.tsv"),
	sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

## Plotting
## ----------

## Get top markers per cluster.

top_markers <- custom_markers[avg_log2FC >= log2(1.5)]
top_markers[, cluster := factor(cluster, levels = c(
	"TA", "early secretory/stem", "early EEC", "EEC", "early goblet", "goblet",
	"PLC", "Enterocyte", "Cancer/misc 1", "Cancer/misc 2"
))]
top_markers <- top_markers[order(cluster, -avg_log2FC)]

top_markers <- top_markers[, head(.SD, 5), by = cluster]
top_markers <- top_markers[["gene"]]

## Heatmap of top markers.

p <- DoHeatmap(
	seurat_integrated, features = top_markers, group.by = "custom_clusters",
	assay = "SCT", slot = "data", angle = 90, cells = select_cells
) +
	scale_fill_viridis_c()

pdf(file.path("results", "custom_clusters", "marker_heatmap.pdf"), width = 16, height = 10)
p; dev.off()

#####################
## LSD1 Expression ##
#####################

## Prepare expression data.

meta_data <- as.data.table(seurat_integrated[[]], keep.rownames = "cell_id")

exp_data <- as.data.table(
        Assays(seurat_integrated, "SCT")@data,
        keep.rownames = "gene"
)[
        gene == "KDM1A"
]
exp_data <- melt(
        exp_data, id.vars = "gene",
        variable.name = "cell_id",
        value.name = "LSD1_SCT_scaled_UMI"
)

LSD1_merged <- merge(meta_data, exp_data, by = "cell_id", all.x = TRUE)

## LSD1 expression by sample.

LSD1_merged[, orig.ident := fct_reorder(
	orig.ident, LSD1_SCT_scaled_UMI,
	median, na.rm = TRUE, .desc = TRUE
)]

sample_colors <- wes_palette("Zissou1", 5, type = "continuous")

p <- ggplot(LSD1_merged, aes(x = orig.ident, y = LSD1_SCT_scaled_UMI)) +
	geom_boxplot(outlier.size = 0.25, aes(fill = orig.ident)) +
	theme_bw() +
	scale_fill_manual(values = sample_colors) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path("results", "gene_plots", "LSD1_sample_boxplot.pdf"), height = 4, width = 6)
p; dev.off()

## LSD1 expression by cluster.

LSD1_subset <- LSD1_merged[orig.ident %in% c("COLON_1", "HT29_EV", "H508_EV")]
LSD1_subset[, custom_clusters := fct_reorder(
	custom_clusters, LSD1_SCT_scaled_UMI,
        median, na.rm = TRUE, .desc = TRUE
)]

cluster_colors <- wes_palette("Zissou1", 10, type = "continuous")

p <- ggplot(LSD1_subset, aes(x = custom_clusters, y = LSD1_SCT_scaled_UMI)) +
        geom_boxplot(outlier.size = 0.25, aes(fill = custom_clusters)) +
        theme_bw() +
        scale_fill_manual(values = cluster_colors) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path("results", "gene_plots", "LSD1_cluster_boxplot.pdf"), height = 4, width = 8)
p; dev.off()

#########################
## Targeted Gene Plots ##
#########################

## Prepare expression data.

target_genes <- c(
	"NEUROG3", "PAX4", "ARX", "RCOR2",
	"NEUROD1", "CHGA", "PYY", "FEV"
)

meta_data <- as.data.table(seurat_integrated[[]], keep.rownames = "cell_id")

exp_data <- as.data.table(
        Assays(seurat_integrated, "SCT")@data,
        keep.rownames = "gene"
)[
        gene %in% target_genes
]
exp_data <- melt(
        exp_data, id.vars = "gene",
        variable.name = "cell_id",
        value.name = "Gene_scaled_UMI"
)

merged <- merge(meta_data, exp_data, by = "cell_id", all.x = TRUE)
merged <- merged[orig.ident %in% c("COLON_1", "HT29_EV", "H508_EV")]

## Plot expression per sample.

sample_colors <- wes_palette("Zissou1", 3, type = "continuous")

p <- ggplot(merged, aes(x = orig.ident, y = Gene_scaled_UMI)) +
	geom_boxplot(outlier.size = 0.25, aes(fill = orig.ident)) +
	theme_bw() +
	scale_fill_manual(values = sample_colors) +
	facet_grid(gene ~ .) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	ylim(0, 3)

pdf(file.path("results", "gene_plots", "targeted_genes_boxplot.pdf"), height = 8, width = 4)
p; dev.off()

#############################
## Differential Expression ##
#############################

options(future.globals.maxSize = 10000 * 1024 ^2)
plan("multiprocess", workers = 6)

Idents(seurat_integrated) <- "custom_clusters"

## Set comparisons.

comparisons <- list(
	"Early_EEC_vs_EEC" = c("early EEC", "EEC"),
	"EEC_vs_Goblet" = c("EEC", "goblet"),
	"Early_EEC_vs_Goblet" = c("early EEC", "goblet")
)

## Run differential expression.

diff_exp <- map(comparisons, function(x) {
	results <- FindMarkers(
		subset(seurat_integrated, subset = orig.ident %in% c("COLON_1", "HT29_EV", "HT29_EV")),
		assay = "SCT", slot = "data", ident.1 = x[1], ident.2 = x[2]
	)
	setDT(results, keep.rownames = "gene")
	return(results)
})

## Format output.

diff_exp <- rbindlist(diff_exp, idcol = "comparison")
diff_exp <- diff_exp[p_val_adj < 0.05 & abs(avg_logFC) > log(1.5)]
diff_exp[, avg_log2FC := log2(exp(avg_logFC))]
diff_exp <- diff_exp[order(comparison, -avg_log2FC)]

## export to table.

if (!dir.exists(file.path("results", "diff_expression"))) {
	dir.create(file.path("results", "diff_expression"))
}

fwrite(
	diff_exp, file.path("results", "diff_expression", "early_EEV_vs_EEC_vs_Goblet.tsv"),
	sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

############################
## Secreted Tumor Factors ##
############################

## Prepare the Data
## ----------

## Load up secreted tumor factors.

sec_tumor <- read_xlsx(file.path("resources", "Secretory_clusters_secreted_peptide_comparison.xlsx"))
setDT(sec_tumor)
setnames(sec_tumor, old = "Tumor-associated", "tumor_associated")

## Get expression values for the secreted tumor factors.

exp_data <- as.data.table(
        Assays(seurat_integrated, "SCT")@data,
        keep.rownames = "gene"
)[
        gene %in% sec_tumor[["Gene"]]
]
exp_data <- melt(
        exp_data, id.vars = "gene",
        variable.name = "cell_id",
        value.name = "norm_counts"
)

## Add back the meta-data to the expression values.

meta_data <- as.data.table(seurat_integrated[[]], keep.rownames = "cell_id")

merged <- merge(meta_data, exp_data, by = "cell_id", all.x = TRUE)
merged <- merged[
	orig.ident %in% c("COLON_1", "HT29_EV", "H508_EV") &
	custom_clusters %in% c("early EEC", "EEC", "goblet")
]

## Merge back in the tumor secratory status info.

setnames(sec_tumor, old = "Gene", new = "gene")
merged <- merge(merged, sec_tumor, all.x = TRUE, by = "gene")


## Plot the Data
## ----------

merged[, c("gene", "custom_clusters") := list(
	factor(gene, levels = sec_tumor[["gene"]]),
	factor(custom_clusters, levels = c("early EEC", "EEC", "goblet"))
)]

p <- ggplot(merged, aes(x = custom_clusters, y = norm_counts, fill = custom_clusters)) +
	geom_violin() +
	theme_minimal() +
	theme(
		strip.text.y = element_text(angle = 0),
		axis.text.x = element_blank(),
		axis.ticks.x = element_blank()
	) +
	facet_grid(gene ~ ., scales = "free") +
	scale_fill_manual(values = wes_palette("Zissou1", 3, type = "continuous"))

if (!dir.exists(file.path("results", "secreted"))) {
        dir.create(file.path("results", "secreted"))
}

ggsave(
	file.path("results", "secreted", "secreted_tumor_factors.pdf"),
	plot = p, device = cairo_pdf, height = 18, width = 4
)

########################
## Differential Genes ##
########################

## Get Differential Genes
## ----------

## Set up comparisons.

compare_clusts <- c("early EEC", "EEC", "goblet")
compare_clusts <- gtools::combinations(3, 2, compare_clusts)
compare_clusts <- as.data.table(as.data.frame(compare_clusts))
setnames(compare_clusts, old = 1:2, new = sprintf("clust_%s", 1:2))

## Subset samples to analyze.

seurat_subset <- subset(
	seurat_integrated, subset =
		orig.ident %in% c("COLON_1", "HT29_EV", "H508_EV") &
		custom_clusters %in% c("early EEC", "EEC", "goblet")
)

Idents(seurat_subset) <- "custom_clusters"

## Differential expression comparisons.

options(future.globals.maxSize = 10000 * 1024 ^2)
plan("multiprocess", workers = 4)

diff_markers <- pmap(compare_clusts, function(clust_1, clust_2) {
	diff_markers <- FindMarkers(
		seurat_subset, assay = "SCT", slot = "data", logfc.threshold = log(1.5),
		min.pct = 0.25, ident.1 = clust_1, ident.2 = clust_2
	)
	return(as.data.table(diff_markers, keep.rownames = "gene"))
})

compare_clusts[, 
	comp_name := str_c(clust_1, clust_2, sep = "_vs_"),
	by = seq_len(nrow(compare_clusts))
]

names(diff_markers) <- compare_clusts[["comp_name"]]
names(diff_markers) <- str_replace(names(diff_markers), " ", "_")

## Format table for output.

diff_markers <- rbindlist(diff_markers, idcol = "comparison")

diff_markers[, avg_log2FC := log2(exp(avg_logFC))]
diff_markers <- diff_markers[abs(avg_log2FC) > log2(1.5) & p_val_adj < 0.05]
setorder(diff_markers, comparison, -avg_log2FC)

## Save table.

fwrite(
	diff_markers, file.path("results", "custom_clusters", "diff_markers_table.tsv"),
	sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

## Go Analysis
## ----------

## Go enrichment.

diff_markers[, change := fifelse(avg_log2FC < 0, "down", "up")]

go_enrichment <- compareCluster(
	gene ~ comparison + change, data = diff_markers, fun = "enrichGO",
	OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP", qvalueCutoff = 0.05
)

simple_go <- simplify(go_enrichment)

## Export table of GO results.

go_table <- as.data.table(go_enrichment)

fwrite(
	go_table,
	file.path("results", "custom_clusters", "enrichment", "diff_markers_GO_table.tsv"),
	sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

simple_go_table <- as.data.table(simple_go)

fwrite(
        simple_go_table,
	file.path("results", "custom_clusters", "enrichment", "diff_markers_simple_GO_table.tsv"),
        sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

## Plot GO results.

p <- dotplot(go_enrichment, x = ~change, showCategory = 10) +
	facet_grid(~comparison)

ggsave(
	file.path("results", "custom_clusters", "enrichment", "diff_markers_GO_dotplot.pdf"),
	plot = p, device = cairo_pdf, height = 8, width = 16
)

p <- dotplot(simple_go, x = ~change, showCategory = 10) +
        facet_grid(~comparison)

ggsave(
        file.path("results", "custom_clusters", "enrichment", "diff_markers_simple_GO_dotplot.pdf"),
        plot = p, device = cairo_pdf, height = 8, width = 16
)

## Pathway Analysis
## ----------

## Change to Entrez IDs.

entrez <- bitr(
	diff_markers[["gene"]], fromType = "SYMBOL", toType = "ENTREZID",
	OrgDb = "org.Hs.eg.db"
)
setDT(entrez)
setnames(entrez, old = "SYMBOL", new = "gene")

entrez <- merge(entrez, diff_markers, by = "gene")

## Reactome pathway.

rp_enrichment <- compareCluster(
	ENTREZID ~ comparison + change, data = entrez, fun = "enrichPathway",
	organism = "human", qvalueCutoff = 0.05, readable = TRUE	
)

## Export table of Reactone DB results.

rp_table <- as.data.table(rp_enrichment)

fwrite(
	rp_table, file.path("results", "custom_clusters", "enrichment", "diff_markers_pathway_table.tsv"),
	sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

## Plot pathway results.

p <- dotplot(rp_enrichment, x = ~change, showCategory = 10) +
        facet_grid(~comparison)

ggsave(
        file.path("results", "custom_clusters", "enrichment", "diff_markers_pathway_dotplot.pdf"),
        plot = p, device = cairo_pdf, height = 8, width = 16
)

##################################
## ATOH1 and NEUROG3 Expression ##
##################################

## Prepare the counts.

counts <- Assays(seurat_integrated, "SCT")@counts
counts <- counts[which(rownames(counts) %in% c("ATOH1", "NEUROG3")), ]
counts <- as.data.table(t(as.matrix(counts)), keep.rownames = "cell_id")
setkey(counts, "cell_id")

## Prepare meta-data.

metadata <- as.data.table(seurat_integrated[[]], keep.rownames = "cell_id")
metadata <- metadata[, .(cell_id, orig.ident, custom_clusters)]
setkey(metadata, "cell_id")

## Merge the meta-data into the counts.

counts <- counts[metadata]

## Annotate expression status.

counts[, status := case_when(
	ATOH1 > 0 & NEUROG3 > 0 ~ "NEUROG3+_ATOH1+",
	ATOH1 > 0 & NEUROG3 == 0 ~ "NEUROG-_ATOH1+",
	ATOH1 == 0 & NEUROG3 > 0 ~ "NEUROG+_ATOH1-",
	ATOH1 == 0 & NEUROG3 == 0 ~ "NEUROG3-_ATOH1-" 
)]

counts <- counts[, .(status_count = .N), by = .(orig.ident, status)]
counts <- counts[
	CJ(orig.ident = unique(counts$orig.ident), status = unique(counts$status)),
	on = .(orig.ident, status)
]
counts[is.na(status_count), status_count := 0]
setorder(counts, orig.ident, status)

counts[, status_frac := status_count / sum(status_count), by = orig.ident]

## Export to table.

fwrite(
	counts, file.path("results", "custom_clusters", "ATOH1_vs_NEUROG3.tsv"),
	sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

##############################
## Selected Gene Expression ##
##############################

select_genes <- c(
	"NEUROG3", "PAX4", "ARX", "RCOR2",
	"NEUROD1", "CHGA", "PYY", "FEV"
)

## Prepare the counts.

counts <- Assays(seurat_integrated, "SCT")@data
counts <- counts[which(rownames(counts) %in% select_genes), ]
counts <- as.data.table(t(as.matrix(counts)), keep.rownames = "cell_id")
setkey(counts, "cell_id")

## Prepare meta-data.

metadata <- as.data.table(seurat_integrated[[]], keep.rownames = "cell_id")
metadata <- metadata[, .(cell_id, orig.ident, custom_clusters)]
setkey(metadata, "cell_id")

## Merge the meta-data into the counts.

counts <- counts[metadata]

## Keep only COLON_1, HT29_EV, and H508_EV.

counts <- counts[orig.ident %in% c("COLON_1", "HT29_EV", "H508_EV")]

## Melting the counts.

counts <- melt(
	counts, id.vars = c("cell_id", "orig.ident", "custom_clusters"),
	measure.vars = select_genes, variable.name = "gene",
	value.name = "Log2_Normalized_Counts"
)

## Getting the fraction of cells per sample that express the gene.

counts[, exp_status := Log2_Normalized_Counts > 0]
counts <- counts[,
	.(perc_exp = (sum(exp_status) / .N) * 100),
	by = .(orig.ident, gene)
]

## Setting the correct factor order for the genes and sample names.

counts[, c("gene", "orig.ident") := list(
	factor(gene, levels = select_genes),
	factor(orig.ident, levels = c("HT29_EV", "H508_EV", "COLON_1"))
)]

## Generating the bar plot.

wes_colors <- as.character(wes_palette("Zissou1", 3, "continuous"))

p <- ggplot(counts, aes(x = orig.ident, y = perc_exp)) +
	geom_col(aes(fill = orig.ident), width = 0.75) +
	theme_minimal() +
	facet_wrap(~ gene, scales = "free", ncol = 1) +
	theme(axis.text.x = element_blank()) +
	scale_fill_manual(values = wes_colors)

pdf(file.path("results", "custom_clusters", "select_genes_barplot.pdf"), height = 12, width = 4)
p; dev.off()
