
## Loaad singularity container.
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

## Dll1 Analysis
## ----------

## Get expression data for Dll1.

expression_data <- as.matrix(Assays(seurat_integrated, "RNA")@counts)
expression_data <- as.data.table(expression_data, keep.rownames = "gene_id")[
	gene_id == "DLL1"
]
expression_data <- melt(
	expression_data, id.vars = "gene_id",
	variable.name = "cell_id", value.name = "UMI_count"
)

## Get metadata for each cell.

meta_data <- as.data.table(seurat_integrated@meta.data, keep.rownames = "cell_id")[,
	.(cell_id, orig.ident)
]

## Merge expression data into metadata.

setkey(expression_data, cell_id)
setkey(meta_data, cell_id)

merged <- merge(expression_data, meta_data)

## Mark Dll1 expression for each cell.

merged[, DLL1_status := ifelse(UMI_count > 0, "positive", "negative")]

## Plot Dll1 status.

if (!dir.exists(file.path("results", "gene_plots"))) {
	dir.create(file.path("results", "gene_plots"))
}

color_palette <- as.character(wes_palette("Zissou1", 2, "continuous"))

p <- ggplot(merged, aes(x = orig.ident, fill = DLL1_status)) +
	geom_bar(position = "fill") +
	scale_fill_manual(values = color_palette) +
	theme_bw() +
	

pdf(file.path("results", "gene_plots", "DLL1_status.pdf"), height = 4, width = 6)
p
dev.off()

## RCOR1 & 2 Expression
## ----------

## Get expression data for RCOR1, RCOR2, and DLL1.

exp_data <- as.matrix(Assays(seurat_integrated, "RNA")@counts)
exp_data <- as.data.table(exp_data, keep.rownames = "gene_id")[
	gene_id %in% c("DLL1", "RCOR1", "RCOR2", "NEUROG3")
]

exp_data <- melt(
	exp_data, id.vars = "gene_id",
	variable.name = "cell_id", value.name = "UMI_count"
)

## Get meta-data.

meta_data <- as.data.table(seurat_integrated@meta.data, keep.rownames = "cell_id")[,
	.(cell_id, orig.ident)
]

## Merge the counts into the meta-data and keep only HT29_EV.

exp_data <- merge(exp_data, meta_data, all.x = TRUE, by = "cell_id")[
	orig.ident == "HT29_EV"
]

## Classify the expression status of the cells.

exp_data <- dcast(exp_data, cell_id + orig.ident ~ gene_id, value.var = "UMI_count")
exp_data <- melt(
	exp_data, measure.vars = c("RCOR1", "RCOR2"),
	variable.name = "RCOR_gene", value.name = "RCOR_exp"
)

exp_data[, DLL1_NEUROG3_status := case_when(
	DLL1 == 0 & NEUROG3 == 0 ~ "DLL1(-)/NEUROG3(-)",
	DLL1 > 0 & NEUROG3 == 0 ~ "DLL1(+)/NEUROG3(-)",
	DLL1 == 0 & NEUROG3 > 0 ~ "DLL1(-)/NEUROG3(+)",
	DLL1 > 0 & NEUROG3 > 0 ~ "DLL1(+)/NEUROG(+)"
)]

exp_data[,
	DLL1_NEUROG3_status := factor(DLL1_NEUROG3_status, levels = c(
		"DLL1(+)/NEUROG(+)", "DLL1(+)/NEUROG3(-)",
		"DLL1(-)/NEUROG3(+)", "DLL1(-)/NEUROG3(-)"
	))
]

## Plot the expression values.

p <- ggplot(exp_data, aes(x = DLL1_NEUROG3_status, y = RCOR_exp, fill = DLL1_NEUROG3_status)) +
	geom_jitter(size = 0.1, color = "grey", width = 0.25) +
	geom_boxplot(outlier.shape = NA, width = 0.25) +
	theme_bw() +
	theme(
		axis.text.x = element_blank(),
		text = element_text(size = 8)
	) +
	scale_fill_manual(values = as.character(wes_palette("Zissou1", 4, "continuous"))) +
	facet_wrap(. ~ RCOR_gene, scales = "free")

DLL1_NEUROG3_file <- file.path("results", "gene_plots", "DLL1_NEUROG3_RCOR.pdf")
pdf(DLL1_NEUROG3_file, height = 3, width = 8); p; dev.off()
