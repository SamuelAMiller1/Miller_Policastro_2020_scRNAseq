
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
