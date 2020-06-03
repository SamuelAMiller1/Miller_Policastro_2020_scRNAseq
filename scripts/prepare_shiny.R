
library("data.table")
library("tidyverse")
library("Seurat")
library("Matrix")

############################
## Prepare Data for Shiny ##
############################

if (!dir.exists(file.path("shiny_app", "datasets"))) {
	dir.create(file.path("shiny_app", "datasets"))
}

## Load seurat object.

seurat_obj <- readRDS(file.path("results", "r_objects", "seurat_expanded.RDS"))

## Prepare Gene Counts
## ----------

## Extract counts.

gene_counts <- as.data.table(
	Assays(seurat_obj, "SCT")@data,
	keep.rownames = "gene"
)

## Make cell ID's the rows, and genes the columns.

gene_counts <- melt(
	gene_counts, id.vars = "gene",
	variable.name = "cell_id",
	value.name = "normalized_counts"
)
gene_counts <- dcast(gene_counts, cell_id ~ gene)

## Save to file.

fwrite(
	gene_counts, file.path("shiny_app", "datasets", "gene_counts.tsv"),
	sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

## Prepare Meta-Data.
## ----------

## Extract meta-data.

meta_data <- as.data.table(seurat_obj[[]], keep.rownames = "cell_id")

## Save to file.

fwrite(
	meta_data, file.path("shiny_app", "datasets", "meta_data.tsv"),
	sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

## Prepare Markers.
## ----------

## Load marker files.

markers <- fread(file.path("results", "custom_clusters", "marker_table.tsv"))

## Change column order.

markers <- markers[,
	.(cluster, gene, pct.1, pct.2, p_val,
	p_val_adj, avg_logFC, avg_log2FC)
]

## Save to file.

fwrite(
	markers, file.path("shiny_app", "datasets", "markers.tsv"),
	sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

## Prepare UMAP
## ----------

## Extract UMAP.

reductions <- as.data.table(Embeddings(seurat_obj, "umap"), keep.rownames = "cell_id")

## Save UMAP embeddings.

fwrite(
	reductions, file.path("shiny_app", "datasets", "umap_embeddings.tsv"),
	sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)
