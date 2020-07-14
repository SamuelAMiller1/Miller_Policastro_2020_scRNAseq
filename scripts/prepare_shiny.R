
library("data.table")
library("tidyverse")
library("Seurat")
library("DBI")
library("RSQLite")


############################
## Prepare Data for Shiny ##
############################

if (!dir.exists(file.path("shiny_app", "datasets"))) {
  dir.create(file.path("shiny_app", "datasets"))
}

## Load seurat object.

seu <- list(
  total = file.path("results", "r_objects", "seurat_expanded.RDS"),
  spliced = file.path("results", "r_objects", "seurat_integrated_spliced.RDS")
)

seu <- map(seu, readRDS)

## Connect to database.

con <- dbConnect(SQLite(), file.path("shiny_app", "datasets", "miller.sqlite"))

## Prepare Gene Counts
## ----------

## Extract counts.

counts <- imap(seu, function(x, y) {
  counts <- x[["SCT"]]@counts
  counts <- as.data.table(t(as.matrix(counts)), keep.rownames = "cell_id")

  counts <- melt(
    counts, id.vars = "cell_id", variable.name = "gene",
    value.name = "exp"
  )

  return(counts)
})

## Save to database.

iwalk(counts, function(x, y) {
  copy_to(
    con, x, str_c(y, "_counts"), temporary = FALSE, overwrite = TRUE,
    indexes = list("gene")
  )
})

## Prepare Meta-Data.
## ----------

metadata <- imap(seu, function(x, y) {
  meta_data <- as.data.table(x[[]], keep.rownames = "cell_id")
  if (y == "total") {
    meta_data[, seurat_clusters := custom_clusters]
  } else if (y == "spliced") {
    meta_data[, seurat_clusters := integrated_snn_res.0.6]
  }
  return(meta_data)
})

## Save to database.

iwalk(metadata, function(x, y) {
  copy_to(
    con, x, str_c(y, "_metadata"), temporary = FALSE, overwrite = TRUE,
    indexes = list("seurat_clusters", "orig.ident")
  )
})

## Prepare Markers.
## ----------

library("future")

options(future.globals.maxSize = 10000 * 1024 ^2)
plan("multiprocess", workers = 3)

## Find markers.

markers <- imap(seu, function(x, y) {
  if (y == "total") {
    Idents(x) <- "custom_clusters"
  } else if (y == "spliced") {
    Idents(x) <- "integrated_snn_res.0.6"
  }

  markers <- FindAllMarkers(
    x, assay = "SCT", slot = "data", logfc.threshold = log(1.5),
    min.pct = 0.25, return.thresh = 0.05
  )

  setDT(markers)
  markers[, avg_log2FC := log2(exp(avg_logFC))]
  markers <- markers[p_val_adj < 0.05]
  markers <- markers[order(cluster, p_val_adj)]
  markers <- markers[,
    .(cluster, gene, pct.1, pct.2, p_val,
    p_val_adj, avg_logFC, avg_log2FC)
  ]

  return(markers)
})

## Write out file.

fwrite(
  rbindlist(markers, idcol = "experiment"),
  file.path("shiny_app", "datasets", "markers.tsv"),
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

## Save to database.

iwalk(markers, function(x, y) {
  copy_to(
    con, x, str_c(y, "_markers"), temporary = FALSE, overwrite = TRUE,
    indexes = list("cluster", "gene")
  )
})

## Prepare UMAP
## ----------

reductions <- map(seu, function(x) {
  reductions <- as.data.table(
    Embeddings(x, "umap"), keep.rownames = "cell_id",
    key = "cell_id"
  )

  return(reductions)
})

## Save to database.

iwalk(reductions, function(x, y) {
  copy_to(con, x, str_c(y, "_reductions"), temporary = FALSE, overwrite = TRUE)
})

## Prepare Enrichment
## ----------

library("clusterProfiler")
library("ReactomePA")
library("org.Hs.eg.db")

## prepare markers for analysis.

enr_markers <- map(markers, function(x) {
  x[, change := fifelse(avg_log2FC > 0, "up", "down")]
  x[,
        group := str_c("cluster", cluster, change, sep = "_"),
        by = seq_len(nrow(x))
  ]

  markers <- split(x, x[["group"]])
  return(markers)
})

## GO analysis.

go_enrichment <- map(enr_markers, function(x) {
  go_enrichment <- map(x, function(y) {
    enriched <- enrichGO(
      gene = y[["gene"]], OrgDb = "org.Hs.eg.db",
      keyType = "SYMBOL", ont = "BP"
    )
    enriched <- as.data.table(enriched)
    enriched <- enriched[p.adjust < 0.05]
    return(enriched)
  })

  go_enrichment <- rbindlist(go_enrichment, idcol = "group")
  return(go_enrichment)
})

go_enrichment <- rbindlist(go_enrichment, idcol = "experiment")

## Reactome analysis.

reactome_enrichment <- map(enr_markers, function(x) {
  enriched <- map(x, function(y) {
    entrez <- bitr(
      y[["gene"]], fromType = "SYMBOL", toType = "ENTREZID",
      OrgDb = "org.Hs.eg.db"
    )
    entrez <- entrez[["ENTREZID"]]

    enriched <- enrichPathway(gene = entrez, organism = "human", readable = TRUE)

    enriched <- as.data.table(enriched)
    enriched <- enriched[p.adjust < 0.05]
    return(enriched)
  })

  enriched <- rbindlist(enriched, idcol = "group")
  return(enriched)
})

reactome_enrichment <- rbindlist(reactome_enrichment, idcol = "experiment")

## Save file of enrichment.

fwrite(
  rbindlist(list(go_enrichment, reactome_enrichment), idcol = "database"),
  file.path("shiny_app", "datasets", "enriched.tsv"),
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

## Save to database.

enriched <- rbindlist(list(go_enrichment, reactome_enrichment), idcol = "database")

iwalk(enriched, function(x, y) {
  copy_to(
    con, x, str_c(y, "_enriched"), temporary = FALSE, overwrite = TRUE,
    indexes = list("database", "group")
  )
})

## Sample Sheet
## ----------

## Samples.

samples <- imap(seu, ~data.table(
  experiment = .y, samples = unique(.x$orig.ident)
))
samples <- rbindlist(samples)

copy_to(con, samples, "samples", temporary = FALSE, overwrite = TRUE)

## Turn of db connection.

dbDisconnect(con)

