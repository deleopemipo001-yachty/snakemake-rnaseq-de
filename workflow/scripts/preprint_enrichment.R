set.seed(1)
suppressPackageStartupMessages({
  library(dplyr)
  library(msigdbr)
})
source(snakemake@input[["ora_functions"]])

read_table <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}
preprint <- read_table(snakemake@input[["preprint"]])
preprint_degs <- read_table(snakemake@input[["preprint_degs"]])
deseq <- read_table(snakemake@input[["deseq"]])
deseq_degs <- read_table(snakemake@input[["deseq_degs"]])
annotation <- read_table(snakemake@input[["annotation"]])
padj_threshold <- as.numeric(snakemake@params[["padj"]])
minimum_set_size <- as.integer(snakemake@params[["min_size"]])
maximum_set_size <- as.integer(snakemake@params[["max_size"]])

preprint_ids <- if ("gene_id" %in% names(preprint)) preprint$gene_id else preprint$gene
preprint_de_ids <- if ("gene_id" %in% names(preprint_degs)) {
  preprint_degs$gene_id
} else {
  preprint_degs$gene
}
preprint_universe <- unique(normalise_gene_id(preprint_ids))
preprint_universe <- preprint_universe[
  !is.na(preprint_universe) & preprint_universe != ""
]
deseq_universe <- unique(normalise_gene_id(deseq$gene_id))
deseq_universe <- deseq_universe[!is.na(deseq_universe) & deseq_universe != ""]

msig <- msigdbr::msigdbr(species = "Homo sapiens")
annotation <- annotation[
  !is.na(annotation$gene_id) & annotation$gene_id != "",
  , drop = FALSE
]
annotation$gene_id <- normalise_gene_id(annotation$gene_id)
annotation_names <- annotation[, c("gene_id", "gene_name"), drop = FALSE]
annotation_names <- annotation_names[!duplicated(annotation_names$gene_id), , drop = FALSE]
if ("ensembl_gene" %in% names(msig) &&
    any(!is.na(msig$ensembl_gene) & msig$ensembl_gene != "")) {
  msig$analysis_gene_id <- normalise_gene_id(msig$ensembl_gene)
} else {
  symbol_map <- annotation[
    !is.na(annotation$gene_name) & annotation$gene_name != "",
    c("gene_name", "gene_id"), drop = FALSE
  ]
  msig$analysis_gene_id <- symbol_map$gene_id[
    match(msig$gene_symbol, symbol_map$gene_name)
  ]
}

collections <- list(
  KEGG = msig[
    msig$gs_cat == "C2" &
      msig$gs_subcat %in% c("CP:KEGG", "CP:KEGG_LEGACY"),
    , drop = FALSE
  ],
  HALLMARK = msig[msig$gs_cat == "H", , drop = FALSE]
)
run_direction <- function(collection_data, universe, query) {
  gene_sets <- make_ora_gene_sets(collection_data, universe)
  run_ora(
    query, gene_sets, universe, annotation_names,
    minimum_set_size, maximum_set_size
  )
}
pathway_padj <- function(table, pathway) {
  values <- table$p.adjust[table$Description == pathway]
  if (!length(values)) NA_real_ else min(values, na.rm = TRUE)
}

results <- list()
for (collection_name in names(collections)) {
  collection_data <- collections[[collection_name]]
  for (direction_name in c("up", "down")) {
    preprint_query <- if (direction_name == "up") {
      preprint_de_ids[preprint_degs$logFC > 0]
    } else {
      preprint_de_ids[preprint_degs$logFC < 0]
    }
    deseq_query <- deseq_degs$gene_id[deseq_degs$direction == direction_name]
    results[[paste(collection_name, direction_name, "preprint", sep = "_")]] <-
      run_direction(collection_data, preprint_universe, preprint_query)
    results[[paste(collection_name, direction_name, "DESeq2", sep = "_")]] <-
      run_direction(collection_data, deseq_universe, deseq_query)
  }
}

agreement_rows <- list()
for (collection_name in names(collections)) {
  for (direction_name in c("up", "down")) {
    preprint_result <- results[[paste(collection_name, direction_name, "preprint", sep = "_")]]
    deseq_result <- results[[paste(collection_name, direction_name, "DESeq2", sep = "_")]]
    pathways <- unique(c(preprint_result$Description, deseq_result$Description))
    if (!length(pathways)) next
    padj_preprint <- vapply(
      pathways, function(pathway) pathway_padj(preprint_result, pathway), numeric(1)
    )
    padj_deseq2 <- vapply(
      pathways, function(pathway) pathway_padj(deseq_result, pathway), numeric(1)
    )
    preprint_significant <- !is.na(padj_preprint) & padj_preprint < padj_threshold
    deseq_significant <- !is.na(padj_deseq2) & padj_deseq2 < padj_threshold
    pathway_status <- ifelse(
      preprint_significant & deseq_significant, "both",
      ifelse(preprint_significant, "preprint only",
             ifelse(deseq_significant, "DESeq2 only", "neither"))
    )
    agreement_rows[[length(agreement_rows) + 1L]] <- data.frame(
      collection = collection_name,
      direction = direction_name,
      pathway = pathways,
      padj_preprint = padj_preprint,
      padj_deseq2 = padj_deseq2,
      status = pathway_status,
      stringsAsFactors = FALSE
    )
  }
}
if (length(agreement_rows)) {
  agreement <- do.call(rbind, agreement_rows)
  agreement$collection_order <- match(agreement$collection, c("KEGG", "HALLMARK"))
  agreement$direction_order <- match(agreement$direction, c("up", "down"))
  agreement$sort_padj <- pmin(
    agreement$padj_preprint, agreement$padj_deseq2, na.rm = TRUE
  )
  agreement$sort_padj[!is.finite(agreement$sort_padj)] <- Inf
  agreement <- agreement[
    order(agreement$collection_order, agreement$direction_order,
          agreement$sort_padj, agreement$pathway),
    c("collection", "direction", "pathway", "padj_preprint", "padj_deseq2", "status"),
    drop = FALSE
  ]
} else {
  agreement <- data.frame(
    collection = character(), direction = character(), pathway = character(),
    padj_preprint = numeric(), padj_deseq2 = numeric(), status = character(),
    stringsAsFactors = FALSE
  )
}
write.table(
  agreement, snakemake@output[["agreement"]], sep = "\t",
  quote = FALSE, row.names = FALSE, na = "NA"
)

normalize_pathway_name <- function(pathway) {
  normalized <- trimws(tolower(as.character(pathway)))
  normalized <- sub("^kegg_", "", normalized)
  normalized <- gsub("[_-]", " ", normalized)
  gsub("\\s+", " ", trimws(normalized))
}
original_input <- snakemake@input[["original"]]
original_path <- if (length(original_input)) as.character(original_input[[1]]) else ""
original <- if (nzchar(original_path) && file.exists(original_path)) {
  tryCatch(
    read.csv(original_path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(error) data.frame()
  )
} else {
  data.frame()
}
find_column <- function(data, candidates) {
  found <- candidates[candidates %in% names(data)]
  if (length(found)) found[[1]] else NA_character_
}
original_name_col <- find_column(
  original, c("Description", "description", "pathway", "Pathway", "ID", "ID.1")
)
if (is.na(original_name_col) && ncol(original) > 0L) {
  original_name_col <- names(original)[[1]]
}
original_padj_col <- find_column(
  original, c("p.adjust", "padj", "adjusted_pvalue", "qvalue")
)
original_names <- if (!is.na(original_name_col)) {
  as.character(original[[original_name_col]])
} else {
  character()
}
original_padj <- if (!is.na(original_padj_col)) {
  suppressWarnings(as.numeric(original[[original_padj_col]]))
} else {
  rep(NA_real_, length(original_names))
}

# msigdbr 7.5.1 exposes the legacy collection under CP:KEGG on some builds and
# CP:KEGG_LEGACY on others. These are the same bundled 2011 legacy pathways.
kegg_names <- unique(as.character(collections$KEGG$gs_name))
kegg_normalized <- normalize_pathway_name(kegg_names)
original_normalized <- normalize_pathway_name(original_names)
matched_index <- match(original_normalized, kegg_normalized)
configured_aliases <- unlist(snakemake@params[["kegg_aliases"]], use.names = TRUE)
if (length(configured_aliases)) {
  alias_keys <- normalize_pathway_name(names(configured_aliases))
  alias_match <- match(original_normalized, alias_keys)
  has_alias <- !is.na(alias_match)
  if (any(has_alias)) {
    matched_index[has_alias] <- match(
      unname(configured_aliases[alias_match[has_alias]]), kegg_names
    )
  }
}
msigdb_names <- kegg_names[matched_index]
match_status <- ifelse(
  is.na(matched_index), "not in MSigDB KEGG legacy set", "matched"
)
get_result_padj <- function(direction, analysis, pathway) {
  result <- results[[paste("KEGG", direction, analysis, sep = "_")]]
  pathway_padj(result, pathway)
}
original_comparison <- data.frame(
  original_pathway = original_names,
  original_padj = original_padj,
  msigdb_name = msigdb_names,
  ours_preprint_up_padj = vapply(seq_along(original_names), function(index) {
    if (is.na(msigdb_names[[index]])) NA_real_ else
      get_result_padj("up", "preprint", msigdb_names[[index]])
  }, numeric(1)),
  ours_preprint_down_padj = vapply(seq_along(original_names), function(index) {
    if (is.na(msigdb_names[[index]])) NA_real_ else
      get_result_padj("down", "preprint", msigdb_names[[index]])
  }, numeric(1)),
  deseq2_up_padj = vapply(seq_along(original_names), function(index) {
    if (is.na(msigdb_names[[index]])) NA_real_ else
      get_result_padj("up", "DESeq2", msigdb_names[[index]])
  }, numeric(1)),
  deseq2_down_padj = vapply(seq_along(original_names), function(index) {
    if (is.na(msigdb_names[[index]])) NA_real_ else
      get_result_padj("down", "DESeq2", msigdb_names[[index]])
  }, numeric(1)),
  match_status = match_status,
  stringsAsFactors = FALSE
)
comparison_table <- capture.output(write.table(
  original_comparison, sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = TRUE, na = "NA"
))
writeLines(
  c(
    "# Original used online enrichKEGG (KEGG REST, June 2026) on combined up+down DEGs with the default universe; this comparison is descriptive.",
    "# MSigDB 7.5.1 KEGG is the 2011 legacy KEGG release; pathways added to KEGG later cannot be matched.",
    comparison_table
  ),
  snakemake@output[["original"]]
)
writeLines(capture.output(sessionInfo()), snakemake@output[["session_info"]])
