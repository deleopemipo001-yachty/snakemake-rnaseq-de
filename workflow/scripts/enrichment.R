set.seed(1)
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(msigdbr)
})
source(snakemake@input[["ora_functions"]])

all_genes <- read.delim(
  snakemake@input[["all_genes"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
degs <- read.delim(
  snakemake@input[["degs"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
annotation <- read.delim(
  snakemake@input[["annotation"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
padj_threshold <- as.numeric(snakemake@params[["padj"]])
minimum_set_size <- as.integer(snakemake@params[["min_size"]])
maximum_set_size <- as.integer(snakemake@params[["max_size"]])

if (!"gene_id" %in% colnames(all_genes)) {
  stop("The DE results table must contain gene_id.", call. = FALSE)
}
if (!all(c("gene_id", "direction") %in% colnames(degs))) {
  stop("The DEG table must contain gene_id and direction.", call. = FALSE)
}
if (!all(c("gene_id", "gene_name") %in% colnames(annotation))) {
  stop("Gene annotation must contain gene_id and gene_name.", call. = FALSE)
}

universe <- unique(normalise_gene_id(all_genes$gene_id))
universe <- universe[!is.na(universe) & universe != ""]
if (length(universe) == 0) {
  stop("The DE results contain no tested genes for the enrichment universe.", call. = FALSE)
}

msig <- msigdbr::msigdbr(species = "Homo sapiens")
required_msig_columns <- c("gs_name", "gs_cat", "gs_subcat", "gene_symbol")
if (!all(required_msig_columns %in% colnames(msig))) {
  stop("The pinned msigdbr data is missing expected gene-set columns.", call. = FALSE)
}

package_version <- as.character(utils::packageVersion("msigdbr"))
if ("db_version" %in% colnames(msig)) {
  database_version <- unique(as.character(msig$db_version))
  database_version <- database_version[!is.na(database_version) & database_version != ""]
  database_version <- paste(database_version, collapse = ", ")
} else {
  database_version <- paste0(package_version, " (bundled with msigdbr)")
}
if (database_version == "") {
  database_version <- paste0(package_version, " (bundled with msigdbr)")
}

annotation <- annotation[
  !is.na(annotation$gene_id) & annotation$gene_id != "",
  ,
  drop = FALSE
]
annotation$gene_id <- normalise_gene_id(annotation$gene_id)
annotation_names <- annotation[, c("gene_id", "gene_name"), drop = FALSE]
annotation_names <- annotation_names[!duplicated(annotation_names$gene_id), , drop = FALSE]

use_ensembl_ids <- "ensembl_gene" %in% colnames(msig) &&
  any(!is.na(msig$ensembl_gene) & msig$ensembl_gene != "")
if (use_ensembl_ids) {
  msig$analysis_gene_id <- normalise_gene_id(msig$ensembl_gene)
} else {
  symbol_map <- annotation[
    !is.na(annotation$gene_name) & annotation$gene_name != "",
    c("gene_name", "gene_id"),
    drop = FALSE
  ]
  msig$analysis_gene_id <- symbol_map$gene_id[
    match(msig$gene_symbol, symbol_map$gene_name)
  ]
}

empty_table <- function() {
  data.frame(
    ID = character(),
    Description = character(),
    Long_description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    p.adjust = numeric(),
    geneID = character(),
    Count = integer(),
    stringsAsFactors = FALSE
  )
}

save_plot_pair <- function(plot, png_path, pdf_path) {
  ggplot2::ggsave(
    filename = png_path,
    plot = plot,
    width = 10,
    height = 7,
    dpi = 300,
    bg = "white"
  )
  ggplot2::ggsave(
    filename = pdf_path,
    plot = plot,
    width = 10,
    height = 7,
    bg = "white"
  )
}

run_one <- function(collection, direction) {
  category_data <- if (collection == "kegg") {
    msig[
      msig$gs_cat == "C2" &
        msig$gs_subcat %in% c("CP:KEGG", "CP:KEGG_LEGACY"),
      ,
      drop = FALSE
    ]
  } else {
    msig[msig$gs_cat == "H", , drop = FALSE]
  }
  gene_sets <- make_ora_gene_sets(category_data, universe)
  query <- degs$gene_id[degs$direction == direction]
  tested_terms <- run_ora(
    query,
    gene_sets,
    universe,
    annotation_names,
    minimum_set_size,
    maximum_set_size
  )

  significant_results <- empty_table()
  if (nrow(tested_terms) > 0) {
    tested_terms$significant <- !is.na(tested_terms$p.adjust) & tested_terms$p.adjust < padj_threshold
    significant_results <- tested_terms[
      tested_terms$significant,
      c("ID", "Description", "Long_description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "geneID", "Count", "significant"),
      drop = FALSE
    ]
    tested_terms$significant <- !is.na(tested_terms$p.adjust) & tested_terms$p.adjust < padj_threshold
  } else {
    tested_terms$significant <- logical()
  }

  table_path <- snakemake@output[[paste0(collection, "_", direction)]]
  write.table(
    tested_terms,
    table_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )

  plot_name <- paste0(collection, "_", direction)
  png_path <- snakemake@output[[paste0(plot_name, "_png")]]
  pdf_path <- snakemake@output[[paste0(plot_name, "_pdf")]]
  collection_title <- if (collection == "kegg") "KEGG" else "Hallmark"
  direction_title <- if (direction == "up") "up-regulated" else "down-regulated"
  title <- paste(collection_title, direction_title, "pathways")
  if (nrow(significant_results) == 0) {
    plot <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No significant terms", size = 6) +
      ggplot2::xlim(-1, 1) +
      ggplot2::ylim(-1, 1) +
      ggplot2::labs(title = title) +
      ggplot2::theme_void()
  } else {
    plot_data <- head(significant_results, 15)
    ratio_parts <- strsplit(plot_data$GeneRatio, "/", fixed = TRUE)
    plot_data$GeneRatioValue <- vapply(
      ratio_parts,
      function(ratio) as.numeric(ratio[[1]]) / as.numeric(ratio[[2]]),
      numeric(1)
    )
    plot_data$Description <- factor(
      plot_data$Description,
      levels = rev(plot_data$Description)
    )
    plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = GeneRatioValue,
        y = Description,
        size = Count,
        colour = p.adjust
      )
    ) +
      ggplot2::geom_point() +
      ggplot2::labs(
        title = title,
        x = "GeneRatio",
        y = NULL,
        size = "Count",
        colour = "Adjusted p-value"
      ) +
      ggplot2::theme_bw()
  }
  save_plot_pair(plot, png_path, pdf_path)
}

dir.create(dirname(snakemake@output[["kegg_up"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(snakemake@output[["kegg_up_png"]]), recursive = TRUE, showWarnings = FALSE)
for (collection in c("kegg", "hallmark")) {
  for (direction in c("up", "down")) {
    run_one(collection, direction)
  }
}

versions <- data.frame(
  item = c("msigdbr_package_version", "gene_set_database_version"),
  version = c(package_version, database_version),
  stringsAsFactors = FALSE
)
write.table(
  versions,
  snakemake@output[["versions"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
writeLines(capture.output(sessionInfo()), snakemake@output[["session_info"]])
