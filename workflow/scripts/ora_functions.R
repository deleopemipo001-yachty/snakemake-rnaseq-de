# Shared hypergeometric ORA implementation used by the DESeq2 and preprint rules.
normalise_gene_id <- function(gene_id) {
  sub("\\.[0-9]+$", "", as.character(gene_id))
}

run_ora <- function(query, gene_sets, universe, annotation_names,
                    minimum_set_size, maximum_set_size) {
  if (nrow(gene_sets) > 0) {
    sizes <- vapply(
      split(gene_sets$gene, gene_sets$ID),
      function(genes) length(unique(genes)),
      integer(1)
    )
    eligible_ids <- names(sizes)[
      sizes >= minimum_set_size & sizes <= maximum_set_size
    ]
    gene_sets <- gene_sets[gene_sets$ID %in% eligible_ids, , drop = FALSE]
  }
  collection_universe <- intersect(universe, unique(gene_sets$gene))
  query <- unique(normalise_gene_id(query))
  query <- intersect(query[!is.na(query) & query != ""], collection_universe)
  number_tested <- length(query)
  number_universe <- length(collection_universe)
  if (number_tested == 0 || nrow(gene_sets) == 0) {
    return(data.frame(
      ID = character(), Description = character(),
      Long_description = character(), GeneRatio = character(),
      BgRatio = character(), pvalue = numeric(), geneID = character(),
      Count = integer(), p.adjust = numeric(), stringsAsFactors = FALSE
    ))
  }
  term_ids <- unique(gene_sets$ID)
  rows <- lapply(term_ids, function(term_id) {
    term <- gene_sets[gene_sets$ID == term_id, , drop = FALSE]
    term_genes <- unique(term$gene)
    set_size <- length(term_genes)
    overlap <- intersect(query, term_genes)
    overlap_count <- length(overlap)
    if (overlap_count == 0) return(NULL)
    pvalue <- stats::phyper(
      overlap_count - 1, set_size, number_universe - set_size,
      number_tested, lower.tail = FALSE
    )
    gene_names <- annotation_names$gene_name[
      match(overlap, annotation_names$gene_id)
    ]
    missing_names <- is.na(gene_names) | gene_names == ""
    gene_names[missing_names] <- overlap[missing_names]
    data.frame(
      ID = term_id,
      Description = term$Description[[1]],
      Long_description = term$Long_description[[1]],
      GeneRatio = paste0(overlap_count, "/", number_tested),
      BgRatio = paste0(set_size, "/", number_universe),
      pvalue = pvalue,
      geneID = paste(gene_names, collapse = "/"),
      Count = as.integer(overlap_count),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    out <- data.frame(
      ID = character(), Description = character(),
      Long_description = character(), GeneRatio = character(),
      BgRatio = character(), pvalue = numeric(), geneID = character(),
      Count = integer(), stringsAsFactors = FALSE
    )
    out$p.adjust <- numeric()
    return(out)
  }
  tested_terms <- do.call(rbind, rows)
  tested_terms$p.adjust <- stats::p.adjust(tested_terms$pvalue, method = "BH")
  tested_terms <- tested_terms[
    order(tested_terms$p.adjust, tested_terms$pvalue), , drop = FALSE
  ]
  tested_terms
}

make_ora_gene_sets <- function(data, universe) {
  term_id <- if ("gs_id" %in% colnames(data)) {
    as.character(data$gs_id)
  } else {
    as.character(data$gs_name)
  }
  description <- if ("gs_description" %in% colnames(data)) {
    as.character(data$gs_description)
  } else {
    as.character(data$gs_name)
  }
  long_description <- description
  description <- as.character(data$gs_name)
  missing_long_description <- is.na(long_description) | long_description == ""
  long_description[missing_long_description] <- description[missing_long_description]
  mapping <- data.frame(
    ID = term_id,
    Description = description,
    Long_description = long_description,
    gene = as.character(data$analysis_gene_id),
    stringsAsFactors = FALSE
  )
  mapping <- mapping[
    !is.na(mapping$ID) & mapping$ID != "" &
      !is.na(mapping$gene) & mapping$gene != "" &
      mapping$gene %in% universe,
    , drop = FALSE
  ]
  dplyr::distinct(mapping)
}
