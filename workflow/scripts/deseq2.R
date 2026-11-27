set.seed(1)
suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
})

counts_path <- snakemake@input[["counts"]]
samples_path <- snakemake@input[["samples"]]
annotation_path <- snakemake@input[["annotation"]]
patient_column <- snakemake@params[["patient_column"]]
condition_column <- snakemake@params[["condition_column"]]
reference <- snakemake@params[["reference"]]
min_count <- as.integer(snakemake@params[["min_count"]])
min_samples <- as.integer(snakemake@params[["min_samples"]])
padj_threshold <- as.numeric(snakemake@params[["padj"]])
lfc_threshold <- as.numeric(snakemake@params[["lfc"]])

counts_table <- read.delim(counts_path, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(counts_table) < 2 || colnames(counts_table)[[1]] != "gene_id") {
  stop("The count matrix must start with a gene_id column and contain sample columns.", call. = FALSE)
}
if (anyDuplicated(counts_table$gene_id)) {
  stop("The count matrix contains duplicate gene_id values.", call. = FALSE)
}
if (anyDuplicated(colnames(counts_table)[-1])) {
  stop("The count matrix contains duplicate sample column names.", call. = FALSE)
}

sample_table <- read.delim(samples_path, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c("sample", patient_column, condition_column)
missing_columns <- setdiff(required_columns, colnames(sample_table))
if (length(missing_columns) > 0) {
  stop(
    paste("The DE sample sheet is missing required columns:", paste(missing_columns, collapse = ", ")),
    call. = FALSE
  )
}
if (anyDuplicated(sample_table$sample)) {
  stop("The DE sample sheet contains duplicate sample names.", call. = FALSE)
}

matrix_samples <- colnames(counts_table)[-1]
sheet_samples <- as.character(sample_table$sample)
missing_samples <- setdiff(sheet_samples, matrix_samples)
extra_samples <- setdiff(matrix_samples, sheet_samples)
if (length(missing_samples) > 0 || length(extra_samples) > 0) {
  stop(
    paste0(
      "Count matrix/sample sheet mismatch. Missing matrix columns: ",
      if (length(missing_samples)) paste(missing_samples, collapse = ", ") else "none",
      "; extra matrix columns: ",
      if (length(extra_samples)) paste(extra_samples, collapse = ", ") else "none",
      ". Match sample names exactly."
    ),
    call. = FALSE
  )
}

count_values <- as.matrix(counts_table[, sheet_samples, drop = FALSE])
storage.mode(count_values) <- "numeric"
if (anyNA(count_values) || any(!is.finite(count_values)) || any(count_values < 0)) {
  stop("The count matrix must contain finite, non-negative counts.", call. = FALSE)
}
if (any(count_values != floor(count_values))) {
  stop("DESeq2 requires integer counts, but the matrix contains non-integer values.", call. = FALSE)
}
storage.mode(count_values) <- "integer"
rownames(count_values) <- counts_table$gene_id

rownames(sample_table) <- as.character(sample_table$sample)
sample_table[[patient_column]] <- factor(sample_table[[patient_column]])
sample_table[[condition_column]] <- factor(sample_table[[condition_column]])
if (!reference %in% levels(sample_table[[condition_column]])) {
  stop(paste0("Reference condition '", reference, "' is absent from the DE sample sheet."), call. = FALSE)
}
if (!all(c("normal", "tumor") %in% levels(sample_table[[condition_column]]))) {
  stop("The DE sample sheet must contain both 'normal' and 'tumor' conditions.", call. = FALSE)
}
sample_table[[condition_column]] <- relevel(sample_table[[condition_column]], ref = reference)

design_formula <- stats::reformulate(c(patient_column, condition_column))
dds <- DESeqDataSetFromMatrix(
  countData = count_values,
  colData = sample_table,
  design = design_formula
)
keep <- rowSums(counts(dds) >= min_count) >= min_samples
dds <- dds[keep, ]
if (nrow(dds) == 0) {
  stop("No genes passed the configured minimum-count filter.", call. = FALSE)
}
dds <- DESeq(dds, quiet = TRUE)

non_reference_levels <- setdiff(levels(dds[[condition_column]]), reference)
if (length(non_reference_levels) != 1) {
  stop(
    paste0(
      "Expected exactly one non-reference condition level, found: ",
      paste(non_reference_levels, collapse = ", ")
    ),
    call. = FALSE
  )
}
coefficient <- paste0(
  condition_column,
  "_",
  non_reference_levels[[1]],
  "_vs_",
  reference
)
if (!coefficient %in% resultsNames(dds)) {
  stop(
    paste0(
      "Could not find tumor-versus-reference coefficient '",
      coefficient,
      "'. Available: ",
      paste(resultsNames(dds), collapse = ", ")
    ),
    call. = FALSE
  )
}
results_table <- results(dds, name = coefficient, independentFiltering = FALSE)
shrunken <- lfcShrink(dds, coef = coefficient, type = "apeglm")

annotation <- read.delim(annotation_path, check.names = FALSE, stringsAsFactors = FALSE)
annotation_columns <- c("gene_id", "gene_name", "gene_biotype")
if (!all(annotation_columns %in% colnames(annotation))) {
  stop("Gene annotation must contain gene_id, gene_name, and gene_biotype columns.", call. = FALSE)
}
if (anyDuplicated(annotation$gene_id)) {
  stop("Gene annotation contains duplicate gene_id values.", call. = FALSE)
}
annotation_index <- match(rownames(results_table), annotation$gene_id)
gene_names <- annotation$gene_name[annotation_index]
gene_biotypes <- annotation$gene_biotype[annotation_index]
missing_names <- is.na(gene_names) | gene_names == ""
gene_names[missing_names] <- rownames(results_table)[missing_names]

unshrunken_lfc <- results_table$log2FoldChange
significant <- !is.na(results_table$padj) &
  results_table$padj < padj_threshold &
  !is.na(unshrunken_lfc) &
  abs(unshrunken_lfc) >= lfc_threshold
direction <- ifelse(
  significant & unshrunken_lfc > 0,
  "up",
  ifelse(significant, "down", "not_significant")
)
all_genes <- data.frame(
  gene_id = rownames(results_table),
  gene_name = gene_names,
  gene_biotype = gene_biotypes,
  baseMean = results_table$baseMean,
  log2FoldChange = unshrunken_lfc,
  lfcSE = results_table$lfcSE,
  log2FoldChange_shrunk = shrunken$log2FoldChange[match(rownames(results_table), rownames(shrunken))],
  pvalue = results_table$pvalue,
  padj = results_table$padj,
  significant = significant,
  direction = direction,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
degs <- all_genes[all_genes$significant, , drop = FALSE]
degs <- degs[order(degs$padj, na.last = TRUE), , drop = FALSE]

output_paths <- snakemake@output
dir.create(dirname(output_paths[["all_genes"]]), recursive = TRUE, showWarnings = FALSE)
write.table(all_genes, output_paths[["all_genes"]], sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(degs, output_paths[["degs"]], sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
summary <- data.frame(
  genes_tested = nrow(all_genes),
  degs = nrow(degs),
  up = sum(degs$direction == "up"),
  down = sum(degs$direction == "down"),
  min_count = min_count,
  min_samples = min_samples,
  padj_threshold = padj_threshold,
  lfc_threshold = lfc_threshold
)
write.table(summary, output_paths[["summary"]], sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(dds, output_paths[["dds"]])
saveRDS(vst(dds, blind = FALSE), output_paths[["vst"]])
writeLines(capture.output(sessionInfo()), output_paths[["session_info"]])
