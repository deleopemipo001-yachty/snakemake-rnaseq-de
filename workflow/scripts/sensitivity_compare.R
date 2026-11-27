set.seed(1)
suppressPackageStartupMessages(library(ggplot2))

full_all <- read.delim(
  snakemake@input[["full_all_genes"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
full_degs <- read.delim(
  snakemake@input[["full_degs"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
sensitivity_all <- read.delim(
  snakemake@input[["sensitivity_all_genes"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
sensitivity_degs <- read.delim(
  snakemake@input[["sensitivity_degs"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_result_columns <- c("gene_id", "log2FoldChange", "direction")
if (!all(required_result_columns %in% colnames(full_all)) ||
    !all(required_result_columns %in% colnames(sensitivity_all))) {
  stop("Both DE result tables must contain gene_id, log2FoldChange, and direction.", call. = FALSE)
}
if (!"gene_id" %in% colnames(full_degs) || !"gene_id" %in% colnames(sensitivity_degs)) {
  stop("Both DEG tables must contain gene_id.", call. = FALSE)
}

full_ids <- unique(as.character(full_degs$gene_id))
sensitivity_ids <- unique(as.character(sensitivity_degs$gene_id))
shared_degs <- intersect(full_ids, sensitivity_ids)
shared_full <- full_degs[match(shared_degs, full_degs$gene_id), , drop = FALSE]
shared_sensitivity <- sensitivity_degs[match(shared_degs, sensitivity_degs$gene_id), , drop = FALSE]
direction_agreement <- if (length(shared_degs) == 0) {
  NA_real_
} else {
  100 * mean(shared_full$direction == shared_sensitivity$direction)
}
retained_percent <- if (length(full_ids) == 0) NA_real_ else 100 * length(shared_degs) / length(full_ids)

common_ids <- intersect(as.character(full_all$gene_id), as.character(sensitivity_all$gene_id))
full_common <- full_all[match(common_ids, full_all$gene_id), , drop = FALSE]
sensitivity_common <- sensitivity_all[match(common_ids, sensitivity_all$gene_id), , drop = FALSE]
complete <- is.finite(full_common$log2FoldChange) & is.finite(sensitivity_common$log2FoldChange)
pearson <- if (sum(complete) < 2) {
  NA_real_
} else {
  stats::cor(
    full_common$log2FoldChange[complete],
    sensitivity_common$log2FoldChange[complete],
    method = "pearson"
  )
}
spearman <- if (sum(complete) < 2) {
  NA_real_
} else {
  stats::cor(
    full_common$log2FoldChange[complete],
    sensitivity_common$log2FoldChange[complete],
    method = "spearman"
  )
}

comparison <- data.frame(
  full_degs = length(full_ids),
  sensitivity_degs = length(sensitivity_ids),
  shared_degs = length(shared_degs),
  full_degs_retained_percent = retained_percent,
  shared_direction_agreement_percent = direction_agreement,
  log2FoldChange_pearson = pearson,
  log2FoldChange_spearman = spearman,
  genes_correlated = sum(complete),
  stringsAsFactors = FALSE
)
comparison_path <- snakemake@output[["comparison"]]
dir.create(dirname(comparison_path), recursive = TRUE, showWarnings = FALSE)
write.table(comparison, comparison_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
writeLines(capture.output(sessionInfo()), snakemake@output[["session_info"]])

full_significant <- unique(as.character(full_degs$gene_id))
sensitivity_significant <- unique(as.character(sensitivity_degs$gene_id))
status <- ifelse(
  common_ids %in% full_significant & common_ids %in% sensitivity_significant,
  "DEG in both",
  ifelse(
    common_ids %in% full_significant,
    "Full analysis only",
    ifelse(common_ids %in% sensitivity_significant, "Sensitivity only", "Not significant in either")
  )
)
plot_data <- data.frame(
  full_log2FoldChange = full_common$log2FoldChange,
  sensitivity_log2FoldChange = sensitivity_common$log2FoldChange,
  status = factor(
    status,
    levels = c("Not significant in either", "Full analysis only", "Sensitivity only", "DEG in both")
  ),
  stringsAsFactors = FALSE
)
plot_data <- plot_data[complete, , drop = FALSE]
excluded <- paste(snakemake@params[["excluded_patients"]], collapse = ", ")
plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = full_log2FoldChange, y = sensitivity_log2FoldChange, colour = status)
) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey80") +
  ggplot2::geom_vline(xintercept = 0, colour = "grey80") +
  ggplot2::geom_point(alpha = 0.7, size = 1.5) +
  ggplot2::labs(
    title = paste("DE sensitivity excluding", excluded),
    x = "Full analysis log2FoldChange",
    y = "Sensitivity analysis log2FoldChange",
    colour = "DEG status"
  ) +
  ggplot2::theme_bw()

png_path <- snakemake@output[["plot_png"]]
pdf_path <- snakemake@output[["plot_pdf"]]
dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(png_path, plot, width = 8, height = 7, dpi = 300, bg = "white")
ggplot2::ggsave(pdf_path, plot, width = 8, height = 7, bg = "white")
