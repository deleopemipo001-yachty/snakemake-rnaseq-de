set.seed(1)
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

dds <- readRDS(snakemake@input[["dds"]])
vst_object <- readRDS(snakemake@input[["vst"]])
all_genes <- read.delim(
  snakemake@input[["all_genes"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
sample_table <- read.delim(
  snakemake@input[["samples"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)

patient_column <- snakemake@params[["patient_column"]]
condition_column <- snakemake@params[["condition_column"]]
reference <- snakemake@params[["reference"]]
padj_threshold <- as.numeric(snakemake@params[["padj"]])
lfc_threshold <- as.numeric(snakemake@params[["lfc"]])

required_sample_columns <- c("sample", patient_column, condition_column)
missing_columns <- setdiff(required_sample_columns, colnames(sample_table))
if (length(missing_columns) > 0) {
  stop(
    paste(
      "The DE sample sheet is missing required columns:",
      paste(missing_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}
if (anyDuplicated(sample_table$sample)) {
  stop("The DE sample sheet contains duplicate sample names.", call. = FALSE)
}

sample_names <- as.character(sample_table$sample)
check_sample_names <- function(observed, source_name) {
  missing <- setdiff(sample_names, observed)
  extra <- setdiff(observed, sample_names)
  if (length(missing) > 0 || length(extra) > 0) {
    stop(
      paste0(
        source_name,
        " sample names do not match the DE sample sheet. Missing: ",
        if (length(missing)) paste(missing, collapse = ", ") else "none",
        "; extra: ",
        if (length(extra)) paste(extra, collapse = ", ") else "none"
      ),
      call. = FALSE
    )
  }
}

check_sample_names(colnames(dds), "dds.rds")
check_sample_names(colnames(vst_object), "vst.rds")
if (!all(c("gene_id", "gene_name", "baseMean", "log2FoldChange", "log2FoldChange_shrunk", "padj") %in% colnames(all_genes))) {
  stop("The DE results table is missing columns required to make figures.", call. = FALSE)
}

sample_table <- sample_table[match(sample_names, sample_table$sample), , drop = FALSE]
metadata <- data.frame(
  sample = sample_names,
  condition = factor(sample_table[[condition_column]]),
  patient = factor(sample_table[[patient_column]]),
  stringsAsFactors = FALSE
)
if (!reference %in% levels(metadata$condition)) {
  stop(paste0("Reference condition '", reference, "' is absent from the DE sample sheet."), call. = FALSE)
}

vst_matrix <- assay(vst_object)[, sample_names, drop = FALSE]
gene_variances <- apply(vst_matrix, 1, stats::var)
gene_variances[!is.finite(gene_variances)] <- -Inf
top_variable_count <- min(500L, nrow(vst_matrix))
if (top_variable_count < 2) {
  stop("At least two genes are required to make PCA figures.", call. = FALSE)
}
top_variable_ids <- names(sort(gene_variances, decreasing = TRUE))[seq_len(top_variable_count)]
top_variable_matrix <- vst_matrix[top_variable_ids, , drop = FALSE]

save_ggplot_pair <- function(plot, png_path, pdf_path, width = 9, height = 7) {
  ggplot2::ggsave(
    filename = png_path,
    plot = plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
  ggplot2::ggsave(
    filename = pdf_path,
    plot = plot,
    width = width,
    height = height,
    device = grDevices::cairo_pdf,
    bg = "white"
  )
}

save_heatmap_pair <- function(plot, png_path, pdf_path, width = 10, height = 8) {
  grDevices::png(filename = png_path, width = width, height = height, units = "in", res = 300)
  grid::grid.newpage()
  grid::grid.draw(plot$gtable)
  grDevices::dev.off()

  grDevices::cairo_pdf(filename = pdf_path, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(plot$gtable)
  grDevices::dev.off()
}

make_pca_plot <- function(expression_matrix, title) {
  pca <- stats::prcomp(t(expression_matrix), center = TRUE, scale. = FALSE)
  variance_percent <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  coordinates <- data.frame(
    sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )
  coordinates <- merge(coordinates, metadata, by = "sample", sort = FALSE)
  coordinates <- coordinates[match(sample_names, coordinates$sample), , drop = FALSE]

  ggplot(coordinates, aes(x = PC1, y = PC2, color = condition, shape = patient, label = sample)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(max.overlaps = Inf, show.legend = FALSE) +
    labs(
      title = title,
      x = sprintf("PC1 (%.1f%% variance)", variance_percent[[1]]),
      y = sprintf("PC2 (%.1f%% variance)", variance_percent[[2]]),
      color = "Condition",
      shape = "Patient"
    ) +
    theme_bw()
}

pca <- make_pca_plot(top_variable_matrix, "PCA of VST values (top 500 variable genes)")
save_ggplot_pair(
  pca,
  snakemake@output[["pca_png"]],
  snakemake@output[["pca_pdf"]]
)

patient_adjusted <- limma::removeBatchEffect(
  top_variable_matrix,
  batch = metadata$patient,
  design = stats::model.matrix(~ condition, data = metadata)
)
pca_patient_adjusted <- make_pca_plot(
  patient_adjusted,
  "PCA after patient effect removal for visualization only"
)
save_ggplot_pair(
  pca_patient_adjusted,
  snakemake@output[["pca_patient_adjusted_png"]],
  snakemake@output[["pca_patient_adjusted_pdf"]]
)

sample_distances <- as.matrix(stats::dist(t(vst_matrix), method = "euclidean"))
rownames(sample_distances) <- sample_names
colnames(sample_distances) <- sample_names
sample_annotations <- data.frame(
  condition = metadata$condition,
  patient = metadata$patient,
  row.names = sample_names
)
distance_heatmap <- pheatmap::pheatmap(
  sample_distances,
  annotation_row = sample_annotations,
  annotation_col = sample_annotations,
  main = "Euclidean distances between samples (VST)",
  silent = TRUE
)
save_heatmap_pair(
  distance_heatmap,
  snakemake@output[["sample_distances_png"]],
  snakemake@output[["sample_distances_pdf"]],
  width = 8,
  height = 7
)

all_genes$significant <- !is.na(all_genes$padj) &
  all_genes$padj < padj_threshold &
  !is.na(all_genes$log2FoldChange) &
  abs(all_genes$log2FoldChange) >= lfc_threshold
all_genes$de_group <- ifelse(
  all_genes$significant & all_genes$log2FoldChange > 0,
  "up",
  ifelse(all_genes$significant, "down", "not significant")
)
de_colors <- c(
  up = "#D73027",
  down = "#4575B4",
  "not significant" = "grey75"
)

ma_plot <- ggplot(
  all_genes,
  aes(x = log10(baseMean + 1), y = log2FoldChange_shrunk, color = de_group)
) +
  geom_point(alpha = 0.75, size = 1.2, na.rm = TRUE) +
  scale_color_manual(values = de_colors, name = "DE status") +
  labs(
    title = "MA plot: mean expression and shrunken log2 fold change",
    x = "log10(baseMean + 1)",
    y = "Shrunken log2 fold change"
  ) +
  theme_bw()
save_ggplot_pair(
  ma_plot,
  snakemake@output[["ma_plot_png"]],
  snakemake@output[["ma_plot_pdf"]]
)

all_genes$volcano_y <- -log10(pmax(all_genes$padj, .Machine$double.xmin))
top_labels <- all_genes[!is.na(all_genes$padj), , drop = FALSE]
top_labels <- top_labels[order(top_labels$padj, top_labels$gene_id), , drop = FALSE]
top_labels <- head(top_labels, 15)
missing_label <- is.na(top_labels$gene_name) | top_labels$gene_name == ""
top_labels$gene_name[missing_label] <- top_labels$gene_id[missing_label]

volcano <- ggplot(
  all_genes,
  aes(x = log2FoldChange, y = volcano_y, color = de_group)
) +
  geom_point(alpha = 0.75, size = 1.2, na.rm = TRUE) +
  geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = top_labels,
    aes(label = gene_name),
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = de_colors, name = "DE status") +
  labs(
    title = "Volcano plot: unshrunken log2 fold change",
    x = "Unshrunken log2 fold change",
    y = "-log10(adjusted p-value)"
  ) +
  theme_bw()
save_ggplot_pair(
  volcano,
  snakemake@output[["volcano_png"]],
  snakemake@output[["volcano_pdf"]]
)

significant_genes <- all_genes[all_genes$significant & !is.na(all_genes$padj), , drop = FALSE]
significant_genes <- significant_genes[
  order(significant_genes$padj, significant_genes$gene_id),
  ,
  drop = FALSE
]
if (nrow(significant_genes) == 0) {
  stop("No significant DEGs are available for the top-gene heatmap.", call. = FALSE)
}
top_degs <- head(significant_genes, 50)
missing_vst_genes <- setdiff(top_degs$gene_id, rownames(vst_matrix))
if (length(missing_vst_genes) > 0) {
  stop(
    paste("Top DEGs are missing from vst.rds:", paste(missing_vst_genes, collapse = ", ")),
    call. = FALSE
  )
}
heatmap_values <- vst_matrix[top_degs$gene_id, sample_names, drop = FALSE]
heatmap_zscores <- t(scale(t(heatmap_values)))
heatmap_zscores[!is.finite(heatmap_zscores)] <- 0
gene_labels <- top_degs$gene_name
missing_label <- is.na(gene_labels) | gene_labels == ""
gene_labels[missing_label] <- top_degs$gene_id[missing_label]
rownames(heatmap_zscores) <- make.unique(gene_labels)

top_gene_heatmap <- pheatmap::pheatmap(
  heatmap_zscores,
  annotation_col = sample_annotations,
  main = "Top significant DEGs (VST, row z-scores)",
  show_rownames = TRUE,
  fontsize_row = 6,
  silent = TRUE
)
save_heatmap_pair(
  top_gene_heatmap,
  snakemake@output[["heatmap_top50_png"]],
  snakemake@output[["heatmap_top50_pdf"]],
  width = 10,
  height = 11
)
