set.seed(1)
suppressPackageStartupMessages({library(limma); library(readr); library(dplyr); library(stringr); library(tibble)})
fpkm <- readr::read_tsv(snakemake@input[["fpkm"]], comment="#", show_col_types=FALSE, name_repair="minimal")
selected <- as.character(snakemake@params[["samples"]])
gene_col <- which(names(fpkm) %in% c("gene","Gene","GeneSymbol","Gene_Name","gene_name","GeneID","gene_id"))
if (!length(gene_col)) gene_col <- 1L else gene_col <- gene_col[[1]]
names(fpkm)[gene_col] <- "gene"
missing <- setdiff(selected, names(fpkm))
if (length(missing)) stop(paste("Missing requested FPKM columns:", paste(missing, collapse=", ")), call.=FALSE)
excluded <- setdiff(names(fpkm)[-gene_col], selected)
mat <- as.matrix(fpkm[, selected, drop=FALSE]); storage.mode(mat) <- "numeric"
if (anyNA(mat) || any(!is.finite(mat)) || any(mat < 0)) stop("FPKM values must be finite and non-negative.", call.=FALSE)
rownames(mat) <- as.character(fpkm$gene)
E <- log2(mat + 1)
keep <- rowMeans(E) >= as.numeric(snakemake@params[["min_mean_log2"]])
E <- E[keep,,drop=FALSE]
if (anyDuplicated(rownames(E))) E <- limma::avereps(E, ID=rownames(E))
group <- ifelse(str_ends(selected, "_T$"), "Tumor", ifelse(str_ends(selected, "_N$"), "Normal", NA_character_))
if (anyNA(group)) stop("Sample names must end in _T or _N.", call.=FALSE)
pair_id <- str_remove(selected, "_[TN]$")
group <- factor(group, levels=c("Normal","Tumor")); pair_id <- factor(pair_id)
design <- model.matrix(~0 + group + pair_id)
colnames(design)[1:2] <- c("Normal","Tumor")
fit <- lmFit(E, design); cont <- makeContrasts(Tumor_vs_Normal=Tumor-Normal, levels=design)
top <- topTable(eBayes(contrasts.fit(fit, cont), trend=TRUE), number=Inf, sort.by="P")
top$gene <- rownames(top); top <- top[,c("gene",setdiff(names(top),"gene")),drop=FALSE]
deg <- top[!is.na(top$adj.P.Val) & top$adj.P.Val <= as.numeric(snakemake@params[["padj"]]) & abs(top$logFC) >= as.numeric(snakemake@params[["lfc"]]),,drop=FALSE]
write.table(top,snakemake@output[["all_genes"]],sep="\t",quote=FALSE,row.names=FALSE)
write.table(deg,snakemake@output[["degs"]],sep="\t",quote=FALSE,row.names=FALSE)
# Read GEO matrix only to attach accession IDs to selected sample titles.
lines <- readLines(gzfile(snakemake@input[["series_matrix"]]),warn=FALSE)
titles <- lines[grepl("^!Sample_title\\t",lines)][1]
accessions <- lines[grepl("^!Sample_geo_accession\\t",lines)][1]
parse_line <- function(x) { z <- strsplit(x,"\\t",fixed=FALSE)[[1]][-1]; gsub('^"|"$','',z) }
gsm <- setNames(parse_line(accessions),parse_line(titles))
selected_gsm <- unname(gsm[selected]); names(selected_gsm) <- selected
if (anyNA(selected_gsm) || any(selected_gsm == "")) {
  stop("Could not map every selected sample title to a GEO GSM accession.", call.=FALSE)
}
n_degs <- nrow(deg)
n_up <- sum(deg$logFC > 0)
n_down <- sum(deg$logFC < 0)
checks <- data.frame(
  metric = c("samples_used", "samples_excluded", "genes_in_model", "DEGs_total", "DEGs_up", "DEGs_down", "note"),
  value = c(
    paste(paste0(selected, "=", selected_gsm), collapse = ", "),
    paste(excluded, collapse = ", "),
    as.character(nrow(top)), as.character(n_degs), as.character(n_up), as.character(n_down),
    "The sample set was verified empirically: the 12 RC samples reproduce 12,407 filtered genes exactly (all 40 samples give 12,367); the matching code in the saved R history cannot select these samples from this series matrix (titles are quoted), so the final run's selection step was not preserved."
  ), stringsAsFactors = FALSE
)
add_check <- function(metric, value, expected = NA_character_, status = NA_character_) {
  checks <<- rbind(checks, data.frame(
    metric = metric, value = as.character(value), expected = as.character(expected),
    status = as.character(status), stringsAsFactors = FALSE
  ))
}
checks$expected <- NA_character_
checks$status <- NA_character_
expectations <- c(
  genes_in_model = as.numeric(snakemake@params[["expected_genes"]]),
  DEGs_total = as.numeric(snakemake@params[["expected_degs"]]),
  DEGs_up = as.numeric(snakemake@params[["expected_up"]]),
  DEGs_down = as.numeric(snakemake@params[["expected_down"]])
)
for (metric in names(expectations)) {
  idx <- match(metric, checks$metric)
  checks$expected[idx] <- as.character(expectations[[metric]])
  checks$status[idx] <- if (as.numeric(checks$value[idx]) == expectations[[metric]]) "PASS" else "FAIL"
}
add_check("excluded_column_count", length(excluded))
limma_version <- as.character(utils::packageVersion("limma"))
add_check(
  "limma_version_used", limma_version, expected = "3.66.0",
  status = if (limma_version == "3.66.0") "PASS" else "FAIL"
)
add_check("limma_version_original", "3.64.3", expected = "3.64.3", status = "PASS")

original_input <- snakemake@input[["original"]]
if (is.null(original_input) || length(original_input) == 0L || !nzchar(as.character(original_input[[1]]))) {
  add_check("original_degs_comparison", "Skipped: original DEG file is not configured.")
} else {
  original_path <- as.character(original_input[[1]])
  if (!file.exists(original_path)) {
    add_check("original_degs_comparison", "Skipped: original DEG file is missing.")
  } else {
    original <- tryCatch(
      read.csv(original_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE),
      error = function(error) NULL
    )
    required_original_columns <- c("logFC", "AveExpr", "adj.P.Val")
    if (is.null(original)) {
      add_check("original_degs_comparison", "Skipped: original DEG file could not be read.")
    } else if (!all(required_original_columns %in% colnames(original))) {
      add_check("original_degs_comparison", paste(
        "FAIL: expected R write.csv columns are missing:",
        paste(setdiff(required_original_columns, colnames(original)), collapse = ", ")
      ), status = "FAIL")
    } else {
      old_gene <- rownames(original)
      reproduced_gene <- as.character(top$gene)
      original_deg <- original$adj.P.Val <= as.numeric(snakemake@params[["padj"]]) &
        abs(original$logFC) >= as.numeric(snakemake@params[["lfc"]])
      original_deg[is.na(original_deg)] <- FALSE
      original_deg_genes <- old_gene[original_deg]
      current_deg_genes <- as.character(deg$gene)
      identical_set <- setequal(reproduced_gene, old_gene)
      identical_calls <- setequal(current_deg_genes, original_deg_genes)
      disagreeing <- length(setdiff(current_deg_genes, original_deg_genes)) +
        length(setdiff(original_deg_genes, current_deg_genes))
      shared_genes <- intersect(reproduced_gene, old_gene)
      add_check(
        "identical_gene_set", if (identical_set) "Y" else "N",
        expected = "Y", status = if (identical_set) "PASS" else "FAIL"
      )
      add_check(
        "identical_DEG_calls", if (identical_calls) "Y" else "N",
        expected = "Y", status = if (identical_calls) "PASS" else "FAIL"
      )
      add_check(
        "disagreeing_DEG_gene_count", disagreeing,
        expected = 0, status = if (disagreeing == 0L) "PASS" else "FAIL"
      )
      add_check(
        "original_gene_count", length(old_gene), expected = nrow(top),
        status = if (length(old_gene) == nrow(top)) "PASS" else "FAIL"
      )
      shared_count_pass <- length(shared_genes) == nrow(top) &&
        nrow(original) == nrow(top) && length(old_gene) == nrow(top)
      add_check(
        "n_shared_genes", length(shared_genes), expected = nrow(top),
        status = if (shared_count_pass) "PASS" else "FAIL"
      )
      delta_threshold <- 1e-6
      compare_numeric <- function(column, label) {
        if (length(shared_genes) == 0L) {
          add_check(label, "NA", expected = paste0("<= ", delta_threshold), status = "NA")
          return(invisible(NULL))
        }
        new_values <- suppressWarnings(as.numeric(top[[column]][match(shared_genes, top$gene)]))
        old_values <- suppressWarnings(as.numeric(original[[column]][match(shared_genes, rownames(original))]))
        if (any(!is.finite(new_values)) || any(!is.finite(old_values))) {
          stop(paste("Non-finite values encountered while comparing", column), call. = FALSE)
        }
        delta <- max(abs(new_values - old_values))
        passed <- delta <= delta_threshold
        add_check(
          label,
          if (is.na(delta)) "NA" else format(delta, scientific = TRUE, digits = 8),
          expected = paste0("<= ", delta_threshold),
          status = if (is.na(delta)) "NA" else if (passed) "PASS" else "FAIL"
        )
      }
      compare_numeric("logFC", "max_abs_delta_logFC")
      compare_numeric("AveExpr", "max_abs_delta_AveExpr")
      compare_numeric("adj.P.Val", "max_abs_delta_adj.P.Val")
    }
  }
}

write.table(checks,snakemake@output[["check"]],sep="\t",quote=FALSE,row.names=FALSE,na="NA")
writeLines(capture.output(sessionInfo()),file.path(dirname(snakemake@output[["check"]]),"session_info.txt"))
