gtf_path <- snakemake@input[["gtf"]]
output_path <- snakemake@output[[1]]

connection <- if (grepl("\\.gz$", gtf_path, ignore.case = TRUE)) {
  gzfile(gtf_path, open = "rt")
} else {
  file(gtf_path, open = "rt")
}

get_attribute <- function(attributes, key) {
  pattern <- paste0("(^|; )", key, ' "([^"]+)"')
  matches <- regexec(pattern, attributes, perl = TRUE)
  values <- regmatches(attributes, matches)
  vapply(values, function(value) {
    if (length(value) >= 3) value[[3]] else NA_character_
  }, character(1))
}

gene_id <- character()
gene_name <- character()
gene_biotype <- character()
repeat {
  lines <- readLines(connection, n = 10000, warn = FALSE)
  if (length(lines) == 0) break
  fields <- strsplit(lines, "\t", fixed = TRUE)
  is_gene <- vapply(fields, function(row) length(row) >= 9 && row[[3]] == "gene", logical(1))
  if (!any(is_gene)) next
  attributes <- vapply(fields[is_gene], function(row) row[[9]], character(1))
  ids <- get_attribute(attributes, "gene_id")
  names <- get_attribute(attributes, "gene_name")
  biotypes <- get_attribute(attributes, "gene_biotype")
  missing_biotype <- is.na(biotypes)
  if (any(missing_biotype)) {
    biotypes[missing_biotype] <- get_attribute(attributes[missing_biotype], "gene_type")
  }
  gene_id <- c(gene_id, ids)
  gene_name <- c(gene_name, names)
  gene_biotype <- c(gene_biotype, biotypes)
}
close(connection)

annotation <- data.frame(
  gene_id = gene_id,
  gene_name = gene_name,
  gene_biotype = gene_biotype,
  stringsAsFactors = FALSE
)
annotation <- annotation[!is.na(annotation$gene_id), , drop = FALSE]
annotation <- annotation[!duplicated(annotation$gene_id), , drop = FALSE]
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.table(annotation, output_path, sep = "\t", quote = FALSE, row.names = FALSE)
