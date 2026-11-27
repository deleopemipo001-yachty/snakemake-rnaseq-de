set.seed(1)
suppressPackageStartupMessages({library(ggplot2)})
readtab <- function(p) read.delim(p, check.names=FALSE, stringsAsFactors=FALSE)
p <- readtab(snakemake@input[["preprint"]]); pd <- readtab(snakemake@input[["preprint_degs"]]); d <- readtab(snakemake@input[["deseq"]]); dd <- readtab(snakemake@input[["deseq_degs"]]); annotation <- readtab(snakemake@input[["annotation"]])
clean <- function(x) sub("\\.[0-9]+$", "", as.character(x))
if (!"gene_id" %in% names(d)) stop("DESeq2 all_genes.tsv must contain gene_id.",call.=FALSE)
pid <- if("gene_id" %in% names(p)) p$gene_id else p$gene
p$ens <- clean(pid); d$ens <- clean(d$gene_id); pd$ens <- clean(if("gene_id"%in%names(pd))pd$gene_id else pd$gene); dd$ens <- clean(dd$gene_id)
shared <- intersect(p$ens, d$ens)
preprint_index <- match(shared, p$ens)
deseq_index <- match(shared, d$ens)
pfc <- as.numeric(p$logFC[preprint_index])
dfc <- as.numeric(d$log2FoldChange[deseq_index])
pdeg <- pd$ens; ddeg <- dd$ens
if (!"gene_id" %in% names(annotation)) stop("Gene annotation table must contain gene_id.", call.=FALSE)
annotation_ids <- unique(clean(annotation$gene_id))
annotation_ids <- annotation_ids[!is.na(annotation_ids) & annotation_ids != ""]
tested_deseq_ids <- unique(d$ens)
deseq_deg_ids <- unique(ddeg)
make_breakdown <- function(set_name, genes) {
  genes <- unique(genes[!is.na(genes) & genes != ""])
  in_annotation <- genes %in% annotation_ids
  tested_by_deseq <- genes %in% tested_deseq_ids
  called_by_deseq <- genes %in% deseq_deg_ids
  category_genes <- list(
    absent_from_Ensembl_116_annotation = genes[!in_annotation],
    present_but_removed_by_DESeq2_count_filter = genes[in_annotation & !tested_by_deseq],
    present_and_tested_but_not_called_DEG_by_DESeq2 = genes[
      in_annotation & tested_by_deseq & !called_by_deseq
    ]
  )
  rows <- lapply(names(category_genes), function(category) {
    category_values <- category_genes[[category]]
    data.frame(
      gene_set = set_name,
      classification = category,
      count = length(category_values),
      genes = paste(category_values, collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  rows[[length(rows) + 1L]] <- data.frame(
    gene_set = set_name,
    classification = "total",
    count = length(genes),
    genes = paste(genes, collapse = ","),
    stringsAsFactors = FALSE
  )
  do.call(rbind, rows)
}
unmatched_preprint_ids <- setdiff(unique(p$ens), unique(d$ens))
preprint_only_deg_ids <- setdiff(unique(pdeg), unique(ddeg))
unmatched_breakdown <- rbind(
  make_breakdown("unmatched_preprint_genes", unmatched_preprint_ids),
  make_breakdown("preprint_only_DEGs", preprint_only_deg_ids)
)
write.table(
  unmatched_breakdown, snakemake@output[["unmatched_breakdown"]],
  sep="\t", quote=FALSE, row.names=FALSE, na="NA"
)
allids <- union(p$ens,d$ens); status <- data.frame(gene_id=allids,preprint=ifelse(allids%in%pdeg,"DEG","not_DEG"),deseq2=ifelse(allids%in%ddeg,"DEG","not_DEG"))
write.table(status,snakemake@output[["overlap"]],sep="\t",quote=FALSE,row.names=FALSE)
calc <- function(universe) { a <- intersect(pdeg,universe); b <- intersect(ddeg,universe); ov <- intersect(a,b); c(preprint_genes=length(intersect(p$ens,universe)),deseq2_genes=length(intersect(d$ens,universe)),shared_tested=length(intersect(p$ens,d$ens)),preprint_DEGs=length(a),deseq2_DEGs=length(b),overlap=length(ov),pct_preprint_recovered=if(length(a))100*length(ov)/length(a) else NA,jaccard=if(length(union(a,b)))length(ov)/length(union(a,b)) else NA,direction_agreement=NA,pearson=NA,spearman=NA) }
full <- calc(allids)
preprint_universe <- intersect(unique(p$ens), unique(d$ens))
in12 <- calc(preprint_universe)
if(length(shared)>1) { x <- pfc; y <- dfc; pear <- suppressWarnings(cor(x,y,use="complete.obs",method="pearson")); spear <- suppressWarnings(cor(x,y,use="complete.obs",method="spearman")); common_deg <- intersect(intersect(pdeg,ddeg),shared); agree <- if(length(common_deg))mean(sign(p$logFC[match(common_deg,p$ens)])==sign(d$log2FoldChange[match(common_deg,d$ens)])) else NA; full[c("direction_agreement","pearson","spearman")] <- c(agree,pear,spear); in12[c("direction_agreement","pearson","spearman")] <- c(agree,pear,spear) }
sumtab <- data.frame(universe=c("all_union","preprint_filtered_universe"),rbind(full,in12),matched=length(shared),unmatched_preprint=length(setdiff(p$ens,d$ens)),unmatched_deseq2=length(setdiff(d$ens,p$ens)),check.names=FALSE)
write.table(sumtab,snakemake@output[["summary"]],sep="\t",quote=FALSE,row.names=FALSE,na="NA")
plotdata <- data.frame(limma=pfc,deseq=dfc,status=ifelse(shared%in%intersect(pdeg,ddeg),"both",ifelse(shared%in%pdeg,"preprint",ifelse(shared%in%ddeg,"DESeq2","neither"))))
r <- suppressWarnings(cor(plotdata$limma,plotdata$deseq,use="complete.obs"))
sc <- ggplot(plotdata,aes(limma,deseq,colour=status))+geom_point(alpha=.7)+geom_abline(slope=1,intercept=0,linetype=2)+labs(x="limma logFC",y="DESeq2 unshrunken log2FC",subtitle=paste("Pearson r =",signif(r,3)))+theme_bw()
venn_circle <- function(center_x, center_y = 0, radius = 1, npoints = 200L) {
  if (length(npoints) != 1L || is.na(npoints) || npoints < 4L) {
    stop("npoints must be one integer of at least 4.", call. = FALSE)
  }
  point_count <- as.integer(npoints)
  angle <- seq(0, 2 * base::pi, length.out = point_count)
  data.frame(
    x = center_x + radius * cos(angle),
    y = center_y + radius * sin(angle)
  )
}
preprint_circle <- venn_circle(center_x = 0, radius = 1, npoints = 200L)
deseq_circle <- venn_circle(center_x = 1, radius = 1, npoints = 200L)
preprint_only <- length(setdiff(pdeg, ddeg))
shared_degs <- length(intersect(pdeg, ddeg))
deseq_only <- length(setdiff(ddeg, pdeg))
ven <- ggplot() +
  geom_path(data = preprint_circle, aes(x, y), colour = "#4477AA", linewidth = 1.2) +
  geom_path(data = deseq_circle, aes(x, y), colour = "#CC6677", linewidth = 1.2) +
  annotate("text", x = -0.48, y = 0, label = preprint_only, size = 6) +
  annotate("text", x = 0.5, y = 0, label = shared_degs, size = 6) +
  annotate("text", x = 1.48, y = 0, label = deseq_only, size = 6) +
  annotate("text", x = -0.55, y = 1.12, label = "Preprint DEGs", colour = "#4477AA") +
  annotate("text", x = 1.55, y = 1.12, label = "DESeq2 DEGs", colour = "#CC6677") +
  coord_equal(xlim = c(-1.2, 2.2), ylim = c(-1.2, 1.45), clip = "off") +
  theme_void()
for (ext in c("png","pdf")) { ggsave(snakemake@output[[paste0("scatter_",ext)]],sc,width=7,height=6,dpi=300); ggsave(snakemake@output[[paste0("venn_",ext)]],ven,width=7,height=6,dpi=300) }
writeLines(capture.output(sessionInfo()),file.path(dirname(snakemake@output[["summary"]]),"compare_session_info.txt"))
