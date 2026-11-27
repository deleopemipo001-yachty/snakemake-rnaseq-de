# Reproduce the GEO FPKM limma analysis from the explicitly selected samples.
rule preprint_limma:
    input:
        fpkm=config["preprint"]["fpkm"],
        series_matrix=config["preprint"]["series_matrix"],
        original=optional_preprint_file(config["preprint"].get("original_degs", "")),
    output:
        all_genes="results/preprint_comparison/preprint_all_genes.tsv",
        degs="results/preprint_comparison/preprint_degs.tsv",
        check="results/preprint_comparison/reproduction_check.tsv",
    log:
        "logs/preprint/preprint_limma.log",
    conda:
        "../envs/limma.yaml"
    threads: 1
    resources:
        mem_mb=6000,
    params:
        samples=config["preprint"]["samples"],
        min_mean_log2=config["preprint"]["min_mean_log2"],
        padj=config["preprint"]["padj"],
        lfc=config["preprint"]["lfc"],
        expected_genes=config["preprint"]["expected_genes"],
        expected_degs=config["preprint"]["expected_degs"],
        expected_up=config["preprint"]["expected_up"],
        expected_down=config["preprint"]["expected_down"],
    script:
        "../scripts/preprint_limma.R"


# Compare reproduced limma results with the paired DESeq2 analysis.
rule preprint_compare:
    input:
        preprint="results/preprint_comparison/preprint_all_genes.tsv",
        preprint_degs="results/preprint_comparison/preprint_degs.tsv",
        deseq="results/deseq2/all_genes.tsv",
        deseq_degs="results/deseq2/degs.tsv",
        annotation="results/deseq2/gene_annotation.tsv",
    output:
        summary="results/preprint_comparison/summary.tsv",
        overlap="results/preprint_comparison/overlap.tsv",
        unmatched_breakdown="results/preprint_comparison/unmatched_breakdown.tsv",
        scatter_png="results/figures/preprint_vs_deseq2_lfc.png",
        scatter_pdf="results/figures/preprint_vs_deseq2_lfc.pdf",
        venn_png="results/figures/preprint_vs_deseq2_venn.png",
        venn_pdf="results/figures/preprint_vs_deseq2_venn.pdf",
    log:
        "logs/preprint/preprint_compare.log",
    conda:
        "../envs/limma.yaml"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "../scripts/preprint_compare.R"


# Run the shared offline ORA method for preprint and DESeq2 directions.
rule preprint_enrichment:
    input:
        preprint="results/preprint_comparison/preprint_all_genes.tsv",
        preprint_degs="results/preprint_comparison/preprint_degs.tsv",
        deseq="results/deseq2/all_genes.tsv",
        deseq_degs="results/deseq2/degs.tsv",
        annotation="results/deseq2/gene_annotation.tsv",
        ora_functions=ora_functions_path,
        original=optional_preprint_file(config["preprint"].get("original_kegg", "")),
    output:
        agreement="results/preprint_comparison/pathway_agreement.tsv",
        original="results/preprint_comparison/original_kegg_vs_ours.tsv",
        session_info="results/preprint_comparison/enrichment_session_info.txt",
    log:
        "logs/preprint/preprint_enrichment.log",
    conda:
        "../envs/enrichment.yaml"
    threads: 1
    resources:
        mem_mb=5000,
    params:
        padj=config["enrichment"]["padj"],
        min_size=config["enrichment"]["min_size"],
        max_size=config["enrichment"]["max_size"],
        kegg_aliases=config["preprint"]["kegg_aliases"],
    script:
        "../scripts/preprint_enrichment.R"


# Explicit target for the complete reproduction and comparison branch.
rule preprint_comparison:
    input:
        "results/preprint_comparison/reproduction_check.tsv",
        "results/preprint_comparison/summary.tsv",
        "results/preprint_comparison/overlap.tsv",
        "results/preprint_comparison/unmatched_breakdown.tsv",
        "results/preprint_comparison/pathway_agreement.tsv",
        "results/preprint_comparison/original_kegg_vs_ours.tsv",
        "results/figures/preprint_vs_deseq2_lfc.png",
        "results/figures/preprint_vs_deseq2_lfc.pdf",
        "results/figures/preprint_vs_deseq2_venn.png",
        "results/figures/preprint_vs_deseq2_venn.pdf",
