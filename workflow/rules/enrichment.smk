# Test up- and down-regulated gene lists against the pinned offline gene sets.
rule pathway_enrichment:
    input:
        all_genes="results/deseq2/all_genes.tsv",
        degs="results/deseq2/degs.tsv",
        annotation="results/deseq2/gene_annotation.tsv",
        ora_functions=ora_functions_path,
    output:
        kegg_up="results/enrichment/kegg_up.tsv",
        kegg_down="results/enrichment/kegg_down.tsv",
        hallmark_up="results/enrichment/hallmark_up.tsv",
        hallmark_down="results/enrichment/hallmark_down.tsv",
        kegg_up_png="results/figures/enrichment_kegg_up.png",
        kegg_up_pdf="results/figures/enrichment_kegg_up.pdf",
        kegg_down_png="results/figures/enrichment_kegg_down.png",
        kegg_down_pdf="results/figures/enrichment_kegg_down.pdf",
        hallmark_up_png="results/figures/enrichment_hallmark_up.png",
        hallmark_up_pdf="results/figures/enrichment_hallmark_up.pdf",
        hallmark_down_png="results/figures/enrichment_hallmark_down.png",
        hallmark_down_pdf="results/figures/enrichment_hallmark_down.pdf",
        versions="results/enrichment/versions.txt",
        session_info="results/enrichment/session_info.txt",
    log:
        "logs/enrichment/pathway_enrichment.log",
    conda:
        "../envs/enrichment.yaml"
    threads: 1
    resources:
        mem_mb=4000,
    params:
        padj=config["enrichment"]["padj"],
        min_size=config["enrichment"]["min_size"],
        max_size=config["enrichment"]["max_size"],
    script:
        "../scripts/enrichment.R"
