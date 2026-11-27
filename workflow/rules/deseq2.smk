# Extract gene identifiers, names, and biotypes from the selected counting GTF.
rule gene_annotation:
    input:
        gtf=gtf_path,
    output:
        "results/deseq2/gene_annotation.tsv",
    log:
        "logs/deseq2/gene_annotation.log",
    conda:
        "../envs/deseq2.yaml"
    threads: 1
    resources:
        mem_mb=1000,
    script:
        "../scripts/gene_annotation.R"


# Fit the paired tumor-versus-normal model and save tables plus plotting objects.
rule deseq2:
    input:
        counts=config["de"]["counts"],
        samples=config["de"]["samples"],
        annotation=config["de"]["gene_annotation"],
    output:
        all_genes="results/deseq2/all_genes.tsv",
        degs="results/deseq2/degs.tsv",
        summary="results/deseq2/summary.tsv",
        dds="results/deseq2/dds.rds",
        vst="results/deseq2/vst.rds",
        session_info="results/deseq2/session_info.txt",
    log:
        "logs/deseq2/deseq2.log",
    conda:
        "../envs/deseq2.yaml"
    threads: 1
    resources:
        mem_mb=9000,
    params:
        patient_column=config["de"]["patient_column"],
        condition_column=config["de"]["condition_column"],
        reference=config["de"]["reference"],
        min_count=config["de"]["min_count"],
        min_samples=config["de"]["min_samples"],
        padj=config["de"]["padj"],
        lfc=config["de"]["lfc"],
    script:
        "../scripts/deseq2.R"


# Make diagnostic plots that show sample relationships and the DE results.
rule de_figures:
    input:
        dds="results/deseq2/dds.rds",
        vst="results/deseq2/vst.rds",
        all_genes="results/deseq2/all_genes.tsv",
        samples=config["de"]["samples"],
    output:
        pca_png="results/figures/pca.png",
        pca_pdf="results/figures/pca.pdf",
        pca_patient_adjusted_png="results/figures/pca_patient_adjusted.png",
        pca_patient_adjusted_pdf="results/figures/pca_patient_adjusted.pdf",
        sample_distances_png="results/figures/sample_distances.png",
        sample_distances_pdf="results/figures/sample_distances.pdf",
        ma_plot_png="results/figures/ma_plot.png",
        ma_plot_pdf="results/figures/ma_plot.pdf",
        volcano_png="results/figures/volcano.png",
        volcano_pdf="results/figures/volcano.pdf",
        heatmap_top50_png="results/figures/heatmap_top50.png",
        heatmap_top50_pdf="results/figures/heatmap_top50.pdf",
    log:
        "logs/deseq2/de_figures.log",
    conda:
        "../envs/deseq2.yaml"
    threads: 1
    resources:
        mem_mb=3000,
    params:
        patient_column=config["de"]["patient_column"],
        condition_column=config["de"]["condition_column"],
        reference=config["de"]["reference"],
        padj=config["de"]["padj"],
        lfc=config["de"]["lfc"],
    script:
        "../scripts/de_figures.R"
