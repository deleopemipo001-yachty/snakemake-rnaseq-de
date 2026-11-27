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
