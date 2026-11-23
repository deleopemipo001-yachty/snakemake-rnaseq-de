# Count read pairs from all sorted sample BAMs using the configured strandedness.
rule featurecounts:
    input:
        bams=expand("results/bam/{sample}.bam", sample=samples.index),
        gtf=gtf_path,
    output:
        counts="results/counts/counts.txt",
        summary="results/counts/counts.txt.summary",
    log:
        "logs/featurecounts/counts.log",
    conda:
        "../envs/subread.yaml"
    threads: 4
    resources:
        mem_mb=9000,
    params:
        strandedness=require_featurecounts_strandedness,
    shell:
        "featureCounts -p --countReadPairs -t exon -g gene_id "
        "-T {threads} -s {params.strandedness} -a {input.gtf:q} "
        "-o {output.counts:q} {input.bams:q} > {log:q} 2>&1"


# Convert featureCounts output into a gene-by-sample matrix for downstream analysis.
rule counts_matrix:
    input:
        "results/counts/counts.txt",
    output:
        "results/counts/counts_matrix.tsv",
    log:
        "logs/featurecounts/counts_matrix.log",
    conda:
        "../envs/subread.yaml"
    threads: 1
    resources:
        mem_mb=1000,
    params:
        sample_names=samples.index.tolist(),
    script:
        "../scripts/counts_matrix.py"


# Compare assignment rates in all supported strandedness modes for one sample.
rule strandedness_check:
    input:
        bam="results/bam/{sample}.bam",
        gtf=gtf_path,
    output:
        "results/qc/strandedness/{sample}.tsv",
    log:
        "logs/strandedness/{sample}.log",
    conda:
        "../envs/subread.yaml"
    threads: 4
    resources:
        mem_mb=6000,
    script:
        "../scripts/strandedness_check.py"
