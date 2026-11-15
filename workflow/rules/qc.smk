# Check each raw mate so read quality issues are visible before trimming.
rule fastqc:
    input:
        fastq=lambda wildcards: sample_fastq(wildcards.sample, wildcards.read),
    output:
        html="results/fastqc/{sample}_{read}_fastqc.html",
        zip="results/fastqc/{sample}_{read}_fastqc.zip",
    log:
        "logs/fastqc/{sample}_{read}.log",
    threads: 1
    resources:
        mem_mb=1000,
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc --threads {threads} --outdir results/fastqc {input.fastq} > {log} 2>&1"


# Trim paired reads with fastp's defaults and automatic paired-end adapter detection.
rule fastp:
    input:
        r1=lambda wildcards: sample_fastq(wildcards.sample, "R1"),
        r2=lambda wildcards: sample_fastq(wildcards.sample, "R2"),
    output:
        r1="results/trimmed/{sample}_R1.fastq.gz",
        r2="results/trimmed/{sample}_R2.fastq.gz",
        html="results/fastp/{sample}.html",
        json="results/fastp/{sample}.json",
    log:
        "logs/fastp/{sample}.log",
    threads: config["fastp"]["threads"]
    resources:
        mem_mb=3000,
    conda:
        "../envs/fastp.yaml"
    shell:
        "fastp --detect_adapter_for_pe --thread {threads} "
        "--in1 {input.r1} --in2 {input.r2} "
        "--out1 {output.r1} --out2 {output.r2} "
        "--html {output.html} --json {output.json} > {log} 2>&1"


# Combine per-read FastQC and per-sample fastp reports into one summary page.
rule multiqc:
    input:
        expand("results/fastqc/{sample}_{read}_fastqc.zip", sample=samples.index, read=["R1", "R2"]),
        expand("results/fastp/{sample}.{ext}", sample=samples.index, ext=["html", "json"]),
    output:
        report="results/multiqc/multiqc_report.html",
    log:
        "logs/multiqc/multiqc.log",
    threads: 1
    resources:
        mem_mb=2000,
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc results/fastqc results/fastp --outdir results/multiqc > {log} 2>&1"
