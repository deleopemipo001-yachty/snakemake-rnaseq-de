# Check each raw mate so read quality issues are visible before trimming.
rule fastqc:
    input:
        fastq=lambda wildcards: sample_fastq(wildcards.sample, wildcards.read),
    output:
        html="results/fastqc/{sample}_{read}_fastqc.html",
        zip="results/fastqc/{sample}_{read}_fastqc.zip",
    log:
        "logs/fastqc/{sample}_{read}.log",
    priority: 0
    conda:
        "../envs/fastqc.yaml"
    threads: 1
    resources:
        mem_mb=1000,
    shell:
        "tmpdir=$(mktemp -d results/fastqc/.fastqc.XXXXXX); "
        "trap 'rm -rf \"$tmpdir\"' EXIT; "
        'fastqc --threads {threads} --outdir "$tmpdir" {input.fastq} > {log} 2>&1 && '
        "find \"$tmpdir\" -name '*_fastqc.html' -exec mv {{}} {output.html} ';' && "
        "find \"$tmpdir\" -name '*_fastqc.zip' -exec mv {{}} {output.zip} ';'"


# Trim paired reads with fastp's defaults and automatic paired-end adapter detection.
rule fastp:
    input:
        r1=lambda wildcards: sample_fastq(wildcards.sample, "R1"),
        r2=lambda wildcards: sample_fastq(wildcards.sample, "R2"),
    output:
        r1=temp("results/trimmed/{sample}_R1.fastq.gz"),
        r2=temp("results/trimmed/{sample}_R2.fastq.gz"),
        html="results/fastp/{sample}.html",
        json="results/fastp/{sample}.json",
    log:
        "logs/fastp/{sample}.log",
    priority: 5
    conda:
        "../envs/fastp.yaml"
    threads: config["fastp"]["threads"]
    resources:
        mem_mb=3000,
    shell:
        "fastp --detect_adapter_for_pe --thread {threads} "
        "--in1 {input.r1} --in2 {input.r2} "
        "--out1 {output.r1} --out2 {output.r2} "
        "--html {output.html} --json {output.json} > {log} 2>&1"


# Build the sample-name replacements that make run-based QC names readable.
rule multiqc_sample_name_map:
    input:
        config["sample_sheet"],
    output:
        "results/multiqc/sample_name_map.tsv",
    log:
        "logs/multiqc/sample_name_map.log",
    conda:
        "../envs/multiqc.yaml"
    threads: 1
    resources:
        mem_mb=1000,
    script:
        "../scripts/make_multiqc_name_map.py"


# Combine per-read FastQC and per-sample fastp reports into one summary page.
rule multiqc:
    input:
        name_map="results/multiqc/sample_name_map.tsv",
        reports=(
            expand(
                "results/fastqc/{sample}_{read}_fastqc.zip",
                sample=samples.index,
                read=["R1", "R2"],
            )
            + expand(
                "results/fastp/{sample}.{ext}",
                sample=samples.index,
                ext=["html", "json"],
            )
            + expand("logs/hisat2/{sample}.summary.txt", sample=samples.index)
            + expand("results/qc/flagstat/{sample}.flagstat", sample=samples.index)
            + ["results/counts/counts.txt.summary"]
        ),
    output:
        report="results/multiqc/multiqc_report.html",
    log:
        "logs/multiqc/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    threads: 1
    resources:
        mem_mb=2000,
    shell:
        "tmpdir=$(mktemp -d results/multiqc/.multiqc.XXXXXX); "
        "trap 'rm -rf \"$tmpdir\"' EXIT; "
        'multiqc --replace-names {input.name_map:q} --outdir "$tmpdir" '
        "--filename multiqc_report.html {input.reports:q} > {log:q} 2>&1 && "
        'mv "$tmpdir/multiqc_report.html" {output.report:q}'
