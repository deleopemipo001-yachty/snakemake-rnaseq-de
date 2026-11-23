# Align each trimmed read pair, sort it by position, and build a BAM index.
rule align:
    input:
        r1="results/trimmed/{sample}_R1.fastq.gz",
        r2="results/trimmed/{sample}_R2.fastq.gz",
        index=index_directory,
        chrom_names="results/qc/chrom_names.txt",
    output:
        bam="results/bam/{sample}.bam",
        bai="results/bam/{sample}.bam.bai",
        summary="logs/hisat2/{sample}.summary.txt",
    log:
        align="logs/hisat2/{sample}.align.log",
        sort="logs/hisat2/{sample}.sort.log",
    priority: 10
    conda:
        "../envs/alignment.yaml"
    threads: 4
    resources:
        mem_mb=9000,
    params:
        index_prefix=index_prefix,
    shell:
        "if [ ! -f {params.index_prefix:q}.1.ht2 ] && "
        "[ ! -f {params.index_prefix:q}.1.ht2l ]; then "
        "echo 'HISAT2 index files are missing for prefix {params.index_prefix}; "
        "expected .1.ht2 or .1.ht2l.' >&2; exit 1; fi && "
        "hisat2 -p {threads} --new-summary --summary-file {output.summary:q} "
        "-x {params.index_prefix:q} -1 {input.r1:q} -2 {input.r2:q} "
        "2> {log.align:q} | "
        "samtools sort -@ 2 -m 512M -o {output.bam:q} - "
        "2> {log.sort:q} && "
        "samtools index -@ 2 {output.bam:q} 2>> {log.sort:q}"


# Summarize mapping statistics for each sorted BAM so MultiQC can report them.
rule flagstat:
    input:
        bam="results/bam/{sample}.bam",
    output:
        "results/qc/flagstat/{sample}.flagstat",
    log:
        "logs/flagstat/{sample}.log",
    conda:
        "../envs/alignment.yaml"
    threads: 2
    resources:
        mem_mb=1000,
    shell:
        "samtools flagstat -@ {threads} {input.bam:q} " "> {output:q} 2> {log:q}"
