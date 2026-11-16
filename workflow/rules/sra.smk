import os


# Download each SRA run into its accession directory and preserve the source file.
rule prefetch_sra:
    output:
        sra=protected("data/sra/{run}/{run}.sra"),
    log:
        "logs/sra/prefetch/{run}.log",
    conda:
        "../envs/sra-tools.yaml"
    threads: 1
    resources:
        mem_mb=1000,
    params:
        output_dir=lambda wildcards, output: os.path.dirname(output.sra),
    shell:
        "mkdir -p {params.output_dir:q} && "
        "for attempt in 1 2 3; do "
        "if prefetch --max-size 10G --output-directory {params.output_dir:q} "
        "{wildcards.run:q} >> {log:q} 2>&1; then break; fi; "
        'if [ "$attempt" -eq 3 ]; then '
        "echo 'prefetch failed after 3 attempts; see this log for details.' >> {log:q}; "
        'exit 1; fi; sleep "$attempt"; done && '
        "test -s {output.sra:q}"


# Validate each downloaded SRA file before it is converted to FASTQ.
rule validate_sra:
    input:
        sra=rules.prefetch_sra.output.sra,
    output:
        "results/qc/sra_validate/{run}.txt",
    log:
        "logs/sra/validate/{run}.log",
    conda:
        "../envs/sra-tools.yaml"
    threads: 1
    resources:
        mem_mb=2000,
    script:
        "../scripts/validate_sra.py"


# Convert a validated SRA run to paired FASTQs and reject mismatched mate counts.
rule fasterq_dump:
    input:
        sra=rules.prefetch_sra.output.sra,
        validation="results/qc/sra_validate/{run}.txt",
    output:
        r1=temp("data/fastq/{run}_1.fastq"),
        r2=temp("data/fastq/{run}_2.fastq"),
    log:
        "logs/sra/fasterq_dump/{run}.log",
    conda:
        "../envs/sra-tools.yaml"
    threads: 4
    resources:
        mem_mb=8000,
    script:
        "../scripts/fasterq_dump.py"
