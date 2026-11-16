import os

ref_config = config["ref"]
index_prefix = os.path.join("resources", ref_config["hisat2_index_prefix"])
index_directory = os.path.dirname(index_prefix)
if ref_config["mode"] == "download":
    gtf_path = os.path.join(
        "resources",
        f"Homo_sapiens.GRCh38.{ref_config['ensembl_release']}.gtf.gz",
    )
else:
    gtf_path = ref_config["gtf"]


if ref_config["mode"] == "download":

    # Download and extract the configured HISAT2 archive; keep it if any step fails.
    rule download_hisat2_index:
        output:
            directory(index_directory),
        log:
            "logs/reference/hisat2_index.log",
        conda:
            "../envs/curl.yaml"
        threads: 1
        resources:
            mem_mb=1000,
        params:
            url=ref_config["hisat2_index_url"],
            archive="resources/grch38_genome.tar.gz",
        shell:
            "mkdir -p resources && "
            "curl --fail --location --retry 5 --retry-all-errors --retry-delay 2 "
            "--continue-at - {params.url:q} --output {params.archive:q} > {log} 2>&1 && "
            "tar -xzf {params.archive:q} -C resources >> {log} 2>&1 && "
            "rm {params.archive:q}"

    # Download the configured Ensembl annotation and retain its compressed source file.
    rule download_gtf:
        output:
            gtf=gtf_path,
        log:
            "logs/reference/gtf.log",
        conda:
            "../envs/curl.yaml"
        threads: 1
        resources:
            mem_mb=1000,
        params:
            url=ref_config["gtf_url"],
        shell:
            "mkdir -p resources && "
            "curl --fail --location --retry 5 --retry-all-errors --retry-delay 2 "
            "--continue-at - {params.url:q} --output {output.gtf:q} > {log} 2>&1"

else:

    # Build a HISAT2 index from the FASTA supplied for local mode.
    rule build_hisat2_index:
        input:
            fasta=ref_config["fasta"],
        output:
            directory(index_directory),
        log:
            "logs/reference/hisat2_index.log",
        conda:
            "../envs/hisat2.yaml"
        threads: 4
        resources:
            mem_mb=8000,
        params:
            prefix=index_prefix,
        shell:
            "mkdir -p {output:q} && "
            "hisat2-build --large-index --threads {threads} "
            "{input.fasta:q} {params.prefix:q} > {log} 2>&1"


# Compare the index and annotation names before downstream rules use the reference.
rule check_chrom_names:
    input:
        index_dir=index_directory,
        gtf=gtf_path,
    output:
        "results/qc/chrom_names.txt",
    log:
        "logs/reference/chrom_names.log",
    conda:
        "../envs/hisat2.yaml"
    threads: 1
    resources:
        mem_mb=1000,
    params:
        prefix=index_prefix,
        mode=ref_config["mode"],
    script:
        "../scripts/check_chrom_names.py"


# Make the selected index, annotation, and name-check report an easy target.
rule references:
    input:
        index=index_directory,
        gtf=gtf_path,
        chrom_names="results/qc/chrom_names.txt",
    log:
        "logs/reference/references.log",
    conda:
        "../envs/hisat2.yaml"
    threads: 1
    resources:
        mem_mb=100,
    shell:
        "touch {log}"
