# Paired-end RNA-seq workflow

A reusable Snakemake workflow for paired-end bulk RNA-seq. The current workflow skeleton performs raw-read quality checks with FastQC, paired-end trimming with fastp, and report aggregation with MultiQC.

## Configure and run

Edit `config/config.yaml` and `config/samples.tsv` for your input data. The sample sheet records each sample's patient, condition, sequencing run, and paired FASTQ paths. Conditions must be `tumor` or `normal`.

Run the workflow with:

```bash
snakemake --cores 4 --sdm conda
```

The test configuration and tiny FASTQ files are in `.test/`. Check the workflow plan with:

```bash
snakemake --directory .test --dry-run
```

The final report is written to `results/multiqc/multiqc_report.html` within the selected working directory.
