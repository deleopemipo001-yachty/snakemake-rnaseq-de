# AGENTS.md — rules for this repository

## Project
A reproducible, containerized Snakemake pipeline for paired-end bulk RNA-Seq:
raw reads → QC/trimming → HISAT2 alignment → featureCounts → DESeq2 (paired design)
→ pathway enrichment → report. First test case: GSE142279 (colorectal tumor vs matched normal).
It must be reusable: anyone runs it on their own data by editing config files only.

## Hardware limits
Runs on a laptop: 4 CPU threads, ~11 GB RAM (WSL2). Never request more than 4 threads
or more than 10000 mem_mb for any rule.

## Structure (official Snakemake workflow template)
- workflow/Snakefile — entry point, includes rule files
- workflow/rules/*.smk — rules grouped by stage (common.smk, qc.smk, ...)
- workflow/envs/*.yaml — ONE conda environment per tool
- workflow/scripts/ — Python/R scripts
- workflow/schemas/ — schemas for config and sample sheet
- config/config.yaml and config/samples.tsv — all user settings
- .test/ — tiny test dataset and its config; the whole workflow must run on it
  with `snakemake --directory .test`

## Rules you must follow
1. Every rule has: log:, threads:, resources: mem_mb, and a short comment explaining
   what it does and why, in plain language.
2. Every conda environment pins exact tool versions. Channels: conda-forge, bioconda.
3. No hardcoded paths or sample names. Everything comes from config and samples.tsv,
   validated against schemas with snakemake.utils.validate.
4. Large intermediate files (raw FASTQ, unsorted BAM) use temp().
5. Plain shell or script rules. Do NOT use Snakemake wrappers.
6. Do NOT make scientific decisions. If something affects results (library strandedness,
   reference genome/annotation versions, thresholds, DESeq2 design, enrichment settings),
   do not pick a value: add a clearly marked config key with a TODO and tell me in your reply.
7. Never download large files or run the full dataset. Test only on .test/.
8. Never commit data, results/, resources/, logs/ or .snakemake/ (keep .gitignore updated).
   The tiny files in .test/data ARE committed.
9. Work on one task at a time. Don't change files unrelated to the task.
10. Do NOT commit or push. I review, run and commit everything myself.
11. Before finishing, make sure `snakemake --lint` passes and the .test dry run works.
    End every reply with: files changed, and the exact commands I should run to verify.