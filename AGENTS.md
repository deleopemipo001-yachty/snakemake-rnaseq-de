# AGENTS.md — rules for this repository

## Project
A reproducible, containerized Snakemake pipeline for paired-end bulk RNA-Seq:
raw reads → QC/trimming → HISAT2 alignment → featureCounts → DESeq2 (paired design)
→ pathway enrichment → report. First test case: GSE142279 (colorectal tumor vs matched normal).
7. Never run large downloads or the full dataset yourself. You may WRITE rules that download
   reference files, but only I run them. Test only on .test/.

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
12. Run `snakefmt workflow/` before finishing every task, and make sure `snakefmt --check workflow/` passes.
13. temp() policy: delete only what can be rebuilt WITHOUT internet. Trimmed FASTQs and
    FASTQs converted from .sra -> temp(). Downloaded files (.sra, reference index, GTF) and
    the user's own input FASTQs -> never temp(); downloaded .sra files -> protected().
14. Reference files: URLs and versions live only in config. Downloads use curl with retries and
    resume (-C -) so an interrupted download continues instead of restarting. .test/ must never
    download the real reference; it uses a tiny local reference instead.
15. After editing any YAML file (config, schema, envs, CI), format it with Prettier 3.8.4 (conda run -n prettier npx --yes prettier@3.8.4 --write <file>) before finishing.
    
## Confirmed scientific decisions (do not change without being told)
- Library strandedness: 0 (unstranded), measured on RC1_T and RC1_N.
- Reference: HISAT2 prebuilt GRCh38 index (Ensembl-based) + Ensembl 116 GTF.
- DESeq2 design: ~ patient + condition; reference level "normal"; contrast tumor vs normal.
- Low-count filter: >= 10 reads in >= 6 samples (min_samples is a config value).
- Significant DEG: padj (BH) < 0.05 and |log2FoldChange| >= 1, using the unshrunken log2FC.
- Shrunken log2FC (apeglm) is reported alongside, and used for ranking and the MA plot.
- Gene names and biotypes come from the same Ensembl GTF used for counting.
- Enrichment: offline, pinned gene sets (msigdbr: KEGG + Hallmark), up and down tested
  separately, universe = all genes tested.

## Rules for R scripts
- Match count-matrix columns to the sample sheet BY NAME, never by position; stop with a
  clear error if any sample is missing or extra.
- Set a fixed random seed. No internet access at run time.
- Write sessionInfo() to the results folder.
- Read all thresholds from config; never hard-code them.
- Save every figure as both PNG (300 dpi) and PDF.