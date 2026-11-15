# Configuration

Edit `config.yaml` to point to a tab-separated sample sheet and set the fastp thread count (1–4). fastp uses default filtering and trimming settings, with paired-end adapter detection enabled by the workflow.

The sample sheet must contain `sample`, `patient`, `condition`, `run`, `fq1`, and `fq2`. `condition` accepts `tumor` or `normal`; FASTQ paths are resolved from the workflow's working directory.
