import pandas as pd
from snakemake.exceptions import WorkflowError
from snakemake.utils import validate

# Validate settings and metadata before any workflow rules are scheduled.
validate(config, "../schemas/config.schema.yaml")
samples = pd.read_csv(config["sample_sheet"], sep="\t", dtype=str)
for column in ["fq1", "fq2"]:
    if column not in samples.columns:
        samples[column] = ""
    samples[column] = samples[column].fillna("").str.strip()
validate(samples, "../schemas/samples.schema.yaml")
de_samples = pd.read_csv(config["de"]["samples"], sep="\t", dtype=str)
for column in ["fq1", "fq2"]:
    if column not in de_samples.columns:
        de_samples[column] = ""
    de_samples[column] = de_samples[column].fillna("").str.strip()
validate(de_samples, "../schemas/samples.schema.yaml")
partial_fastq_paths = samples["fq1"].ne("") != samples["fq2"].ne("")
if partial_fastq_paths.any():
    invalid_samples = ", ".join(samples.loc[partial_fastq_paths, "sample"])
    raise WorkflowError(
        "Each sample must provide both fq1 and fq2 paths, or leave both empty "
        f"to download reads from its SRA run accession. Invalid samples: {invalid_samples}"
    )
samples = samples.set_index("sample", drop=False).sort_index()
sra_runs = sorted(samples.loc[samples["fq1"].eq(""), "run"].unique())


def sample_fastq(sample, read):
    """Return a local FASTQ path or the FASTQ generated from the sample's SRA run."""
    row = samples.loc[sample]
    if row["fq1"] and row["fq2"]:
        column = "fq1" if read == "R1" else "fq2"
        return row[column]
    mate = "1" if read == "R1" else "2"
    return f"data/fastq/{row['run']}_{mate}.fastq"


def require_featurecounts_strandedness(wildcards):
    strandedness = config["featurecounts"]["strandedness"]
    if strandedness is None:
        raise WorkflowError(
            "featureCounts was requested, but featurecounts.strandedness is null. "
            "Run strandedness_check for a sample, review its assigned-pair "
            "percentages, then set featurecounts.strandedness to 0, 1, or 2."
        )
    return strandedness
