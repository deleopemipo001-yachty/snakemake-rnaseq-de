import pandas as pd
from snakemake.exceptions import WorkflowError
from snakemake.utils import validate

# Validate settings and metadata before any workflow rules are scheduled.
validate(config, "../schemas/config.schema.yaml")
samples = pd.read_csv(config["sample_sheet"], sep="\t", dtype=str)
validate(samples, "../schemas/samples.schema.yaml")
samples = samples.set_index("sample", drop=False).sort_index()


def sample_fastq(sample, read):
    """Return the configured raw FASTQ path for a sample and mate."""
    column = "fq1" if read == "R1" else "fq2"
    return samples.loc[sample, column]


def require_featurecounts_strandedness(wildcards):
    strandedness = config["featurecounts"]["strandedness"]
    if strandedness is None:
        raise WorkflowError(
            "featureCounts was requested, but featurecounts.strandedness is null. "
            "Run strandedness_check for a sample, review its assigned-pair "
            "percentages, then set featurecounts.strandedness to 0, 1, or 2."
        )
    return strandedness
