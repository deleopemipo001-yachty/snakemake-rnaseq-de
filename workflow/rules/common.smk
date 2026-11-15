import pandas as pd
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
