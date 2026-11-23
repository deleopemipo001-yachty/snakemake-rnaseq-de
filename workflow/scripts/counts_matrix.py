import csv
from pathlib import Path


sample_names = list(snakemake.params.sample_names)
with open(snakemake.input[0], newline="") as counts_file:
    reader = csv.reader(counts_file, delimiter="\t")
    header = next((row for row in reader if row and row[0] == "Geneid"), None)
    if header is None:
        raise ValueError("featureCounts output has no Geneid header row")

    metadata_columns = 6
    bam_columns = {
        Path(column).stem: index
        for index, column in enumerate(header[metadata_columns:], start=metadata_columns)
    }
    missing_samples = [sample for sample in sample_names if sample not in bam_columns]
    if missing_samples:
        raise ValueError(
            "featureCounts output is missing BAM columns for samples: "
            + ", ".join(missing_samples)
        )

    with open(snakemake.output[0], "w", newline="") as matrix_file:
        writer = csv.writer(matrix_file, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", *sample_names])
        for row in reader:
            if row:
                writer.writerow([row[0], *(row[bam_columns[sample]] for sample in sample_names)])
