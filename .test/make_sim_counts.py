#!/usr/bin/env python3
"""Create deterministic paired tumor/normal test counts from the real Ensembl GTF."""

import argparse
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd


SEED = 20261127
ATTR_PATTERN = re.compile(r'([A-Za-z0-9_]+) "([^"]+)"')


def read_protein_coding_genes(gtf_path):
    opener = gzip.open if gtf_path.suffix == ".gz" else open
    genes = {}
    with opener(gtf_path, "rt") as gtf:
        for line in gtf:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attributes = dict(ATTR_PATTERN.findall(fields[8]))
            biotype = attributes.get("gene_biotype", attributes.get("gene_type"))
            if biotype != "protein_coding":
                continue
            gene_id = attributes.get("gene_id")
            if gene_id:
                genes[gene_id] = attributes.get("gene_name", gene_id)
    return sorted(genes.items())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    project_root = Path(__file__).resolve().parents[1]
    parser.add_argument(
        "--gtf",
        type=Path,
        default=project_root / "resources/Homo_sapiens.GRCh38.116.gtf.gz",
        help="Ensembl GTF used for counting (default: downloaded GRCh38.116 GTF)",
    )
    parser.add_argument("--outdir", type=Path, default=Path(__file__).parent / "data")
    args = parser.parse_args()

    genes = read_protein_coding_genes(args.gtf)
    if len(genes) < 2000:
        raise SystemExit(
            f"Expected at least 2000 protein-coding genes in {args.gtf}; found {len(genes)}"
        )
    genes = genes[:2000]
    gene_ids = [gene_id for gene_id, _ in genes]
    gene_names = [gene_name for _, gene_name in genes]

    samples = []
    for patient_number in range(1, 4):
        patient = f"SIM{patient_number}"
        for condition, suffix in (("normal", "N"), ("tumor", "T")):
            samples.append(
                {
                    "sample": f"{patient}_{suffix}",
                    "patient": patient,
                    "condition": condition,
                    "run": f"{patient}_{suffix}",
                    "fq1": "",
                    "fq2": "",
                }
            )
    sample_table = pd.DataFrame(samples)

    rng = np.random.default_rng(SEED)
    gene_count = len(gene_ids)
    baseline_means = rng.lognormal(mean=4.0, sigma=1.0, size=gene_count)
    dispersions = np.full(gene_count, 0.15)
    patient_names = [f"SIM{i}" for i in range(1, 4)]
    patient_effects = dict(zip(patient_names, rng.lognormal(0, 0.25, len(patient_names))))
    changed_indices = rng.choice(gene_count, size=round(gene_count * 0.10), replace=False)
    half = len(changed_indices) // 2
    fold_changes = np.ones(gene_count)
    fold_changes[changed_indices[:half]] = 2.5
    fold_changes[changed_indices[half:]] = 0.4

    count_matrix = np.zeros((gene_count, len(sample_table)), dtype=np.int64)
    for sample_index, sample in enumerate(samples):
        tumor_factor = fold_changes if sample["condition"] == "tumor" else 1.0
        means = baseline_means * patient_effects[sample["patient"]] * tumor_factor
        size = 1.0 / dispersions
        probabilities = size / (size + means)
        count_matrix[:, sample_index] = rng.negative_binomial(size, probabilities)

    args.outdir.mkdir(parents=True, exist_ok=True)
    counts = pd.DataFrame(count_matrix, columns=sample_table["sample"])
    counts.insert(0, "gene_id", gene_ids)
    counts.to_csv(args.outdir / "sim_counts_matrix.tsv", sep="\t", index=False)
    sample_table.to_csv(args.outdir / "sim_samples.tsv", sep="\t", index=False)
    annotation = pd.DataFrame(
        {
            "gene_id": gene_ids,
            "gene_name": gene_names,
            "gene_biotype": "protein_coding",
        }
    )
    annotation.to_csv(args.outdir / "sim_gene_annotation.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
