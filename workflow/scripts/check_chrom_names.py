import gzip
import subprocess
from pathlib import Path


log_path = Path(snakemake.log[0])
inspection = subprocess.run(
    ["hisat2-inspect", "-n", snakemake.params.prefix],
    capture_output=True,
    text=True,
)
if inspection.returncode:
    message = (
        f"Could not list chromosome names from HISAT2 index "
        f"{snakemake.params.prefix}: {inspection.stderr.strip()}"
    )
    log_path.write_text(message + "\n")
    raise Exception(message)

index_names = {
    line.split(maxsplit=1)[0]
    for line in inspection.stdout.splitlines()
    if line.split(maxsplit=1)
}
open_gtf = gzip.open if snakemake.input.gtf.endswith(".gz") else open
gtf_names = set()
with open_gtf(snakemake.input.gtf, "rt") as gtf_file:
    for line in gtf_file:
        if line and not line.startswith("#"):
            gtf_names.add(line.split("\t", 1)[0])

expected = [str(number) for number in range(1, 23)] + ["X", "Y"]
errors = []
for chromosome in expected:
    aliases = (chromosome, f"chr{chromosome}")
    if not any(name in index_names for name in aliases):
        errors.append(f"{chromosome} missing from HISAT2 index")
    if not any(name in gtf_names for name in aliases):
        errors.append(f"{chromosome} missing from GTF")
    if not index_names.intersection(gtf_names).intersection(aliases):
        errors.append(
            f"{chromosome} has different chromosome names between index and GTF"
        )

missing_gtf = (
    sorted(gtf_names - index_names) if snakemake.params.mode == "local" else []
)
if missing_gtf:
    errors.append(
        "GTF chromosomes missing from local HISAT2 index: " + ", ".join(missing_gtf)
    )

if errors:
    message = "Reference chromosome-name check failed: " + "; ".join(errors)
    log_path.write_text(message + "\n")
    raise Exception(message)

report = [
    "HISAT2 index chromosome names:",
    *sorted(index_names),
    "",
    "GTF chromosome names:",
    *sorted(gtf_names),
]
Path(snakemake.output[0]).write_text("\n".join(report) + "\n")
log_path.write_text(
    f"Chromosome names agree for chromosomes 1-22, X, and Y. "
    f"Mode: {snakemake.params.mode}.\n"
)
