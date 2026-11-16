import subprocess
import tempfile
from pathlib import Path


output_r1 = Path(snakemake.output.r1)
output_r2 = Path(snakemake.output.r2)
output_dir = output_r1.parent
output_dir.mkdir(parents=True, exist_ok=True)
temp_root = Path("data/tmp/fasterq-dump")
temp_root.mkdir(parents=True, exist_ok=True)

# Remove incomplete files from a prior failed conversion before regenerating them.
for output in (output_r1, output_r2):
    output.unlink(missing_ok=True)

log_lines = []
with tempfile.TemporaryDirectory(
    prefix=f"{snakemake.wildcards.run}-", dir=temp_root
) as temp_dir:
    command = [
        "fasterq-dump",
        "--split-3",
        "--threads",
        str(snakemake.threads),
        "--temp",
        temp_dir,
        "--outdir",
        str(output_dir),
        snakemake.input.sra,
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    log_lines.extend(
        [f"Command: {' '.join(command)}", result.stdout, result.stderr]
    )
    if result.returncode:
        Path(snakemake.log[0]).write_text("\n".join(log_lines))
        raise Exception(
            f"fasterq-dump failed for run {snakemake.wildcards.run}; "
            f"see {snakemake.log[0]} for details."
        )


def read_count(fastq_path):
    with fastq_path.open("rb") as fastq_file:
        line_count = sum(1 for _ in fastq_file)
    if line_count % 4:
        raise Exception(f"FASTQ file has an incomplete record: {fastq_path}")
    return line_count // 4


if not output_r1.is_file() or not output_r2.is_file():
    raise Exception(
        f"fasterq-dump did not produce both paired FASTQs for run "
        f"{snakemake.wildcards.run}; expected {output_r1} and {output_r2}."
    )

reads_r1 = read_count(output_r1)
reads_r2 = read_count(output_r2)
log_lines.append(f"R1 read pairs: {reads_r1}")
log_lines.append(f"R2 read pairs: {reads_r2}")
if reads_r1 != reads_r2:
    output_r1.unlink(missing_ok=True)
    output_r2.unlink(missing_ok=True)
    Path(snakemake.log[0]).write_text("\n".join(log_lines))
    raise Exception(
        f"Read-pair count mismatch for run {snakemake.wildcards.run}: "
        f"R1 has {reads_r1} reads and R2 has {reads_r2} reads."
    )

Path(snakemake.log[0]).write_text("\n".join(log_lines))
