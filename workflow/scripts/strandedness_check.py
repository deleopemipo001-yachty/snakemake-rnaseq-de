import subprocess
import tempfile
from pathlib import Path


log_lines = []
rows = ["strandedness\tassigned_read_pairs\tpercent_assigned"]

with tempfile.TemporaryDirectory(prefix="featurecounts-strandedness-") as temp_dir:
    temp_path = Path(temp_dir)
    for strandedness in (0, 1, 2):
        output_prefix = temp_path / f"counts_s{strandedness}.txt"
        command = [
            "featureCounts",
            "-p",
            "--countReadPairs",
            "-t",
            "exon",
            "-g",
            "gene_id",
            "-T",
            str(snakemake.threads),
            "-s",
            str(strandedness),
            "-a",
            snakemake.input.gtf,
            "-o",
            str(output_prefix),
            snakemake.input.bam,
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        log_lines.extend(
            [
                f"Command: {' '.join(command)}",
                result.stdout,
                result.stderr,
            ]
        )
        if result.returncode:
            Path(snakemake.log[0]).write_text("\n".join(log_lines))
            raise Exception(
                f"featureCounts failed for strandedness -s {strandedness}: "
                f"{result.stderr.strip()}"
            )

        assigned = 0
        total = 0
        summary_path = Path(f"{output_prefix}.summary")
        with summary_path.open() as summary:
            next(summary, None)
            for line in summary:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 2 and fields[-1].isdigit():
                    count = int(fields[-1])
                    total += count
                    if fields[0] == "Assigned":
                        assigned = count

        percentage = (100 * assigned / total) if total else 0.0
        rows.append(f"{strandedness}\t{assigned}\t{percentage:.2f}")

Path(snakemake.output[0]).write_text("\n".join(rows) + "\n")
Path(snakemake.log[0]).write_text("\n".join(log_lines))
