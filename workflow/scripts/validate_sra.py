import subprocess
from pathlib import Path


command = ["vdb-validate", snakemake.input.sra]
result = subprocess.run(command, capture_output=True, text=True)
details = result.stdout + result.stderr
Path(snakemake.output[0]).write_text(details)
Path(snakemake.log[0]).write_text(
    f"Command: {' '.join(command)}\nExit code: {result.returncode}\n{details}"
)
if result.returncode:
    raise Exception(
        f"SRA validation failed for run {snakemake.wildcards.run}; "
        f"conversion was stopped. See {snakemake.log[0]} and {snakemake.output[0]}."
    )
