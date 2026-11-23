import csv


with open(snakemake.input[0], newline="") as sample_sheet:
    rows = csv.DictReader(sample_sheet, delimiter="\t")
    with open(snakemake.output[0], "w", newline="") as name_map:
        writer = csv.writer(name_map, delimiter="\t", lineterminator="\n")
        for row in rows:
            run = row["run"].strip()
            sample = row["sample"].strip()
            writer.writerow([run, sample])
            writer.writerow([f"{run}_1", f"{sample}_R1"])
            writer.writerow([f"{run}_2", f"{sample}_R2"])
