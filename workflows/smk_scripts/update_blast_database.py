import subprocess
from pathlib import Path


def main():
    fasta = Path(str(snakemake.input.fasta))  # type: ignore[name-defined]
    updated_db = Path(str(snakemake.output.updated_db))  # type: ignore[name-defined]

    # Derive project-specific blast_db directory from output path
    blast_db_dir = updated_db.parent
    blast_db_dir.mkdir(parents=True, exist_ok=True)

    combined_fasta = blast_db_dir / "combined_sequences.fasta"

    log_path = Path(str(snakemake.log[0]))  # type: ignore[name-defined]
    log_path.parent.mkdir(parents=True, exist_ok=True)

    # Base path without suffix for makeblastdb -out (it appends .nhr, .nin, .nsq)
    db_base = updated_db.with_suffix("")

    with log_path.open("w") as log:
        # Append new sequences to combined FASTA (replicates original behaviour)
        proc_cat = subprocess.run(
            ["bash", "-c", f"cat {fasta} >> {combined_fasta}"],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
        if proc_cat.returncode != 0:
            raise RuntimeError(f"Failed to append {fasta} to {combined_fasta}")

        # Rebuild BLAST database
        cmd = [
            "makeblastdb",
            "-in",
            str(combined_fasta),
            "-dbtype",
            "nucl",
            "-out",
            str(db_base),
            "-title",
            "Updated_Repeat_Database",
        ]
        proc_db = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=False)

    if proc_db.returncode != 0:
        raise RuntimeError(f"makeblastdb failed with exit code {proc_db.returncode}")


if __name__ == "__main__":
    main()

