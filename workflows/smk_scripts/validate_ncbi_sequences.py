import subprocess
from pathlib import Path


def main():
    fasta = Path(str(snakemake.input.fasta))  # type: ignore[name-defined]
    metadata = Path(str(snakemake.input.metadata))  # type: ignore[name-defined]
    validation_report = Path(str(snakemake.output.validation_report))  # type: ignore[name-defined]

    validation_report.parent.mkdir(parents=True, exist_ok=True)

    log_path = Path(str(snakemake.log[0]))  # type: ignore[name-defined]
    log_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "python3",
        "post_tarean/validate_ncbi_sequences.py",
        "--fasta",
        str(fasta),
        "--metadata",
        str(metadata),
        "--output",
        str(validation_report),
    ]

    with log_path.open("w") as log:
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=False)

    if proc.returncode != 0:
        raise RuntimeError(f"validate_ncbi_sequences failed with exit code {proc.returncode}")


if __name__ == "__main__":
    main()

