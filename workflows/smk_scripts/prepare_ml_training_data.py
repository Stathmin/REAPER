import subprocess
from pathlib import Path


def main():
    fasta_files = [str(p) for p in snakemake.input.fasta_files]  # type: ignore[name-defined]
    training_data = Path(str(snakemake.output.training_data))  # type: ignore[name-defined]
    labels = Path(str(snakemake.output.labels))  # type: ignore[name-defined]

    training_data.parent.mkdir(parents=True, exist_ok=True)
    labels.parent.mkdir(parents=True, exist_ok=True)

    log_path = Path(str(snakemake.log[0]))  # type: ignore[name-defined]
    log_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "python3",
        "ml/prepare_training_data.py",
        "--input-files",
        *fasta_files,
        "--output-fasta",
        str(training_data),
        "--output-labels",
        str(labels),
    ]

    with log_path.open("w") as log:
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=False)

    if proc.returncode != 0:
        raise RuntimeError(f"prepare_training_data failed with exit code {proc.returncode}")


if __name__ == "__main__":
    main()

