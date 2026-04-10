import subprocess
from pathlib import Path


def main():
    metadata_files = [str(p) for p in snakemake.input.metadata_files]  # type: ignore[name-defined]
    summary = Path(str(snakemake.output.summary))  # type: ignore[name-defined]

    summary.parent.mkdir(parents=True, exist_ok=True)

    log_path = Path(str(snakemake.log[0]))  # type: ignore[name-defined]
    log_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "python3",
        "post_tarean/create_taxonomy_summary.py",
        "--input-files",
        *metadata_files,
        "--output",
        str(summary),
    ]

    with log_path.open("w") as log:
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=False)

    if proc.returncode != 0:
        raise RuntimeError(f"create_taxonomy_summary failed with exit code {proc.returncode}")


if __name__ == "__main__":
    main()

