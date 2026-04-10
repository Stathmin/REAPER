import subprocess
from pathlib import Path


def main():
    metadata = Path(str(snakemake.input.metadata))  # type: ignore[name-defined]
    freshness_report = Path(str(snakemake.output.freshness_report))  # type: ignore[name-defined]
    max_age_days = int(snakemake.params.max_age_days)  # type: ignore[name-defined]

    freshness_report.parent.mkdir(parents=True, exist_ok=True)

    log_path = Path(str(snakemake.log[0]))  # type: ignore[name-defined]
    log_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "python3",
        "post_tarean/check_ncbi_freshness.py",
        "--metadata",
        str(metadata),
        "--max-age-days",
        str(max_age_days),
        "--output",
        str(freshness_report),
    ]

    with log_path.open("w") as log:
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=False)

    if proc.returncode != 0:
        raise RuntimeError(f"check_ncbi_freshness failed with exit code {proc.returncode}")


if __name__ == "__main__":
    main()

