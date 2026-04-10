from pathlib import Path


def main():
    fasta_files = list(snakemake.input.fasta_files)  # type: ignore[name-defined]
    validation_reports = list(snakemake.input.validation_reports)  # type: ignore[name-defined]
    taxonomy_summary = Path(str(snakemake.input.taxonomy_summary))  # type: ignore[name-defined]
    ml_data = Path(str(snakemake.input.ml_data))  # type: ignore[name-defined]

    completion_marker = Path(str(snakemake.output.completion_marker))  # type: ignore[name-defined]
    completion_marker.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        f"NCBI data gathering completed for project {snakemake.wildcards.project}",  # type: ignore[attr-defined]
        f"Generated {len(fasta_files)} FASTA files",
        f"Generated {len(validation_reports)} validation reports",
        f"Taxonomy summary: {taxonomy_summary}",
        f"ML training data: {ml_data}",
    ]

    completion_marker.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

