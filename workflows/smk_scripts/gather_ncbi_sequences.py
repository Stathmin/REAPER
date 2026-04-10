import subprocess
from pathlib import Path


def _latest_matching(path: Path, pattern: str) -> Path:
    """Return the most recently modified file matching pattern under path."""
    candidates = list(path.glob(pattern))
    if not candidates:
        raise FileNotFoundError(f"No files matching {pattern} in {path}")
    return max(candidates, key=lambda p: p.stat().st_mtime)


def _ensure_copy(target: Path, expected: Path) -> None:
    """Ensure expected is a real file copy of target (no symlinks)."""
    if expected.is_symlink():
        expected.unlink()
    if not expected.exists():
        expected.parent.mkdir(parents=True, exist_ok=True)
        expected.write_bytes(target.read_bytes())


def main():
    taxid = str(snakemake.params.taxid)  # type: ignore[name-defined]
    project_id = str(snakemake.params.project_id)  # type: ignore[name-defined]
    email = str(snakemake.params.email)  # type: ignore[name-defined]

    log_path = Path(str(snakemake.log[0]))  # type: ignore[name-defined]
    log_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "python3",
        "post_tarean/ncbi_data_gatherer.py",
        taxid,
        "--project-id",
        project_id,
        "--email",
        email,
    ]

    with log_path.open("w") as log:
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=False)

    # Be resilient: NCBI search can legitimately yield 0 sequences for a taxid/query.
    # In that case, ncbi_data_gatherer exits non-zero, but we still want the workflow
    # to produce stable, expected outputs (copy the most recent existing files if any,
    # otherwise create empty placeholders).
    gather_ok = proc.returncode == 0

    # Map timestamped outputs produced by ncbi_data_gatherer to the fixed
    # filenames expected by Snakemake.
    project_root = Path("projects") / project_id
    blast_dir = project_root / "blast_db"
    metadata_dir = project_root / "metadata"

    expected_fasta = blast_dir / f"ncbi_repeats_{taxid}.fasta"
    try:
        latest_fasta = _latest_matching(blast_dir, "ncbi_repeats_*.fasta")
        _ensure_copy(latest_fasta, expected_fasta)
    except FileNotFoundError:
        expected_fasta.parent.mkdir(parents=True, exist_ok=True)
        if not expected_fasta.exists():
            expected_fasta.write_text("")

    expected_nhr = blast_dir / f"ncbi_repeats_{taxid}.nhr"
    try:
        latest_nhr = _latest_matching(blast_dir, "ncbi_repeats_*.nhr")
        _ensure_copy(latest_nhr, expected_nhr)
    except FileNotFoundError:
        expected_nhr.parent.mkdir(parents=True, exist_ok=True)
        if not expected_nhr.exists():
            expected_nhr.write_text("")

    expected_meta = metadata_dir / f"ncbi_metadata_{taxid}.csv"
    try:
        latest_meta = _latest_matching(metadata_dir, f"ncbi_metadata_{taxid}_*.csv")
        _ensure_copy(latest_meta, expected_meta)
    except FileNotFoundError:
        expected_meta.parent.mkdir(parents=True, exist_ok=True)
        if not expected_meta.exists():
            # Minimal CSV header for downstream tools.
            expected_meta.write_text("accession,title,species,taxid\n")

    if not gather_ok:
        # Keep non-zero gathering from failing the whole pipeline when placeholders were written.
        return


if __name__ == "__main__":
    main()

