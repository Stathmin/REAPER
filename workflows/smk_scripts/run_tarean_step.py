import os
import argparse
import shutil
import hashlib
import json
import subprocess
import time
from pathlib import Path
from typing import Optional


def _yaml_scalar(v):
    if v is None:
        return "null"
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return str(v)
    s = str(v)
    # Quote scalars that may confuse YAML parsing.
    if s == "" or any(ch in s for ch in [":", "#", "{", "}", "[", "]", ",", "\n", "\r", "\t", "\"", "'"]) or s.strip() != s:
        s = s.replace("\\", "\\\\").replace("\"", "\\\"")
        return f"\"{s}\""
    return s


def _write_yaml_atomic(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    lines = []
    for k in sorted(data.keys()):
        v = data[k]
        if isinstance(v, dict):
            lines.append(f"{k}:")
            for kk in sorted(v.keys()):
                lines.append(f"  {kk}: {_yaml_scalar(v[kk])}")
        else:
            lines.append(f"{k}: {_yaml_scalar(v)}")
    tmp.write_text("\n".join(lines) + "\n")
    tmp.replace(path)


def _parse_yaml_scalar(raw: str):
    s = raw.strip()
    if s == "null":
        return None
    if s == "true":
        return True
    if s == "false":
        return False
    if len(s) >= 2 and s[0] == "\"" and s[-1] == "\"":
        # Minimal unescape matching _yaml_scalar
        body = s[1:-1]
        return body.replace("\\\"", "\"").replace("\\\\", "\\")
    # Best-effort number parsing
    try:
        if "." in s:
            return float(s)
        return int(s)
    except Exception:
        return s


def _read_yaml_simple(path: Path) -> dict:
    """Read YAML written by _write_yaml_atomic (1-level nested dicts only)."""
    out: dict = {}
    cur_key = None
    try:
        lines = path.read_text(errors="replace").splitlines()
    except Exception:
        return out
    for line in lines:
        if not line.strip():
            continue
        if line.startswith("  ") and cur_key is not None:
            # nested
            if ":" not in line:
                continue
            k, v = line.strip().split(":", 1)
            if cur_key not in out or not isinstance(out.get(cur_key), dict):
                out[cur_key] = {}
            out[cur_key][k.strip()] = _parse_yaml_scalar(v)
            continue
        # top-level
        cur_key = None
        if ":" not in line:
            continue
        k, v = line.split(":", 1)
        k = k.strip()
        if v.strip() == "":
            out[k] = {}
            cur_key = k
        else:
            out[k] = _parse_yaml_scalar(v)
    return out


def _signature_cache_key(sig: dict) -> str:
    """Stable cache key for deciding whether seqclust can be skipped."""
    key_payload = {
        # These sections capture semantically relevant inputs + resolved params.
        "paths": sig.get("paths") or {},
        "read_preparation": sig.get("read_preparation") or {},
        "read_numbers": sig.get("read_numbers") or {},
        "resolved_params": sig.get("resolved_params") or {},
    }
    blob = json.dumps(key_payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()
    return hashlib.md5(blob).hexdigest()


def _count_fasta_headers_best_effort(path: Path, max_bytes: int = 200 * 1024 * 1024):
    """Count FASTA records by '>' headers when file is not too large.

    Returns (records_or_none, bytes, counted_bool).
    """
    try:
        size = path.stat().st_size
    except FileNotFoundError:
        return None, None, False
    if size > max_bytes:
        return None, size, False
    n = 0
    with path.open("rb") as f:
        for line in f:
            if line.startswith(b">"):
                n += 1
    return n, size, True


def _append_seqclust_summary(log_path: Path, exit_code: int) -> None:
    try:
        data = log_path.read_text(errors="replace")
    except Exception:
        return

    def _count(substr: str) -> int:
        return data.count(substr)

    def _extract_error_context(lines: list[str], max_lines: int = 40) -> list[str]:
        """Extract a high-signal error excerpt near the end of the log."""
        if not lines:
            return []

        # Prefer the last full Python traceback block if present.
        for i in range(len(lines) - 1, -1, -1):
            if lines[i].startswith("Traceback (most recent call last):"):
                start = i
                break
        else:
            start = None

        if start is not None:
            excerpt = lines[start : start + max_lines]
            return excerpt

        # Otherwise, fall back to the last REvalError / Error / RuntimeError-ish line window.
        markers = ("REvalError", "RuntimeError", "ERROR ", "Error in ")
        for i in range(len(lines) - 1, -1, -1):
            if any(m in lines[i] for m in markers):
                start = max(0, i - 10)
                return lines[start : i + 1]
        return lines[-min(max_lines, len(lines)) :]

    # Heuristic markers (kept compatible with old wrap_seqclust).
    rserv_started = _count("Rserv started")
    db_created = _count("Building a new DB")
    sequences_loaded = _count("sequences in")  # noisy but useful
    blast_finished = _count("all to all blast finished")
    clustering_done = _count("louvain clustering")
    assembly_done = _count("assembling")
    ltr_detection = _count("detecting LTR")
    warnings = _count("Warning:")
    errors = _count("Error:")

    # Tail for quick context
    all_lines = data.splitlines()
    tail_lines = all_lines[-5:]
    err_ctx = _extract_error_context(all_lines)

    lines = []
    lines.append("")
    lines.append("SEQCLUST_EXECUTION_SUMMARY")
    lines.append(f"  exit_code: {exit_code}")
    lines.append("  progress_indicators:")
    if rserv_started:
        lines.append("    - rserv_started: true")
    if db_created:
        lines.append("    - db_created: true")
    if sequences_loaded:
        lines.append("    - sequences_loaded: true")
    if blast_finished:
        lines.append("    - blast_finished: true")
    if clustering_done:
        lines.append("    - clustering_done: true")
    if assembly_done:
        lines.append("    - assembly_done: true")
    if ltr_detection:
        lines.append("    - ltr_detection_attempted: true")
    if warnings:
        lines.append(f"  warnings: {warnings}")
    if errors:
        lines.append(f"  errors: {errors}")
    if err_ctx:
        lines.append("  error_context:")
        for el in err_ctx:
            # Keep YAML-ish indentation stable even for empty lines.
            lines.append(f"    {el}" if el else "    ")
    lines.append("  tail:")
    for tl in tail_lines:
        lines.append(f"    {tl}")
    lines.append("")

    try:
        with log_path.open("a") as f:
            f.write("\n".join(lines) + "\n")
    except Exception:
        return


def run_seqclust(
    *,
    project: str,
    sample: str,
    prepared_path: Path,
    tarean_dir: Path,
    tarean_done: Path,
    tarean_log: Path,
    seqclust_log: Path,
    signature_ripe: Optional[Path],
    threads: str,
    assembly_min: str,
    mincl_percent: str,
    min_lcov: str,
    merge_threshold: str,
    r_value: str,
    options: str,
    paired: bool,
    automatic_filtering: bool,
    tarean_mode: bool,
    keep_names: bool,
    cleanup: bool,
    domain_search: Optional[str],
    prefix_length: Optional[int],
    rserv_port: int,
    cleanup_after_prepare: bool,
    reads_per_assembly=None,
    sample_prefix=None,
    pythonhashseed=None,
    temp_dir=None,
    proportions_json=None,
) -> None:
    seqclust_log.parent.mkdir(parents=True, exist_ok=True)
    tarean_dir.mkdir(parents=True, exist_ok=True)
    ORIGINAL_DIR = Path.cwd()

    LOG = seqclust_log
    with LOG.open("a") as log:
        log.write(f"{time.ctime()}: Starting TAREAN analysis for {sample}\n")
        log.write(f"{time.ctime()}: Using {threads} threads\n")
        log.write(f"{time.ctime()}: Prepared input={prepared_path}\n")
        log.write(f"{time.ctime()}: TAREAN dir={tarean_dir}\n")
        log.write(f"{time.ctime()}: seqclust cleanup={cleanup_after_prepare}\n")

    def _rmdir_if_empty(p: Path) -> None:
        try:
            p.rmdir()
        except Exception:
            return

    def _cleanup_success_tmp() -> None:
        if not cleanup_after_prepare:
            return
        # Always clean run-local scratch in the tarean_dir.
        for scratch in (tarean_dir / "tmp", tarean_dir / "Rserv"):
            try:
                if scratch.exists():
                    shutil.rmtree(scratch, ignore_errors=True)
            except Exception:
                pass
        # Drop empty tmp bases (best-effort, non-recursive).
        # - iterative prep often uses `<tarean_dir_parent>/tmp`
        _rmdir_if_empty(tarean_dir.parent / "tmp")
        # - pipeline tmp bases requested: projects/<p>/samples/<s>/tmp and projects/<p>/comparative/<a>/tmp
        td = str(tarean_dir)
        if "/samples/" in td and "/comparative/" not in td:
            _rmdir_if_empty(Path(f"projects/{project}/samples/{sample}/tmp"))
        if "/comparative/" in td:
            _rmdir_if_empty(Path(f"projects/{project}/comparative/{sample}/tmp"))

    cluster_table_path = tarean_dir / "CLUSTER_TABLE.csv"
    # Lock short-circuit: if a ripe signature exists AND matches the current effective inputs/params,
    # skip expensive seqclust even if Snakemake invoked us due to provenance triggers.
    # (This preserves caching for code-only changes, but avoids incorrect reuse under config/param changes.)
    if tarean_done.exists() and signature_ripe is not None and signature_ripe.exists():
        if cluster_table_path.exists() and cluster_table_path.stat().st_size > 0:
            prior = _read_yaml_simple(signature_ripe)
            prior_key = str(prior.get("cache_key") or "")
            # Compute current key from current effective inputs/params.
            prepared_records, prepared_bytes, prepared_counted = _count_fasta_headers_best_effort(prepared_path)
            cur_sig = {
                "paths": {
                    "prepared_input": str(prepared_path),
                    "tarean_dir": str(tarean_dir),
                },
                "read_preparation": {
                    "reads_per_assembly_planned": reads_per_assembly,
                    "sample_prefix": sample_prefix,
                    "pythonhashseed": pythonhashseed,
                    "temp_dir": temp_dir,
                    "comparative_proportions_json": proportions_json,
                },
                "read_numbers": {
                    "prepared_fasta_bytes": prepared_bytes,
                    "prepared_fasta_records": prepared_records,
                    "prepared_fasta_records_counted": prepared_counted,
                },
                "resolved_params": {
                    "assembly_min": assembly_min,
                    "mincl_percent": mincl_percent,
                    "min_lcov": min_lcov,
                    "merge_threshold": merge_threshold,
                    "r_value_kb": r_value,
                    "options": options,
                    "paired": paired,
                    "automatic_filtering": automatic_filtering,
                    "tarean_mode": tarean_mode,
                    "keep_names": keep_names,
                    "cleanup": cleanup,
                    "domain_search": domain_search,
                    "threads": threads,
                },
            }
            if prefix_length is not None:
                cur_sig["resolved_params"]["prefix_length"] = prefix_length
            cur_key = _signature_cache_key(cur_sig)

            if prior_key and prior_key == cur_key:
                with LOG.open("a") as log:
                    log.write(
                        f"{time.ctime()}: LOCKED: ripe signature matches cache_key={cur_key}; "
                        f"skipping seqclust for {sample}\n"
                    )
                _cleanup_success_tmp()
                if "/samples/" in str(tarean_done) and "/comparative/" not in str(tarean_done):
                    try:
                        from project_manager import ProjectManager
                        ProjectManager().update_sample_tarean_status(project, sample, "complete")
                    except Exception:
                        pass
                return
            else:
                with LOG.open("a") as log:
                    log.write(
                        f"{time.ctime()}: UNLOCKED: ripe signature present but cache_key mismatch "
                        f"(prior={prior_key or 'missing'} current={cur_key}); will run seqclust\n"
                    )

    signature_draft = tarean_dir / "seqclust_signature.draft.yaml"
    if signature_ripe is None:
        signature_ripe = tarean_dir / "seqclust_signature.ripe.yaml"

    prepared_records, prepared_bytes, prepared_counted = _count_fasta_headers_best_effort(prepared_path)

    # Write draft signature at the start of the run.
    sig = {
        "status": "draft",
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "project": str(project),
        "sample": str(sample),
        "paths": {
            "prepared_input": str(prepared_path),
            "tarean_dir": str(tarean_dir),
        },
        "read_preparation": {
            "reads_per_assembly_planned": reads_per_assembly,
            "sample_prefix": sample_prefix,
            "pythonhashseed": pythonhashseed,
            "temp_dir": temp_dir,
            "comparative_proportions_json": proportions_json,
        },
        "read_numbers": {
            "prepared_fasta_bytes": prepared_bytes,
            "prepared_fasta_records": prepared_records,
            "prepared_fasta_records_counted": prepared_counted,
        },
        "resolved_params": {
            "assembly_min": assembly_min,
            "mincl_percent": mincl_percent,
            "min_lcov": min_lcov,
            "merge_threshold": merge_threshold,
            "r_value_kb": r_value,
            "options": options,
            "paired": paired,
            "automatic_filtering": automatic_filtering,
            "tarean_mode": tarean_mode,
            "keep_names": keep_names,
            "cleanup": cleanup,
            "domain_search": domain_search,
            "threads": threads,
        },
    }
    # comparative-only param (if present in params)
    if prefix_length is not None:
        sig["resolved_params"]["prefix_length"] = prefix_length
    sig["cache_key"] = _signature_cache_key(sig)
    _write_yaml_atomic(signature_draft, sig)

    # Clean up any existing Rserve processes for this sample (port-based to avoid pkill self-match)
    try:
        procs = subprocess.run(
            ["lsof", "-t", f"-i:{rserv_port}"],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
        )
        pids = [p.strip() for p in procs.stdout.splitlines() if p.strip()]
        for pid in pids:
            subprocess.run(["kill", pid], check=False)
    except FileNotFoundError:
        # lsof not available; skip best-effort cleanup
        pass
    # We keep behaviour simple here; if needed we could actually kill PIDs similar to bash script.
    time.sleep(2)

    # Set Rserv directory to project directory instead of /tmp
    os.environ["R_USER_LIBS_USER"] = str((Path.cwd() / tarean_dir / "Rserv").resolve())
    tarean_dir_rserv = tarean_dir / "Rserv"
    tarean_dir_rserv.mkdir(parents=True, exist_ok=True)
    # Ensure Rserve uses a run-local workdir (otherwise it defaults to /tmp/Rserv).
    os.environ["RSERVE_WORKDIR"] = str(tarean_dir_rserv.resolve())

    # Set unique Rserv port to avoid conflicts
    os.environ["RSERVE_PORT"] = str(rserv_port)

    # Runtime patching hook (no edits to external `repex_tarean/` needed).
    # This ensures `lib.r2py.create_connection()` starts Rserve with our workdir+port.
    shims_dir = (ORIGINAL_DIR / "shims").resolve()
    os.environ["PYTHONPATH"] = (
        f"{shims_dir}{os.pathsep}{os.environ.get('PYTHONPATH', '')}".rstrip(os.pathsep)
    )
    os.environ["REPORTR_PATCH_RSERVE"] = "1"

    # Per-run tmpdir to avoid /tmp collisions when running multiple seqclust jobs.
    # Put it under the run directory to keep job-local artifacts contained.
    per_run_tmp = (tarean_dir / "tmp").resolve()
    per_run_tmp.mkdir(parents=True, exist_ok=True)
    os.environ["TMPDIR"] = str(per_run_tmp)
    os.environ["TMP"] = str(per_run_tmp)
    os.environ["TEMP"] = str(per_run_tmp)

    # Stable python hashing (used across read-prep scripts and for determinism).
    if pythonhashseed is not None:
        os.environ["PYTHONHASHSEED"] = str(pythonhashseed)

    # Audit log for deterministic runtime isolation.
    with LOG.open("a") as log:
        log.write(f"{time.ctime()}: Rserv path (R_USER_LIBS_USER)={os.environ.get('R_USER_LIBS_USER')}\n")
        log.write(f"{time.ctime()}: RSERVE_PORT={os.environ.get('RSERVE_PORT')}\n")
        log.write(f"{time.ctime()}: TMPDIR/TMP/TEMP={per_run_tmp}\n")

    with LOG.open("a") as log:
        log.write(f"{time.ctime()}: Starting seqclust process (no external timeout wrapper)\n")
        log.write(
            f"{time.ctime()}: Command (logical): seqclust -v {tarean_dir} "
            f"{'-p ' if paired else ''}"
            f"{'-A ' if automatic_filtering else ''}"
            f"{'-t ' if tarean_mode else ''}"
            f"{'-k ' if keep_names else ''}"
            f"{'-C ' if cleanup else ''}"
            f"-a {assembly_min} -c {threads} -m {mincl_percent} -o {min_lcov} -M {merge_threshold} "
            f"-r {r_value} -opt {options} "
            f"{('-P ' + str(prefix_length) + ' ') if prefix_length is not None else ''}"
            f"{('-D ' + str(domain_search) + ' ') if domain_search else ''}"
            f"{prepared_path}\n"
        )

    # Run the in-repo binary built from `repex_tarean/` (compiled by installer).
    seqclust_bin = ORIGINAL_DIR / "repex_tarean" / "seqclust"
    if not seqclust_bin.exists():
        raise FileNotFoundError(
            f"seqclust not found at {seqclust_bin}. "
            "Run `python3 install_reportr.py` to build repex_tarean/seqclust."
        )
    seqclust_cwd = (ORIGINAL_DIR / "repex_tarean").resolve()
    # IMPORTANT: Only seqclust (and its Python/R stack, e.g. pyRserve) should run in
    # the `repeatexplorer` conda env. The Snakemake wrapper itself runs in `reportr`.
    conda_exe = os.environ.get("CONDA_EXE", "conda")
    cmd = [
        conda_exe,
        "run",
        "-n",
        "repeatexplorer",
        str(seqclust_bin),
        "-v",
        str(ORIGINAL_DIR / tarean_dir),
        "-a",
        assembly_min,
        "-c",
        threads,
        "-m",
        mincl_percent,
        "-o",
        min_lcov,
        "-M",
        merge_threshold,
        "-r",
        r_value,
        "-opt",
        options,
        str(ORIGINAL_DIR / prepared_path),
    ]
    if paired:
        cmd.insert(cmd.index("-a"), "-p")
    if automatic_filtering:
        cmd.insert(cmd.index("-a"), "-A")
    if tarean_mode:
        cmd.insert(cmd.index("-a"), "-t")
    if keep_names:
        cmd.insert(cmd.index("-a"), "-k")
    if cleanup:
        cmd.insert(cmd.index("-a"), "-C")
    if domain_search:
        # Append domain search choice (e.g. DIAMOND)
        cmd.insert(cmd.index(str(ORIGINAL_DIR / prepared_path)), "-D")
        cmd.insert(cmd.index(str(ORIGINAL_DIR / prepared_path)), str(domain_search))
    if prefix_length is not None:
        cmd.insert(cmd.index(str(ORIGINAL_DIR / prepared_path)), "-P")
        cmd.insert(cmd.index(str(ORIGINAL_DIR / prepared_path)), str(prefix_length))
    if cleanup_after_prepare:
        cmd.append("--cleanup")

    with LOG.open("a") as log:
        proc = subprocess.run(
            cmd,
            cwd=seqclust_cwd,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )

    # Clean up Rserve process for this sample (port-based best-effort)
    try:
        procs = subprocess.run(
            ["lsof", "-t", f"-i:{rserv_port}"],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
        )
        pids = [p.strip() for p in procs.stdout.splitlines() if p.strip()]
        for pid in pids:
            subprocess.run(["kill", pid], check=False)
    except FileNotFoundError:
        pass

    tarean_exit_code = proc.returncode
    _append_seqclust_summary(LOG, tarean_exit_code)

    core_dir = tarean_dir
    # Minimal set of outputs expected when seqclust runs with cleanup (-C).
    # Note: sequences.db/hitsort.db are typically removed by cleanup.
    required_core = [
        core_dir / "CLUSTER_TABLE.csv",
    ]
    core_ok = all(f.exists() and f.stat().st_size > 0 for f in required_core)
    # TAREAN consensus FASTAs may be absent for samples without satellites,
    # so treat them as a "nice to have" rather than a hard requirement.

    if tarean_exit_code == 0 or core_ok:
        with LOG.open("a") as log:
            if tarean_exit_code == 0:
                log.write(f"{time.ctime()}: TAREAN analysis completed successfully\n")
                log.write("TAREAN_SUCCESS: exit_code=0\n")
                tarean_log.write_text(f"{time.ctime()}: TAREAN analysis completed successfully\n")
                tarean_done.write_text(f"TAREAN_COMPLETE: {time.ctime()}\n")
            else:
                log.write(
                    f"{time.ctime()}: TAREAN warning: late-stage failure (exit code {tarean_exit_code}) "
                    "but core outputs are present; marking sample as usable\n"
                )
                log.write(
                    f"TAREAN_WARNING: exit_code={tarean_exit_code} (core outputs present, see above for details)\n"
                )
                tarean_log.write_text(
                    f"{time.ctime()}: TAREAN analysis completed with warnings for {sample}\n"
                )
                tarean_done.write_text(f"TAREAN_COMPLETE_WITH_WARNINGS: {time.ctime()}\n")
        with LOG.open("a") as log:
            log.write(f"{time.ctime()}: Analysis completed for {sample}\n")

        # Success-only cleanup of run-local scratch directories that seqclust may leave behind.
        # Covers both solo and comparative runs because both call this same runner.
        _cleanup_success_tmp()
        # Promote to ripe signature.
        sig["status"] = "ripe"
        sig["completed_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
        _write_yaml_atomic(signature_ripe, sig)
        # Update sample metadata.yaml with tarean_status (solo runs only).
        if "/samples/" in str(tarean_done) and "/comparative/" not in str(tarean_done):
            try:
                from project_manager import ProjectManager
                status = "complete" if tarean_exit_code == 0 else "complete_with_warnings"
                ProjectManager().update_sample_tarean_status(project, sample, status)
            except Exception:  # non-fatal
                pass
        return

    with LOG.open("a") as log:
        log.write(f"{time.ctime()}: TAREAN analysis failed with exit code {tarean_exit_code}\n")
        log.write(
            f"TAREAN_ERROR: exit_code={tarean_exit_code} (core outputs missing, see above for details)\n"
        )
    raise RuntimeError(f"seqclust failed with exit code {tarean_exit_code}")


if __name__ == "__main__":
    # Snakemake executes this script by injecting a global `snakemake` object.
    # For best-practice testing (and for non-Snakemake callers), we also provide
    # a small CLI that maps 1:1 to the same implementation.
    if "snakemake" in globals():
        project = snakemake.params.project  # type: ignore[name-defined]
        sample = snakemake.params.sample  # type: ignore[name-defined]
        assembly_min = str(snakemake.params.assembly_min)  # type: ignore[name-defined]
        mincl_percent = str(snakemake.params.mincl_percent)  # type: ignore[name-defined]
        min_lcov = str(snakemake.params.min_lcov)  # type: ignore[name-defined]
        merge_threshold = str(snakemake.params.merge_threshold)  # type: ignore[name-defined]
        r_value = str(snakemake.params.r_value)  # type: ignore[name-defined]
        options = str(snakemake.params.options)  # type: ignore[name-defined]
        paired = bool(snakemake.params.paired)  # type: ignore[name-defined]
        automatic_filtering = bool(snakemake.params.automatic_filtering)  # type: ignore[name-defined]
        tarean_mode = bool(snakemake.params.tarean_mode)  # type: ignore[name-defined]
        keep_names = bool(snakemake.params.keep_names)  # type: ignore[name-defined]
        cleanup = bool(snakemake.params.cleanup)  # type: ignore[name-defined]
        domain_search = getattr(snakemake.params, "domain_search", None)  # type: ignore[name-defined]
        try:
            prefix_length = int(getattr(snakemake.params, "prefix_length"))  # type: ignore[name-defined]
        except Exception:
            prefix_length = None
        # IMPORTANT: use Snakemake-scaled threads, not raw config values.
        threads = str(snakemake.threads)  # type: ignore[name-defined]
        tarean_dir = Path(str(snakemake.params.tarean_dir))  # type: ignore[name-defined]
        rserv_port = int(snakemake.params.rserv_port)  # type: ignore[name-defined]
        cleanup_after_prepare = bool(snakemake.params.cleanup_after_prepare)  # type: ignore[name-defined]

        # Support both solo and comparative rules:
        # - solo uses named input 'prepared'
        # - comparative uses a different key; fall back to first input item
        try:
            prepared_path = Path(str(snakemake.input.prepared))  # type: ignore[name-defined]
        except Exception:
            prepared_path = Path(str(snakemake.input[0]))  # type: ignore[name-defined]
        try:
            tarean_done = Path(str(snakemake.output.tarean_done))  # type: ignore[name-defined]
        except Exception:
            # comparative rule uses a different token name
            tarean_done = Path(str(snakemake.output.comparative_complete))  # type: ignore[name-defined]
        tarean_log = Path(str(snakemake.output.tarean_log))  # type: ignore[name-defined]
        signature_ripe = None
        try:
            signature_ripe = Path(str(snakemake.output.seqclust_signature_ripe))  # type: ignore[name-defined]
        except Exception:
            signature_ripe = None

        seqclust_log = Path(str(snakemake.log[0]))  # type: ignore[name-defined]

        run_seqclust(
            project=str(project),
            sample=str(sample),
            prepared_path=prepared_path,
            tarean_dir=tarean_dir,
            tarean_done=tarean_done,
            tarean_log=tarean_log,
            seqclust_log=seqclust_log,
            signature_ripe=signature_ripe,
            threads=threads,
            assembly_min=assembly_min,
            mincl_percent=mincl_percent,
            min_lcov=min_lcov,
            merge_threshold=merge_threshold,
            r_value=r_value,
            options=options,
            paired=paired,
            automatic_filtering=automatic_filtering,
            tarean_mode=tarean_mode,
            keep_names=keep_names,
            cleanup=cleanup,
            domain_search=domain_search,
            prefix_length=prefix_length,
            rserv_port=rserv_port,
            cleanup_after_prepare=cleanup_after_prepare,
            reads_per_assembly=getattr(snakemake.params, "reads_per_assembly", None),  # type: ignore[name-defined]
            sample_prefix=getattr(snakemake.params, "sample_prefix", None),  # type: ignore[name-defined]
            pythonhashseed=getattr(snakemake.params, "pythonhashseed", None),  # type: ignore[name-defined]
            temp_dir=getattr(snakemake.params, "temp_dir", None),  # type: ignore[name-defined]
            proportions_json=getattr(snakemake.params, "proportions_json", None),  # type: ignore[name-defined]
        )
    else:
        p = argparse.ArgumentParser(description="Run seqclust with RepOrtR runtime isolation")
        p.add_argument("--project", required=True)
        p.add_argument("--sample", required=True)
        p.add_argument("--prepared", required=True, type=Path)
        p.add_argument("--tarean-dir", required=True, type=Path)
        p.add_argument("--tarean-done", required=True, type=Path)
        p.add_argument("--tarean-log", required=True, type=Path)
        p.add_argument("--seqclust-log", required=True, type=Path)
        p.add_argument("--threads", required=True)
        p.add_argument("--assembly-min", required=True)
        p.add_argument("--mincl-percent", required=True)
        p.add_argument("--min-lcov", required=True)
        p.add_argument("--merge-threshold", default="0")
        p.add_argument("--r-value-kb", required=True)
        p.add_argument("--options", default="ILLUMINA_SENSITIVE_BLASTPLUS")
        p.add_argument("--paired", action="store_true")
        p.add_argument("--automatic-filtering", action="store_true")
        p.add_argument("--tarean-mode", action="store_true")
        p.add_argument("--keep-names", action="store_true")
        p.add_argument("--cleanup", action="store_true")
        p.add_argument("--domain-search", default=None)
        p.add_argument("--prefix-length", type=int, default=None)
        p.add_argument("--rserv-port", type=int, required=True)
        p.add_argument("--cleanup-after-prepare", action="store_true")
        args = p.parse_args()
        if int(args.rserv_port) == 0:
            import socket

            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            s.bind(("", 0))
            s.listen(1)
            args.rserv_port = s.getsockname()[1]
            s.close()

        run_seqclust(
            project=str(args.project),
            sample=str(args.sample),
            prepared_path=args.prepared,
            tarean_dir=args.tarean_dir,
            tarean_done=args.tarean_done,
            tarean_log=args.tarean_log,
            seqclust_log=args.seqclust_log,
            signature_ripe=None,
            threads=str(args.threads),
            assembly_min=str(args.assembly_min),
            mincl_percent=str(args.mincl_percent),
            min_lcov=str(args.min_lcov),
            merge_threshold=str(args.merge_threshold),
            r_value=str(args.r_value_kb),
            options=str(args.options),
            paired=bool(args.paired),
            automatic_filtering=bool(args.automatic_filtering),
            tarean_mode=bool(args.tarean_mode),
            keep_names=bool(args.keep_names),
            cleanup=bool(args.cleanup),
            domain_search=args.domain_search,
            prefix_length=args.prefix_length,
            rserv_port=int(args.rserv_port),
            cleanup_after_prepare=bool(args.cleanup_after_prepare),
        )

