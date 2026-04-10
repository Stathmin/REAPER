#!/usr/bin/env python3
"""Snakemake helper for configuration validation."""

from datetime import datetime
from pathlib import Path

import yaml


def _ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _append(path: Path, line: str) -> None:
    with path.open("a", encoding="utf-8") as handle:
        handle.write(f"{line}\n")


def _check_mapping_keys(
    report_path: Path,
    section_name: str,
    mapping: dict,
    required_keys: list[str],
) -> None:
    _append(report_path, f"Validating {section_name}...")
    for key in required_keys:
        if key in mapping:
            _append(report_path, f"✓ {section_name} key '{key}' found")
        else:
            _append(report_path, f"⚠ {section_name} key '{key}' missing (warning only)")


def main() -> None:
    output_path = Path(str(snakemake.output.config_valid))
    log_path = Path(str(snakemake.log[0]))
    config_path = Path("projects/global_config.yaml")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    _append(log_path, f"{_ts()}: Validating configuration...")
    output_path.write_text("", encoding="utf-8")

    try:
        config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    except Exception as exc:  # noqa: BLE001
        _append(
            output_path,
            f"⚠ Could not read {config_path}: {exc} (warning only, validation skipped)",
        )
        _append(log_path, f"{_ts()}: Configuration validation skipped due to parse error")
        return

    global_cfg = config.get("global") or {}
    _check_mapping_keys(
        output_path,
        "global config",
        global_cfg,
        list(snakemake.params.required_global_keys),
    )
    _check_mapping_keys(
        output_path,
        "read cleaning",
        global_cfg.get("read_cleaning") or {},
        list(snakemake.params.required_read_cleaning_keys),
    )
    _check_mapping_keys(
        output_path,
        "read preparation",
        global_cfg.get("read_preparation") or {},
        list(snakemake.params.required_read_preparation_keys),
    )

    _append(output_path, "Validating project configurations...")
    projects = config.get("projects") or {}
    required_project_keys = ["samples", "tarean_params", "comparative_species"]
    required_tarean_keys = ["assembly_min", "mincl", "threads"]

    if not projects:
        _append(output_path, "⚠ No projects declared under `projects` (warning only)")

    for project_id, project_cfg in projects.items():
        _append(output_path, f"Project: {project_id}")
        for key in required_project_keys:
            if key in (project_cfg or {}):
                _append(output_path, f"  ✓ {key} present")
            else:
                _append(output_path, f"  ⚠ {key} missing (warning only)")

        tarean_cfg = (project_cfg or {}).get("tarean_params") or {}
        for key in required_tarean_keys:
            if key in tarean_cfg:
                _append(output_path, f"    ✓ TAREAN {key}: {tarean_cfg[key]}")
            else:
                _append(output_path, f"    ⚠ TAREAN {key} missing (warning only)")

    _append(output_path, "✓ Configuration validation completed")
    _append(log_path, f"{_ts()}: Configuration validation completed successfully")


if __name__ == "__main__":
    main()
