#!/usr/bin/env python3
"""
Smoke tests for the refactored post-TAREAN pipeline.

These tests focus on the high-level orchestration behaviour of
`RepeatAnalyzer.run_full_pipeline`, verifying that:
- The new step-based orchestrator is invoked without errors.
- BLAST can be toggled via configuration / the `report_only` flag.

They deliberately stub out heavy BLAST and I/O work so they run quickly.
"""

import os
import tempfile
from types import SimpleNamespace

import pandas as pd

from post_tarean.config import PipelineConfig
from post_tarean import pipeline as post_pipeline


def _make_test_analyzer(tmpdir: str) -> post_pipeline.RepeatAnalyzer:
    """Create a lightweight RepeatAnalyzer instance without hitting real config."""
    # Construct a bare instance without running __init__
    analyzer = post_pipeline.RepeatAnalyzer.__new__(post_pipeline.RepeatAnalyzer)

    cfg = PipelineConfig()
    cfg.analysis.output_dir = tmpdir
    cfg.analysis.enabled_steps = ["blast", "summary", "quality_gating"]

    analyzer.config = cfg
    analyzer.config_manager = SimpleNamespace(project_id="test_project")

    return analyzer


def test_run_full_pipeline_with_blast(monkeypatch):
    """Pipeline should run BLAST and summary steps when both are enabled."""
    with tempfile.TemporaryDirectory() as tmpdir:
        analyzer = _make_test_analyzer(tmpdir)

        # Stub path resolution to avoid real filesystem layout.
        fake_tarean_path = os.path.join(tmpdir, "tarean")
        os.makedirs(fake_tarean_path, exist_ok=True)

        def fake_resolve_tarean_path(_cm, _cfg, _index):
            return fake_tarean_path, [fake_tarean_path]

        monkeypatch.setattr(post_pipeline, "resolve_tarean_path", fake_resolve_tarean_path)

        calls = {"blast": 0, "summary": 0}

        def fake_run_blast(subjects):
            calls["blast"] += 1
            # Minimal DataFrame with required columns
            return pd.DataFrame(
                {
                    "qseqid": [],
                    "evalue": [],
                    "length": [],
                    "qlength": [],
                    "task": [],
                    "sseqid": [],
                    "pident": [],
                }
            )

        def fake_generate_summary(index, path, blast_results=None):
            calls["summary"] += 1
            assert index == "test_sample"
            assert path == fake_tarean_path
            assert blast_results is not None
            # Minimal columns for quality gating + summary merge paths.
            return pd.DataFrame(
                {
                    "Cluster": ["test_sample_CL1"],
                    "size, %": [1.2],
                    "TAREAN_annotation": ["Putative satellites (high confidence)"],
                    "best_hit": ["some_hit"],
                    "best_evalue": [1e-50],
                    "best_pident": [98.0],
                    "coverage": [0.95],
                }
            )

        # Stub heavy report writers to no-ops
        monkeypatch.setattr(analyzer, "run_blast_analysis", fake_run_blast)
        monkeypatch.setattr(analyzer, "generate_summary_report", fake_generate_summary)
        monkeypatch.setattr(analyzer, "create_excel_report", lambda *a, **k: None)
        monkeypatch.setattr(analyzer, "create_word_report", lambda *a, **k: None)

        ok = post_pipeline.RepeatAnalyzer.run_full_pipeline(
            analyzer, "test_sample", output_dir=tmpdir, report_only=False
        )

        assert ok is True
        assert calls["blast"] == 1
        assert calls["summary"] == 1
        assert os.path.exists(os.path.join(tmpdir, "test_sample_confidence.json"))


def test_run_full_pipeline_report_only_skips_blast(monkeypatch):
    """When report_only is set, BLAST should be skipped even if enabled in config."""
    with tempfile.TemporaryDirectory() as tmpdir:
        analyzer = _make_test_analyzer(tmpdir)

        fake_tarean_path = os.path.join(tmpdir, "tarean")
        os.makedirs(fake_tarean_path, exist_ok=True)

        def fake_resolve_tarean_path(_cm, _cfg, _index):
            return fake_tarean_path, [fake_tarean_path]

        monkeypatch.setattr(post_pipeline, "resolve_tarean_path", fake_resolve_tarean_path)

        calls = {"blast": 0, "summary": 0}

        def fake_run_blast(subjects):
            calls["blast"] += 1
            return pd.DataFrame()

        def fake_generate_summary(index, path, blast_results=None):
            calls["summary"] += 1
            # In report-only mode we expect BLAST results to be None
            assert blast_results is None
            return pd.DataFrame(
                {
                    "Cluster": ["test_sample_CL1"],
                    "size, %": [0.2],
                    "TAREAN_annotation": ["Unclassified repeat (No evidence)"],
                    "best_hit": ["No significant hits"],
                    "best_evalue": [float("inf")],
                    "best_pident": [0.0],
                    "coverage": [0.0],
                }
            )

        monkeypatch.setattr(analyzer, "run_blast_analysis", fake_run_blast)
        monkeypatch.setattr(analyzer, "generate_summary_report", fake_generate_summary)
        monkeypatch.setattr(analyzer, "create_excel_report", lambda *a, **k: None)
        monkeypatch.setattr(analyzer, "create_word_report", lambda *a, **k: None)

        ok = post_pipeline.RepeatAnalyzer.run_full_pipeline(
            analyzer, "test_sample", output_dir=tmpdir, report_only=True
        )

        assert ok is True
        assert calls["blast"] == 0
        assert calls["summary"] == 1
        assert os.path.exists(os.path.join(tmpdir, "test_sample_confidence.json"))

