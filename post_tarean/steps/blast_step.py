"""
BLAST step for the post-TAREAN pipeline.

This step is intentionally thin and delegates the actual BLAST invocation to
the existing ``RepeatAnalyzer.run_blast_analysis`` method so that behaviour is
preserved while giving us a composable unit for orchestration.
"""

from typing import Any, Dict

from .base import PostTareanContext, PostTareanStep


class BlastStep(PostTareanStep):
    """Run BLAST for the given subject and store results in the shared state."""

    name = "blast"

    def run(self, analyzer: Any, context: PostTareanContext, state: Dict[str, Any]) -> Dict[str, Any]:
        # Incremental caching: if enabled globally and a previous BLAST result
        # file exists for this subject, reuse it (no BLAST recomputation).
        import yaml
        from pathlib import Path
        import pandas as pd

        import os

        global_cfg_path = os.environ.get("REPORTR_GLOBAL_CONFIG", "projects/global_config.yaml")
        global_cfg = yaml.safe_load(open(global_cfg_path, "r")) or {}
        incremental_enabled = bool(
            (global_cfg.get("global") or {}).get("incremental_processing", False)
        )

        blast_results_path = Path(context.output_dir) / f"{context.subject}_blast_results.csv"
        blast_results_path.parent.mkdir(parents=True, exist_ok=True)

        blast_results: Any = None
        if incremental_enabled and blast_results_path.exists():
            try:
                blast_results = pd.read_csv(blast_results_path)
            except Exception:
                blast_results = None

        if blast_results is None:
            # Reuse configuration stored on the analyzer; the method already applies
            # defaults for DBs, tasks, and threads.
            blast_results = analyzer.run_blast_analysis([context.subject])

            # Persist for downstream rules + incremental reuse.
            try:
                blast_results.to_csv(blast_results_path, index=False)
            except Exception:
                # Don't fail the pipeline if writing fails; reports may still
                # be generated from in-memory results.
                pass

        else:
            # Keep behaviour transparent in logs.
            try:
                analyzer.logger.info(
                    f"Incremental BLAST: reusing existing results from {blast_results_path}"
                )
            except Exception:
                pass

        state = dict(state)
        state["blast_results"] = blast_results
        return state

