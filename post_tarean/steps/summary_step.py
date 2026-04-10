"""
Summary/report step for the post-TAREAN pipeline.

This step builds on top of the existing ``RepeatAnalyzer`` methods to generate
merged summary tables and, depending on configuration, write Excel/Word/CSV
reports for a single subject.
"""

import os
from typing import Any, Dict

from .base import PostTareanContext, PostTareanStep


class SummaryStep(PostTareanStep):
    """Generate summary tables and reports for a post-TAREAN subject."""

    name = "summary"

    def run(self, analyzer: Any, context: PostTareanContext, state: Dict[str, Any]) -> Dict[str, Any]:
        blast_results = state.get("blast_results")

        # Always ensure the output directory exists.
        os.makedirs(context.output_dir, exist_ok=True)

        # Use the existing high-level method to construct the merged DataFrame.
        summary_data = analyzer.generate_summary_report(
            context.subject,
            str(context.tarean_path),
            blast_results=blast_results,
        )

        # Derive the base name for all report artefacts from the analyzer config.
        analysis_cfg = analyzer.config.analysis
        base_name = os.path.join(
            str(context.output_dir),
            f"{context.subject}_{analysis_cfg.output_prefix}",
        )

        # Excel report
        if getattr(analysis_cfg, "create_excel", False):
            analyzer.create_excel_report(summary_data, f"{base_name}.xlsx")

        # Word report
        if getattr(analysis_cfg, "create_word", False):
            analyzer.create_word_report(context.subject, summary_data, f"{base_name}.docx")

        # CSV summary
        if getattr(analysis_cfg, "create_csv", False):
            summary_data.to_csv(f"{base_name}.csv", index=False)

        state = dict(state)
        state["summary_data"] = summary_data
        state["base_name"] = base_name
        return state

