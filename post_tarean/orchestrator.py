"""
Lightweight orchestrator for post-TAREAN steps.

This module wires together the existing ``RepeatAnalyzer`` with a sequence of
small, composable steps defined in ``post_tarean.steps``.
"""

from dataclasses import dataclass
from typing import Any, Dict, Iterable, Mapping

from post_tarean.steps.base import PostTareanContext, PostTareanStep
from post_tarean.steps.blast_step import BlastStep
from post_tarean.steps.summary_step import SummaryStep
from post_tarean.steps.quality_gating_step import QualityGatingStep


@dataclass
class PostTareanOrchestrator:
    """
    Orchestrate the execution of a sequence of post-TAREAN steps.

    The orchestrator is intentionally thin: it does not contain domain logic,
    only the glue that runs steps in order and passes a shared state dict.
    """

    analyzer: Any
    context: PostTareanContext
    steps: Mapping[str, PostTareanStep]

    def run(self, enabled_steps: Iterable[str]) -> Dict[str, Any]:
        """
        Run the selected steps in order, returning the final shared state.
        """
        state: Dict[str, Any] = {}
        for name in enabled_steps:
            step = self.steps.get(name)
            if step is None:
                continue
            state = step.run(self.analyzer, self.context, state)
        return state


def build_default_steps() -> Dict[str, PostTareanStep]:
    """
    Construct the default step registry.

    The returned dictionary does not encode ordering; that is decided by the
    caller via the ``enabled_steps`` sequence passed to ``PostTareanOrchestrator.run``.
    """
    steps: Dict[str, PostTareanStep] = {
        "blast": BlastStep(),
        "summary": SummaryStep(),
        "quality_gating": QualityGatingStep(),
    }
    return steps


