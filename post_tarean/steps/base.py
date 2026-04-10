"""
Base abstractions for post-TAREAN steps.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Protocol


@dataclass
class PostTareanContext:
    """
    Minimal, serializable context describing a single post-TAREAN run.

    Attributes
    ----------
    project_id:
        Logical project identifier (matches keys in ``projects/global_config.yaml``).
    subject:
        Sample name or comparative analysis identifier.
    mode:
        Either ``\"sample\"`` or ``\"comparative\"``; controls how paths are
        interpreted and which rules depend on the outputs.
    tarean_path:
        Filesystem path to the TAREAN/RepeatExplorer output directory for this
        subject.
    output_dir:
        Directory where post-TAREAN artifacts (BLAST results, summaries,
        reports) should be written.
    """

    project_id: str
    subject: str
    mode: str
    tarean_path: Path
    output_dir: Path


class PostTareanStep(Protocol):
    """
    Simple protocol for a post-TAREAN processing step.

    Each step receives the shared ``RepeatAnalyzer`` instance, a context object,
    and a mutable state dictionary that can be used to pass intermediate
    results between steps (e.g., BLAST DataFrames).
    """

    name: str

    def run(self, analyzer: Any, context: PostTareanContext, state: Dict[str, Any]) -> Dict[str, Any]:
        ...

