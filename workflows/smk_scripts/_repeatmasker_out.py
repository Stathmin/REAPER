from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List


# RepeatMasker repeat field can contain our tagged consensus IDs.
# Historical tag shapes:
# - Old: SMPL=KA1|ITER=1|RANK=2|ORIG=...
# - Current: SMPL=KA1|ORG=...|GENOMES=...|ITER=1|RANK=2|ORIG=...
TAG_RE = re.compile(
    r"SMPL=(?P<smpl>[^|]+)\|"
    r"(?:(?:ORG=[^|]*)\|(?:GENOMES=[^|]*)\|)?"
    r"ITER=(?P<iter>[0-9]+)\|RANK=(?P<rank>[0-9]+)\|ORIG=(?P<orig>.+)"
)


@dataclass(frozen=True)
class RMHit:
    score: float
    pdiv: float
    query: str
    qbegin: int
    qend: int
    strand: str
    repeat: str


def iter_hits(lines: Iterable[str]) -> Iterator[RMHit]:
    for ln in lines:
        if not ln.strip():
            continue
        if ln.lstrip().startswith("SW") or ln.lstrip().startswith("score"):
            continue
        parts = ln.replace("(", "").replace(")", "").split()
        if len(parts) < 10:
            continue
        try:
            score = float(parts[0])
            pdiv = float(parts[1])
            query = parts[4]
            qbegin = int(parts[5])
            qend = int(parts[6])
            strand = parts[8]
            repeat = parts[9]
        except Exception:
            continue
        yield RMHit(score=score, pdiv=pdiv, query=query, qbegin=qbegin, qend=qend, strand=strand, repeat=repeat)


def parse_repeat_tags(rep: str) -> Dict[str, str]:
    m = TAG_RE.fullmatch(rep)
    if not m:
        return {"smpl": "", "iter": "", "rank": "", "orig": rep}
    d = m.groupdict()
    return {"smpl": d["smpl"], "iter": d["iter"], "rank": d["rank"], "orig": d["orig"]}


def safe_quantile(values: List[float], q: float) -> float:
    if not values:
        return float("nan")
    s = sorted(values)
    if len(s) == 1:
        return s[0]
    pos = (len(s) - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return s[lo]
    frac = pos - lo
    return s[lo] * (1 - frac) + s[hi] * frac


def infer_prefix(query_id: str, allowed: Iterable[str]) -> str:
    # Comparative headers are like KA1read123_f; allow only configured prefixes.
    for p in allowed:
        if query_id.startswith(p):
            return p
    return ""

