#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


class UnionFind:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}
        self.rank: dict[str, int] = {}

    def add(self, x: str) -> None:
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0

    def find(self, x: str) -> str:
        self.add(x)
        p = self.parent[x]
        if p != x:
            self.parent[x] = self.find(p)
        return self.parent[x]

    def union(self, a: str, b: str) -> None:
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


def main() -> int:
    ap = argparse.ArgumentParser(description="Cluster repeats by connected components from LocalDB edges.")
    ap.add_argument("--edges-tsv", required=True, help="Filtered edges table with qseqid,sseqid")
    ap.add_argument("--nodes-fasta", required=True, help="Combined query FASTA (defines node universe)")
    ap.add_argument("--out-clusters-tsv", required=True)
    ap.add_argument("--out-cluster-sizes-tsv", required=True)
    args = ap.parse_args()

    edges_p = Path(args.edges_tsv)
    out_clusters = Path(args.out_clusters_tsv)
    out_sizes = Path(args.out_cluster_sizes_tsv)
    out_clusters.parent.mkdir(parents=True, exist_ok=True)
    out_sizes.parent.mkdir(parents=True, exist_ok=True)

    # Universe from FASTA headers (one line starting with >ID)
    nodes: list[str] = []
    with Path(args.nodes_fasta).open() as f:
        for ln in f:
            if ln.startswith(">"):
                nodes.append(ln[1:].strip().split()[0])

    uf = UnionFind()
    for n in nodes:
        uf.add(n)

    if edges_p.exists() and edges_p.stat().st_size:
        df = pd.read_csv(edges_p, sep="\t")
        if not df.empty:
            for _, r in df.iterrows():
                a = str(r.get("qseqid", "")).strip()
                b = str(r.get("sseqid", "")).strip()
                if not a or not b:
                    continue
                if a == b:
                    continue
                uf.union(a, b)

    # Map components to stable cluster IDs by sorting component keys.
    comp_to_members: dict[str, list[str]] = {}
    for n in sorted(nodes):
        root = uf.find(n)
        comp_to_members.setdefault(root, []).append(n)

    clusters_rows: list[dict] = []
    sizes_rows: list[dict] = []
    for i, (root, members) in enumerate(sorted(comp_to_members.items(), key=lambda kv: (-len(kv[1]), kv[0])), start=1):
        cid = f"CL{i:05d}"
        sizes_rows.append({"cluster_id": cid, "n_nodes": len(members)})
        for m in members:
            clusters_rows.append({"repeat": m, "cluster_id": cid})

    pd.DataFrame(clusters_rows).to_csv(out_clusters, sep="\t", index=False)
    pd.DataFrame(sizes_rows).to_csv(out_sizes, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

