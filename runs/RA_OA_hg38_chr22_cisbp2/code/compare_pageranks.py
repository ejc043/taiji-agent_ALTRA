#!/usr/bin/env python3
"""Compare this run's GeneRanks.tsv against a reference Taiji run.

Reports Spearman + Pearson correlation, top-N overlap, RMSE/MAE, and the
top worst-error TFs. Auto-detects motif-ID format (HOCOMOCO `*_HUMAN.*` vs
cisBP plain gene symbols) and normalizes to gene symbols on both sides
before joining.

Usage:
    python compare_pageranks.py \\
        --mine-dir   ../Output/Partial \\
        --ref-dir    /path/to/reference/taiji_results \\
        --pairs      RA_11=ra999 OA_02=oa1316
"""
from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

HOCOMOCO_RE = re.compile(r"^([^_]+)_HUMAN(?:\.|$)")


def to_symbol(motif_id: str) -> str:
    """HOCOMOCO `BCL6_HUMAN.H11MO.0.A` -> `BCL6`; cisBP `BCL6` -> `BCL6`."""
    m = HOCOMOCO_RE.match(motif_id)
    return m.group(1) if m else motif_id


def load_ranks(path: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    with path.open() as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            sym_i = next(i for i, h in enumerate(header) if h.lower() in ("tf", "gene", "geneid"))
        except StopIteration:
            sym_i = 0
        try:
            score_i = next(i for i, h in enumerate(header)
                           if "rank" in h.lower() or "score" in h.lower())
        except StopIteration:
            score_i = 1
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(sym_i, score_i):
                continue
            sym = to_symbol(parts[sym_i])
            try:
                out[sym] = float(parts[score_i])
            except ValueError:
                continue
    return out


def spearman(xs: list[float], ys: list[float]) -> float:
    def ranks(vs: list[float]) -> list[float]:
        order = sorted(range(len(vs)), key=lambda i: vs[i])
        r = [0.0] * len(vs)
        i = 0
        while i < len(vs):
            j = i
            while j + 1 < len(vs) and vs[order[j + 1]] == vs[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    rx, ry = ranks(xs), ranks(ys)
    return pearson(rx, ry)


def pearson(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    return num / (dx * dy) if dx and dy else float("nan")


def compare_one(mine: dict[str, float], ref: dict[str, float], label: str) -> None:
    shared = sorted(set(mine) & set(ref))
    only_mine = set(mine) - set(ref)
    only_ref = set(ref) - set(mine)
    print(f"\n=== {label} ===")
    print(f"  mine:   {len(mine):>5} TFs")
    print(f"  ref:    {len(ref):>5} TFs")
    print(f"  shared: {len(shared):>5}   only_mine={len(only_mine)}  only_ref={len(only_ref)}")
    if not shared:
        print("  (no shared TFs — check motif-ID format)")
        return

    xs = [mine[s] for s in shared]
    ys = [ref[s] for s in shared]
    print(f"  spearman = {spearman(xs, ys):.4f}")
    print(f"  pearson  = {pearson(xs, ys):.4f}")

    diffs = [(s, mine[s] - ref[s]) for s in shared]
    sq = [d * d for _, d in diffs]
    abs_d = [abs(d) for _, d in diffs]
    rmse = math.sqrt(sum(sq) / len(sq))
    mae = sum(abs_d) / len(abs_d)
    rng = max(ys) - min(ys) if ys else 1.0
    print(f"  RMSE     = {rmse:.3e}")
    print(f"  MAE      = {mae:.3e}")
    print(f"  NRMSE    = {(rmse / rng * 100):.2f}% (of ref score range {rng:.3e})")

    diffs.sort(key=lambda t: -abs(t[1]))
    print("  top-10 |err| TFs:")
    for s, d in diffs[:10]:
        print(f"    {s:<14} mine={mine[s]:.3e}  ref={ref[s]:.3e}  err={d:+.3e}")

    for k in (10, 25, 50, 100):
        top_mine = sorted(mine, key=lambda s: -mine[s])[:k]
        top_ref = sorted(ref, key=lambda s: -ref[s])[:k]
        overlap = len(set(top_mine) & set(top_ref))
        print(f"  top-{k:<3} overlap: {overlap}/{k}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--mine-dir", required=True, type=Path,
                   help="dir with <sample>_output/GeneRanks.tsv (typically ../Output/Partial)")
    p.add_argument("--ref-dir", required=True, type=Path,
                   help="dir with <ref_sample>/<ref_sample>_taiji_outputs/GeneRanks.tsv")
    p.add_argument("--pairs", nargs="+", required=True,
                   help="space-separated pairs MINE=REF (e.g. RA_11=ra999 OA_02=oa1316)")
    args = p.parse_args()

    for pair in args.pairs:
        mine_name, ref_name = pair.split("=", 1)
        mine_path = args.mine_dir / f"{mine_name}_output" / "GeneRanks.tsv"
        ref_path = args.ref_dir / ref_name / f"{ref_name}_taiji_outputs" / "GeneRanks.tsv"
        if not mine_path.is_file():
            print(f"MISSING mine: {mine_path}")
            continue
        if not ref_path.is_file():
            print(f"MISSING ref:  {ref_path}")
            continue
        compare_one(load_ranks(mine_path), load_ranks(ref_path), f"{mine_name} vs {ref_name}")


if __name__ == "__main__":
    main()
