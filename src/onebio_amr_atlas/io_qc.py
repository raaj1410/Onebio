from __future__ import annotations
import gzip, statistics
from pathlib import Path
from typing import Dict, Any

def _open_maybe_gz(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")

def fastq_sample_stats(fq: Path, max_reads: int = 50000) -> Dict[str, Any]:
    lens = []
    qs = []
    n_rates = []
    reads = 0
    with _open_maybe_gz(fq) as f:
        while True:
            h = f.readline()
            if not h:
                break
            seq = f.readline().strip()
            plus = f.readline()
            qual = f.readline().strip()
            if not qual:
                break
            reads += 1
            if reads > max_reads:
                break
            L = len(seq)
            lens.append(L)
            if qual:
                q = sum((ord(c)-33) for c in qual) / max(1, len(qual))
                qs.append(q)
            if L:
                n_rates.append(sum(1 for c in seq.upper() if c == "N") / L)
    if not lens:
        return {"reads_sampled": 0}
    def n50(arr):
        total = sum(arr); acc = 0
        for x in sorted(arr, reverse=True):
            acc += x
            if acc >= total/2:
                return x
        return arr[-1]
    return {
        "reads_sampled": len(lens),
        "len_mean": sum(lens)/len(lens),
        "len_median": statistics.median(lens),
        "len_min": min(lens),
        "len_max": max(lens),
        "len_n50": n50(lens),
        "q_mean": sum(qs)/len(qs) if qs else None,
        "q_median": statistics.median(qs) if qs else None,
        "n_rate_mean": sum(n_rates)/len(n_rates) if n_rates else None,
        "n_rate_median": statistics.median(n_rates) if n_rates else None,
    }
