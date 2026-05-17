from __future__ import annotations
from pathlib import Path
from typing import Dict, Any, Iterator, List, Tuple
import gzip

def _open_gff(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")

def parse_gff_attributes(attr: str) -> Dict[str, str]:
    d: Dict[str,str] = {}
    for part in attr.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            k,v = part.split("=", 1)
            d[k] = v
        elif " " in part:
            bits = part.split(" ", 1)
            d[bits[0]] = bits[1].strip('"')
    return d

def iter_gff_rows(gff: Path) -> Iterator[Dict[str, Any]]:
    with _open_gff(gff) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = cols
            yield {
                "seqid": seqid,
                "source": source,
                "type": ftype,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "attrs": parse_gff_attributes(attrs),
            }

def build_locus_table(gff: Path, prefer_type: str = "CDS") -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for r in iter_gff_rows(gff):
        if r["type"] != prefer_type:
            continue
        a = r["attrs"]
        locus = a.get("locus_tag") or a.get("locus") or a.get("ID") or a.get("Name")
        if not locus:
            continue
        gene = a.get("gene") or a.get("gene_name") or a.get("Name")
        product = a.get("product") or a.get("Note")
        rows.append({
            "locus_tag": locus,
            "gene": gene,
            "product": product,
            "seqid": r["seqid"],
            "start": r["start"],
            "end": r["end"],
            "strand": r["strand"],
            "length": abs(r["end"] - r["start"]) + 1,
        })
    if rows:
        return rows
    # fallback to gene features
    for r in iter_gff_rows(gff):
        if r["type"] != "gene":
            continue
        a = r["attrs"]
        locus = a.get("locus_tag") or a.get("ID") or a.get("Name")
        if not locus:
            continue
        gene = a.get("gene") or a.get("Name")
        rows.append({
            "locus_tag": locus,
            "gene": gene,
            "product": None,
            "seqid": r["seqid"],
            "start": r["start"],
            "end": r["end"],
            "strand": r["strand"],
            "length": abs(r["end"] - r["start"]) + 1,
        })
    return rows

def detect_featurecounts_keys(gff: Path) -> Tuple[str, str]:
    seen = 0
    has_locus = 0
    for r in iter_gff_rows(gff):
        if r["type"] == "CDS":
            seen += 1
            if "locus_tag" in r["attrs"]:
                has_locus += 1
            if seen >= 2000:
                break
    if seen == 0:
        return ("gene", "ID")
    if has_locus > 0:
        return ("CDS", "locus_tag")
    return ("CDS", "ID")
