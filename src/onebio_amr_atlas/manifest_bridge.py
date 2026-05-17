from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Iterable

import pandas as pd

from .utils import ensure_dir

RUN_RE = re.compile(r"(?:SRR|ERR|DRR)\d+", re.IGNORECASE)
BS_RE = re.compile(r"(?:SAMN|SAMEA|SAMD)\d+", re.IGNORECASE)
ASM_RE = re.compile(r"(?:GCA|GCF)_\d+\.\d+", re.IGNORECASE)

MANIFEST_COLUMNS = [
    "biosample_accession",
    "bioproject_accession",
    "taxon_id",
    "strain",
    "genome_id",
    "assembly_accession_best",
    "evidence_tier",
    "rnaseq_run_count",
    "run_list",
    "ena_fastq_ftp",
    "ena_fastq_md5",
    "access_method",
]


def _read_table(path: Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(path, dtype=str).fillna("")
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False).fillna("")


def _token(text: object, rx: re.Pattern[str]) -> str:
    m = rx.search(str(text or ""))
    return m.group(0).upper() if m else ""


def _runs(text: object) -> list[str]:
    return sorted(set(m.group(0).upper() for m in RUN_RE.finditer(str(text or ""))))


def _truthy_series(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y", "t"})


def _split_antibiotics(value: object) -> Iterable[str]:
    for x in re.split(r"[;,]+", str(value or "")):
        x = " ".join(x.strip().split())
        if x and x.lower() not in {"nan", "none", "null"}:
            yield x


def checker_csv_to_manifest(
    checker_table: Path,
    out_dir: Path,
    only_final_hits: bool = True,
    only_rnaseq_available: bool = True,
    max_rows: int | None = None,
) -> tuple[Path, Path, Path, Path]:
    """Convert RNA-seq AMR checker output into strict Module A/B/C inputs.

    The checker is discovery. The manifest is execution. Keeping that boundary
    explicit prevents the usual pipeline soup where nobody knows which CSV is
    supposed to feed which script. Civilization inches forward.
    """
    out_dir = ensure_dir(Path(out_dir))
    df = _read_table(Path(checker_table))
    original_rows = len(df)

    for c in [
        "biosample_accession", "assembly_accession_best", "assembly_accession",
        "rnaseq_run_ids", "run_list", "final_conservative_hit", "rnaseq_available",
        "antibiotics_resistant", "antibiotics_susceptible", "antibiotics_intermediate",
    ]:
        if c not in df.columns:
            df[c] = ""

    df["biosample_accession"] = df["biosample_accession"].map(lambda x: _token(x, BS_RE))

    # Assembly fallback: checker output may use assembly_accession or assembly_accession_best.
    df["assembly_accession_best"] = df["assembly_accession_best"].where(
        df["assembly_accession_best"].astype(str).str.strip().ne(""),
        df["assembly_accession"],
    )
    df["assembly_accession_best"] = df["assembly_accession_best"].map(lambda x: _token(x, ASM_RE))

    if only_final_hits and "final_conservative_hit" in df.columns:
        df = df[_truthy_series(df["final_conservative_hit"])].copy()
    if only_rnaseq_available and "rnaseq_available" in df.columns:
        df = df[_truthy_series(df["rnaseq_available"])].copy()

    df["run_list"] = df.apply(
        lambda r: ";".join(_runs(r.get("rnaseq_run_ids", "")) or _runs(r.get("run_list", ""))),
        axis=1,
    )
    df["rnaseq_run_count"] = df["run_list"].map(lambda x: len(_runs(x)))

    df = df[(df["biosample_accession"] != "") & (df["assembly_accession_best"] != "") & (df["run_list"] != "")].copy()
    df = df.drop_duplicates(subset=["biosample_accession", "assembly_accession_best"], keep="first")

    if max_rows is not None and int(max_rows) > 0:
        df = df.head(int(max_rows)).copy()

    for c in MANIFEST_COLUMNS:
        if c not in df.columns:
            df[c] = ""

    # Preserve useful defaults.
    df["evidence_tier"] = df["evidence_tier"].where(df["evidence_tier"].astype(str).str.strip().ne(""), "Strong same-BioSample RNA-seq + AST + genome")
    df["access_method"] = df["access_method"].where(df["access_method"].astype(str).str.strip().ne(""), "FASTQ resolver: local -> ENA manifest -> ENA run lookup -> SRA Toolkit")

    manifest = df[MANIFEST_COLUMNS].copy()
    manifest_path = out_dir / "manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)

    ast_rows = []
    for _, r in df.iterrows():
        bs = str(r.get("biosample_accession", ""))
        for col, pheno in [
            ("antibiotics_resistant", "R"),
            ("antibiotics_susceptible", "S"),
            ("antibiotics_intermediate", "I"),
        ]:
            for ab in _split_antibiotics(r.get(col, "")):
                ast_rows.append({"biosample_accession": bs, "antibiotic": ab, "phenotype": pheno})
    ast = pd.DataFrame(ast_rows).drop_duplicates() if ast_rows else pd.DataFrame(columns=["biosample_accession", "antibiotic", "phenotype"])
    ast_path = out_dir / "ast_long.tsv"
    ast.to_csv(ast_path, sep="\t", index=False)

    selected_path = out_dir / "checker_selected_rows.csv"
    df.to_csv(selected_path, index=False)

    qc = {
        "input_table": str(checker_table),
        "input_rows": int(original_rows),
        "selected_rows": int(len(df)),
        "manifest_rows": int(len(manifest)),
        "ast_rows": int(len(ast)),
        "only_final_hits": bool(only_final_hits),
        "only_rnaseq_available": bool(only_rnaseq_available),
        "missing_after_filter": {
            "biosample_accession": int((df["biosample_accession"] == "").sum()) if not df.empty else 0,
            "assembly_accession_best": int((df["assembly_accession_best"] == "").sum()) if not df.empty else 0,
            "run_list": int((df["run_list"] == "").sum()) if not df.empty else 0,
        },
    }
    qc_path = out_dir / "manifest_qc.json"
    qc_path.write_text(json.dumps(qc, indent=2), encoding="utf-8")

    return manifest_path, ast_path, selected_path, qc_path
