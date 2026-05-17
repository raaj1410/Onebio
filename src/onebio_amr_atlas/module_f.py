from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple, Optional
import re
import json

import pandas as pd

from .utils import ensure_dir, log, PipelineError


def _pick_col(cols: List[str], keys: List[str]) -> Optional[str]:
    cols_l = [c.lower() for c in cols]
    for k in keys:
        for i, c in enumerate(cols_l):
            if k in c:
                return cols[i]
    return None


def _esc(s: object) -> str:
    return str(s).replace('"', "'")


def _sid(s: object, max_len: int = 60) -> str:
    x = re.sub(r"[^A-Za-z0-9_]+", "_", str(s))
    if re.match(r"^[0-9]", x):
        x = "n_" + x
    return x[:max_len]


def _read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _load_ast(ast_path: Path) -> pd.DataFrame:
    if not ast_path.exists():
        return pd.DataFrame(columns=["biosample_accession", "antibiotic", "phenotype"])
    ast = _read_tsv(ast_path)
    # normalise minimal expected cols
    if "biosample_accession" not in ast.columns:
        bs_col = _pick_col(list(ast.columns), ["biosample_accession", "biosample", "sample"])  # type: ignore
        if bs_col:
            ast = ast.rename(columns={bs_col: "biosample_accession"})
    if "antibiotic" not in ast.columns:
        ab_col = _pick_col(list(ast.columns), ["antibiotic", "drug", "antimicrobial"])  # type: ignore
        if ab_col:
            ast = ast.rename(columns={ab_col: "antibiotic"})
    if "phenotype" not in ast.columns:
        ph_col = _pick_col(list(ast.columns), ["phenotype", "sir", "interpret", "result"])  # type: ignore
        if ph_col:
            ast = ast.rename(columns={ph_col: "phenotype"})
    keep = [c for c in ["biosample_accession", "antibiotic", "phenotype", "mic", "unit"] if c in ast.columns]
    if keep:
        ast = ast[keep].copy()
    return ast


def _ast_items_for_bs(ast: pd.DataFrame, bs: str) -> List[Tuple[str, str]]:
    if ast.empty:
        return []
    sub = ast[ast["biosample_accession"].astype(str).str.strip() == bs].copy()
    if sub.empty:
        return []

    # collapse per antibiotic
    items: List[Tuple[str, str]] = []
    mic_col = "mic" if "mic" in sub.columns else None
    unit_col = "unit" if "unit" in sub.columns else None

    for ab, grp in sub.groupby(sub["antibiotic"].astype(str).str.strip(), sort=True):
        ab = (ab or "UNKNOWN_DRUG").strip() or "UNKNOWN_DRUG"
        labels: List[str] = []
        for _, rr in grp.iterrows():
            bits: List[str] = []
            ph = str(rr.get("phenotype", "") or "").strip()
            if ph:
                bits.append(f"S/I/R: {ph}")
            if mic_col:
                mic = str(rr.get(mic_col, "") or "").strip()
                if mic:
                    u = str(rr.get(unit_col, "") or "").strip() if unit_col else ""
                    bits.append(f"MIC: {mic}{(' ' + u) if u else ''}")
            if bits:
                labels.append(" | ".join(bits))
        # unique
        labels = list(dict.fromkeys([x for x in labels if x]))
        items.append((ab, "\\n".join(labels) if labels else "AST_RESULT"))
    return items


def _dot_for_bs(
    *,
    bs: str,
    asm: str,
    amr_df: pd.DataFrame,
    locus_df: pd.DataFrame,
    tpm_row: Optional[pd.Series],
    ast_items: List[Tuple[str, str]],
    top_n: Optional[int],
    title_suffix: str,
) -> str:
    """Generate a DOT graph for one BioSample.

    - If top_n is provided, we select top_n loci by TPM (when available), else first top_n.
    - If top_n is None, include all unique loci.
    """

    df = amr_df.copy()
    if "locus_tag" not in df.columns:
        raise PipelineError(f"amr_hits.tsv for {bs} missing locus_tag column")

    # enrich with locus meta
    if not locus_df.empty and "locus_tag" in locus_df.columns:
        cols = [c for c in ["locus_tag", "gene", "product"] if c in locus_df.columns]
        locus_small = locus_df[cols].drop_duplicates("locus_tag")
        df = df.merge(locus_small, on="locus_tag", how="left")

    df = df.drop_duplicates("locus_tag")

    # add tpm
    has_tpm = tpm_row is not None
    if has_tpm:
        # tpm_row is indexed by locus_tag
        def get_tpm(lt: str) -> float:
            try:
                v = tpm_row.get(lt)  # type: ignore
            except Exception:
                v = None
            try:
                return float(v) if v is not None else 0.0
            except Exception:
                return 0.0
        df["__tpm"] = df["locus_tag"].astype(str).map(get_tpm)
    else:
        df["__tpm"] = 0.0

    # choose rows
    if top_n is not None:
        if has_tpm:
            df = df.sort_values("__tpm", ascending=False)
        df = df.head(int(top_n))

    # DOT header
    dot: List[str] = []
    dot.append("digraph KG {")
    dot.append("  rankdir=TB; overlap=false; splines=polyline;")
    dot.append('  graph [bgcolor="white", ranksep=0.7, nodesep=0.35];')
    dot.append('  node  [shape=box, style="rounded,filled", fillcolor="#111827", fontcolor="white", fontname="Helvetica", fontsize=12];')
    dot.append('  edge  [color="#9ca3af", fontcolor="#6b7280", fontname="Helvetica", fontsize=10, arrowsize=0.7];')

    dot.append(f'  BS  [label="BioSample\\n{_esc(bs)}"];')
    asm_label = f"Assembly\\n{_esc(asm)}" if asm else "Assembly"
    dot.append(f'  ASM [label="{asm_label}"];')
    dot.append('  BS -> ASM [label="HAS_ASSEMBLY"];')

    # AST nodes
    ab_nodes: List[str] = []
    for i, (ab, lbl) in enumerate(ast_items):
        aid = f"AB{i}"
        ab_nodes.append(aid)
        dot.append(f'  {aid} [label="Antibiotic\\n{_esc(ab)}", fillcolor="#0f766e"];')
        dot.append(f'  BS -> {aid} [label="AST_RESULT\\n{_esc(lbl)}"];')

    # Gene + determinant nodes
    gene_nodes: List[str] = []
    det_nodes: List[str] = []

    # ensure expected cols exist
    c_det = "amr_gene" if "amr_gene" in df.columns else ("element_symbol" if "element_symbol" in df.columns else None)
    c_cls = "drug_class" if "drug_class" in df.columns else ("class" if "class" in df.columns else None)

    for i, r in df.reset_index(drop=True).iterrows():
        lt = str(r.get("locus_tag", "") or "").strip()
        if not lt or lt.lower() == "nan":
            continue
        g = f"G{i}"
        d = f"D{i}"
        gene_nodes.append(g)
        det_nodes.append(d)

        gene = str(r.get("gene", "") or "").strip()
        prod = str(r.get("product", "") or "").strip()
        det = str(r.get(c_det, "") or "").strip() if c_det else ""
        det = det if det else "UNKNOWN_DET"
        dcls = str(r.get(c_cls, "") or "").strip() if c_cls else ""
        dcls = dcls if dcls else "UNKNOWN"

        g_label = f"Gene/Locus\\nlocus_tag: { _esc(lt) }"
        if gene and gene.lower() != "nan":
            g_label += f"\\ngene: {_esc(gene)}"
        if prod and prod.lower() != "nan":
            short = prod[:55] + ("…" if len(prod) > 55 else "")
            g_label += f"\\n{_esc(short)}"

        d_label = f"AMR Determinant\\n{_esc(det)}\\nclass: {_esc(dcls)}"

        dot.append(f'  {g} [label="{g_label}"];')
        dot.append(f'  {d} [label="{d_label}"];')

        dot.append(f'  ASM -> {g} [label="HAS_GENE"];')

        tpm = float(r.get("__tpm", 0.0) or 0.0)
        if has_tpm:
            dot.append(f'  BS -> {g} [label="EXPRESSES\\nTPM: {tpm:.2f}"];')
        else:
            dot.append(f'  BS -> {g} [label="EXPRESSES\\nTPM: NA"];')

        dot.append(f'  {g} -> {d} [label="ANNOTATED_AS"];')

    # ranks to keep slide-like style
    top_rank = ["ASM"] + ab_nodes
    if top_rank:
        dot.append("  { rank=same; " + " ".join(top_rank) + " }")
    if gene_nodes:
        dot.append("  { rank=same; " + " ".join(gene_nodes) + " }")
    if det_nodes:
        dot.append("  { rank=same; " + " ".join(det_nodes) + " }")

    # footer
    dot.append("}")
    return "\n".join(dot)


def build_kg_dot_all(
    manifest_tsv: Path,
    genomes_dir: Path,
    features_dir: Path,
    out_dir: Path,
    *,
    top_n: int = 10,
    make_full: bool = True,
) -> Path:
    """Module F: Generate BioSample knowledge-graph DOT files.

    Outputs:
      out_dir/top{top_n}/kg_<BS>_top{top_n}_with_ast.dot
      out_dir/full/kg_<BS>_full_with_ast.dot (optional)
      out_dir/kg_index.tsv
    """

    ensure_dir(out_dir)

    man = _read_tsv(manifest_tsv)
    if "biosample_accession" not in man.columns:
        raise PipelineError("manifest.tsv missing biosample_accession")

    ast_path = Path(manifest_tsv).parent / "ast_long.tsv"
    ast = _load_ast(ast_path)

    X_expr_path = Path(features_dir) / "X_expr_tpm.tsv"
    if X_expr_path.exists():
        X_expr = pd.read_csv(X_expr_path, sep="\t", index_col=0)
    else:
        # allow KG generation without expression, but topN selection will be arbitrary
        log("WARN: Missing X_expr_tpm.tsv. Expression edges will show TPM: NA and topN selection is arbitrary.")
        X_expr = pd.DataFrame()

    out_top = ensure_dir(Path(out_dir) / f"top{int(top_n)}")
    out_full = ensure_dir(Path(out_dir) / "full") if make_full else None

    rows: List[Dict[str, str]] = []

    for _, r in man.iterrows():
        bs = str(r.get("biosample_accession", "") or "").strip()
        if not bs or bs.lower() == "nan":
            continue

        asm = str(r.get("assembly_accession_best", "") or "").strip()
        bs_dir = Path(genomes_dir) / bs
        amr_path = bs_dir / "amr_hits.tsv"
        locus_path = bs_dir / "locus_map.tsv"

        if not amr_path.exists():
            log(f"WARN: {bs} missing amr_hits.tsv, skipping KG")
            continue
        if not locus_path.exists():
            log(f"WARN: {bs} missing locus_map.tsv, KG will have limited gene labels")

        try:
            amr_df = _read_tsv(amr_path)
            locus_df = _read_tsv(locus_path) if locus_path.exists() else pd.DataFrame()
        except Exception as e:
            log(f"WARN: {bs} failed reading AMR/locus TSVs: {e}")
            continue

        tpm_row = None
        if not X_expr.empty and bs in X_expr.index.astype(str):
            try:
                tpm_row = X_expr.loc[bs]
            except Exception:
                tpm_row = None

        ast_items = _ast_items_for_bs(ast, bs)

        # topN
        top_dot = _dot_for_bs(
            bs=bs,
            asm=asm,
            amr_df=amr_df,
            locus_df=locus_df,
            tpm_row=tpm_row,
            ast_items=ast_items,
            top_n=int(top_n),
            title_suffix=f"top{int(top_n)}",
        )
        top_path = out_top / f"kg_{bs}_top{int(top_n)}_with_ast.dot"
        top_path.write_text(top_dot, encoding="utf-8")

        full_path = ""
        if make_full and out_full is not None:
            full_dot = _dot_for_bs(
                bs=bs,
                asm=asm,
                amr_df=amr_df,
                locus_df=locus_df,
                tpm_row=tpm_row,
                ast_items=ast_items,
                top_n=None,
                title_suffix="full",
            )
            fp = out_full / f"kg_{bs}_full_with_ast.dot"
            fp.write_text(full_dot, encoding="utf-8")
            full_path = str(fp)

        rows.append({
            "biosample_accession": bs,
            "assembly_accession_best": asm,
            "amr_hits_rows": str(len(amr_df)),
            "ast_rows": str(len(ast_items)),
            "top_dot": str(top_path),
            "full_dot": full_path,
            "has_expr": "1" if tpm_row is not None else "0",
        })

    idx = Path(out_dir) / "kg_index.tsv"
    pd.DataFrame(rows).to_csv(idx, sep="\t", index=False)
    log(f"KG index: {idx}")
    log(f"KG top{int(top_n)} dir: {out_top}")
    if make_full and out_full is not None:
        log(f"KG full dir: {out_full}")
    return idx
