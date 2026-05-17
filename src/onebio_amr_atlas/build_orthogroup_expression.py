from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

import pandas as pd


PANAROO_META_COLS = {
    "Gene",
    "Non-unique Gene name",
    "Annotation",
    "No. isolates",
    "No. sequences",
    "Avg sequences per isolate",
    "Genome Fragment",
    "Order within Fragment",
    "Accessory Fragment",
    "Accessory Order with Fragment",
    "QC",
    "Min group size nuc",
    "Max group size nuc",
    "Avg group size nuc",
}


# -----------------------------
# Generic helpers
# -----------------------------

def clean_token(x: str) -> str:
    x = str(x or "").strip().strip('"').strip("'")
    if not x or x.lower() == "nan":
        return ""
    x = re.sub(r"(_len|_stop|_pseudo)$", "", x)
    return x


def split_gene_cell(cell: str) -> List[str]:
    raw = str(cell or "").strip()
    if not raw or raw.lower() == "nan":
        return []
    out = []
    for part in raw.split(";"):
        part = clean_token(part)
        if part:
            out.append(part)
    return out


def parse_gff_attributes(attr: str) -> Dict[str, str]:
    d: Dict[str, str] = {}
    for part in str(attr or "").strip().split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
        elif " " in part:
            k, v = part.split(" ", 1)
            d[k] = v.strip('"')
    return d


def iter_gff_rows(gff_path: Path):
    with gff_path.open("rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            yield {
                "seqid": cols[0],
                "source": cols[1],
                "type": cols[2],
                "start": cols[3],
                "end": cols[4],
                "strand": cols[6],
                "attrs": parse_gff_attributes(cols[8]),
            }


def canonical_feature_id(attrs: Dict[str, str]) -> str:
    """
    Must match the same precedence used in your genome-processing code:
    locus_tag -> locus -> ID -> Name
    """
    for key in ["locus_tag", "locus", "ID", "Name"]:
        val = clean_token(attrs.get(key, ""))
        if val:
            return val
    return ""


def alias_variants(raw: str) -> Set[str]:
    """
    Generate robust aliases for matching Panaroo IDs back to the canonical
    feature key used in locus_map / TPM tables.
    """
    raw = clean_token(raw)
    if not raw:
        return set()

    vals: Set[str] = {raw}

    # Strip common feature prefixes used in GFF/Panaroo IDs
    queue = [raw]
    prefixes = [
        "cds-",
        "gene-",
        "rna-",
        "trna-",
        "rrna-",
        "ncrna-",
        "transcript-",
        "id-",
    ]
    while queue:
        x = queue.pop()
        for pref in prefixes:
            if x.lower().startswith(pref):
                y = x[len(pref):]
                if y and y not in vals:
                    vals.add(y)
                    queue.append(y)

    # Strip trailing version suffixes like .1 .2
    new_vals = set()
    for x in vals:
        new_vals.add(x)
        x2 = re.sub(r"\.\d+$", "", x)
        if x2:
            new_vals.add(x2)
    vals = new_vals

    # If Dbxref-like token appears, keep both full and right-hand side
    new_vals = set(vals)
    for x in vals:
        if ":" in x:
            rhs = x.split(":", 1)[1]
            rhs = clean_token(rhs)
            if rhs:
                new_vals.add(rhs)
                new_vals.add(re.sub(r"\.\d+$", "", rhs))
    vals = new_vals

    return {v for v in vals if v and v.lower() != "nan"}


# -----------------------------
# Sample discovery
# -----------------------------

def find_biosamples(genomes_dir: Path) -> List[str]:
    return sorted([p.name for p in genomes_dir.iterdir() if p.is_dir()])


def load_locus_map(genomes_dir: Path, bs: str) -> pd.DataFrame:
    p = genomes_dir / bs / "locus_map.tsv"
    if not p.exists():
        return pd.DataFrame(columns=["locus_tag"])
    df = pd.read_csv(p, sep="\t", dtype=str).fillna("")
    if "locus_tag" not in df.columns:
        return pd.DataFrame(columns=["locus_tag"])
    return df


# -----------------------------
# Alias map construction
# -----------------------------

def build_sample_alias_map(genomes_dir: Path, bs: str) -> Tuple[Dict[str, str], Dict[str, int]]:
    """
    Build a sample-specific alias -> canonical_feature_id map.

    Canonical feature ID is the exact key used in TPM/locus_map.
    Aliases are gathered from raw GFF attributes so Panaroo IDs can map back.
    """
    gff = genomes_dir / bs / "annotation.gff"
    locus_df = load_locus_map(genomes_dir, bs)

    canonical_loci: Set[str] = set(locus_df["locus_tag"].astype(str).map(clean_token).tolist())
    canonical_loci.discard("")

    alias_to_targets: Dict[str, Set[str]] = defaultdict(set)

    # Always seed with the canonical keys themselves
    for loc in canonical_loci:
        for a in alias_variants(loc):
            alias_to_targets[a].add(loc)

    # If GFF is present, harvest richer aliases
    if gff.exists():
        for row in iter_gff_rows(gff):
            if row["type"] != "CDS":
                continue
            attrs = row["attrs"]
            canon = canonical_feature_id(attrs)
            canon = clean_token(canon)
            if not canon:
                continue

            # Prefer the locus_map key if available
            if canonical_loci and canon not in canonical_loci:
                # Sometimes canonical GFF choice differs from locus_map key.
                # Try to connect through direct aliases.
                possible = set()
                for key in ["locus_tag", "locus", "ID", "Name"]:
                    for v in alias_variants(attrs.get(key, "")):
                        if v in canonical_loci:
                            possible.add(v)
                if len(possible) == 1:
                    canon = next(iter(possible))

            # Skip if still not mappable to a real expression key
            if canonical_loci and canon not in canonical_loci:
                continue

            alias_fields = [
                "locus_tag",
                "locus",
                "ID",
                "Name",
                "old_locus_tag",
                "protein_id",
                "transcript_id",
                "Parent",
            ]

            for key in alias_fields:
                for a in alias_variants(attrs.get(key, "")):
                    alias_to_targets[a].add(canon)

            # Parse Dbxref accessions like RefSeq:WP_..., Genbank:ABC..., etc.
            dbxref = str(attrs.get("Dbxref", "") or "")
            if dbxref:
                for item in dbxref.split(","):
                    item = clean_token(item)
                    if not item:
                        continue
                    for a in alias_variants(item):
                        alias_to_targets[a].add(canon)
                    if ":" in item:
                        rhs = item.split(":", 1)[1]
                        for a in alias_variants(rhs):
                            alias_to_targets[a].add(canon)

    # Keep only unique aliases
    resolved: Dict[str, str] = {}
    ambiguous = 0
    for alias, targets in alias_to_targets.items():
        targets = {clean_token(x) for x in targets if clean_token(x)}
        if len(targets) == 1:
            resolved[alias] = next(iter(targets))
        elif len(targets) > 1:
            ambiguous += 1

    stats = {
        "canonical_locus_count": len(canonical_loci),
        "resolved_alias_count": len(resolved),
        "ambiguous_alias_count": ambiguous,
    }
    return resolved, stats


# -----------------------------
# Panaroo -> expression mapping
# -----------------------------

def build_locus_to_orthogroup(
    gpa_csv: Path,
    genomes_dir: Path,
    out_dir: Path,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    gpa = pd.read_csv(gpa_csv, dtype=str).fillna("")
    biosamples = [bs for bs in find_biosamples(genomes_dir) if bs in gpa.columns]
    if not biosamples:
        raise RuntimeError(
            "No Panaroo sample columns matched BioSample folder names. "
            "Check that your Panaroo inputs were named with BioSample IDs."
        )

    map_rows = []
    cov_rows = []
    unmatched_rows = []

    for bs in biosamples:
        alias_map, alias_stats = build_sample_alias_map(genomes_dir, bs)

        total_tokens = 0
        matched_tokens = 0
        unmatched_examples = []

        for _, row in gpa.iterrows():
            orthogroup = clean_token(row.get("Gene", ""))
            if not orthogroup:
                continue

            tokens = split_gene_cell(row.get(bs, ""))
            if not tokens:
                continue

            for tok in tokens:
                total_tokens += 1
                hit = ""
                tried = []
                for a in alias_variants(tok):
                    tried.append(a)
                    if a in alias_map:
                        hit = alias_map[a]
                        break

                if hit:
                    matched_tokens += 1
                    map_rows.append({
                        "biosample_accession": bs,
                        "orthogroup": orthogroup,
                        "panaroo_token": tok,
                        "locus_tag": hit,
                    })
                else:
                    if len(unmatched_examples) < 200:
                        unmatched_examples.append(tok)
                        unmatched_rows.append({
                            "biosample_accession": bs,
                            "orthogroup": orthogroup,
                            "panaroo_token": tok,
                            "normalized_variants": ";".join(sorted(alias_variants(tok))),
                        })

        cov_rows.append({
            "biosample_accession": bs,
            "panaroo_tokens_seen": total_tokens,
            "tokens_matched_to_locus_map": matched_tokens,
            "match_rate": 0.0 if total_tokens == 0 else matched_tokens / total_tokens,
            **alias_stats,
        })

    mapping_df = pd.DataFrame(map_rows).drop_duplicates(
        subset=["biosample_accession", "orthogroup", "locus_tag"]
    )
    coverage_df = pd.DataFrame(cov_rows).sort_values(
        ["match_rate", "biosample_accession"], ascending=[False, True]
    )
    unmatched_df = pd.DataFrame(unmatched_rows)

    mapping_df.to_csv(out_dir / "locus_to_orthogroup.tsv", sep="\t", index=False)
    coverage_df.to_csv(out_dir / "panaroo_locus_match_coverage.tsv", sep="\t", index=False)
    unmatched_df.to_csv(out_dir / "panaroo_unmatched_examples.tsv", sep="\t", index=False)

    return mapping_df, coverage_df, unmatched_df


def aggregate_run_tpm_to_orthogroup(
    run_index_tsv: Path,
    mapping_df: pd.DataFrame,
    aggregate: str = "sum",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    run_index = pd.read_csv(run_index_tsv, sep="\t", dtype=str).fillna("")
    ok = run_index[run_index["status"] == "ok"].copy()
    if ok.empty:
        raise RuntimeError("No successful runs found in rnaseq_run_index.tsv")

    per_sample_map: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    for _, r in mapping_df.iterrows():
        bs = clean_token(r["biosample_accession"])
        loc = clean_token(r["locus_tag"])
        og = clean_token(r["orthogroup"])
        if bs and loc and og:
            per_sample_map[(bs, loc)].add(og)

    # Keep only unambiguous locus -> orthogroup assignments within sample
    final_map = {
        k: next(iter(v))
        for k, v in per_sample_map.items()
        if len(v) == 1
    }

    long_rows = []
    cov_rows = []

    for _, r in ok.iterrows():
        bs = clean_token(r["biosample_accession"])
        run_id = clean_token(r["run_id"])
        tpm_path = Path(str(r.get("tpm_tsv", "")).strip())

        if not tpm_path.exists():
            cov_rows.append({
                "biosample_accession": bs,
                "run_id": run_id,
                "total_loci_in_tpm": 0,
                "loci_mapped_to_orthogroup": 0,
                "mapping_rate": 0.0,
                "status": "missing_tpm",
            })
            continue

        tpm = pd.read_csv(tpm_path, sep="\t", dtype=str).fillna("")
        if "locus_tag" not in tpm.columns or "tpm" not in tpm.columns:
            cov_rows.append({
                "biosample_accession": bs,
                "run_id": run_id,
                "total_loci_in_tpm": 0,
                "loci_mapped_to_orthogroup": 0,
                "mapping_rate": 0.0,
                "status": "bad_tpm_schema",
            })
            continue

        tpm["locus_tag"] = tpm["locus_tag"].astype(str).map(clean_token)
        tpm["tpm"] = pd.to_numeric(tpm["tpm"], errors="coerce").fillna(0.0)
        tpm["orthogroup"] = tpm["locus_tag"].map(lambda x: final_map.get((bs, x), ""))

        total_loci = int(len(tpm))
        mapped_loci = int((tpm["orthogroup"] != "").sum())

        cov_rows.append({
            "biosample_accession": bs,
            "run_id": run_id,
            "total_loci_in_tpm": total_loci,
            "loci_mapped_to_orthogroup": mapped_loci,
            "mapping_rate": 0.0 if total_loci == 0 else mapped_loci / total_loci,
            "status": "ok",
        })

        sub = tpm[tpm["orthogroup"] != ""].copy()
        if sub.empty:
            continue

        if aggregate == "sum":
            agg = sub.groupby("orthogroup", as_index=False)["tpm"].sum()
        elif aggregate == "max":
            agg = sub.groupby("orthogroup", as_index=False)["tpm"].max()
        elif aggregate == "mean":
            agg = sub.groupby("orthogroup", as_index=False)["tpm"].mean()
        else:
            raise ValueError(f"Unsupported aggregate mode: {aggregate}")

        for _, rr in agg.iterrows():
            long_rows.append({
                "biosample_accession": bs,
                "run_id": run_id,
                "orthogroup": clean_token(rr["orthogroup"]),
                "tpm": float(rr["tpm"]),
            })

    long_df = pd.DataFrame(long_rows)
    cov_df = pd.DataFrame(cov_rows).sort_values(
        ["mapping_rate", "biosample_accession", "run_id"],
        ascending=[False, True, True]
    )

    return long_df, cov_df


def collapse_runs_to_biosample(long_df: pd.DataFrame) -> pd.DataFrame:
    if long_df.empty:
        return pd.DataFrame()

    collapsed = (
        long_df.groupby(["biosample_accession", "orthogroup"], as_index=False)["tpm"]
        .median()
    )
    wide = collapsed.pivot_table(
        index="biosample_accession",
        columns="orthogroup",
        values="tpm",
        aggfunc="first",
        fill_value=0.0,
    )
    wide = wide.sort_index(axis=0).sort_index(axis=1)
    return wide


def build_copy_number_matrix(mapping_df: pd.DataFrame) -> pd.DataFrame:
    if mapping_df.empty:
        return pd.DataFrame()

    cn = (
        mapping_df.groupby(["biosample_accession", "orthogroup"])["locus_tag"]
        .nunique()
        .reset_index(name="copy_number")
    )
    wide = cn.pivot_table(
        index="biosample_accession",
        columns="orthogroup",
        values="copy_number",
        aggfunc="first",
        fill_value=0,
    )
    wide = wide.sort_index(axis=0).sort_index(axis=1)
    return wide


def build_orthogroup_amr_summary(
    genomes_dir: Path,
    mapping_df: pd.DataFrame,
) -> pd.DataFrame:
    if mapping_df.empty:
        return pd.DataFrame()

    per_sample_map: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    for _, r in mapping_df.iterrows():
        bs = clean_token(r["biosample_accession"])
        loc = clean_token(r["locus_tag"])
        og = clean_token(r["orthogroup"])
        if bs and loc and og:
            per_sample_map[(bs, loc)].add(og)

    final_map = {
        k: next(iter(v))
        for k, v in per_sample_map.items()
        if len(v) == 1
    }

    rows = []
    for bs in find_biosamples(genomes_dir):
        amr_path = genomes_dir / bs / "amr_hits.tsv"
        if not amr_path.exists():
            continue

        amr = pd.read_csv(amr_path, sep="\t", dtype=str).fillna("")
        if "locus_tag" not in amr.columns:
            continue

        amr["biosample_accession"] = bs
        amr["locus_tag"] = amr["locus_tag"].astype(str).map(clean_token)
        amr["orthogroup"] = amr["locus_tag"].map(lambda x: final_map.get((bs, x), ""))

        sub = amr[amr["orthogroup"] != ""].copy()
        if sub.empty:
            continue

        keep_cols = [c for c in [
            "biosample_accession", "orthogroup", "locus_tag", "amr_gene", "drug_class",
            "subclass", "element_type", "element_subtype", "method",
            "identity", "coverage", "accession"
        ] if c in sub.columns]

        rows.append(sub[keep_cols])

    if not rows:
        return pd.DataFrame()

    return pd.concat(rows, ignore_index=True).drop_duplicates()


def main() -> None:
    p = argparse.ArgumentParser(
        description="Robust Panaroo -> orthogroup expression builder using GFF alias recovery."
    )
    p.add_argument("--panaroo", required=True, type=Path,
                   help="Path to Panaroo gene_presence_absence.csv")
    p.add_argument("--run-index", required=True, type=Path,
                   help="Path to rnaseq_run_index.tsv")
    p.add_argument("--genomes", required=True, type=Path,
                   help="Path to work/genomes")
    p.add_argument("--out", required=True, type=Path,
                   help="Output directory")
    p.add_argument("--aggregate", choices=["sum", "max", "mean"], default="sum",
                   help="How to combine multiple loci from the same orthogroup within a run")
    args = p.parse_args()

    out_dir = args.out.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    print("[1/4] Building robust locus -> orthogroup map")
    mapping_df, pan_cov, _ = build_locus_to_orthogroup(
        gpa_csv=args.panaroo.resolve(),
        genomes_dir=args.genomes.resolve(),
        out_dir=out_dir,
    )

    print("[2/4] Aggregating TPM to orthogroups")
    long_df, expr_cov = aggregate_run_tpm_to_orthogroup(
        run_index_tsv=args.run_index.resolve(),
        mapping_df=mapping_df,
        aggregate=args.aggregate,
    )
    long_df.to_csv(out_dir / "orthogroup_tpm_long.tsv", sep="\t", index=False)
    expr_cov.to_csv(out_dir / "orthogroup_expression_mapping_coverage.tsv", sep="\t", index=False)

    print("[3/4] Building wide matrices")
    X_og = collapse_runs_to_biosample(long_df)
    X_og.to_csv(out_dir / "X_orthogroup_tpm.tsv", sep="\t")

    X_copy = build_copy_number_matrix(mapping_df)
    X_copy.to_csv(out_dir / "X_orthogroup_copy_number.tsv", sep="\t")

    print("[4/4] Linking AMR hits to orthogroups")
    og_amr = build_orthogroup_amr_summary(
        genomes_dir=args.genomes.resolve(),
        mapping_df=mapping_df,
    )
    og_amr.to_csv(out_dir / "orthogroup_amr_hits.tsv", sep="\t", index=False)

    n_bs_pan = int((pan_cov["match_rate"] > 0).sum()) if not pan_cov.empty else 0
    n_bs_expr = int((expr_cov["mapping_rate"] > 0).sum()) if not expr_cov.empty else 0

    print()
    print(f"[OK] locus_to_orthogroup.tsv                  -> {out_dir / 'locus_to_orthogroup.tsv'}")
    print(f"[OK] panaroo_locus_match_coverage.tsv         -> {out_dir / 'panaroo_locus_match_coverage.tsv'}")
    print(f"[OK] panaroo_unmatched_examples.tsv           -> {out_dir / 'panaroo_unmatched_examples.tsv'}")
    print(f"[OK] orthogroup_tpm_long.tsv                  -> {out_dir / 'orthogroup_tpm_long.tsv'}")
    print(f"[OK] orthogroup_expression_mapping_coverage.tsv -> {out_dir / 'orthogroup_expression_mapping_coverage.tsv'}")
    print(f"[OK] X_orthogroup_tpm.tsv                     -> {out_dir / 'X_orthogroup_tpm.tsv'}")
    print(f"[OK] X_orthogroup_copy_number.tsv             -> {out_dir / 'X_orthogroup_copy_number.tsv'}")
    print(f"[OK] orthogroup_amr_hits.tsv                  -> {out_dir / 'orthogroup_amr_hits.tsv'}")
    print()
    print(f"[SUMMARY] BioSamples with Panaroo matches > 0: {n_bs_pan}")
    print(f"[SUMMARY] Runs with expression matches > 0:   {n_bs_expr}")


if __name__ == "__main__":
    main()