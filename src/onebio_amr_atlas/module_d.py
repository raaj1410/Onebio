from __future__ import annotations
from pathlib import Path
from typing import Tuple
import pandas as pd
from .utils import ensure_dir, log, PipelineError

def build_features(manifest_tsv: Path, genomes_dir: Path, rnaseq_dir: Path, out_dir: Path,
                   expr_threshold_tpm: float = 1.0) -> Tuple[Path, Path, Path]:
    ensure_dir(out_dir)
    run_index_path = Path(rnaseq_dir)/"rnaseq_run_index.tsv"
    if not run_index_path.exists():
        raise PipelineError("Missing rnaseq_run_index.tsv. Run Module C first.")

    run_index = pd.read_csv(run_index_path, sep="\t", dtype=str)
    ok = run_index[run_index["status"]=="ok"].copy()
    if ok.empty:
        raise PipelineError("No successful runs with TPM. Run Module C with --full.")

    # long expression
    rows = []
    for _, r in ok.iterrows():
        bs = str(r["biosample_accession"])
        tpm_path = r.get("tpm_tsv","")
        if not tpm_path or str(tpm_path).lower()=="nan":
            continue
        tpm = pd.read_csv(tpm_path, sep="\t", dtype=str)
        tpm["tpm"] = tpm["tpm"].astype(float)
        for _, tr in tpm.iterrows():
            rows.append({"biosample_accession": bs, "locus_tag": str(tr["locus_tag"]), "tpm": float(tr["tpm"])})

    expr_long = pd.DataFrame(rows)
    X_expr = expr_long.pivot_table(index="biosample_accession", columns="locus_tag", values="tpm", aggfunc="median", fill_value=0.0)
    X_expr_path = Path(out_dir)/"X_expr_tpm.tsv"
    X_expr.to_csv(X_expr_path, sep="\t")

    # Genome matrix (presence/absence) from locus_map.tsv
    man = pd.read_csv(manifest_tsv, sep="\t", dtype=str)
    all_loci = set()
    loci_by_bs = {}
    for _, r in man.iterrows():
        bs = str(r.get("biosample_accession",""))
        locus_map = Path(genomes_dir)/bs/"locus_map.tsv"
        if not locus_map.exists():
            continue
        loci = pd.read_csv(locus_map, sep="\t", dtype=str)["locus_tag"].dropna().astype(str).tolist()
        loci_by_bs[bs] = set(loci)
        all_loci.update(loci)

    genome_rows = []
    for bs, present in loci_by_bs.items():
        genome_rows.append({"biosample_accession": bs, **{l: (1 if l in present else 0) for l in all_loci}})
    X_genome = pd.DataFrame(genome_rows).set_index("biosample_accession")
    X_genome_path = Path(out_dir)/"X_genome.tsv"
    X_genome.to_csv(X_genome_path, sep="\t")

    # Fused: join expression + AST numeric labels
    ast_path = Path(manifest_tsv).parent/"ast_long.tsv"
    if ast_path.exists():
        ast = pd.read_csv(ast_path, sep="\t", dtype=str)
        ast_wide = ast.pivot_table(index="biosample_accession", columns="antibiotic", values="phenotype", aggfunc="first")
        def enc(x):
            if x=="R": return 1.0
            if x=="S": return 0.0
            if x=="I": return 0.5
            return float("nan")
        ast_num = ast_wide.apply(lambda col: col.map(enc))
    else:
        ast_num = pd.DataFrame(index=X_expr.index)

    fused = X_expr.copy()
    for col in ast_num.columns:
        fused[f"AST::{col}"] = ast_num[col]
    fused_path = Path(out_dir)/"X_fused.tsv"
    fused.to_csv(fused_path, sep="\t")

    log(f"X_genome: {X_genome_path}")
    log(f"X_expr_tpm: {X_expr_path}")
    log(f"X_fused: {fused_path}")
    return X_genome_path, X_expr_path, fused_path
