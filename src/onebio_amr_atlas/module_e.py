
from __future__ import annotations
from pathlib import Path
import sqlite3
import pandas as pd
from .utils import ensure_dir, log, PipelineError
from .gff import build_locus_table


EXPECTED_AMR_COLS = [
    "locus_tag", "amr_gene", "drug_class", "subclass",
    "element_type", "element_subtype", "method",
    "identity", "coverage", "accession",
]


def _read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _load_locus_rows(genomes_dir: Path, biosample_accession: str):
    bs = str(biosample_accession)
    locus_map = Path(genomes_dir) / bs / "locus_map.tsv"
    gff = Path(genomes_dir) / bs / "annotation.gff"

    rows = []
    if locus_map.exists() and locus_map.stat().st_size > 0:
        locus_df = _read_tsv(locus_map)
        for col in ["locus_tag", "gene", "product", "seqid", "start", "end", "strand", "length"]:
            if col not in locus_df.columns:
                locus_df[col] = ""
        locus_df = locus_df.drop_duplicates(subset=["locus_tag"])
        for _, r in locus_df.iterrows():
            rows.append((
                bs,
                str(r.get("locus_tag", "")),
                str(r.get("gene", "")),
                str(r.get("product", "")),
                str(r.get("seqid", "")),
                int(pd.to_numeric(r.get("start", 0), errors="coerce") if str(r.get("start", "")).strip() else 0),
                int(pd.to_numeric(r.get("end", 0), errors="coerce") if str(r.get("end", "")).strip() else 0),
                str(r.get("strand", "")),
                int(pd.to_numeric(r.get("length", 0), errors="coerce") if str(r.get("length", "")).strip() else 0),
            ))
        return rows

    if gff.exists():
        locus = build_locus_table(gff, prefer_type="CDS")
        for lr in locus:
            rows.append((
                bs,
                str(lr.get("locus_tag", "")),
                str(lr.get("gene", "")),
                str(lr.get("product", "")),
                str(lr.get("seqid", "")),
                int(lr.get("start") or 0),
                int(lr.get("end") or 0),
                str(lr.get("strand", "")),
                int(lr.get("length") or 0),
            ))
    return rows


def build_atlas(manifest_tsv: Path, genomes_dir: Path, rnaseq_dir: Path, features_dir: Path, out_dir: Path) -> Path:
    ensure_dir(out_dir)
    db_path = Path(out_dir) / "atlas.sqlite"
    if db_path.exists():
        db_path.unlink()

    X_expr_path = Path(features_dir) / "X_expr_tpm.tsv"
    if not X_expr_path.exists():
        raise PipelineError("Missing X_expr_tpm.tsv. Run Module D first.")
    X_expr = pd.read_csv(X_expr_path, sep="\t", index_col=0)

    man = pd.read_csv(manifest_tsv, sep="\t", dtype=str).fillna("")
    ast_path = Path(manifest_tsv).parent / "ast_long.tsv"
    ast = _read_tsv(ast_path) if ast_path.exists() else pd.DataFrame(columns=["biosample_accession", "antibiotic", "phenotype"])

    con = sqlite3.connect(str(db_path))
    cur = con.cursor()
    cur.executescript("""
    PRAGMA journal_mode=WAL;

    CREATE TABLE biosample (
        biosample_accession TEXT PRIMARY KEY,
        bioproject_accession TEXT,
        taxon_id TEXT,
        strain TEXT,
        genome_id TEXT,
        assembly_accession_best TEXT,
        evidence_tier TEXT
    );

    CREATE TABLE locus (
        biosample_accession TEXT,
        locus_tag TEXT,
        gene TEXT,
        product TEXT,
        seqid TEXT,
        start INT,
        end INT,
        strand TEXT,
        length INT,
        PRIMARY KEY (biosample_accession, locus_tag)
    );

    CREATE TABLE expression (
        biosample_accession TEXT,
        locus_tag TEXT,
        tpm REAL,
        PRIMARY KEY (biosample_accession, locus_tag)
    );

    CREATE TABLE ast (
        biosample_accession TEXT,
        antibiotic TEXT,
        phenotype TEXT
    );

    CREATE TABLE amr_hits (
        biosample_accession TEXT,
        locus_tag TEXT,
        amr_gene TEXT,
        drug_class TEXT,
        subclass TEXT,
        element_type TEXT,
        element_subtype TEXT,
        method TEXT,
        identity REAL,
        coverage REAL,
        accession TEXT
    );

    CREATE INDEX idx_locus_bs ON locus(biosample_accession);
    CREATE INDEX idx_expr_bs ON expression(biosample_accession);
    CREATE INDEX idx_expr_locus ON expression(locus_tag);
    CREATE INDEX idx_ast_bs ON ast(biosample_accession);
    CREATE INDEX idx_ast_ab ON ast(antibiotic);
    CREATE INDEX idx_amr_bs ON amr_hits(biosample_accession);
    CREATE INDEX idx_amr_locus ON amr_hits(locus_tag);
    CREATE INDEX idx_amr_gene ON amr_hits(amr_gene);
    CREATE INDEX idx_amr_class ON amr_hits(drug_class);
    """)

    bios_rows = []
    for _, r in man.iterrows():
        bios_rows.append((
            str(r.get("biosample_accession", "")),
            str(r.get("bioproject_accession", "")),
            str(r.get("taxon_id", "")),
            str(r.get("strain", "")),
            str(r.get("genome_id", "")),
            str(r.get("assembly_accession_best", "")),
            str(r.get("evidence_tier", "")),
        ))
    cur.executemany("INSERT INTO biosample VALUES (?,?,?,?,?,?,?)", bios_rows)

    locus_rows = []
    expr_rows = []
    bs_expr_list = X_expr.index.astype(str).tolist()

    for bs in bs_expr_list:
        locus_rows.extend(_load_locus_rows(genomes_dir, bs))
        row = X_expr.loc[bs]
        for locus_tag, tpm in row.items():
            expr_rows.append((bs, str(locus_tag), float(tpm)))

    if locus_rows:
        cur.executemany("INSERT OR REPLACE INTO locus VALUES (?,?,?,?,?,?,?,?,?)", locus_rows)
    if expr_rows:
        cur.executemany("INSERT OR REPLACE INTO expression VALUES (?,?,?)", expr_rows)

    if not ast.empty:
        for col in ["biosample_accession", "antibiotic", "phenotype"]:
            if col not in ast.columns:
                ast[col] = ""
        ast_rows = [
            (str(a), str(b), str(p))
            for a, b, p in ast[["biosample_accession", "antibiotic", "phenotype"]].values.tolist()
        ]
        cur.executemany("INSERT INTO ast VALUES (?,?,?)", ast_rows)

    amr_rows = []
    for bs in man["biosample_accession"].astype(str).dropna().tolist():
        amr_path = Path(genomes_dir) / bs / "amr_hits.tsv"
        if not amr_path.exists() or amr_path.stat().st_size == 0:
            continue
        try:
            amr = _read_tsv(amr_path)
        except Exception as e:
            log(f"WARNING: failed reading AMR hits for {bs}: {e}")
            continue
        if amr.empty:
            continue

        for col in EXPECTED_AMR_COLS:
            if col not in amr.columns:
                amr[col] = ""

        amr = amr[EXPECTED_AMR_COLS].copy()
        amr = amr.drop_duplicates(subset=["locus_tag", "amr_gene", "accession", "method"])

        amr["identity"] = pd.to_numeric(amr["identity"], errors="coerce")
        amr["coverage"] = pd.to_numeric(amr["coverage"], errors="coerce")

        for r in amr.itertuples(index=False):
            amr_rows.append((
                bs,
                str(r.locus_tag),
                str(r.amr_gene),
                str(r.drug_class),
                str(r.subclass),
                str(r.element_type),
                str(r.element_subtype),
                str(r.method),
                None if pd.isna(r.identity) else float(r.identity),
                None if pd.isna(r.coverage) else float(r.coverage),
                str(r.accession),
            ))

    if amr_rows:
        cur.executemany(
            "INSERT INTO amr_hits VALUES (?,?,?,?,?,?,?,?,?,?,?)",
            amr_rows,
        )

    cur.executescript("""
    DROP VIEW IF EXISTS amr_expression;
    CREATE VIEW amr_expression AS
    SELECT
        a.biosample_accession,
        a.locus_tag,
        l.gene,
        l.product,
        a.amr_gene,
        a.drug_class,
        a.subclass,
        a.element_type,
        a.element_subtype,
        a.method,
        a.identity,
        a.coverage,
        a.accession,
        e.tpm
    FROM amr_hits a
    LEFT JOIN locus l
      ON a.biosample_accession = l.biosample_accession
     AND a.locus_tag = l.locus_tag
    LEFT JOIN expression e
      ON a.biosample_accession = e.biosample_accession
     AND a.locus_tag = e.locus_tag;

    DROP VIEW IF EXISTS amr_summary_by_biosample;
    CREATE VIEW amr_summary_by_biosample AS
    SELECT
        biosample_accession,
        COUNT(*) AS amr_hit_n,
        COUNT(DISTINCT amr_gene) AS amr_gene_n,
        COUNT(DISTINCT drug_class) AS drug_class_n
    FROM amr_hits
    GROUP BY biosample_accession;
    """)

    con.commit()
    con.close()
    log(f"Atlas DB: {db_path}")
    return db_path


def query_gene(db: Path, gene_like: str, limit: int = 20) -> pd.DataFrame:
    con = sqlite3.connect(str(db))
    q = """
    SELECT l.gene, l.locus_tag, e.biosample_accession, e.tpm
    FROM locus l
    JOIN expression e
      ON l.biosample_accession = e.biosample_accession
     AND l.locus_tag = e.locus_tag
    WHERE l.gene LIKE ?
    ORDER BY e.tpm DESC
    LIMIT ?
    """
    df = pd.read_sql_query(q, con, params=[gene_like, limit])
    con.close()
    return df
