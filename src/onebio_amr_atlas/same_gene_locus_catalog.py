#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import sqlite3
from pathlib import Path


def has_table(con: sqlite3.Connection, name: str) -> bool:
    row = con.execute(
        """
        SELECT 1
        FROM sqlite_master
        WHERE type IN ('table', 'view') AND name = ?
        LIMIT 1
        """,
        (name,),
    ).fetchone()
    return row is not None


def table_columns(con: sqlite3.Connection, name: str) -> list[str]:
    rows = con.execute(f"PRAGMA table_info({name})").fetchall()
    return [r[1] for r in rows]


def detect_expression_value_col(con: sqlite3.Connection) -> str:
    cols = set(table_columns(con, "expression"))
    for c in ("tpm", "value", "count"):
        if c in cols:
            return c
    raise RuntimeError(
        "Could not find expression value column in 'expression'. "
        "Expected one of: tpm, value, count"
    )


def stream_query_to_tsv(
    con: sqlite3.Connection,
    query: str,
    params: tuple,
    out_tsv: Path,
    fetch_size: int = 50000,
) -> int:
    cur = con.cursor()
    cur.execute(query, params)
    headers = [d[0] for d in cur.description]

    n = 0
    with out_tsv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(headers)
        while True:
            rows = cur.fetchmany(fetch_size)
            if not rows:
                break
            writer.writerows(rows)
            n += len(rows)
    return n


def write_small_html_from_tsv(tsv_path: Path, html_path: Path, title: str, max_rows: int = 1000) -> None:
    rows = []
    with tsv_path.open("r", encoding="utf-8") as fh:
        for i, line in enumerate(fh):
            rows.append(line.rstrip("\n").split("\t"))
            if i >= max_rows:
                break

    if not rows:
        html = f"<html><body><h1>{title}</h1><p>No rows.</p></body></html>"
        html_path.write_text(html, encoding="utf-8")
        return

    header = rows[0]
    body = rows[1:]

    def esc(x: str) -> str:
        return (
            x.replace("&", "&amp;")
             .replace("<", "&lt;")
             .replace(">", "&gt;")
        )

    parts = [
        "<html><head><meta charset='utf-8'>",
        f"<title>{esc(title)}</title>",
        """
        <style>
        body { font-family: Arial, sans-serif; margin: 24px; }
        h1 { font-size: 22px; }
        table { border-collapse: collapse; width: 100%; }
        th, td { border: 1px solid #ddd; padding: 6px; font-size: 12px; }
        th { background: #f4f4f4; position: sticky; top: 0; }
        tr:nth-child(even) { background: #fafafa; }
        </style>
        """,
        "</head><body>",
        f"<h1>{esc(title)}</h1>",
        "<table><thead><tr>",
    ]
    for h in header:
        parts.append(f"<th>{esc(h)}</th>")
    parts.append("</tr></thead><tbody>")
    for row in body:
        parts.append("<tr>")
        for cell in row:
            parts.append(f"<td>{esc(cell)}</td>")
        parts.append("</tr>")
    parts.append("</tbody></table></body></html>")

    html_path.write_text("".join(parts), encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Fast SQLite-first catalog for same-gene / different-loci / different-expression"
    )
    ap.add_argument("--db", required=True, help="Path to atlas.sqlite")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--expr-threshold", type=float, default=1.0, help="Threshold to count a locus as expressed")
    ap.add_argument("--min-loci", type=int, default=2, help="Keep genes with at least this many unique locus tags")
    ap.add_argument("--min-biosamples", type=int, default=2, help="Keep genes seen in at least this many BioSamples")
    ap.add_argument("--min-delta-expr", type=float, default=5.0, help="Keep genes with at least this max-min expression spread")
    args = ap.parse_args()

    db = Path(args.db)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    con = sqlite3.connect(str(db))
    con.execute("PRAGMA temp_store = MEMORY;")
    con.execute("PRAGMA cache_size = -200000;")   # ~200 MB cache
    con.execute("PRAGMA synchronous = OFF;")
    con.execute("PRAGMA journal_mode = OFF;")
    con.execute("PRAGMA automatic_index = ON;")

    if not has_table(con, "locus"):
        raise RuntimeError("Missing table/view: locus")
    if not has_table(con, "expression"):
        raise RuntimeError("Missing table/view: expression")

    locus_cols = set(table_columns(con, "locus"))
    expr_cols = set(table_columns(con, "expression"))

    needed_locus = {"biosample_accession", "locus_tag", "gene"}
    if not needed_locus.issubset(locus_cols):
        raise RuntimeError(
            f"locus is missing required columns: {sorted(needed_locus - locus_cols)}"
        )
    if not {"biosample_accession", "locus_tag"}.issubset(expr_cols):
        raise RuntimeError("expression must contain biosample_accession and locus_tag")

    expr_value_col = detect_expression_value_col(con)

    print("[1/6] Creating helpful indexes")
    con.executescript(
        """
        CREATE INDEX IF NOT EXISTS idx_locus_bs_lt
          ON locus (biosample_accession, locus_tag);

        CREATE INDEX IF NOT EXISTS idx_locus_gene
          ON locus (gene);

        CREATE INDEX IF NOT EXISTS idx_expr_bs_lt
          ON expression (biosample_accession, locus_tag);
        """
    )
    if has_table(con, "amr_hits"):
        con.executescript(
            """
            CREATE INDEX IF NOT EXISTS idx_amr_bs_lt
              ON amr_hits (biosample_accession, locus_tag);

            CREATE INDEX IF NOT EXISTS idx_amr_gene
              ON amr_hits (amr_gene);
            """
        )

    print("[2/6] Creating temp locus-expression table")
    optional_cols = []
    for c in ["product", "seqid", "start", "end", "strand", "length"]:
        if c in locus_cols:
            optional_cols.append(f"l.{c}")

    temp_sql = f"""
    DROP TABLE IF EXISTS tmp_locus_expr;

    CREATE TEMP TABLE tmp_locus_expr AS
    SELECT
        l.biosample_accession,
        l.locus_tag,
        TRIM(COALESCE(l.gene, '')) AS gene
        {"," if optional_cols else ""} {" , ".join(optional_cols)}
        , COALESCE(e.{expr_value_col}, 0.0) AS expr_value
    FROM locus l
    LEFT JOIN expression e
      ON l.biosample_accession = e.biosample_accession
     AND l.locus_tag = e.locus_tag
    WHERE TRIM(COALESCE(l.gene, '')) <> '';
    """
    con.executescript(temp_sql)
    con.executescript(
        """
        CREATE INDEX idx_tmp_gene ON tmp_locus_expr (gene);
        CREATE INDEX idx_tmp_bs_gene ON tmp_locus_expr (biosample_accession, gene);
        CREATE INDEX idx_tmp_bs_lt ON tmp_locus_expr (biosample_accession, locus_tag);
        """
    )

    print("[3/6] Building gene summary")
    con.execute("DROP TABLE IF EXISTS tmp_sample_gene;")
    con.execute(
        """
        CREATE TEMP TABLE tmp_sample_gene AS
        SELECT
            biosample_accession,
            gene,
            COUNT(DISTINCT locus_tag) AS copy_count_in_sample,
            SUM(expr_value) AS sample_gene_total_expr,
            MAX(expr_value) AS sample_gene_max_expr,
            MIN(expr_value) AS sample_gene_min_expr
        FROM tmp_locus_expr
        GROUP BY biosample_accession, gene
        """
    )
    con.executescript(
        """
        CREATE INDEX idx_tmp_sg_bs_gene ON tmp_sample_gene (biosample_accession, gene);
        CREATE INDEX idx_tmp_sg_gene ON tmp_sample_gene (gene);
        """
    )

    gene_summary_sql = """
    WITH gene_stats AS (
        SELECT
            le.gene AS gene,
            COUNT(DISTINCT le.biosample_accession) AS biosample_n,
            COUNT(DISTINCT le.locus_tag) AS unique_locus_tags,
            COUNT(*) AS rows_total,
            SUM(CASE WHEN le.expr_value >= ? THEN 1 ELSE 0 END) AS rows_expr_positive,
            COUNT(DISTINCT CASE WHEN sg.copy_count_in_sample > 1 THEN le.biosample_accession END) AS samples_with_multicopy,
            MAX(sg.copy_count_in_sample) AS max_copy_count_in_any_sample,
            AVG(le.expr_value) AS mean_expr,
            MIN(le.expr_value) AS min_expr,
            MAX(le.expr_value) AS max_expr,
            MAX(le.expr_value) - MIN(le.expr_value) AS delta_max_min_expr,
            CASE
                WHEN AVG(le.expr_value) > 0
                THEN sqrt(AVG(le.expr_value * le.expr_value) - AVG(le.expr_value) * AVG(le.expr_value))
                ELSE NULL
            END AS std_expr,
            CASE
                WHEN AVG(le.expr_value) > 0
                THEN (
                    sqrt(AVG(le.expr_value * le.expr_value) - AVG(le.expr_value) * AVG(le.expr_value))
                    / AVG(le.expr_value)
                )
                ELSE NULL
            END AS cv_expr
        FROM tmp_locus_expr le
        JOIN tmp_sample_gene sg
          ON le.biosample_accession = sg.biosample_accession
         AND le.gene = sg.gene
        GROUP BY le.gene
    )
    SELECT *
    FROM gene_stats
    ORDER BY
        samples_with_multicopy DESC,
        delta_max_min_expr DESC,
        biosample_n DESC,
        unique_locus_tags DESC,
        gene
    """

    n_gene_summary = stream_query_to_tsv(
        con,
        gene_summary_sql,
        (args.expr_threshold,),
        outdir / "gene_summary.tsv",
    )

    print("[4/6] Building qualifying-gene catalog")
    con.execute("DROP TABLE IF EXISTS tmp_qualifying_genes;")
    con.execute(
        """
        CREATE TEMP TABLE tmp_qualifying_genes AS
        SELECT gene
        FROM (
            WITH gene_stats AS (
                SELECT
                    le.gene AS gene,
                    COUNT(DISTINCT le.biosample_accession) AS biosample_n,
                    COUNT(DISTINCT le.locus_tag) AS unique_locus_tags,
                    MAX(le.expr_value) - MIN(le.expr_value) AS delta_max_min_expr
                FROM tmp_locus_expr le
                GROUP BY le.gene
            )
            SELECT *
            FROM gene_stats
        )
        WHERE unique_locus_tags >= ?
          AND biosample_n >= ?
          AND delta_max_min_expr >= ?
        """,
        (args.min_loci, args.min_biosamples, args.min_delta_expr),
    )
    con.executescript(
        """
        CREATE INDEX idx_tmp_qg_gene ON tmp_qualifying_genes (gene);
        """
    )

    common_detail_sql = """
    SELECT
        le.gene,
        le.biosample_accession,
        le.locus_tag,
        le.expr_value,
        sg.copy_count_in_sample,
        sg.sample_gene_total_expr,
        CASE
            WHEN sg.sample_gene_total_expr > 0
            THEN le.expr_value * 1.0 / sg.sample_gene_total_expr
            ELSE 0.0
        END AS sample_fraction_of_gene_expr,
        ROW_NUMBER() OVER (
            PARTITION BY le.biosample_accession, le.gene
            ORDER BY le.expr_value DESC, le.locus_tag
        ) AS rank_within_sample_gene
        {extra_cols}
    FROM tmp_locus_expr le
    JOIN tmp_sample_gene sg
      ON le.biosample_accession = sg.biosample_accession
     AND le.gene = sg.gene
    JOIN tmp_qualifying_genes qg
      ON le.gene = qg.gene
    {where_clause}
    ORDER BY
        le.gene,
        le.biosample_accession,
        rank_within_sample_gene,
        le.locus_tag
    """

    extra_cols = ""
    for c in ["product", "seqid", "start", "end", "strand", "length"]:
        if c in locus_cols:
            extra_cols += f", le.{c}"

    n_full_catalog = stream_query_to_tsv(
        con,
        common_detail_sql.format(extra_cols=extra_cols, where_clause=""),
        (),
        outdir / "same_gene_different_loci_catalog.tsv",
    )

    print("[5/6] Building within-sample multicopy catalog")
    n_multicopy = stream_query_to_tsv(
        con,
        common_detail_sql.format(
            extra_cols=extra_cols,
            where_clause="WHERE sg.copy_count_in_sample > 1"
        ),
        (),
        outdir / "within_sample_multicopy_catalog.tsv",
    )

    if has_table(con, "amr_hits") and "amr_gene" in table_columns(con, "amr_hits"):
        print("[6/6] Building AMR-family catalogs")
        amr_summary_sql = """
        SELECT
            a.amr_gene,
            COUNT(DISTINCT le.biosample_accession) AS biosample_n,
            COUNT(DISTINCT le.locus_tag) AS unique_locus_tags,
            COUNT(*) AS rows_total,
            SUM(CASE WHEN le.expr_value >= ? THEN 1 ELSE 0 END) AS rows_expr_positive,
            AVG(le.expr_value) AS mean_expr,
            MIN(le.expr_value) AS min_expr,
            MAX(le.expr_value) AS max_expr,
            MAX(le.expr_value) - MIN(le.expr_value) AS delta_max_min_expr
        FROM tmp_locus_expr le
        JOIN amr_hits a
          ON le.biosample_accession = a.biosample_accession
         AND le.locus_tag = a.locus_tag
        WHERE TRIM(COALESCE(a.amr_gene, '')) <> ''
        GROUP BY a.amr_gene
        ORDER BY
            delta_max_min_expr DESC,
            biosample_n DESC,
            unique_locus_tags DESC,
            a.amr_gene
        """
        stream_query_to_tsv(
            con,
            amr_summary_sql,
            (args.expr_threshold,),
            outdir / "amr_family_summary.tsv",
        )

        amr_detail_sql = f"""
        SELECT
            a.amr_gene,
            le.gene,
            le.biosample_accession,
            le.locus_tag,
            le.expr_value
            {extra_cols}
        FROM tmp_locus_expr le
        JOIN amr_hits a
          ON le.biosample_accession = a.biosample_accession
         AND le.locus_tag = a.locus_tag
        WHERE TRIM(COALESCE(a.amr_gene, '')) <> ''
        ORDER BY
            a.amr_gene,
            le.biosample_accession,
            le.expr_value DESC,
            le.locus_tag
        """
        stream_query_to_tsv(
            con,
            amr_detail_sql,
            (),
            outdir / "amr_family_locus_catalog.tsv",
        )
    else:
        print("[6/6] No amr_hits/amr_gene found, skipping AMR-family outputs")

    con.close()

    print("[final] Writing quick HTML previews")
    write_small_html_from_tsv(
        outdir / "gene_summary.tsv",
        outdir / "gene_summary_top1000.html",
        "Gene summary (top 1000 rows)",
        max_rows=1000,
    )
    write_small_html_from_tsv(
        outdir / "within_sample_multicopy_catalog.tsv",
        outdir / "within_sample_multicopy_top1000.html",
        "Within-sample multicopy catalog (top 1000 rows)",
        max_rows=1000,
    )

    summary_lines = [
        f"gene_summary_rows\t{n_gene_summary}",
        f"same_gene_different_loci_catalog_rows\t{n_full_catalog}",
        f"within_sample_multicopy_catalog_rows\t{n_multicopy}",
        f"expr_threshold\t{args.expr_threshold}",
        f"min_loci\t{args.min_loci}",
        f"min_biosamples\t{args.min_biosamples}",
        f"min_delta_expr\t{args.min_delta_expr}",
    ]
    (outdir / "summary.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    print("\n[OK] Wrote:")
    print(f"  {outdir / 'gene_summary.tsv'}")
    print(f"  {outdir / 'same_gene_different_loci_catalog.tsv'}")
    print(f"  {outdir / 'within_sample_multicopy_catalog.tsv'}")
    if (outdir / "amr_family_summary.tsv").exists():
        print(f"  {outdir / 'amr_family_summary.tsv'}")
        print(f"  {outdir / 'amr_family_locus_catalog.tsv'}")
    print(f"  {outdir / 'summary.txt'}")
    print(f"  {outdir / 'gene_summary_top1000.html'}")
    print(f"  {outdir / 'within_sample_multicopy_top1000.html'}")


if __name__ == "__main__":
    main()