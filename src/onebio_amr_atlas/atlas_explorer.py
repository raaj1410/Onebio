import os
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st


# ============================================================
# Helpers
# ============================================================
def norm_text(x):
    if pd.isna(x):
        return ""
    return str(x).strip()


@st.cache_data(show_spinner=False)
def run_query(db_path: str, sql: str, params: tuple = ()) -> pd.DataFrame:
    con = sqlite3.connect(db_path)
    try:
        df = pd.read_sql_query(sql, con, params=params)
    finally:
        con.close()
    return df
@st.cache_data(show_spinner=False)
def search_genes(db_path: str, query: str, limit: int = 100) -> pd.DataFrame:
    con = sqlite3.connect(db_path)
    try:
        df = pd.read_sql_query(
            """
            SELECT gene, COUNT(*) AS locus_rows
            FROM locus
            WHERE gene IS NOT NULL
              AND TRIM(gene) <> ''
              AND gene LIKE ?
            GROUP BY gene
            ORDER BY locus_rows DESC, gene
            LIMIT ?
            """,
            con,
            params=(f"%{query}%", limit)
        )
    finally:
        con.close()
    return df


@st.cache_data(show_spinner=False)
def get_multicopy_gene_summary_sql(
    db_path: str,
    gene_space: str,
    antibiotic: str,
    only_multicopy: bool,
    top_n: int = 100
) -> pd.DataFrame:
    con = sqlite3.connect(db_path)
    try:
        if gene_space == "AMR genes only":
            if antibiotic == "ALL":
                sql = f"""
                WITH per_sample AS (
                    SELECT
                        h.amr_gene AS gene_name,
                        h.biosample_accession,
                        COUNT(DISTINCT h.locus_tag) AS locus_copies,
                        SUM(COALESCE(e.tpm, 0)) AS total_gene_tpm
                    FROM amr_hits h
                    LEFT JOIN expression e
                      ON h.biosample_accession = e.biosample_accession
                     AND h.locus_tag = e.locus_tag
                    WHERE h.amr_gene IS NOT NULL
                      AND TRIM(h.amr_gene) <> ''
                    GROUP BY h.amr_gene, h.biosample_accession
                )
                SELECT
                    gene_name,
                    COUNT(DISTINCT biosample_accession) AS carrier_biosamples,
                    SUM(CASE WHEN locus_copies > 1 THEN 1 ELSE 0 END) AS multicopy_biosamples,
                    MAX(locus_copies) AS max_locus_copies,
                    AVG(locus_copies) AS mean_locus_copies,
                    AVG(total_gene_tpm) AS mean_total_gene_tpm,
                    AVG(CASE WHEN locus_copies > 1 THEN total_gene_tpm END) AS mean_total_gene_tpm_multicopy,
                    MAX(total_gene_tpm) AS max_total_gene_tpm
                FROM per_sample
                {"WHERE locus_copies > 1" if only_multicopy else ""}
                GROUP BY gene_name
                ORDER BY multicopy_biosamples DESC, carrier_biosamples DESC, max_locus_copies DESC, mean_total_gene_tpm DESC
                LIMIT ?
                """
                params = (top_n,)
            else:
                sql = f"""
                WITH eligible_samples AS (
                    SELECT DISTINCT biosample_accession
                    FROM ast
                    WHERE antibiotic = ?
                ),
                per_sample AS (
                    SELECT
                        h.amr_gene AS gene_name,
                        h.biosample_accession,
                        COUNT(DISTINCT h.locus_tag) AS locus_copies,
                        SUM(COALESCE(e.tpm, 0)) AS total_gene_tpm
                    FROM amr_hits h
                    JOIN eligible_samples s
                      ON h.biosample_accession = s.biosample_accession
                    LEFT JOIN expression e
                      ON h.biosample_accession = e.biosample_accession
                     AND h.locus_tag = e.locus_tag
                    WHERE h.amr_gene IS NOT NULL
                      AND TRIM(h.amr_gene) <> ''
                    GROUP BY h.amr_gene, h.biosample_accession
                )
                SELECT
                    gene_name,
                    COUNT(DISTINCT biosample_accession) AS carrier_biosamples,
                    SUM(CASE WHEN locus_copies > 1 THEN 1 ELSE 0 END) AS multicopy_biosamples,
                    MAX(locus_copies) AS max_locus_copies,
                    AVG(locus_copies) AS mean_locus_copies,
                    AVG(total_gene_tpm) AS mean_total_gene_tpm,
                    AVG(CASE WHEN locus_copies > 1 THEN total_gene_tpm END) AS mean_total_gene_tpm_multicopy,
                    MAX(total_gene_tpm) AS max_total_gene_tpm
                FROM per_sample
                {"WHERE locus_copies > 1" if only_multicopy else ""}
                GROUP BY gene_name
                ORDER BY multicopy_biosamples DESC, carrier_biosamples DESC, max_locus_copies DESC, mean_total_gene_tpm DESC
                LIMIT ?
                """
                params = (antibiotic, top_n)

        else:
            if antibiotic == "ALL":
                sql = f"""
                WITH per_sample AS (
                    SELECT
                        l.gene AS gene_name,
                        l.biosample_accession,
                        COUNT(DISTINCT l.locus_tag) AS locus_copies,
                        SUM(COALESCE(e.tpm, 0)) AS total_gene_tpm
                    FROM locus l
                    LEFT JOIN expression e
                      ON l.biosample_accession = e.biosample_accession
                     AND l.locus_tag = e.locus_tag
                    WHERE l.gene IS NOT NULL
                      AND TRIM(l.gene) <> ''
                    GROUP BY l.gene, l.biosample_accession
                )
                SELECT
                    gene_name,
                    COUNT(DISTINCT biosample_accession) AS carrier_biosamples,
                    SUM(CASE WHEN locus_copies > 1 THEN 1 ELSE 0 END) AS multicopy_biosamples,
                    MAX(locus_copies) AS max_locus_copies,
                    AVG(locus_copies) AS mean_locus_copies,
                    AVG(total_gene_tpm) AS mean_total_gene_tpm,
                    AVG(CASE WHEN locus_copies > 1 THEN total_gene_tpm END) AS mean_total_gene_tpm_multicopy,
                    MAX(total_gene_tpm) AS max_total_gene_tpm
                FROM per_sample
                {"WHERE locus_copies > 1" if only_multicopy else ""}
                GROUP BY gene_name
                ORDER BY multicopy_biosamples DESC, carrier_biosamples DESC, max_locus_copies DESC, mean_total_gene_tpm DESC
                LIMIT ?
                """
                params = (top_n,)
            else:
                sql = f"""
                WITH eligible_samples AS (
                    SELECT DISTINCT biosample_accession
                    FROM ast
                    WHERE antibiotic = ?
                ),
                per_sample AS (
                    SELECT
                        l.gene AS gene_name,
                        l.biosample_accession,
                        COUNT(DISTINCT l.locus_tag) AS locus_copies,
                        SUM(COALESCE(e.tpm, 0)) AS total_gene_tpm
                    FROM locus l
                    JOIN eligible_samples s
                      ON l.biosample_accession = s.biosample_accession
                    LEFT JOIN expression e
                      ON l.biosample_accession = e.biosample_accession
                     AND l.locus_tag = e.locus_tag
                    WHERE l.gene IS NOT NULL
                      AND TRIM(l.gene) <> ''
                    GROUP BY l.gene, l.biosample_accession
                )
                SELECT
                    gene_name,
                    COUNT(DISTINCT biosample_accession) AS carrier_biosamples,
                    SUM(CASE WHEN locus_copies > 1 THEN 1 ELSE 0 END) AS multicopy_biosamples,
                    MAX(locus_copies) AS max_locus_copies,
                    AVG(locus_copies) AS mean_locus_copies,
                    AVG(total_gene_tpm) AS mean_total_gene_tpm,
                    AVG(CASE WHEN locus_copies > 1 THEN total_gene_tpm END) AS mean_total_gene_tpm_multicopy,
                    MAX(total_gene_tpm) AS max_total_gene_tpm
                FROM per_sample
                {"WHERE locus_copies > 1" if only_multicopy else ""}
                GROUP BY gene_name
                ORDER BY multicopy_biosamples DESC, carrier_biosamples DESC, max_locus_copies DESC, mean_total_gene_tpm DESC
                LIMIT ?
                """
                params = (antibiotic, top_n)

        df = pd.read_sql_query(sql, con, params=params)
    finally:
        con.close()

    return df


@st.cache_data(show_spinner=False)
def get_db_stats(db_path: str):
    stats = {}
    queries = {
        "biosamples": "SELECT COUNT(*) AS n FROM biosample",
        "bioprojects": "SELECT COUNT(DISTINCT bioproject_accession) AS n FROM biosample",
        "ast_rows": "SELECT COUNT(*) AS n FROM ast",
        "antibiotics": "SELECT COUNT(DISTINCT antibiotic) AS n FROM ast",
        "amr_hits": "SELECT COUNT(*) AS n FROM amr_hits",
        "amr_genes": "SELECT COUNT(DISTINCT amr_gene) AS n FROM amr_hits",
        "expression_rows": "SELECT COUNT(*) AS n FROM expression",
        "loci": "SELECT COUNT(*) AS n FROM locus",
    }
    for k, q in queries.items():
        stats[k] = int(run_query(db_path, q).iloc[0, 0])
    return stats


@st.cache_data(show_spinner=False)
def get_distinct_options(db_path: str):
    options = {}
    option_queries = {
        "bioprojects": """
            SELECT DISTINCT bioproject_accession
            FROM biosample
            WHERE bioproject_accession IS NOT NULL AND bioproject_accession <> ''
            ORDER BY bioproject_accession
        """,
        "evidence_tiers": """
            SELECT DISTINCT evidence_tier
            FROM biosample
            WHERE evidence_tier IS NOT NULL AND evidence_tier <> ''
            ORDER BY evidence_tier
        """,
        "antibiotics": """
            SELECT DISTINCT antibiotic
            FROM ast
            WHERE antibiotic IS NOT NULL AND antibiotic <> ''
            ORDER BY antibiotic
        """,
        "phenotypes": """
            SELECT DISTINCT phenotype
            FROM ast
            WHERE phenotype IS NOT NULL AND phenotype <> ''
            ORDER BY phenotype
        """,
        "drug_classes": """
            SELECT DISTINCT drug_class
            FROM amr_hits
            WHERE drug_class IS NOT NULL AND drug_class <> ''
            ORDER BY drug_class
        """,
        "subclasses": """
            SELECT DISTINCT subclass
            FROM amr_hits
            WHERE subclass IS NOT NULL AND subclass <> ''
            ORDER BY subclass
        """,
        "amr_genes": """
            SELECT DISTINCT amr_gene
            FROM amr_hits
            WHERE amr_gene IS NOT NULL AND amr_gene <> ''
            ORDER BY amr_gene
        """,
	"all_genes": """
    	    SELECT DISTINCT gene
            FROM locus
            WHERE gene IS NOT NULL
            AND TRIM(gene) <> ''
            ORDER BY gene
        """,
    }

    for key, sql in option_queries.items():
        df = run_query(db_path, sql)
        options[key] = df.iloc[:, 0].dropna().astype(str).tolist()

    return options


def build_filter_conditions(
    selected_projects,
    selected_evidence,
    selected_antibiotics,
    selected_phenotypes,
    selected_drug_classes,
    selected_subclasses,
    selected_amr_genes,
    identity_min,
    coverage_min,
):
    conditions = []
    params = []

    if selected_projects:
        conditions.append("b.bioproject_accession IN ({})".format(",".join(["?"] * len(selected_projects))))
        params.extend(selected_projects)

    if selected_evidence:
        conditions.append("b.evidence_tier IN ({})".format(",".join(["?"] * len(selected_evidence))))
        params.extend(selected_evidence)

    if selected_antibiotics:
        conditions.append("a.antibiotic IN ({})".format(",".join(["?"] * len(selected_antibiotics))))
        params.extend(selected_antibiotics)

    if selected_phenotypes:
        conditions.append("a.phenotype IN ({})".format(",".join(["?"] * len(selected_phenotypes))))
        params.extend(selected_phenotypes)

    if selected_drug_classes:
        conditions.append("h.drug_class IN ({})".format(",".join(["?"] * len(selected_drug_classes))))
        params.extend(selected_drug_classes)

    if selected_subclasses:
        conditions.append("h.subclass IN ({})".format(",".join(["?"] * len(selected_subclasses))))
        params.extend(selected_subclasses)

    if selected_amr_genes:
        conditions.append("h.amr_gene IN ({})".format(",".join(["?"] * len(selected_amr_genes))))
        params.extend(selected_amr_genes)

    if identity_min > 0:
        conditions.append("(h.identity IS NULL OR h.identity >= ?)")
        params.append(identity_min)

    if coverage_min > 0:
        conditions.append("(h.coverage IS NULL OR h.coverage >= ?)")
        params.append(coverage_min)

    return conditions, tuple(params)


def df_download_button(df: pd.DataFrame, label: str, filename: str):
    st.download_button(
        label=label,
        data=df.to_csv(index=False).encode("utf-8"),
        file_name=filename,
        mime="text/csv"
    )


def show_empty(msg="No data found for the current filters."):
    st.warning(msg)




def render_atlas_explorer(initial_db_path: str = "work/run1/atlas/atlas.sqlite", default_work_dir: str = "work/run1") -> None:
    # ============================================================
    # Sidebar: database + filters
    # ============================================================
    st.sidebar.title("SQLite AMR Atlas")

    default_db = initial_db_path
    db_path = st.sidebar.text_input("SQLite database path", value=default_db)

    if not db_path or not os.path.exists(db_path):
        st.error("SQLite file not found. Update the path in the sidebar.")
        st.stop()

    stats = get_db_stats(db_path)
    options = get_distinct_options(db_path)

    st.sidebar.markdown("### Global filters")
    selected_projects = st.sidebar.multiselect("Bioproject", options["bioprojects"])
    selected_evidence = st.sidebar.multiselect("Evidence tier", options["evidence_tiers"])
    selected_antibiotics = st.sidebar.multiselect("Antibiotic", options["antibiotics"])
    selected_phenotypes = st.sidebar.multiselect("Phenotype", options["phenotypes"])
    selected_drug_classes = st.sidebar.multiselect("Drug class", options["drug_classes"])
    selected_subclasses = st.sidebar.multiselect("Subclass", options["subclasses"])
    selected_amr_genes = st.sidebar.multiselect("AMR gene", options["amr_genes"])

    identity_min = st.sidebar.slider("Minimum identity", min_value=0.0, max_value=100.0, value=0.0, step=1.0)
    coverage_min = st.sidebar.slider("Minimum coverage", min_value=0.0, max_value=100.0, value=0.0, step=1.0)

    row_limit = st.sidebar.slider("Row limit for large tables", min_value=100, max_value=5000, value=1000, step=100)

    conditions, params = build_filter_conditions(
        selected_projects,
        selected_evidence,
        selected_antibiotics,
        selected_phenotypes,
        selected_drug_classes,
        selected_subclasses,
        selected_amr_genes,
        identity_min,
        coverage_min,
    )

    where_clause = ""
    if conditions:
        where_clause = "WHERE " + " AND ".join(conditions)

    st.sidebar.markdown("### What this app shows")
    st.sidebar.markdown(
        """
    - **Overview**: whole atlas at a glance  
    - **Filtered cohort**: one biosample per row  
    - **Sample deep dive**: one selected biosample  
    - **Antibiotic explorer**: AST + AMR context  
    - **AMR gene explorer**: where one gene appears and how strongly it is expressed  
    - **Gene/locus search**: search any locus, gene, or product text
    """
    )

    # ============================================================
    # Page header
    # ============================================================
    st.title("Interactive SQLite AMR Atlas")
    st.caption("Browse biosamples, AST phenotypes, AMR hits, and expression directly from your SQLite atlas.")

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Biosamples", stats["biosamples"])
    c2.metric("Bioprojects", stats["bioprojects"])
    c3.metric("AMR hits", stats["amr_hits"])
    c4.metric("Expression rows", stats["expression_rows"])

    # ============================================================
    # Tabs
    # ============================================================
    tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8, tab9 = st.tabs([
        "1. Overview",
        "2. Filtered cohort",
        "3. Sample deep dive",
        "4. Antibiotic explorer",
        "5. AMR gene explorer",
        "6. Gene / locus search",
        "7. Same gene, different loci",
        "8. Orthogroup explorer",
        "9. Knowledge graph",
    ])

    # ============================================================
    # Tab 1: Overview
    # ============================================================
    with tab1:
        st.subheader("Atlas overview")
        st.markdown(
            """
    This page is for a quick common-sense read of the atlas.

    - **Bioproject chart** shows which studies contribute most samples
    - **Antibiotic chart** shows which drugs have the most AST records
    - **Drug class chart** shows which AMR classes dominate the atlas
    """
        )

        proj_df = run_query(
            db_path,
            """
            SELECT bioproject_accession, COUNT(*) AS n_samples
            FROM biosample
            GROUP BY bioproject_accession
            ORDER BY n_samples DESC
            LIMIT 20
            """
        )
        if not proj_df.empty:
            fig = px.bar(
                proj_df,
                x="bioproject_accession",
                y="n_samples",
                title="Top bioprojects by biosample count"
            )
            st.plotly_chart(fig, use_container_width=True)

        abx_df = run_query(
            db_path,
            """
            SELECT antibiotic, COUNT(*) AS n_rows
            FROM ast
            GROUP BY antibiotic
            ORDER BY n_rows DESC
            LIMIT 25
            """
        )
        if not abx_df.empty:
            fig = px.bar(
                abx_df,
                x="antibiotic",
                y="n_rows",
                title="AST rows per antibiotic"
            )
            st.plotly_chart(fig, use_container_width=True)

        cls_df = run_query(
            db_path,
            """
            SELECT drug_class, COUNT(*) AS n_hits
            FROM amr_hits
            WHERE drug_class IS NOT NULL AND drug_class <> ''
            GROUP BY drug_class
            ORDER BY n_hits DESC
            LIMIT 20
            """
        )
        if not cls_df.empty:
            fig = px.bar(
                cls_df,
                x="drug_class",
                y="n_hits",
                title="AMR hit counts by drug class"
            )
            st.plotly_chart(fig, use_container_width=True)

    # ============================================================
    # Tab 2: Filtered cohort
    # ============================================================
    with tab2:
        st.subheader("Filtered biosample cohort")
        st.markdown(
            """
    Each row here is one **biosample**.

    Useful columns:
    - **amr_hit_n** = number of AMR hit rows in that biosample
    - **amr_gene_n** = number of distinct AMR genes
    - **drug_class_n** = number of distinct AMR classes
    - **ast_drug_n** = number of antibiotics tested for AST
    """
        )

        cohort_sql = f"""
            SELECT
                b.biosample_accession,
                b.bioproject_accession,
                b.strain,
                b.genome_id,
                b.assembly_accession_best,
                b.evidence_tier,
                COALESCE(s.amr_hit_n, 0) AS amr_hit_n,
                COALESCE(s.amr_gene_n, 0) AS amr_gene_n,
                COALESCE(s.drug_class_n, 0) AS drug_class_n,
                COUNT(DISTINCT a.antibiotic) AS ast_drug_n,
                GROUP_CONCAT(DISTINCT a.antibiotic) AS ast_antibiotics,
                GROUP_CONCAT(DISTINCT a.phenotype) AS phenotypes_seen
            FROM biosample b
            LEFT JOIN amr_summary_by_biosample s
                ON b.biosample_accession = s.biosample_accession
            LEFT JOIN ast a
                ON b.biosample_accession = a.biosample_accession
            LEFT JOIN amr_hits h
                ON b.biosample_accession = h.biosample_accession
            {where_clause}
            GROUP BY
                b.biosample_accession, b.bioproject_accession, b.strain,
                b.genome_id, b.assembly_accession_best, b.evidence_tier,
                s.amr_hit_n, s.amr_gene_n, s.drug_class_n
            ORDER BY amr_hit_n DESC, ast_drug_n DESC
            LIMIT {row_limit}
        """

        cohort_df = run_query(db_path, cohort_sql, params)

        if cohort_df.empty:
            show_empty()
        else:
            c1, c2, c3 = st.columns(3)
            c1.metric("Filtered biosamples", len(cohort_df))
            c2.metric("Median AMR genes", int(cohort_df["amr_gene_n"].median()))
            c3.metric("Median AST drugs", int(cohort_df["ast_drug_n"].median()))

            st.dataframe(cohort_df, use_container_width=True, height=450)
            df_download_button(cohort_df, "Download filtered cohort CSV", "filtered_cohort.csv")

    # ============================================================
    # Tab 3: Sample deep dive
    # ============================================================
    with tab3:
        st.subheader("Sample deep dive")

        if "cohort_df" not in locals() or cohort_df.empty:
            show_empty("No filtered biosamples available. Check the filters first.")
        else:
            sample_choice = st.selectbox(
                "Choose biosample",
                cohort_df["biosample_accession"].tolist()
            )

            meta_df = run_query(
                db_path,
                """
                SELECT *
                FROM biosample
                WHERE biosample_accession = ?
                """,
                (sample_choice,)
            )

            ast_sample = run_query(
                db_path,
                """
                SELECT antibiotic, phenotype
                FROM ast
                WHERE biosample_accession = ?
                ORDER BY antibiotic
                """,
                (sample_choice,)
            )

            amr_sample = run_query(
                db_path,
                """
                SELECT
                    locus_tag, amr_gene, drug_class, subclass,
                    element_type, element_subtype, method,
                    identity, coverage, accession
                FROM amr_hits
                WHERE biosample_accession = ?
                ORDER BY drug_class, amr_gene
                """,
                (sample_choice,)
            )

            amr_expr_sample = run_query(
                db_path,
                """
                SELECT
                    locus_tag, gene, product, amr_gene, drug_class,
                    subclass, identity, coverage, tpm
                FROM amr_expression
                WHERE biosample_accession = ?
                ORDER BY tpm DESC
                """,
                (sample_choice,)
            )

            top_expr = run_query(
                db_path,
                """
                SELECT
                    e.locus_tag,
                    l.gene,
                    l.product,
                    e.tpm
                FROM expression e
                LEFT JOIN locus l
                  ON e.biosample_accession = l.biosample_accession
                 AND e.locus_tag = l.locus_tag
                WHERE e.biosample_accession = ?
                ORDER BY e.tpm DESC
                LIMIT 50
                """,
                (sample_choice,)
            )

            m1, m2, m3, m4 = st.columns(4)
            m1.metric("AST rows", len(ast_sample))
            m2.metric("AMR hits", len(amr_sample))
            m3.metric("Expressed AMR rows", len(amr_expr_sample))
            m4.metric("Top locus TPM", round(float(top_expr["tpm"].max()), 2) if not top_expr.empty else 0.0)

            st.markdown("### Sample metadata")
            st.dataframe(meta_df, use_container_width=True)

            subtab1, subtab2, subtab3, subtab4 = st.tabs([
                "AST profile",
                "AMR hits",
                "AMR expression",
                "Top expressed loci",
            ])

            with subtab1:
                if ast_sample.empty:
                    show_empty("No AST rows for this biosample.")
                else:
                    st.dataframe(ast_sample, use_container_width=True)
                    fig = px.bar(
                        ast_sample.groupby("phenotype", as_index=False).size(),
                        x="phenotype",
                        y="size",
                        title="Phenotype counts for this biosample"
                    )
                    st.plotly_chart(fig, use_container_width=True)

            with subtab2:
                if amr_sample.empty:
                    show_empty("No AMR hits for this biosample.")
                else:
                    st.dataframe(amr_sample, use_container_width=True)
                    df_download_button(amr_sample, "Download AMR hits CSV", f"{sample_choice}_amr_hits.csv")

            with subtab3:
                if amr_expr_sample.empty:
                    show_empty("No AMR expression rows for this biosample.")
                else:
                    st.dataframe(amr_expr_sample, use_container_width=True)

                    fig = px.bar(
                        amr_expr_sample.head(20),
                        x="amr_gene",
                        y="tpm",
                        color="drug_class",
                        title="Top 20 expressed AMR-linked entries"
                    )
                    st.plotly_chart(fig, use_container_width=True)

                    fig2 = px.scatter(
                        amr_expr_sample,
                        x="identity",
                        y="tpm",
                        color="drug_class",
                        hover_data=["amr_gene", "locus_tag", "gene", "product"],
                        title="AMR expression vs identity"
                    )
                    st.plotly_chart(fig2, use_container_width=True)

            with subtab4:
                if top_expr.empty:
                    show_empty("No expression rows for this biosample.")
                else:
                    st.dataframe(top_expr, use_container_width=True)
                    fig = px.bar(
                        top_expr.head(20),
                        x="locus_tag",
                        y="tpm",
                        hover_data=["gene", "product"],
                        title="Top 20 expressed loci"
                    )
                    st.plotly_chart(fig, use_container_width=True)

    # ============================================================
    # Tab 4: Antibiotic explorer
    # ============================================================
    with tab4:
        st.subheader("Antibiotic explorer")

        abx_choice = st.selectbox("Choose antibiotic", options["antibiotics"], key="antibiotic_tab_choice")

        st.markdown(
            """
    This view answers:
    - how many samples are **R** or **S** for this antibiotic
    - which AMR genes are most commonly present in those samples
    - what AMR expression looks like for those samples
    """
        )

        abx_ast = run_query(
            db_path,
            """
            SELECT
                a.biosample_accession,
                b.bioproject_accession,
                b.strain,
                a.antibiotic,
                a.phenotype
            FROM ast a
            LEFT JOIN biosample b
              ON a.biosample_accession = b.biosample_accession
            WHERE a.antibiotic = ?
            ORDER BY a.phenotype, b.bioproject_accession, a.biosample_accession
            """,
            (abx_choice,)
        )

        if abx_ast.empty:
            show_empty("No AST data for this antibiotic.")
        else:
            c1, c2, c3 = st.columns(3)
            c1.metric("Rows", len(abx_ast))
            c2.metric("Resistant", int((abx_ast["phenotype"] == "R").sum()))
            c3.metric("Susceptible", int((abx_ast["phenotype"] == "S").sum()))

            fig = px.histogram(
                abx_ast,
                x="phenotype",
                color="phenotype",
                title=f"Phenotype distribution for {abx_choice}"
            )
            st.plotly_chart(fig, use_container_width=True)

            proj_breakdown = (
                abx_ast.groupby(["bioproject_accession", "phenotype"], as_index=False)
                .size()
                .rename(columns={"size": "n"})
            )
            if not proj_breakdown.empty:
                fig2 = px.bar(
                    proj_breakdown,
                    x="bioproject_accession",
                    y="n",
                    color="phenotype",
                    title=f"{abx_choice}: phenotype counts by bioproject"
                )
                st.plotly_chart(fig2, use_container_width=True)

            # Join this antibiotic's samples to AMR and expression
            abx_join = run_query(
                db_path,
                """
                SELECT
                    a.biosample_accession,
                    a.phenotype,
                    h.amr_gene,
                    h.drug_class,
                    h.subclass,
                    COALESCE(x.tpm, 0) AS tpm
                FROM ast a
                LEFT JOIN amr_hits h
                  ON a.biosample_accession = h.biosample_accession
                LEFT JOIN amr_expression x
                  ON h.biosample_accession = x.biosample_accession
                 AND h.locus_tag = x.locus_tag
                WHERE a.antibiotic = ?
                  AND h.amr_gene IS NOT NULL
                  AND h.amr_gene <> ''
                """,
                (abx_choice,)
            )

            if not abx_join.empty:
                st.markdown("### Top AMR genes in these samples")
                gene_counts = (
                    abx_join.groupby(["phenotype", "amr_gene"], as_index=False)
                    .agg(
                        carriers=("biosample_accession", "nunique"),
                        mean_tpm=("tpm", "mean")
                    )
                )
                top_genes = (
                    gene_counts.groupby("amr_gene", as_index=False)["carriers"]
                    .sum()
                    .sort_values("carriers", ascending=False)
                    .head(20)["amr_gene"]
                    .tolist()
                )
                gene_counts_top = gene_counts[gene_counts["amr_gene"].isin(top_genes)]

                fig3 = px.bar(
                    gene_counts_top,
                    x="amr_gene",
                    y="carriers",
                    color="phenotype",
                    barmode="group",
                    hover_data=["mean_tpm"],
                    title=f"Top AMR genes among samples tested for {abx_choice}"
                )
                st.plotly_chart(fig3, use_container_width=True)

                tpm_gene = (
                    abx_join.groupby(["phenotype", "amr_gene"], as_index=False)
                    .agg(median_tpm=("tpm", "median"))
                )
                tpm_gene = tpm_gene[tpm_gene["amr_gene"].isin(top_genes)]
                fig4 = px.bar(
                    tpm_gene,
                    x="amr_gene",
                    y="median_tpm",
                    color="phenotype",
                    barmode="group",
                    title=f"Median AMR TPM by phenotype for {abx_choice}"
                )
                st.plotly_chart(fig4, use_container_width=True)

                st.dataframe(gene_counts.sort_values(["carriers", "mean_tpm"], ascending=[False, False]), use_container_width=True)
                df_download_button(gene_counts, "Download antibiotic gene summary CSV", f"{abx_choice}_gene_summary.csv")

    # ============================================================
    # Tab 5: AMR gene explorer
    # ============================================================
    with tab5:
        st.subheader("AMR gene explorer")

        gene_choice = st.selectbox("Choose AMR gene", options["amr_genes"], key="amr_gene_tab_choice")

        st.markdown(
            """
    This page shows where one AMR gene appears:
    - which biosamples carry it
    - which studies it appears in
    - what its TPM looks like
    - which antibiotics and phenotypes are linked to those carriers
    """
        )

        gene_df = run_query(
            db_path,
            """
            SELECT
                h.biosample_accession,
                b.bioproject_accession,
                b.strain,
                h.locus_tag,
                h.amr_gene,
                h.drug_class,
                h.subclass,
                h.identity,
                h.coverage,
                COALESCE(x.tpm, 0) AS tpm
            FROM amr_hits h
            LEFT JOIN biosample b
              ON h.biosample_accession = b.biosample_accession
            LEFT JOIN amr_expression x
              ON h.biosample_accession = x.biosample_accession
             AND h.locus_tag = x.locus_tag
            WHERE h.amr_gene = ?
            ORDER BY tpm DESC
            LIMIT 5000
            """,
            (gene_choice,)
        )

        if gene_df.empty:
            show_empty("No records for this AMR gene.")
        else:
            c1, c2, c3 = st.columns(3)
            c1.metric("Carrier rows", len(gene_df))
            c2.metric("Carrier biosamples", gene_df["biosample_accession"].nunique())
            c3.metric("Median TPM", round(float(gene_df["tpm"].median()), 2))

            proj_gene = (
                gene_df.groupby("bioproject_accession", as_index=False)
                .agg(carriers=("biosample_accession", "nunique"))
                .sort_values("carriers", ascending=False)
                .head(20)
            )
            fig = px.bar(
                proj_gene,
                x="bioproject_accession",
                y="carriers",
                title=f"Carriers of {gene_choice} by bioproject"
            )
            st.plotly_chart(fig, use_container_width=True)

            fig2 = px.histogram(
                gene_df,
                x="tpm",
                nbins=40,
                title=f"TPM distribution for {gene_choice}"
            )
            st.plotly_chart(fig2, use_container_width=True)

            gene_ast = run_query(
                db_path,
                """
                SELECT
                    h.biosample_accession,
                    a.antibiotic,
                    a.phenotype,
                    COALESCE(x.tpm, 0) AS tpm
                FROM amr_hits h
                LEFT JOIN ast a
                  ON h.biosample_accession = a.biosample_accession
                LEFT JOIN amr_expression x
                  ON h.biosample_accession = x.biosample_accession
                 AND h.locus_tag = x.locus_tag
                WHERE h.amr_gene = ?
                  AND a.antibiotic IS NOT NULL
                """,
                (gene_choice,)
            )

            if not gene_ast.empty:
                st.markdown("### Antibiotic / phenotype context for this gene")
                gene_ast_sum = (
                    gene_ast.groupby(["antibiotic", "phenotype"], as_index=False)
                    .agg(
                        biosamples=("biosample_accession", "nunique"),
                        median_tpm=("tpm", "median")
                    )
                )
                fig3 = px.bar(
                    gene_ast_sum,
                    x="antibiotic",
                    y="biosamples",
                    color="phenotype",
                    barmode="group",
                    hover_data=["median_tpm"],
                    title=f"{gene_choice}: carriers across antibiotics"
                )
                st.plotly_chart(fig3, use_container_width=True)

                st.dataframe(gene_ast_sum.sort_values(["biosamples", "median_tpm"], ascending=[False, False]), use_container_width=True)

            st.markdown("### Raw carrier rows")
            st.dataframe(gene_df, use_container_width=True)
            df_download_button(gene_df, "Download gene carrier table CSV", f"{gene_choice}_carriers.csv")

    # ============================================================
    # Tab 6: Gene/locus search
    # ============================================================
    with tab6:
        st.subheader("Gene / locus search")

        search_mode = st.radio(
            "Search by",
            ["locus_tag", "gene", "product"],
            horizontal=True
        )
        search_term = st.text_input("Enter search text")

        if search_term:
            if search_mode == "locus_tag":
                sql = """
                    SELECT
                        l.biosample_accession,
                        b.bioproject_accession,
                        l.locus_tag,
                        l.gene,
                        l.product,
                        l.seqid,
                        l.start,
                        l.end,
                        l.strand,
                        l.length,
                        COALESCE(e.tpm, 0) AS tpm
                    FROM locus l
                    LEFT JOIN biosample b
                      ON l.biosample_accession = b.biosample_accession
                    LEFT JOIN expression e
                      ON l.biosample_accession = e.biosample_accession
                     AND l.locus_tag = e.locus_tag
                    WHERE l.locus_tag LIKE ?
                    ORDER BY tpm DESC
                    LIMIT ?
                """
            elif search_mode == "gene":
                sql = """
                    SELECT
                        l.biosample_accession,
                        b.bioproject_accession,
                        l.locus_tag,
                        l.gene,
                        l.product,
                        l.seqid,
                        l.start,
                        l.end,
                        l.strand,
                        l.length,
                        COALESCE(e.tpm, 0) AS tpm
                    FROM locus l
                    LEFT JOIN biosample b
                      ON l.biosample_accession = b.biosample_accession
                    LEFT JOIN expression e
                      ON l.biosample_accession = e.biosample_accession
                     AND l.locus_tag = e.locus_tag
                    WHERE l.gene LIKE ?
                    ORDER BY tpm DESC
                    LIMIT ?
                """
            else:
                sql = """
                    SELECT
                        l.biosample_accession,
                        b.bioproject_accession,
                        l.locus_tag,
                        l.gene,
                        l.product,
                        l.seqid,
                        l.start,
                        l.end,
                        l.strand,
                        l.length,
                        COALESCE(e.tpm, 0) AS tpm
                    FROM locus l
                    LEFT JOIN biosample b
                      ON l.biosample_accession = b.biosample_accession
                    LEFT JOIN expression e
                      ON l.biosample_accession = e.biosample_accession
                     AND l.locus_tag = e.locus_tag
                    WHERE l.product LIKE ?
                    ORDER BY tpm DESC
                    LIMIT ?
                """

            search_df = run_query(db_path, sql, (f"%{search_term}%", row_limit))

            if search_df.empty:
                show_empty("No matching rows found.")
            else:
                st.dataframe(search_df, use_container_width=True)
                df_download_button(search_df, "Download search results CSV", f"{search_mode}_search_results.csv")

                fig = px.scatter(
                    search_df,
                    x="length",
                    y="tpm",
                    hover_data=["biosample_accession", "locus_tag", "gene", "product"],
                    title=f"Search results: {search_mode} contains '{search_term}'"
                )
                st.plotly_chart(fig, use_container_width=True)
    # ============================================================
    # Tab 7: Same Gene, Different loci
    # ============================================================
    with tab7:
        st.subheader("Same gene, different loci")
        st.markdown(
            """
    This section explores whether the **same gene** appears in **multiple loci** within a biosample,
    and whether those copies show different expression.

    ### What the key numbers mean
    - **Gene instance** = one `biosample_accession + locus_tag`
    - **Locus copies in sample** = number of distinct loci carrying that gene inside one biosample
    - **Total gene TPM per sample** = sum of TPM across all copies of that gene in that biosample
    - **Multi-copy biosample** = a biosample where the same gene appears in more than one locus
    """
        )

        # ------------------------------------------------------------
        # ROW 1: main controls
        # ------------------------------------------------------------
        r1c1, r1c2, r1c3, r1c4 = st.columns([1.2, 1.2, 1.0, 1.0])

        gene_space = r1c1.radio(
            "Gene source",
            ["AMR genes only", "All annotated genes"],
            horizontal=False,
            key="tab7_gene_space"
        )

        antibiotic_choices = ["ALL"] + options["antibiotics"]
        tab7_abx = r1c2.selectbox(
            "Choose antibiotic",
            antibiotic_choices,
            key="tab7_antibiotic_choice"
        )

        only_multicopy = r1c3.checkbox(
            "Only multicopy records",
            value=False,
            key="tab7_only_multicopy"
        )

        top_n_genes = r1c4.slider(
            "Top genes",
            min_value=20,
            max_value=500,
            value=100,
            step=20,
            key="tab7_topn_genes"
        )

        # ------------------------------------------------------------
        # ROW 2: gene selection controls
        # ------------------------------------------------------------
        tab7_gene = "ALL"
        gene_hits = pd.DataFrame()

        if gene_space == "AMR genes only":
            st.markdown("### Gene selection")
            amr_choices = ["ALL"] + options["amr_genes"]
            tab7_gene = st.selectbox(
                "Choose AMR gene or keep ALL",
                amr_choices,
                key="tab7_gene_choice_amr"
            )
            st.caption("Use ALL to rank genes by multicopy behaviour, or choose one gene for detailed copy-level analysis.")

        else:
            st.markdown("### Annotated gene search")
            s1, s2 = st.columns([1.4, 1.2])

            gene_search = s1.text_input(
                "Search annotated gene name",
                value="",
                placeholder="Type a gene name like oqxA, acrB, rpoB...",
                key="tab7_gene_search"
            )

            if gene_search.strip():
                gene_hits = search_genes(db_path, gene_search.strip(), limit=100)

                if gene_hits.empty:
                    s2.warning("No matching genes found.")
                    tab7_gene = "ALL"
                else:
                    tab7_gene = s2.selectbox(
                        "Select one matching gene or keep ALL",
                        ["ALL"] + gene_hits["gene"].tolist(),
                        key="tab7_gene_choice_all"
                    )

                    st.caption("Matching annotated genes")
                    st.dataframe(gene_hits, use_container_width=True, height=180)
            else:
                tab7_gene = "ALL"
                st.info("Leave search blank to summarise all annotated genes, or type a gene name to inspect one specific gene.")

        # ------------------------------------------------------------
        # FAST SUMMARY MODE: Gene = ALL
        # ------------------------------------------------------------
        if tab7_gene == "ALL":
            gene_summary = get_multicopy_gene_summary_sql(
                db_path=db_path,
                gene_space=gene_space,
                antibiotic=tab7_abx,
                only_multicopy=only_multicopy,
                top_n=top_n_genes
            )

            if gene_summary.empty:
                show_empty("No genes found for the current selection.")
            else:
                m1, m2, m3, m4 = st.columns(4)
                m1.metric("Genes shown", gene_summary["gene_name"].nunique())
                m2.metric("Total carrier biosamples", int(gene_summary["carrier_biosamples"].sum()))
                m3.metric("Genes with multicopy", int((gene_summary["multicopy_biosamples"] > 0).sum()))
                m4.metric("Highest copy count", int(gene_summary["max_locus_copies"].max()))

                fig1 = px.bar(
                    gene_summary.sort_values(
                        ["multicopy_biosamples", "carrier_biosamples", "max_locus_copies", "mean_total_gene_tpm"],
                        ascending=True
                    ),
                    x="multicopy_biosamples",
                    y="gene_name",
                    orientation="h",
                    hover_data=[
                        "carrier_biosamples",
                        "max_locus_copies",
                        "mean_locus_copies",
                        "mean_total_gene_tpm",
                        "mean_total_gene_tpm_multicopy",
                        "max_total_gene_tpm"
                    ],
                    title="Genes ranked by multicopy biosamples"
                )
                st.plotly_chart(fig1, use_container_width=True)

                fig2 = px.scatter(
                    gene_summary,
                    x="carrier_biosamples",
                    y="multicopy_biosamples",
                    size="mean_total_gene_tpm",
                    hover_name="gene_name",
                    hover_data=[
                        "max_locus_copies",
                        "mean_locus_copies",
                        "mean_total_gene_tpm_multicopy",
                        "max_total_gene_tpm"
                    ],
                    title="Common genes vs multicopy behaviour"
                )
                st.plotly_chart(fig2, use_container_width=True)

                st.markdown("### Gene multicopy summary")
                st.dataframe(gene_summary, use_container_width=True)
                df_download_button(
                    gene_summary,
                    "Download gene multicopy summary CSV",
                    f"tab7_{gene_space.replace(' ', '_')}_{tab7_abx}_gene_summary.csv"
                )

        # ------------------------------------------------------------
        # DETAILED MODE: one selected gene
        # ------------------------------------------------------------
        else:
            params = []

            if gene_space == "AMR genes only":
                if tab7_abx == "ALL":
                    sql = """
                        WITH ast_summary AS (
                            SELECT
                                biosample_accession,
                                GROUP_CONCAT(DISTINCT antibiotic || ':' || phenotype) AS ast_profile
                            FROM ast
                            GROUP BY biosample_accession
                        )
                        SELECT
                            b.biosample_accession,
                            b.bioproject_accession,
                            b.strain,
                            asts.ast_profile,
                            h.amr_gene AS gene_name,
                            h.locus_tag,
                            l.gene AS locus_gene,
                            l.product,
                            l.seqid,
                            l.start,
                            l.end,
                            l.strand,
                            l.length,
                            h.drug_class,
                            h.subclass,
                            h.identity,
                            h.coverage,
                            COALESCE(e.tpm, 0) AS tpm
                        FROM amr_hits h
                        JOIN biosample b
                          ON h.biosample_accession = b.biosample_accession
                        LEFT JOIN locus l
                          ON h.biosample_accession = l.biosample_accession
                         AND h.locus_tag = l.locus_tag
                        LEFT JOIN expression e
                          ON h.biosample_accession = e.biosample_accession
                         AND h.locus_tag = e.locus_tag
                        LEFT JOIN ast_summary asts
                          ON h.biosample_accession = asts.biosample_accession
                        WHERE h.amr_gene = ?
                    """
                    params.append(tab7_gene)

                else:
                    sql = """
                        SELECT
                            a.antibiotic,
                            a.phenotype,
                            b.biosample_accession,
                            b.bioproject_accession,
                            b.strain,
                            h.amr_gene AS gene_name,
                            h.locus_tag,
                            l.gene AS locus_gene,
                            l.product,
                            l.seqid,
                            l.start,
                            l.end,
                            l.strand,
                            l.length,
                            h.drug_class,
                            h.subclass,
                            h.identity,
                            h.coverage,
                            COALESCE(e.tpm, 0) AS tpm
                        FROM ast a
                        JOIN amr_hits h
                          ON a.biosample_accession = h.biosample_accession
                        JOIN biosample b
                          ON h.biosample_accession = b.biosample_accession
                        LEFT JOIN locus l
                          ON h.biosample_accession = l.biosample_accession
                         AND h.locus_tag = l.locus_tag
                        LEFT JOIN expression e
                          ON h.biosample_accession = e.biosample_accession
                         AND h.locus_tag = e.locus_tag
                        WHERE a.antibiotic = ?
                          AND h.amr_gene = ?
                    """
                    params.extend([tab7_abx, tab7_gene])

            else:
                if tab7_abx == "ALL":
                    sql = """
                        WITH ast_summary AS (
                            SELECT
                                biosample_accession,
                                GROUP_CONCAT(DISTINCT antibiotic || ':' || phenotype) AS ast_profile
                            FROM ast
                            GROUP BY biosample_accession
                        )
                        SELECT
                            b.biosample_accession,
                            b.bioproject_accession,
                            b.strain,
                            asts.ast_profile,
                            l.gene AS gene_name,
                            l.locus_tag,
                            l.gene AS locus_gene,
                            l.product,
                            l.seqid,
                            l.start,
                            l.end,
                            l.strand,
                            l.length,
                            NULL AS drug_class,
                            NULL AS subclass,
                            NULL AS identity,
                            NULL AS coverage,
                            COALESCE(e.tpm, 0) AS tpm
                        FROM locus l
                        JOIN biosample b
                          ON l.biosample_accession = b.biosample_accession
                        LEFT JOIN expression e
                          ON l.biosample_accession = e.biosample_accession
                         AND l.locus_tag = e.locus_tag
                        LEFT JOIN ast_summary asts
                          ON l.biosample_accession = asts.biosample_accession
                        WHERE l.gene = ?
                    """
                    params.append(tab7_gene)

                else:
                    sql = """
                        SELECT
                            a.antibiotic,
                            a.phenotype,
                            b.biosample_accession,
                            b.bioproject_accession,
                            b.strain,
                            l.gene AS gene_name,
                            l.locus_tag,
                            l.gene AS locus_gene,
                            l.product,
                            l.seqid,
                            l.start,
                            l.end,
                            l.strand,
                            l.length,
                            NULL AS drug_class,
                            NULL AS subclass,
                            NULL AS identity,
                            NULL AS coverage,
                            COALESCE(e.tpm, 0) AS tpm
                        FROM ast a
                        JOIN locus l
                          ON a.biosample_accession = l.biosample_accession
                        JOIN biosample b
                          ON l.biosample_accession = b.biosample_accession
                        LEFT JOIN expression e
                          ON l.biosample_accession = e.biosample_accession
                         AND l.locus_tag = e.locus_tag
                        WHERE a.antibiotic = ?
                          AND l.gene = ?
                    """
                    params.extend([tab7_abx, tab7_gene])

            gene_loci_df = run_query(db_path, sql, tuple(params))

            if gene_loci_df.empty:
                show_empty("No rows found for the current selection.")
            else:
                gene_loci_df["instance_id"] = (
                    gene_loci_df["biosample_accession"].astype(str)
                    + " | "
                    + gene_loci_df["locus_tag"].astype(str)
                )

                if tab7_abx == "ALL":
                    group_cols = [
                        "gene_name", "biosample_accession", "bioproject_accession",
                        "strain", "ast_profile"
                    ]
                else:
                    group_cols = [
                        "gene_name", "antibiotic", "phenotype",
                        "biosample_accession", "bioproject_accession", "strain"
                    ]

                sample_copy_df = (
                    gene_loci_df.groupby(group_cols, as_index=False)
                    .agg(
                        locus_copies=("locus_tag", "nunique"),
                        total_gene_tpm=("tpm", "sum"),
                        mean_copy_tpm=("tpm", "mean"),
                        median_copy_tpm=("tpm", "median"),
                        max_copy_tpm=("tpm", "max"),
                        min_copy_tpm=("tpm", "min")
                    )
                )

                sample_copy_df["copy_status"] = np.where(
                    sample_copy_df["locus_copies"] > 1,
                    "Multi-copy",
                    "Single-copy"
                )

                if only_multicopy:
                    sample_copy_df = sample_copy_df[sample_copy_df["locus_copies"] > 1].copy()
                    keep_pairs = sample_copy_df[["gene_name", "biosample_accession"]].drop_duplicates()
                    gene_loci_df = gene_loci_df.merge(
                        keep_pairs,
                        on=["gene_name", "biosample_accession"],
                        how="inner"
                    )

                if sample_copy_df.empty:
                    show_empty("No rows remain after applying the multi-copy filter.")
                else:
                    m1, m2, m3, m4 = st.columns(4)
                    m1.metric("Gene instances", gene_loci_df["instance_id"].nunique())
                    m2.metric("Carrier biosamples", sample_copy_df["biosample_accession"].nunique())
                    m3.metric("Multi-copy biosamples", int((sample_copy_df["locus_copies"] > 1).sum()))
                    m4.metric("Median total gene TPM", round(float(sample_copy_df["total_gene_tpm"].median()), 2))

                    if tab7_abx == "ALL":
                        st.markdown("### Per-biosample copy summary")
                        st.dataframe(
                            sample_copy_df.sort_values(["locus_copies", "total_gene_tpm"], ascending=[False, False]),
                            use_container_width=True
                        )

                        fig1 = px.box(
                            sample_copy_df,
                            x="copy_status",
                            y="total_gene_tpm",
                            points="all",
                            hover_data=["biosample_accession", "bioproject_accession", "locus_copies", "ast_profile"],
                            title=f"Total {tab7_gene} TPM by copy status across all antibiotics"
                        )
                        st.plotly_chart(fig1, use_container_width=True)

                        copy_dist = (
                            sample_copy_df.groupby("locus_copies", as_index=False)
                            .size()
                            .rename(columns={"size": "biosample_count"})
                        )
                        fig2 = px.bar(
                            copy_dist,
                            x="locus_copies",
                            y="biosample_count",
                            title=f"How many copies of {tab7_gene} occur per biosample?"
                        )
                        st.plotly_chart(fig2, use_container_width=True)

                    else:
                        phenotype_summary = (
                            sample_copy_df.groupby("phenotype", as_index=False)
                            .agg(
                                biosamples=("biosample_accession", "nunique"),
                                mean_locus_copies=("locus_copies", "mean"),
                                multi_copy_biosamples=("copy_status", lambda s: int((s == "Multi-copy").sum())),
                                mean_total_gene_tpm=("total_gene_tpm", "mean"),
                                median_total_gene_tpm=("total_gene_tpm", "median"),
                                max_total_gene_tpm=("total_gene_tpm", "max")
                            )
                        )

                        st.markdown("### Phenotype-level summary")
                        st.dataframe(phenotype_summary, use_container_width=True)

                        p1, p2 = st.columns(2)

                        with p1:
                            fig1 = px.box(
                                sample_copy_df,
                                x="phenotype",
                                y="total_gene_tpm",
                                color="copy_status",
                                points="all",
                                hover_data=["biosample_accession", "bioproject_accession", "locus_copies", "max_copy_tpm"],
                                title=f"Total {tab7_gene} TPM per biosample by phenotype"
                            )
                            st.plotly_chart(fig1, use_container_width=True)

                        with p2:
                            copy_dist = (
                                sample_copy_df.groupby(["phenotype", "locus_copies"], as_index=False)
                                .size()
                                .rename(columns={"size": "biosample_count"})
                            )
                            fig2 = px.bar(
                                copy_dist,
                                x="locus_copies",
                                y="biosample_count",
                                color="phenotype",
                                barmode="group",
                                title=f"How many locus copies of {tab7_gene} occur per biosample?"
                            )
                            st.plotly_chart(fig2, use_container_width=True)

                    st.markdown("### Copy-level expression")
                    color_col = "phenotype" if tab7_abx != "ALL" else "bioproject_accession"

                    hover_cols = [
                        "biosample_accession", "locus_tag", "locus_gene", "product",
                        "seqid", "start", "end", "strand"
                    ]
                    if "identity" in gene_loci_df.columns:
                        hover_cols += ["identity", "coverage"]

                    fig3 = px.strip(
                        gene_loci_df,
                        x=color_col,
                        y="tpm",
                        color=color_col,
                        hover_data=hover_cols,
                        title=f"Each dot = one locus instance carrying {tab7_gene}"
                    )
                    st.plotly_chart(fig3, use_container_width=True)

                    instance_summary = (
                        gene_loci_df.groupby(["instance_id"], as_index=False)
                        .agg(
                            biosample_accession=("biosample_accession", "first"),
                            bioproject_accession=("bioproject_accession", "first"),
                            strain=("strain", "first"),
                            locus_tag=("locus_tag", "first"),
                            locus_gene=("locus_gene", "first"),
                            product=("product", "first"),
                            seqid=("seqid", "first"),
                            start=("start", "first"),
                            end=("end", "first"),
                            strand=("strand", "first"),
                            length=("length", "first"),
                            tpm=("tpm", "first")
                        )
                        .sort_values("tpm", ascending=False)
                    )

                    st.markdown("### Top locus instances by TPM")
                    top_instances = instance_summary.head(30).copy()
                    fig4 = px.bar(
                        top_instances.sort_values("tpm", ascending=True),
                        x="tpm",
                        y="instance_id",
                        orientation="h",
                        hover_data=["biosample_accession", "bioproject_accession", "locus_tag", "locus_gene", "product"],
                        title=f"Top 30 locus instances carrying {tab7_gene}"
                    )
                    st.plotly_chart(fig4, use_container_width=True)

                    st.markdown("### Per-biosample copy summary")
                    st.dataframe(
                        sample_copy_df.sort_values(["locus_copies", "total_gene_tpm"], ascending=[False, False]),
                        use_container_width=True
                    )
                    df_download_button(
                        sample_copy_df,
                        f"Download per-biosample copy summary for {tab7_gene}",
                        f"tab7_{tab7_abx}_{tab7_gene}_per_biosample_copy_summary.csv"
                    )

                    # ------------------------------------------------------------
                    # NEW: copy-number to locus drilldown
                    # ------------------------------------------------------------
                    st.markdown("### Drilldown: inspect the loci behind a copy number")

                    drill_c1, drill_c2 = st.columns([1.0, 1.8])

                    available_copy_numbers = sorted(sample_copy_df["locus_copies"].dropna().unique().tolist())

                    selected_copy_number = drill_c1.selectbox(
                        "Choose copy number",
                        available_copy_numbers,
                        key=f"tab7_copy_number_{tab7_gene}_{tab7_abx}"
                    )

                    copy_filtered_samples = sample_copy_df[
                        sample_copy_df["locus_copies"] == selected_copy_number
                    ].copy()

                    if "phenotype" in copy_filtered_samples.columns:
                        copy_filtered_samples["sample_label"] = (
                            copy_filtered_samples["biosample_accession"].astype(str)
                            + " | phenotype=" + copy_filtered_samples["phenotype"].astype(str)
                            + " | total_tpm=" + copy_filtered_samples["total_gene_tpm"].round(2).astype(str)
                        )
                    else:
                        copy_filtered_samples["sample_label"] = (
                            copy_filtered_samples["biosample_accession"].astype(str)
                            + " | total_tpm=" + copy_filtered_samples["total_gene_tpm"].round(2).astype(str)
                        )

                    selected_sample_label = drill_c2.selectbox(
                        "Choose biosample to inspect",
                        copy_filtered_samples["sample_label"].tolist(),
                        key=f"tab7_sample_label_{tab7_gene}_{tab7_abx}"
                    )

                    selected_sample_row = copy_filtered_samples[
                        copy_filtered_samples["sample_label"] == selected_sample_label
                    ].iloc[0]

                    selected_biosample = selected_sample_row["biosample_accession"]

                    locus_drill_df = gene_loci_df[
                        gene_loci_df["biosample_accession"] == selected_biosample
                    ].copy()

                    locus_drill_df = locus_drill_df.sort_values(
                        ["seqid", "start", "end", "tpm"],
                        ascending=[True, True, True, False]
                    )

                    st.markdown("#### Selected biosample summary")
                    s1, s2, s3, s4 = st.columns(4)
                    s1.metric("Biosample", str(selected_biosample))
                    s2.metric("Copy count", int(selected_sample_row["locus_copies"]))
                    s3.metric("Total gene TPM", round(float(selected_sample_row["total_gene_tpm"]), 2))
                    s4.metric("Max copy TPM", round(float(selected_sample_row["max_copy_tpm"]), 2))

                    drill_show_cols = [
                        "biosample_accession", "locus_tag", "seqid", "start", "end", "strand",
                        "locus_gene", "product", "tpm"
                    ]

                    for extra_col in ["drug_class", "subclass", "identity", "coverage"]:
                        if extra_col in locus_drill_df.columns:
                            drill_show_cols.append(extra_col)

                    drill_show_cols = [c for c in drill_show_cols if c in locus_drill_df.columns]

                    st.markdown("#### Locus instances for this biosample")
                    st.dataframe(
                        locus_drill_df[drill_show_cols],
                        use_container_width=True
                    )

                    df_download_button(
                        locus_drill_df[drill_show_cols],
                        f"Download loci for {selected_biosample}",
                        f"tab7_{tab7_gene}_{selected_biosample}_loci.csv"
                    )

                    if {"seqid", "start", "end", "locus_tag", "tpm"}.issubset(locus_drill_df.columns):
                        plot_df = locus_drill_df.copy()
                        plot_df["midpoint"] = (plot_df["start"] + plot_df["end"]) / 2
                        plot_df["span"] = plot_df["end"] - plot_df["start"] + 1

                        fig_drill = px.scatter(
                            plot_df,
                            x="midpoint",
                            y="seqid",
                            size="span",
                            color="tpm",
                            hover_data=["locus_tag", "start", "end", "strand", "product"],
                            title=f"Locus positions for {tab7_gene} in {selected_biosample}"
                        )
                        st.plotly_chart(fig_drill, use_container_width=True)

                    st.markdown("### Raw locus-instance table")
                    st.dataframe(
                        gene_loci_df.sort_values(["tpm"], ascending=False),
                        use_container_width=True
                    )
                    df_download_button(
                        gene_loci_df,
                        f"Download raw locus-instance table for {tab7_gene}",
                        f"tab7_{tab7_abx}_{tab7_gene}_raw_locus_instances.csv"
                    )

    # ============================================================
    # Tab 8: Orthogroup explorer
    # ============================================================
    with tab8:
        st.subheader("Orthogroup explorer")
        st.markdown(
            """
    This tab reads outputs from the Panaroo orthogroup-expression builder.

    Expected files inside the selected orthogroup output folder:
    - `orthogroup_tpm_long.tsv`
    - `X_orthogroup_tpm.tsv`
    - `X_orthogroup_copy_number.tsv`
    - `orthogroup_amr_hits.tsv`
    - `orthogroup_expression_mapping_coverage.tsv`
    """
        )

        default_og_dir = str(Path(default_work_dir) / "orthogroups")
        og_dir_txt = st.text_input("Orthogroup output folder", value=default_og_dir, key="tab8_og_dir")
        og_dir = Path(og_dir_txt)

        if not og_dir.exists():
            st.warning("Orthogroup output folder not found yet. Run the orthogroup builder after Panaroo produces gene_presence_absence.csv.")
            st.code(
                f"onebio-atlas orthogroup --panaroo <gene_presence_absence.csv> --run-index {Path(default_work_dir) / 'rnaseq' / 'rnaseq_run_index.tsv'} --genomes {Path(default_work_dir) / 'genomes'} --out {og_dir}",
                language="bash",
            )
        else:
            cov_path = og_dir / "orthogroup_expression_mapping_coverage.tsv"
            long_path = og_dir / "orthogroup_tpm_long.tsv"
            copy_path = og_dir / "X_orthogroup_copy_number.tsv"
            amr_path = og_dir / "orthogroup_amr_hits.tsv"

            if cov_path.exists():
                cov_df = pd.read_csv(cov_path, sep="\t")
                st.markdown("### Mapping coverage")
                c1, c2, c3 = st.columns(3)
                c1.metric("Runs evaluated", len(cov_df))
                c2.metric("Median mapping rate", round(float(cov_df.get("mapping_rate", pd.Series([0])).median()), 3))
                c3.metric("Mapped runs", int((cov_df.get("mapping_rate", pd.Series(dtype=float)) > 0).sum()))
                st.dataframe(cov_df, use_container_width=True)

            if long_path.exists():
                og_long = pd.read_csv(long_path, sep="\t")
                st.markdown("### Orthogroup expression summary")
                if not og_long.empty:
                    og_summary = (
                        og_long.groupby("orthogroup", as_index=False)
                        .agg(
                            biosamples=("biosample_accession", "nunique"),
                            runs=("run_id", "nunique"),
                            mean_tpm=("tpm", "mean"),
                            median_tpm=("tpm", "median"),
                            max_tpm=("tpm", "max"),
                        )
                        .sort_values(["max_tpm", "biosamples"], ascending=[False, False])
                    )
                    topn = st.slider("Top orthogroups", 10, 200, 50, 10, key="tab8_topn")
                    st.dataframe(og_summary.head(topn), use_container_width=True)
                    fig = px.bar(
                        og_summary.head(topn).sort_values("max_tpm", ascending=True),
                        x="max_tpm",
                        y="orthogroup",
                        orientation="h",
                        hover_data=["biosamples", "runs", "mean_tpm", "median_tpm"],
                        title="Top orthogroups by maximum TPM",
                    )
                    st.plotly_chart(fig, use_container_width=True)
                    df_download_button(og_summary, "Download orthogroup expression summary", "orthogroup_expression_summary.csv")

                    chosen_og = st.selectbox("Inspect one orthogroup", og_summary["orthogroup"].head(1000).tolist(), key="tab8_og_choice")
                    detail = og_long[og_long["orthogroup"] == chosen_og].copy().sort_values("tpm", ascending=False)
                    st.markdown("### Selected orthogroup detail")
                    st.dataframe(detail, use_container_width=True)
                    fig2 = px.strip(
                        detail,
                        x="biosample_accession",
                        y="tpm",
                        hover_data=["run_id", "orthogroup"],
                        title=f"TPM distribution for {chosen_og}",
                    )
                    st.plotly_chart(fig2, use_container_width=True)

            if copy_path.exists():
                copy_df = pd.read_csv(copy_path, sep="\t", index_col=0)
                st.markdown("### Orthogroup copy-number matrix")
                st.dataframe(copy_df.head(200), use_container_width=True)

            if amr_path.exists():
                og_amr = pd.read_csv(amr_path, sep="\t")
                st.markdown("### AMR-linked orthogroups")
                if og_amr.empty:
                    st.info("No AMR hits linked to orthogroups.")
                else:
                    st.dataframe(og_amr, use_container_width=True)
                    amr_sum = (
                        og_amr.groupby(["orthogroup", "amr_gene"], as_index=False)
                        .agg(biosamples=("biosample_accession", "nunique"), loci=("locus_tag", "nunique"))
                        .sort_values(["biosamples", "loci"], ascending=[False, False])
                    )
                    fig3 = px.bar(
                        amr_sum.head(50).sort_values("biosamples", ascending=True),
                        x="biosamples",
                        y="orthogroup",
                        color="amr_gene",
                        orientation="h",
                        hover_data=["loci"],
                        title="AMR-linked orthogroups by carrier biosamples",
                    )
                    st.plotly_chart(fig3, use_container_width=True)
                    df_download_button(amr_sum, "Download AMR-linked orthogroup summary", "orthogroup_amr_summary.csv")


    # ============================================================
    # Tab 9: Knowledge graph visualisation
    # ============================================================
    with tab9:
        st.subheader("Knowledge graph visualisation")
        st.markdown(
            """
    This tab visualises DOT files produced by Module F. For large graphs, use the top-N graph first unless you want the browser to behave like it has seen the face of God.
    """
        )

        default_kg_dir = str(Path(default_work_dir) / "kg")
        kg_dir_txt = st.text_input("Knowledge graph folder", value=default_kg_dir, key="tab9_kg_dir")
        kg_dir = Path(kg_dir_txt)

        dot_files = []
        if kg_dir.exists():
            dot_files = sorted(kg_dir.glob("**/*.dot"))

        if not dot_files:
            st.warning("No DOT files found yet. Run the knowledge graph module first.")
            st.code(
                f"onebio-atlas kg --manifest <manifest.tsv> --genomes {Path(default_work_dir) / 'genomes'} --features {Path(default_work_dir) / 'features'} --out {kg_dir} --top-n 10 --no-full",
                language="bash",
            )
        else:
            selected_dot = st.selectbox(
                "Choose graph",
                dot_files,
                format_func=lambda p: str(Path(p).relative_to(kg_dir)) if str(p).startswith(str(kg_dir)) else str(p),
                key="tab9_dot_choice",
            )
            dot_text = Path(selected_dot).read_text(encoding="utf-8", errors="replace")
            c1, c2 = st.columns([1, 1])
            c1.metric("DOT files", len(dot_files))
            c2.metric("Selected size KB", round(Path(selected_dot).stat().st_size / 1024, 1))
            st.download_button(
                "Download selected DOT",
                data=dot_text.encode("utf-8"),
                file_name=Path(selected_dot).name,
                mime="text/vnd.graphviz",
            )
            try:
                st.graphviz_chart(dot_text, use_container_width=True)
            except Exception as e:
                st.error(f"Graphviz rendering failed: {e}")
                st.code(dot_text[:8000], language="dot")
