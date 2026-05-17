from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Iterable

import pandas as pd

from .utils import ensure_dir, log, PipelineError, require_tools
from .manifest_bridge import checker_csv_to_manifest
from .module_b import run_genome_processing
from .module_c import run_rnaseq_standard
from .module_d import build_features
from .module_e import build_atlas
from .module_f import build_kg_dot_all


def _norm_steps(steps: str | Iterable[str]) -> list[str]:
    if isinstance(steps, str):
        parts = [x.strip().lower() for x in steps.split(",")]
    else:
        parts = [str(x).strip().lower() for x in steps]
    return [x for x in parts if x]


def run_selected_modules(
    manifest: Path,
    work_dir: Path,
    steps: str | Iterable[str] = "genome,rnaseq,features,atlas,kg",
    threads: int = 8,
    max_runs_per_biosample: int = 1,
    fastq_root: Path | None = None,
    amr_mode: str = "auto",
    assembly_preference: str = "genbank",
    force_rnaseq: bool = False,
) -> dict[str, str]:
    """Run selected modules with stable folder names.

    Output layout:
      work_dir/genomes
      work_dir/rnaseq
      work_dir/features
      work_dir/atlas
      work_dir/kg
    """
    work_dir = ensure_dir(Path(work_dir))
    manifest = Path(manifest)
    steps_l = _norm_steps(steps)
    out: dict[str, str] = {"work_dir": str(work_dir), "manifest": str(manifest)}

    genomes = work_dir / "genomes"
    rnaseq = work_dir / "rnaseq"
    features = work_dir / "features"
    atlas = work_dir / "atlas"
    kg = work_dir / "kg"

    if "genome" in steps_l or "genomes" in steps_l or "all" in steps_l:
        log("Module B: genome retrieval + locus map + AMRFinderPlus layer")
        run_genome_processing(manifest, genomes, threads=threads, prefer=assembly_preference, amr=amr_mode)
        out["genomes"] = str(genomes)

    if "rnaseq" in steps_l or "rna-seq" in steps_l or "all" in steps_l:
        log("Module C: FASTQ fallback download + QC + alignment + counts + TPM")
        idx = run_rnaseq_standard(
            manifest,
            genomes,
            rnaseq,
            fastq_root=fastq_root,
            threads=threads,
            max_runs_per_biosample=max_runs_per_biosample,
            force=force_rnaseq,
        )
        out["rnaseq"] = str(rnaseq)
        out["rnaseq_run_index"] = str(idx)

    if "features" in steps_l or "feature" in steps_l or "all" in steps_l:
        log("Module D: BioSample-level genome/expression/fused matrices")
        build_features(manifest, genomes, rnaseq, features)
        out["features"] = str(features)

    if "atlas" in steps_l or "sqlite" in steps_l or "all" in steps_l:
        log("Module E: SQLite atlas")
        db = build_atlas(manifest, genomes, rnaseq, features, atlas)
        out["atlas"] = str(atlas)
        out["atlas_db"] = str(db)

    if "kg" in steps_l or "knowledge-graph" in steps_l or "knowledge_graph" in steps_l or "all" in steps_l:
        log("Module F: knowledge graph DOT files")
        build_kg_dot_all(manifest, genomes, features, kg, top_n=10, make_full=False)
        out["kg"] = str(kg)

    summary_path = work_dir / "pipeline_outputs.json"
    summary_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    out["summary_json"] = str(summary_path)
    return out


def run_checker_subprocess(
    species: str,
    out_root: Path,
    cache_root: Path,
    ncbi_ast_export: Path | None = None,
    disable_ebi: bool = False,
    disable_ncbi: bool = False,
    ebi_release: str = "auto",
) -> int:
    cmd = [
        "python", "-m", "onebio_amr_atlas.rnaseq_amr_checker",
        "--species", species,
        "--out-root", str(out_root),
        "--cache-root", str(cache_root),
        "--ebi-release", ebi_release,
    ]
    if ncbi_ast_export:
        cmd += ["--ncbi-ast-export", str(ncbi_ast_export)]
    if disable_ebi:
        cmd.append("--disable-ebi")
    if disable_ncbi:
        cmd.append("--disable-ncbi")
    return subprocess.call(cmd)


def write_demo_data(out_dir: Path) -> dict[str, str]:
    """Create tiny demo inputs for UI testing without downloading the universe."""
    out_dir = ensure_dir(Path(out_dir))
    demo = ensure_dir(out_dir / "demo")
    checker = demo / "demo_checker_output.csv"
    checker.write_text(
        "species,genome_id,genome_name,strain,taxon_id,biosample_accession,bioproject_accession,assembly_accession_best,final_conservative_hit,rnaseq_available,rnaseq_run_count,rnaseq_run_ids,ena_fastq_ftp,ena_fastq_md5,evidence_tier,antibiotics_resistant,antibiotics_susceptible,antibiotics_intermediate,amr_records_total,amr_sources\n"
        "Klebsiella pneumoniae,DEMO_KP_001,Demo isolate KP001,KP001,573,SAMN00000001,PRJNA000001,GCA_000000001.1,True,True,1,SRR00000001,,,Strong same-BioSample RNA-seq + AST + genome,ampicillin;ciprofloxacin,amikacin,cefepime,3,BV-BRC;NCBI_AST_BROWSER\n",
        encoding="utf-8",
    )
    manifest_dir = ensure_dir(demo / "manifest")
    manifest, ast, selected, qc = checker_csv_to_manifest(checker, manifest_dir, only_final_hits=True, only_rnaseq_available=True)

    # Minimal toy genome assets for local smoke-testing of file layout.
    toy_bs = ensure_dir(demo / "toy_genomes" / "SAMN00000001")
    (toy_bs / "genome.fna").write_text(
        ">contig1\n"
        "ATGAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTTAA\n",
        encoding="utf-8",
    )
    (toy_bs / "annotation.gff").write_text(
        "##gff-version 3\n"
        "contig1\tDemo\tCDS\t1\t57\t.\t+\t0\tID=cds-demo001;locus_tag=demo001;gene=blaDEM;product=demo beta-lactamase\n",
        encoding="utf-8",
    )
    (toy_bs / "locus_map.tsv").write_text(
        "locus_tag\tgene\tproduct\tseqid\tstart\tend\tstrand\tlength\n"
        "demo001\tblaDEM\tdemo beta-lactamase\tcontig1\t1\t57\t+\t57\n",
        encoding="utf-8",
    )
    (toy_bs / "amr_hits.tsv").write_text(
        "locus_tag\tamr_gene\tdrug_class\tsubclass\telement_type\telement_subtype\tmethod\tidentity\tcoverage\taccession\n"
        "demo001\tblaDEM\tbeta-lactam\tpenicillin\tAMR\tbeta-lactamase\tdemo\t100\t100\tDEMO001\n",
        encoding="utf-8",
    )

    # Local FASTQ so Module C can be smoke-tested without touching SRA/ENA.
    import gzip
    fq_dir = ensure_dir(demo / "local_fastq" / "SRR00000001")
    with gzip.open(fq_dir / "SRR00000001.fastq.gz", "wt", encoding="utf-8") as fh:
        fh.write("@demo_read_1\nATGAAACCCGGGTTTAAACCCGGGTTTAAA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n")

    return {
        "checker_csv": str(checker),
        "manifest": str(manifest),
        "ast_long": str(ast),
        "selected_rows": str(selected),
        "manifest_qc": str(qc),
        "toy_genomes": str(demo / "toy_genomes"),
        "local_fastq": str(demo / "local_fastq"),
    }


def run_orthogroup_analysis(
    panaroo_gene_presence_absence: Path,
    run_index: Path,
    genomes: Path,
    out_dir: Path,
    aggregate: str = "sum",
) -> dict[str, str]:
    """Build orthogroup-level TPM/copy-number/AMR outputs from Panaroo results."""
    from .build_orthogroup_expression import (
        build_locus_to_orthogroup,
        aggregate_run_tpm_to_orthogroup,
        collapse_runs_to_biosample,
        build_copy_number_matrix,
        build_orthogroup_amr_summary,
    )

    out_dir = ensure_dir(Path(out_dir))
    mapping_df, pan_cov, unmatched = build_locus_to_orthogroup(
        Path(panaroo_gene_presence_absence), Path(genomes), out_dir
    )
    long_df, expr_cov = aggregate_run_tpm_to_orthogroup(
        Path(run_index), mapping_df, aggregate=aggregate
    )
    long_df.to_csv(out_dir / "orthogroup_tpm_long.tsv", sep="\t", index=False)
    expr_cov.to_csv(out_dir / "orthogroup_expression_mapping_coverage.tsv", sep="\t", index=False)
    collapse_runs_to_biosample(long_df).to_csv(out_dir / "X_orthogroup_tpm.tsv", sep="\t")
    build_copy_number_matrix(mapping_df).to_csv(out_dir / "X_orthogroup_copy_number.tsv", sep="\t")
    og_amr = build_orthogroup_amr_summary(Path(genomes), mapping_df)
    og_amr.to_csv(out_dir / "orthogroup_amr_hits.tsv", sep="\t", index=False)

    summary = {
        "out_dir": str(out_dir),
        "locus_to_orthogroup": str(out_dir / "locus_to_orthogroup.tsv"),
        "panaroo_locus_match_coverage": str(out_dir / "panaroo_locus_match_coverage.tsv"),
        "panaroo_unmatched_examples": str(out_dir / "panaroo_unmatched_examples.tsv"),
        "orthogroup_tpm_long": str(out_dir / "orthogroup_tpm_long.tsv"),
        "orthogroup_expression_mapping_coverage": str(out_dir / "orthogroup_expression_mapping_coverage.tsv"),
        "X_orthogroup_tpm": str(out_dir / "X_orthogroup_tpm.tsv"),
        "X_orthogroup_copy_number": str(out_dir / "X_orthogroup_copy_number.tsv"),
        "orthogroup_amr_hits": str(out_dir / "orthogroup_amr_hits.tsv"),
    }
    (out_dir / "orthogroup_outputs.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def run_same_gene_locus_catalog(
    db: Path,
    out_dir: Path,
    expr_threshold: float = 1.0,
    min_loci: int = 2,
    min_biosamples: int = 2,
    min_delta_expr: float = 5.0,
) -> dict[str, str]:
    """Run the SQLite-first same-gene/different-loci catalog builder."""
    out_dir = ensure_dir(Path(out_dir))
    cmd = [
        "python", "-m", "onebio_amr_atlas.same_gene_locus_catalog",
        "--db", str(db),
        "--outdir", str(out_dir),
        "--expr-threshold", str(expr_threshold),
        "--min-loci", str(min_loci),
        "--min-biosamples", str(min_biosamples),
        "--min-delta-expr", str(min_delta_expr),
    ]
    subprocess.check_call(cmd)
    summary = {"out_dir": str(out_dir)}
    for p in sorted(out_dir.glob("*")):
        if p.is_file():
            summary[p.stem] = str(p)
    (out_dir / "same_gene_outputs.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary
