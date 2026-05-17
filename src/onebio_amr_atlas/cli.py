from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

from .manifest_bridge import checker_csv_to_manifest
from .pipeline import run_selected_modules, run_checker_subprocess, write_demo_data, run_orthogroup_analysis, run_same_gene_locus_catalog
from .module_a import build_manifest
from .module_b import run_genome_processing
from .module_c import run_rnaseq_standard
from .module_d import build_features
from .module_e import build_atlas, query_gene
from .module_f import build_kg_dot_all
from .doctor import run_doctor
from .utils import log, PipelineError


def main() -> None:
    p = argparse.ArgumentParser(
        prog="onebio-atlas",
        description="Single UI/CLI for RNA-seq + genome + AST + AMR atlas construction",
    )
    sp = p.add_subparsers(dest="cmd", required=True)

    p_ui = sp.add_parser("ui", help="Launch browser UI")
    p_ui.add_argument("--port", type=int, default=8501)
    p_ui.add_argument("--host", default="localhost")

    p_demo = sp.add_parser("demo-data", help="Write tiny demo inputs")
    p_demo.add_argument("--out", required=True, type=Path)

    p_disc = sp.add_parser("discover", help="Run the multi-source RNA-seq AMR checker")
    p_disc.add_argument("--species", required=True, help="Comma-separated species list")
    p_disc.add_argument("--out-root", required=True, type=Path)
    p_disc.add_argument("--cache-root", required=True, type=Path)
    p_disc.add_argument("--ncbi-ast-export", type=Path)
    p_disc.add_argument("--disable-ebi", action="store_true")
    p_disc.add_argument("--disable-ncbi", action="store_true")
    p_disc.add_argument("--ebi-release", default="auto")

    p_prep = sp.add_parser("prepare-manifest", help="Convert checker CSV/XLSX into manifest.tsv + ast_long.tsv")
    p_prep.add_argument("--checker-table", required=True, type=Path)
    p_prep.add_argument("--out", required=True, type=Path)
    p_prep.add_argument("--include-non-final", action="store_true", help="Do not filter final_conservative_hit")
    p_prep.add_argument("--include-no-rnaseq", action="store_true", help="Do not filter rnaseq_available")
    p_prep.add_argument("--max-rows", type=int, default=0)

    p_run = sp.add_parser("run", help="Run selected pipeline modules")
    p_run.add_argument("--manifest", required=True, type=Path)
    p_run.add_argument("--work-dir", required=True, type=Path)
    p_run.add_argument("--steps", default="genome,rnaseq,features,atlas,kg", help="Comma list: genome,rnaseq,features,atlas,kg or all")
    p_run.add_argument("--threads", type=int, default=8)
    p_run.add_argument("--max-runs-per-biosample", type=int, default=1)
    p_run.add_argument("--fastq-root", type=Path)
    p_run.add_argument("--amr", choices=["auto", "skip", "required"], default="auto")
    p_run.add_argument("--prefer", choices=["genbank", "refseq"], default="genbank")
    p_run.add_argument("--force-rnaseq", action="store_true")

    p_og = sp.add_parser("orthogroup", help="Build orthogroup TPM/copy-number/AMR outputs from Panaroo gene_presence_absence.csv")
    p_og.add_argument("--panaroo", required=True, type=Path, help="Panaroo gene_presence_absence.csv")
    p_og.add_argument("--run-index", required=True, type=Path, help="rnaseq_run_index.tsv")
    p_og.add_argument("--genomes", required=True, type=Path, help="work/run1/genomes folder")
    p_og.add_argument("--out", required=True, type=Path)
    p_og.add_argument("--aggregate", choices=["sum", "max", "mean"], default="sum")

    p_sg = sp.add_parser("same-gene", help="Build same-gene/different-loci SQLite catalog outputs")
    p_sg.add_argument("--db", required=True, type=Path, help="atlas.sqlite")
    p_sg.add_argument("--out", required=True, type=Path)
    p_sg.add_argument("--expr-threshold", type=float, default=1.0)
    p_sg.add_argument("--min-loci", type=int, default=2)
    p_sg.add_argument("--min-biosamples", type=int, default=2)
    p_sg.add_argument("--min-delta-expr", type=float, default=5.0)

    # Keep the older module-level commands because experienced users will want them.
    pA = sp.add_parser("build-manifest", help="Legacy Module A: XLSX to strict manifest folder")
    pA.add_argument("--xlsx", required=True, type=Path)
    pA.add_argument("--out", required=True, type=Path)

    pB = sp.add_parser("genome", help="Module B: Download assembly+GFF and build locus maps")
    pB.add_argument("--manifest", required=True, type=Path)
    pB.add_argument("--out", required=True, type=Path)
    pB.add_argument("--threads", type=int, default=8)
    pB.add_argument("--prefer", choices=["genbank", "refseq"], default="genbank")
    pB.add_argument("--amr", choices=["auto", "skip", "required"], default="auto")

    pC = sp.add_parser("rnaseq", help="Module C: Download FASTQ, QC, map, count, TPM")
    pC.add_argument("--manifest", required=True, type=Path)
    pC.add_argument("--genomes", required=True, type=Path)
    pC.add_argument("--out", required=True, type=Path)
    pC.add_argument("--threads", type=int, default=8)
    pC.add_argument("--max-runs-per-biosample", type=int, default=999)
    pC.add_argument("--fastq-root", type=Path)
    pC.add_argument("--force", action="store_true")

    pD = sp.add_parser("features", help="Module D: Build feature matrices")
    pD.add_argument("--manifest", required=True, type=Path)
    pD.add_argument("--genomes", required=True, type=Path)
    pD.add_argument("--rnaseq", required=True, type=Path)
    pD.add_argument("--out", required=True, type=Path)
    pD.add_argument("--expr-threshold-tpm", type=float, default=1.0)

    pE = sp.add_parser("atlas", help="Module E: Build SQLite atlas")
    pE.add_argument("--manifest", required=True, type=Path)
    pE.add_argument("--genomes", required=True, type=Path)
    pE.add_argument("--rnaseq", required=True, type=Path)
    pE.add_argument("--features", required=True, type=Path)
    pE.add_argument("--out", required=True, type=Path)

    pF = sp.add_parser("kg", help="Module F: Generate knowledge graph DOT files")
    pF.add_argument("--manifest", required=True, type=Path)
    pF.add_argument("--genomes", required=True, type=Path)
    pF.add_argument("--features", required=True, type=Path)
    pF.add_argument("--out", required=True, type=Path)
    pF.add_argument("--top-n", type=int, default=10)
    pF.add_argument("--no-full", action="store_true")

    pQ = sp.add_parser("query-gene", help="Query gene expression from atlas.sqlite")
    pQ.add_argument("--db", required=True, type=Path)
    pQ.add_argument("--gene-like", required=True)
    pQ.add_argument("--limit", type=int, default=20)

    pDoc = sp.add_parser("doctor", help="Diagnose low mapping for one reference + FASTQs")
    pDoc.add_argument("--ref", required=True, type=Path)
    pDoc.add_argument("--fastq", required=True, type=Path, nargs="+")
    pDoc.add_argument("--out", required=True, type=Path)
    pDoc.add_argument("--threads", type=int, default=8)
    pDoc.add_argument("--full", action="store_true")

    args = p.parse_args()
    try:
        if args.cmd == "ui":
            subprocess.call(["streamlit", "run", str(Path(__file__).with_name("ui.py")), "--server.port", str(args.port), "--server.address", args.host])
        elif args.cmd == "demo-data":
            print(json.dumps(write_demo_data(args.out), indent=2))
        elif args.cmd == "discover":
            raise SystemExit(run_checker_subprocess(args.species, args.out_root, args.cache_root, args.ncbi_ast_export, args.disable_ebi, args.disable_ncbi, args.ebi_release))
        elif args.cmd == "prepare-manifest":
            paths = checker_csv_to_manifest(
                args.checker_table,
                args.out,
                only_final_hits=not args.include_non_final,
                only_rnaseq_available=not args.include_no_rnaseq,
                max_rows=(args.max_rows if args.max_rows > 0 else None),
            )
            print("\n".join(str(x) for x in paths))
        elif args.cmd == "run":
            out = run_selected_modules(
                args.manifest,
                args.work_dir,
                steps=args.steps,
                threads=args.threads,
                max_runs_per_biosample=args.max_runs_per_biosample,
                fastq_root=args.fastq_root,
                amr_mode=args.amr,
                assembly_preference=args.prefer,
                force_rnaseq=args.force_rnaseq,
            )
            print(json.dumps(out, indent=2))
        elif args.cmd == "orthogroup":
            out = run_orthogroup_analysis(args.panaroo, args.run_index, args.genomes, args.out, aggregate=args.aggregate)
            print(json.dumps(out, indent=2))
        elif args.cmd == "same-gene":
            out = run_same_gene_locus_catalog(
                args.db, args.out,
                expr_threshold=args.expr_threshold,
                min_loci=args.min_loci,
                min_biosamples=args.min_biosamples,
                min_delta_expr=args.min_delta_expr,
            )
            print(json.dumps(out, indent=2))
        elif args.cmd == "build-manifest":
            build_manifest(args.xlsx, args.out)
        elif args.cmd == "genome":
            run_genome_processing(args.manifest, args.out, threads=args.threads, prefer=args.prefer, amr=args.amr)
        elif args.cmd == "rnaseq":
            run_rnaseq_standard(args.manifest, args.genomes, args.out, fastq_root=args.fastq_root, threads=args.threads, max_runs_per_biosample=args.max_runs_per_biosample, force=args.force)
        elif args.cmd == "features":
            build_features(args.manifest, args.genomes, args.rnaseq, args.out, expr_threshold_tpm=args.expr_threshold_tpm)
        elif args.cmd == "atlas":
            build_atlas(args.manifest, args.genomes, args.rnaseq, args.features, args.out)
        elif args.cmd == "kg":
            build_kg_dot_all(args.manifest, args.genomes, args.features, args.out, top_n=args.top_n, make_full=(not args.no_full))
        elif args.cmd == "query-gene":
            df = query_gene(args.db, args.gene_like, limit=args.limit)
            print(df.to_string(index=False))
        elif args.cmd == "doctor":
            run_doctor(args.ref, list(args.fastq), args.out, threads=args.threads, full=args.full)
    except PipelineError as e:
        log(f"ERROR: {e}")
        raise SystemExit(2)


def ui_main() -> None:
    subprocess.call([
        "streamlit", "run", str(Path(__file__).with_name("ui.py")),
        "--server.port", "8501",
        "--server.address", "localhost",
    ])


if __name__ == "__main__":
    main()
