from __future__ import annotations

import json
import os
import sqlite3
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st

from onebio_amr_atlas.manifest_bridge import checker_csv_to_manifest
from onebio_amr_atlas.pipeline import run_selected_modules, write_demo_data, run_orthogroup_analysis, run_same_gene_locus_catalog
from onebio_amr_atlas.atlas_explorer import render_atlas_explorer

st.set_page_config(page_title="OneBio AMR Atlas", layout="wide")


def run_live(cmd: list[str], cwd: str | None = None) -> int:
    box = st.empty()
    log_lines: list[str] = []
    p = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
    assert p.stdout is not None
    for line in p.stdout:
        log_lines.append(line.rstrip())
        box.code("\n".join(log_lines[-120:]), language="bash")
    return p.wait()


def read_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".tsv", ".tab"}:
        return pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return pd.read_csv(path, dtype=str).fillna("")


st.title("OneBio AMR Atlas")
st.caption("Genome + RNA-seq + AST + AMR integration without forcing users to memorise twelve scripts. A rare act of mercy.")

with st.sidebar:
    st.header("Workflow")
    mode = st.selectbox(
        "Choose task",
        [
            "1. Discover candidate isolates",
            "2. Build module manifest",
            "3. Run pipeline modules",
            "4. Explore SQLite atlas",
            "5. Orthogroup and locus analyses",
            "6. Demo data",
        ],
    )
    st.divider()
    st.write("Recommended flow: Discovery → Manifest → Modules → Atlas explorer")

if mode.startswith("1"):
    st.subheader("Discover strains with genome + RNA-seq + AST phenotype")
    st.write("This wraps the RNA-seq AMR checker. Keep the species list small for testing unless you enjoy watching progress bars age.")
    species = st.text_area("Species list", "Klebsiella pneumoniae", height=100)
    out_root = Path(st.text_input("Output folder", "outputs/checker"))
    cache_root = Path(st.text_input("Cache folder", "cache/checker"))
    ncbi_ast = st.text_input("Optional NCBI AST Browser CSV/TSV path", "")
    c1, c2, c3 = st.columns(3)
    disable_ebi = c1.checkbox("Disable EBI AMR Portal", value=False)
    disable_ncbi = c2.checkbox("Disable NCBI AST", value=False)
    ebi_release = c3.text_input("EBI release", "auto")
    if st.button("Run discovery", type="primary"):
        cmd = [
            "python", "-m", "onebio_amr_atlas.rnaseq_amr_checker",
            "--species", species,
            "--out-root", str(out_root),
            "--cache-root", str(cache_root),
            "--ebi-release", ebi_release,
        ]
        if ncbi_ast.strip():
            cmd += ["--ncbi-ast-export", ncbi_ast.strip()]
        if disable_ebi:
            cmd.append("--disable-ebi")
        if disable_ncbi:
            cmd.append("--disable-ncbi")
        rc = run_live(cmd)
        st.success("Discovery finished" if rc == 0 else f"Discovery exited with code {rc}")

elif mode.startswith("2"):
    st.subheader("Build manifest for downstream modules")
    st.write("Use the combined rnaseq_true CSV or a species-level checker CSV. Output becomes the clean input to genome/RNA-seq/atlas modules.")
    input_mode = st.radio("Input method", ["File path", "Upload file"], horizontal=True)
    checker_path: Path | None = None
    if input_mode == "File path":
        txt = st.text_input("Checker table path", "outputs/checker/rnaseq_true_combined.csv")
        checker_path = Path(txt) if txt.strip() else None
    else:
        up = st.file_uploader("Upload checker CSV/TSV/XLSX", type=["csv", "tsv", "xlsx", "xls"])
        if up is not None:
            suffix = Path(up.name).suffix
            tmp = Path(tempfile.mkdtemp()) / f"checker{suffix}"
            tmp.write_bytes(up.getvalue())
            checker_path = tmp
    out_dir = Path(st.text_input("Manifest output folder", "work/manifest"))
    c1, c2, c3 = st.columns(3)
    only_final = c1.checkbox("Keep only final conservative hits", value=True)
    only_rnaseq = c2.checkbox("Keep only RNA-seq available", value=True)
    max_rows = c3.number_input("Limit rows for test run", min_value=0, value=0, step=1)
    if st.button("Create manifest", type="primary") and checker_path:
        manifest, ast, selected, qc = checker_csv_to_manifest(
            checker_path,
            out_dir,
            only_final_hits=only_final,
            only_rnaseq_available=only_rnaseq,
            max_rows=(int(max_rows) if int(max_rows) > 0 else None),
        )
        st.success("Manifest created")
        st.json(json.loads(qc.read_text()))
        st.write("Manifest:", str(manifest))
        st.dataframe(pd.read_csv(manifest, sep="\t", dtype=str).fillna(""), use_container_width=True)
        st.write("AST long:", str(ast))

elif mode.startswith("3"):
    st.subheader("Run selected modules")
    manifest = Path(st.text_input("manifest.tsv path", "work/manifest/manifest.tsv"))
    work_dir = Path(st.text_input("Work folder", "work/run1"))
    workflow = st.selectbox(
        "What should be done?",
        [
            "Full pipeline: genome,rnaseq,features,atlas,kg",
            "Genome only",
            "RNA-seq only",
            "Features only",
            "Atlas only",
            "Knowledge graph only",
            "Custom",
        ],
    )
    custom = st.text_input("Custom steps", "genome,rnaseq,features,atlas,kg")
    step_map = {
        "Full pipeline: genome,rnaseq,features,atlas,kg": "genome,rnaseq,features,atlas,kg",
        "Genome only": "genome",
        "RNA-seq only": "rnaseq",
        "Features only": "features",
        "Atlas only": "atlas",
        "Knowledge graph only": "kg",
        "Custom": custom,
    }
    c1, c2, c3, c4 = st.columns(4)
    threads = c1.number_input("Threads", min_value=1, max_value=64, value=8)
    max_runs = c2.number_input("Max runs/BioSample", min_value=1, max_value=999, value=1)
    amr = c3.selectbox("AMRFinderPlus mode", ["auto", "skip", "required"])
    prefer = c4.selectbox("Assembly preference", ["genbank", "refseq"])
    fastq_root_txt = st.text_input("Optional local FASTQ root", "")
    force = st.checkbox("Force RNA-seq rerun", value=False)
    if st.button("Run modules", type="primary"):
        selected_steps = [x.strip() for x in step_map[workflow].split(",") if x.strip()]
        if not selected_steps:
            st.error("No steps selected.")
        else:
            overall = st.progress(0, text="Starting selected modules")
            status_box = st.empty()
            outputs: dict[str, str] = {}
            for i, step in enumerate(selected_steps, start=1):
                label = f"Step {i}/{len(selected_steps)}: {step}"
                status_box.info(label)
                step_bar = st.progress(0, text=f"Running {step}")
                try:
                    step_out = run_selected_modules(
                        manifest=manifest,
                        work_dir=work_dir,
                        steps=step,
                        threads=int(threads),
                        max_runs_per_biosample=int(max_runs),
                        fastq_root=(Path(fastq_root_txt) if fastq_root_txt.strip() else None),
                        amr_mode=amr,
                        assembly_preference=prefer,
                        force_rnaseq=force,
                    )
                    outputs.update(step_out)
                    step_bar.progress(100, text=f"Completed {step}")
                    overall.progress(int(i / len(selected_steps) * 100), text=f"Completed {step}")
                except Exception as e:
                    step_bar.progress(100, text=f"Failed at {step}")
                    st.error(f"Pipeline stopped at step `{step}`: {e}")
                    st.info("Rerun with the same work folder. Completed steps are reused where the module has .done files or existing outputs.")
                    break
            else:
                status_box.success("Selected modules completed")
                st.json(outputs)

elif mode.startswith("4"):
    render_atlas_explorer(initial_db_path="work/run1/atlas/atlas.sqlite", default_work_dir="work/run1")

elif mode.startswith("5"):
    st.subheader("Orthogroup and locus-level analyses")
    analysis = st.selectbox(
        "Choose analysis",
        [
            "Prepare Panaroo inputs",
            "Build orthogroup expression from Panaroo",
            "Build same-gene different-loci catalog",
        ],
    )

    if analysis == "Prepare Panaroo inputs":
        genomes = Path(st.text_input("Genomes folder", "work/run1/genomes", key="pan_genomes"))
        out_dir = Path(st.text_input("Panaroo prep output folder", "work/run1/panaroo", key="pan_out"))
        st.markdown("This creates standardised `.gff` + `.fa` files and a `panaroo_pairs.tsv`. You still run Panaroo separately because Panaroo is not a small toy dependency.")
        if st.button("Prepare Panaroo inputs", type="primary"):
            from onebio_amr_atlas.prepare_panaroo_inputs import safe_link_or_copy
            input_dir = out_dir / "inputs"
            input_dir.mkdir(parents=True, exist_ok=True)
            pair_lines = []
            missing_lines = ["biosample_accession\treason"]
            for bs_dir in sorted(genomes.iterdir() if genomes.exists() else []):
                if not bs_dir.is_dir():
                    continue
                bs = bs_dir.name.strip()
                gff = bs_dir / "annotation.gff"
                fasta = bs_dir / "genome.fna"
                if not gff.exists() or not fasta.exists():
                    missing_lines.append(f"{bs}\tmissing annotation.gff or genome.fna")
                    continue
                out_gff = input_dir / f"{bs}.gff"
                out_fa = input_dir / f"{bs}.fa"
                safe_link_or_copy(gff, out_gff)
                safe_link_or_copy(fasta, out_fa)
                pair_lines.append(f"{out_gff}\t{out_fa}")
            out_dir.mkdir(parents=True, exist_ok=True)
            (out_dir / "panaroo_pairs.tsv").write_text("\n".join(pair_lines) + ("\n" if pair_lines else ""), encoding="utf-8")
            (out_dir / "missing_inputs.tsv").write_text("\n".join(missing_lines) + "\n", encoding="utf-8")
            st.success("Panaroo inputs prepared")
            st.code(f"panaroo -i {input_dir}/*.gff -o {out_dir}/panaroo_out --clean-mode strict --remove-invalid-genes -t 8", language="bash")
            st.write("Pair list:", str(out_dir / "panaroo_pairs.tsv"))

    elif analysis == "Build orthogroup expression from Panaroo":
        panaroo = Path(st.text_input("Panaroo gene_presence_absence.csv", "work/run1/panaroo/panaroo_out/gene_presence_absence.csv"))
        run_index = Path(st.text_input("RNA-seq run index", "work/run1/rnaseq/rnaseq_run_index.tsv"))
        genomes = Path(st.text_input("Genomes folder", "work/run1/genomes"))
        out_dir = Path(st.text_input("Orthogroup output folder", "work/run1/orthogroups"))
        aggregate = st.selectbox("Aggregate multiple loci", ["sum", "max", "mean"], index=0)
        if st.button("Build orthogroup outputs", type="primary"):
            try:
                with st.status("Building orthogroup expression outputs", expanded=True) as status:
                    st.write("Mapping Panaroo tokens back to locus tags")
                    out = run_orthogroup_analysis(panaroo, run_index, genomes, out_dir, aggregate=aggregate)
                    status.update(label="Orthogroup outputs completed", state="complete")
                st.json(out)
            except Exception as e:
                st.error(str(e))

    else:
        db = Path(st.text_input("atlas.sqlite path", "work/run1/atlas/atlas.sqlite"))
        out_dir = Path(st.text_input("Output folder", "work/run1/same_gene_loci"))
        c1, c2, c3, c4 = st.columns(4)
        expr_threshold = c1.number_input("Expression threshold TPM", value=1.0)
        min_loci = c2.number_input("Minimum loci", min_value=1, value=2)
        min_biosamples = c3.number_input("Minimum biosamples", min_value=1, value=2)
        min_delta = c4.number_input("Minimum expression spread", value=5.0)
        if st.button("Build same-gene locus catalog", type="primary"):
            try:
                out = run_same_gene_locus_catalog(db, out_dir, float(expr_threshold), int(min_loci), int(min_biosamples), float(min_delta))
                st.success("Catalog completed")
                st.json(out)
            except Exception as e:
                st.error(str(e))

else:
    st.subheader("Generate demo files")
    out = Path(st.text_input("Demo output folder", "onebio_demo"))
    if st.button("Write demo data", type="primary"):
        paths = write_demo_data(out)
        st.success("Demo files written")
        st.json(paths)
