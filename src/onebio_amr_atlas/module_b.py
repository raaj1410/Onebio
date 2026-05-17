from __future__ import annotations
from pathlib import Path
from typing import Optional, Tuple

import pandas as pd

from .utils import ensure_dir, log, PipelineError, which, run_cmd
from .net import ncbi_download_assembly_and_gff
from .gff import build_locus_table


def _normalise_cols(cols):
    out = []
    for c in cols:
        c2 = str(c).strip().lower()
        c2 = c2.replace("%", "pct_")
        for ch in [" ", "-", "/", "(", ")", "[", "]", ":"]:
            c2 = c2.replace(ch, "_")
        while "__" in c2:
            c2 = c2.replace("__", "_")
        out.append(c2.strip("_"))
    return out


def _pick_col(df: pd.DataFrame, candidates) -> Optional[str]:
    """Pick a column from df using a list of candidate names."""
    cols = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols:
            return cols[cand.lower()]
    norm = {n: c for n, c in zip(_normalise_cols(df.columns), df.columns)}
    for cand in candidates:
        cc = str(cand).strip().lower()
        cc = cc.replace("%", "pct_")
        for ch in [" ", "-", "/", "(", ")", "[", "]", ":"]:
            cc = cc.replace(ch, "_")
        while "__" in cc:
            cc = cc.replace("__", "_")
        cc = cc.strip("_")
        if cc in norm:
            return norm[cc]
    return None


def _run_amrfinder(genome_fna: Path, gff: Path, locus_map_tsv: Path, out_dir: Path,
                   amr_mode: str = "auto") -> Tuple[Optional[Path], Optional[Path]]:
    """Run AMRFinderPlus and create a locus_tag-mapped amr_hits.tsv.

    Returns (amrfinder_raw_tsv, amr_hits_tsv) or (None, None) if skipped in auto mode.
    """
    amrfinder_bin = which("amrfinder")
    if amrfinder_bin is None:
        msg = "AMRFinderPlus (amrfinder) not found in PATH."
        if amr_mode == "required":
            raise PipelineError(msg + " Install NCBI AMRFinderPlus or run with --amr skip.")
        log("WARNING: " + msg + " Skipping AMR layer (auto mode).")
        return None, None

    ensure_dir(out_dir)
    raw_path = out_dir / "amrfinder.tsv"
    hits_path = out_dir / "amr_hits.tsv"

    # AMRFinderPlus v4+ treats --gff as redundant unless a protein FASTA is also provided.
    # If we have a proteins file, use it; otherwise run on nucleotides only.
    # (We still map hits back to locus_tag using contig/start/stop columns + locus_map.tsv.)
    proteins = None
    for cand in ["proteins.faa", "protein.faa", "cds.faa", "genes.faa"]:
        p = out_dir / cand
        if p.exists() and p.stat().st_size > 0:
            proteins = p
            break

    if proteins is not None:
        cmd = [amrfinder_bin, "-n", str(genome_fna), "-p", str(proteins), "-g", str(gff), "--plus", "-o", str(raw_path)]
    else:
        cmd = [amrfinder_bin, "-n", str(genome_fna), "--plus", "-o", str(raw_path)]

    try:
        run_cmd(cmd, log_path=out_dir / "amrfinder.log", check=True)
    except Exception as e:
        # If the tool complains that --gff is redundant, rerun without it.
        msg = str(e).lower()
        if "gff is redundant" in msg:
            if proteins is not None:
                cmd2 = [amrfinder_bin, "-n", str(genome_fna), "-p", str(proteins), "--plus", "-o", str(raw_path)]
            else:
                cmd2 = [amrfinder_bin, "-n", str(genome_fna), "--plus", "-o", str(raw_path)]
            run_cmd(cmd2, log_path=out_dir / "amrfinder.log", check=True)
        else:
            raise

    if not raw_path.exists() or raw_path.stat().st_size == 0:
        log(f"WARNING: AMRFinderPlus produced empty output for {out_dir.name}")
        return raw_path, None

    amr = pd.read_csv(raw_path, sep="\t", dtype=str).fillna("")
    if amr.empty:
        hits_path.write_text(
            "locus_tag\tamr_gene\tdrug_class\tsubclass\telement_type\telement_subtype\tmethod\tidentity\tcoverage\taccession\n",
            encoding="utf-8",
        )
        return raw_path, hits_path

    locus = pd.read_csv(locus_map_tsv, sep="\t", dtype=str).fillna("")
    # coordinate key in locus table
    locus["start"] = pd.to_numeric(locus["start"], errors="coerce").fillna(-1).astype(int)
    locus["end"] = pd.to_numeric(locus["end"], errors="coerce").fillna(-1).astype(int)
    if "strand" not in locus.columns:
        locus["strand"] = "."
    locus["strand"] = locus["strand"].astype(str).where(locus["strand"].isin(["+", "-", "."]), ".")
    locus["_key"] = (
        locus["seqid"].astype(str)
        + "|" + locus[["start", "end"]].min(axis=1).astype(str)
        + "|" + locus[["start", "end"]].max(axis=1).astype(str)
        + "|" + locus["strand"].astype(str)
    )
    key_to_locus = locus.drop_duplicates(subset=["_key"]).set_index("_key")["locus_tag"].to_dict()
    locus_set = set(locus["locus_tag"].astype(str))

    # Build per-contig interval lists for overlap-based mapping (more robust than exact coordinate match).
    # We keep strand-specific lists but can fall back to strand-agnostic mapping.
    locus["_s"] = locus[["start", "end"]].min(axis=1)
    locus["_e"] = locus[["start", "end"]].max(axis=1)
    locus["_s"] = pd.to_numeric(locus["_s"], errors="coerce").fillna(-1).astype(int)
    locus["_e"] = pd.to_numeric(locus["_e"], errors="coerce").fillna(-1).astype(int)

    per_contig = {}
    per_contig_any = {}
    for contig, sub in locus.groupby("seqid", sort=False):
        sub2 = sub[["locus_tag","_s","_e","strand"]].copy()
        sub2 = sub2.sort_values("_s").reset_index(drop=True)
        per_contig_any[str(contig)] = sub2
        for st, sub3 in sub2.groupby("strand", sort=False):
            per_contig[(str(contig), str(st))] = sub3.reset_index(drop=True)

    def _map_by_overlap(contig: str, s: int, e: int, strand: str) -> str:
        """Return best-matching locus_tag by interval overlap; empty string if none."""
        if s < 0 or e < 0:
            return ""
        a = min(s, e); b = max(s, e)
        # first try same strand, then strand-agnostic
        cand = per_contig.get((str(contig), str(strand)))
        if cand is None or cand.empty:
            cand = per_contig_any.get(str(contig))
        if cand is None or cand.empty:
            return ""
        # simple scan over nearby genes (small n ~6k; hits ~100, fast enough)
        best_tag = ""
        best_ol = 0
        for _, r in cand.iterrows():
            gs, ge = int(r["_s"]), int(r["_e"])
            if gs > b + 50:  # break early; genes sorted by start
                break
            if ge < a - 50:
                continue
            ol = min(b, ge) - max(a, gs) + 1
            if ol <= 0:
                continue
            if ol > best_ol:
                best_ol = ol
                best_tag = str(r["locus_tag"])
        # require a reasonable overlap (>= 60% of hit length) to avoid wrong mapping
        hit_len = b - a + 1
        if best_ol >= int(0.60 * hit_len):
            return best_tag
        return ""

    # columns in AMRFinder output
    # (v4+ default headers include: Protein id, Contig id, Start, Stop, Strand, Element symbol, Element name, Type, Subtype, Class, Subclass, Method, % Coverage of reference, % Identity to reference, Closest reference accession, ...)
    c_prot = _pick_col(amr, ["protein id", "protein_id", "protein_identifier", "locus_tag", "locus tag", "id"])
    c_contig = _pick_col(amr, ["contig id", "contig_id", "contig", "sequence id", "sequence_id"])
    c_start = _pick_col(amr, ["start", "begin"])
    c_stop = _pick_col(amr, ["stop", "end"])
    c_strand = _pick_col(amr, ["strand"])

    # gene/determinant identifiers
    c_gene = _pick_col(amr, ["element symbol", "element_symbol", "gene_symbol", "gene symbol", "gene"])
    c_name = _pick_col(amr, ["element name", "element_name", "closest reference name"])
    c_eltype = _pick_col(amr, ["type", "element type", "element_type"])
    c_elsub = _pick_col(amr, ["subtype", "element subtype", "element_subtype"])

    # drug classes
    c_class = _pick_col(amr, ["class", "drug class", "drug_class"])
    c_subclass = _pick_col(amr, ["subclass"])

    c_method = _pick_col(amr, ["method"])
    c_ident = _pick_col(amr, [
        "pct_identity_to_reference",
        "pct_identity_to_reference_sequence",
        "% identity to reference",
        "% identity to reference sequence",
        "percent_identity",
    ])
    c_cov = _pick_col(amr, [
        "pct_coverage_of_reference",
        "pct_coverage_of_reference_sequence",
        "% coverage of reference",
        "% coverage of reference sequence",
        "percent_coverage",
    ])
    c_acc = _pick_col(amr, ["closest reference accession", "accession of closest sequence", "accession", "refseq_accession", "acc"])

    # resolve locus_tag
    locus_tag = None

    # 1) Direct ID mapping (rarely available in -n mode, but keep it)
    if c_prot is not None:
        s = amr[c_prot].astype(str).str.strip()
        if (s.isin(locus_set).mean()) >= 0.2:
            locus_tag = s

    # 2) Exact coordinate mapping: contig|min|max|strand
    idx = amr.index  # all rows

    if locus_tag is None and all(x is not None for x in [c_contig, c_start, c_stop]):
        st = pd.to_numeric(amr.loc[idx, c_start].astype(str).str.replace(r"[^0-9]", "", regex=True), errors="coerce")
        en = pd.to_numeric(amr.loc[idx, c_stop].astype(str).str.replace(r"[^0-9]", "", regex=True), errors="coerce")
        strand = amr.loc[idx, c_strand].astype(str) if c_strand is not None else "."
        strand = strand.where(strand.isin(["+", "-", "."]), ".")
        mn = pd.concat([st, en], axis=1).min(axis=1)
        mx = pd.concat([st, en], axis=1).max(axis=1)
        key = (
            amr.loc[idx, c_contig].astype(str).str.strip()
            + "|" + mn.fillna(-1).astype(int).astype(str)
            + "|" + mx.fillna(-1).astype(int).astype(str)
            + "|" + strand.astype(str)
        )
        locus_tag = key.map(key_to_locus).fillna("")

        # 3) If exact match maps poorly, fall back to overlap mapping (robust to off-by-one / feature-type differences)
        mapped_frac = (locus_tag.astype(str).str.strip() != "").mean()
        if mapped_frac < 0.50:
            # do overlap-based mapping per row
            lt2 = []
            for cont, a, b, stnd in zip(
                amr.loc[idx, c_contig].astype(str).str.strip(),
                mn.fillna(-1).astype(int).tolist(),
                mx.fillna(-1).astype(int).tolist(),
                strand.astype(str).tolist(),
            ):
                lt2.append(_map_by_overlap(cont, int(a), int(b), stnd))
            locus_tag = pd.Series(lt2, index=amr.index).fillna("")

    if locus_tag is None:
        log(f"WARNING: Could not map AMRFinder hits to locus_tag for {out_dir.name}. Writing hits with blank locus_tag.")
        locus_tag = pd.Series([""] * len(amr), index=amr.index)

    def _series_or_blank(df: pd.DataFrame, col: str | None) -> pd.Series:
        """Return df[col] as a string Series if present; otherwise a blank Series of same length."""
        if col is None:
            return pd.Series([""] * len(df), index=df.index)
        return df[col].astype(str)

    out = pd.DataFrame({
        "locus_tag": pd.Series(locus_tag, index=amr.index).astype(str),
        "amr_gene": _series_or_blank(amr, c_gene),
        "drug_class": _series_or_blank(amr, c_class),
        "subclass": _series_or_blank(amr, c_subclass),
        "element_type": _series_or_blank(amr, c_eltype),
        "element_subtype": _series_or_blank(amr, (c_name if c_name is not None else c_elsub)),
        "method": _series_or_blank(amr, c_method),
        "identity": _series_or_blank(amr, c_ident),
        "coverage": _series_or_blank(amr, c_cov),
        "accession": _series_or_blank(amr, c_acc),
    }).fillna("")

    # prefer subtype when gene symbol missing
    out.loc[out["amr_gene"].astype(str).str.strip() == "", "amr_gene"] = out["element_subtype"].astype(str)

    mapped = out["locus_tag"].astype(str).str.strip() != ""

    # Write unmapped hits for debugging (so you can fix mapping without losing information)
    if (~mapped).any():
        dbg = out.loc[~mapped].copy()
        idx = dbg.index
        if c_contig is not None:
            dbg["contig"] = amr.loc[idx, c_contig].astype(str).values
        if c_start is not None:
            dbg["start"] = amr.loc[idx, c_start].astype(str).values
        if c_stop is not None:
            dbg["stop"] = amr.loc[idx, c_stop].astype(str).values
        if c_strand is not None:
            dbg["strand"] = amr.loc[idx, c_strand].astype(str).values
        dbg_path = out_dir / "amr_unmapped.tsv"
        dbg.to_csv(dbg_path, sep="\t", index=False)
        log(f"AMR: unmapped hits written to {dbg_path} (n={int((~mapped).sum())})")

    if mapped.any():
        out = out.loc[mapped].copy()
    else:
        # keep empty with headers
        out = out.iloc[0:0].copy()
    out = out.drop_duplicates(subset=["locus_tag", "amr_gene", "accession", "method"]) \
             .sort_values(["locus_tag", "amr_gene"]) \
             .reset_index(drop=True)

    out.to_csv(hits_path, sep="\t", index=False)
    return raw_path, hits_path

def run_genome_processing(manifest_tsv: Path, out_dir: Path, threads: int = 8, prefer: str = "genbank",
                         amr: str = "auto") -> Path:
    ensure_dir(out_dir)
    man = pd.read_csv(manifest_tsv, sep="\t", dtype=str)

    outputs = []
    for _, r in man.iterrows():
        bs = str(r.get("biosample_accession","")).strip()
        asm = str(r.get("assembly_accession_best","")).strip()
        if not bs or bs.lower() == "nan" or not asm or asm.lower() == "nan":
            continue

        bs_dir = ensure_dir(out_dir / bs)
        ok_flag = bs_dir / ".done"
        if ok_flag.exists():
            genome_fna = bs_dir / "genome.fna"
            gff = bs_dir / "annotation.gff"
            locus_path = bs_dir / "locus_map.tsv"
            amr_raw = bs_dir / "amrfinder.tsv"
            amr_hits = bs_dir / "amr_hits.tsv"

            if amr != "skip" and genome_fna.exists() and gff.exists() and locus_path.exists() and (not amr_hits.exists()):
                try:
                    _run_amrfinder(genome_fna, gff, locus_path, bs_dir, amr_mode=amr)
                except Exception as e:
                    if amr == "required":
                        raise
                    log(f"WARNING: AMR layer failed for {bs}: {e}")

            outputs.append({
                "biosample_accession": bs,
                "assembly_accession": asm,
                "genome_fna": str(genome_fna),
                "gff": str(gff),
                "locus_map": str(locus_path),
                "amrfinder": str(amr_raw) if amr_raw.exists() else "",
                "amr_hits": str(amr_hits) if amr_hits.exists() else "",
            })
            continue

        log(f"Genome: {bs} assembly={asm}")
        genome_fna, gff = ncbi_download_assembly_and_gff(asm, bs_dir, prefer=prefer)

        locus_rows = build_locus_table(gff, prefer_type="CDS")
        if not locus_rows:
            raise PipelineError(f"No loci parsed from GFF for {bs}")
        locus_df = pd.DataFrame(locus_rows).drop_duplicates(subset=["locus_tag"])
        locus_path = bs_dir / "locus_map.tsv"
        locus_df.to_csv(locus_path, sep="\t", index=False)

        amr_raw, amr_hits = (None, None)
        if amr != "skip":
            try:
                amr_raw, amr_hits = _run_amrfinder(genome_fna, gff, locus_path, bs_dir, amr_mode=amr)
            except Exception as e:
                if amr == "required":
                    raise
                log(f"WARNING: AMR layer failed for {bs}: {e}")

        ok_flag.write_text("ok\n", encoding="utf-8")
        outputs.append({
            "biosample_accession": bs,
            "assembly_accession": asm,
            "genome_fna": str(genome_fna),
            "gff": str(gff),
            "locus_map": str(locus_path),
            "amrfinder": str(amr_raw) if amr_raw else "",
            "amr_hits": str(amr_hits) if amr_hits else "",
        })

    index_path = out_dir / "genome_index.tsv"
    pd.DataFrame(outputs).to_csv(index_path, sep="\t", index=False)
    log(f"Genome index: {index_path}")
    return index_path
