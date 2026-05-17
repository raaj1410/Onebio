from __future__ import annotations

from pathlib import Path
from typing import Dict, Any, List, Tuple
import gzip
import io
import json
import re
import shutil
import subprocess
import pandas as pd

from .utils import ensure_dir, log, PipelineError, require_tools, run_cmd, safe_split_semi
from .io_qc import fastq_sample_stats
from .gff import detect_featurecounts_keys
from .net import download_fastqs_resilient


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def _tool_exists(name: str) -> bool:
    return shutil.which(name) is not None


def _parse_flagstat_text(txt: str) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for line in txt.splitlines():
        line = line.strip()
        if "in total" in line:
            m = re.search(r"^(\d+)\s+\+\s+(\d+)\s+in total", line)
            if m:
                out["total"] = int(m.group(1))
        if " mapped (" in line and "secondary" not in line:
            m = re.search(r"^(\d+)\s+\+\s+(\d+)\s+mapped\s+\(([-\d\.]+)%", line)
            if m:
                out["mapped"] = int(m.group(1))
                out["mapped_pct"] = float(m.group(3))
    return out


def _parse_bowtie2_stderr(txt: str) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    m = re.search(r"([0-9.]+)%\s+overall alignment rate", txt)
    if m:
        out["overall_alignment_rate"] = float(m.group(1))
    return out


def _detect_run_technology(platform_hint: str = "", model_hint: str = "", read_qc: Dict[str, Any] | None = None) -> str:
    text = f"{platform_hint} {model_hint}".lower()
    if any(x in text for x in ["oxford", "nanopore", "ont", "promethion", "minion", "gridion"]):
        return "ont"
    if any(x in text for x in ["illumina", "nextseq", "novaseq", "hiseq", "miseq"]):
        return "illumina"

    read_qc = read_qc or {}
    mean_len = read_qc.get("len_mean")
    if mean_len is not None:
        try:
            mean_len = float(mean_len)
            if mean_len >= 300:
                return "ont"
        except Exception:
            pass
    return "illumina"


def _classify_fastqs(fastqs: List[Path]) -> Dict[str, Any]:
    """
    Classify FASTQ inputs into paired or single groups.
    Supports one file, one pair, or multiple lane files.
    """
    fastqs = sorted([Path(x) for x in fastqs])
    out: Dict[str, Any] = {"layout": "single", "r1": [], "r2": [], "single": []}

    if len(fastqs) == 1:
        out["single"] = [str(fastqs[0])]
        return out

    r1, r2, single = [], [], []
    for fq in fastqs:
        name = fq.name
        if re.search(r"(_R?1([_\.\-]|$)|\.1\.f(ast)?q)", name, flags=re.IGNORECASE):
            r1.append(str(fq))
        elif re.search(r"(_R?2([_\.\-]|$)|\.2\.f(ast)?q)", name, flags=re.IGNORECASE):
            r2.append(str(fq))
        else:
            single.append(str(fq))

    if r1 and r2 and len(r1) == len(r2) and not single:
        out["layout"] = "paired"
        out["r1"] = sorted(r1)
        out["r2"] = sorted(r2)
        return out

    # fallback: if exactly 2 files and neither could be classified, treat as a pair
    if len(fastqs) == 2 and len(single) == 2:
        out["layout"] = "paired"
        out["r1"] = [str(fastqs[0])]
        out["r2"] = [str(fastqs[1])]
        return out

    out["single"] = [str(x) for x in fastqs]
    return out


def _merge_fastqs(inputs: List[str], out_fq_gz: Path) -> Path:
    """
    Merge multiple FASTQ or FASTQ.GZ files into one gzipped FASTQ.
    This is mainly for fastp, which takes one R1 and one R2 input path.
    """
    out_fq_gz = Path(out_fq_gz)
    ensure_dir(out_fq_gz.parent)
    if len(inputs) == 1:
        return Path(inputs[0])

    with gzip.open(out_fq_gz, "wb") as oh:
        for p in inputs:
            p = str(p)
            if p.endswith(".gz"):
                with gzip.open(p, "rb") as ih:
                    shutil.copyfileobj(ih, oh)
            else:
                with open(p, "rb") as ih:
                    shutil.copyfileobj(ih, oh)
    return out_fq_gz


def _prepare_fastp_inputs(fastqs: List[Path], work_dir: Path) -> Dict[str, Any]:
    """
    Prepare a single input R1/R2 or single-end FASTQ for fastp.
    Multiple lane files are merged into gzipped temporary inputs.
    """
    fq_info = _classify_fastqs(fastqs)
    prepared: Dict[str, Any] = {"layout": fq_info["layout"], "r1": [], "r2": [], "single": []}
    prep_dir = ensure_dir(work_dir / "_fastp_inputs")

    if fq_info["layout"] == "paired":
        r1 = _merge_fastqs(fq_info["r1"], prep_dir / "merged_R1.fastq.gz")
        r2 = _merge_fastqs(fq_info["r2"], prep_dir / "merged_R2.fastq.gz")
        prepared["r1"] = [str(r1)]
        prepared["r2"] = [str(r2)]
    else:
        single = _merge_fastqs(fq_info["single"], prep_dir / "merged_single.fastq.gz")
        prepared["single"] = [str(single)]
    return prepared


# -----------------------------------------------------------------------------
# QC and trimming
# -----------------------------------------------------------------------------

def run_fastqc(fastqs: List[Path], out_dir: Path, threads: int = 4) -> List[Path]:
    ensure_dir(out_dir)
    require_tools(["fastqc"])
    cmd = ["fastqc", "-t", str(threads), "-o", str(out_dir), *[str(x) for x in fastqs]]
    run_cmd(cmd, log_path=out_dir / "fastqc.log")
    return sorted(out_dir.glob("*_fastqc.html"))


def run_fastp(
    fastqs: List[Path],
    out_dir: Path,
    threads: int = 8,
    min_length: int = 30,
    qualified_q: int = 20,
) -> Tuple[List[Path], Path, Path]:
    """
    fastp for short-read trimming only.
    Returns trimmed FASTQ paths and fastp json/html reports.
    """
    ensure_dir(out_dir)
    require_tools(["fastp"])

    fq_info = _prepare_fastp_inputs(fastqs, out_dir)
    html = out_dir / "fastp.html"
    json_path = out_dir / "fastp.json"

    if fq_info["layout"] == "paired":
        r1_out = out_dir / "trimmed_R1.fastq.gz"
        r2_out = out_dir / "trimmed_R2.fastq.gz"
        cmd = [
            "fastp",
            "-w", str(threads),
            "-i", fq_info["r1"][0],
            "-I", fq_info["r2"][0],
            "-o", str(r1_out),
            "-O", str(r2_out),
            "--detect_adapter_for_pe",
            "--qualified_quality_phred", str(qualified_q),
            "--length_required", str(min_length),
            "--json", str(json_path),
            "--html", str(html),
        ]
        run_cmd(cmd, log_path=out_dir / "fastp.log")
        return [r1_out, r2_out], json_path, html

    single_out = out_dir / "trimmed.fastq.gz"
    cmd = [
        "fastp",
        "-w", str(threads),
        "-i", fq_info["single"][0],
        "-o", str(single_out),
        "--qualified_quality_phred", str(qualified_q),
        "--length_required", str(min_length),
        "--json", str(json_path),
        "--html", str(html),
    ]
    run_cmd(cmd, log_path=out_dir / "fastp.log")
    return [single_out], json_path, html


def run_chopper_ont(
    fastqs: List[Path],
    out_dir: Path,
    min_q: int = 8,
    min_len: int = 200,
    headcrop: int = 0,
    tailcrop: int = 0,
) -> List[Path]:
    """
    ONT trimming / filtering using chopper.
    If chopper is not installed, caller should decide whether to skip or fail.
    """
    ensure_dir(out_dir)
    require_tools(["chopper"])
    trimmed = []

    for fq in fastqs:
        fq = Path(fq)
        out_fq = out_dir / f"{fq.stem}.trimmed.fastq.gz"
        in_cat = ["gunzip", "-c", str(fq)] if str(fq).endswith(".gz") else ["cat", str(fq)]
        p1 = subprocess.Popen(in_cat, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        cmd = ["chopper", "-q", str(min_q), "-l", str(min_len)]
        if headcrop > 0:
            cmd += ["--headcrop", str(headcrop)]
        if tailcrop > 0:
            cmd += ["--tailcrop", str(tailcrop)]

        p2 = subprocess.Popen(cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p1.stdout.close()
        with open(out_fq, "wb") as oh:
            p3 = subprocess.run(["gzip", "-c"], stdin=p2.stdout, stdout=oh, stderr=subprocess.PIPE)
        p2.stdout.close()
        err1 = p1.stderr.read() if p1.stderr else b""
        err2 = p2.stderr.read() if p2.stderr else b""
        rc1 = p1.wait()
        rc2 = p2.wait()
        if rc1 != 0:
            raise PipelineError(f"input streaming failed for {fq.name}: {err1[:2000]}")
        if rc2 != 0:
            raise PipelineError(f"chopper failed for {fq.name}: {err2[:2000]}")
        if p3.returncode != 0:
            raise PipelineError(f"gzip failed for {fq.name}: {p3.stderr[:2000]}")
        trimmed.append(out_fq)

    return trimmed


def run_nanoplot(fastqs: List[Path], out_dir: Path, threads: int = 4) -> None:
    ensure_dir(out_dir)
    require_tools(["NanoPlot"])
    cmd = ["NanoPlot", "--threads", str(threads), "--outdir", str(out_dir), "--fastq", *[str(x) for x in fastqs]]
    run_cmd(cmd, log_path=out_dir / "nanoplot.log")


def parse_fastp_summary(json_path: Path) -> Dict[str, Any]:
    if not json_path.exists():
        return {}
    try:
        data = json.loads(json_path.read_text(encoding="utf-8"))
        before = data.get("summary", {}).get("before_filtering", {})
        after = data.get("summary", {}).get("after_filtering", {})
        return {
            "reads_before": before.get("total_reads"),
            "reads_after": after.get("total_reads"),
            "bases_before": before.get("total_bases"),
            "bases_after": after.get("total_bases"),
            "q20_rate_before": before.get("q20_rate"),
            "q20_rate_after": after.get("q20_rate"),
            "q30_rate_before": before.get("q30_rate"),
            "q30_rate_after": after.get("q30_rate"),
            "gc_before": before.get("gc_content"),
            "gc_after": after.get("gc_content"),
        }
    except Exception:
        return {}


# -----------------------------------------------------------------------------
# Alignment
# -----------------------------------------------------------------------------

def build_bowtie2_index(ref_fna: Path, out_dir: Path) -> Path:
    ensure_dir(out_dir)
    require_tools(["bowtie2-build"])
    prefix = out_dir / "bt2_index" / "genome"
    ensure_dir(prefix.parent)

    expected = Path(str(prefix) + ".1.bt2")
    if not expected.exists():
        run_cmd(["bowtie2-build", str(ref_fna), str(prefix)], log_path=out_dir / "bowtie2_build.log")
    return prefix


def align_short_reads_bowtie2(ref_fna: Path, fastqs: List[Path], out_dir: Path, threads: int = 8) -> Tuple[Path, Dict[str, Any]]:
    ensure_dir(out_dir)
    require_tools(["bowtie2", "bowtie2-build", "samtools"])
    fq_info = _classify_fastqs(fastqs)
    index_prefix = build_bowtie2_index(ref_fna, out_dir)
    bam = out_dir / "aligned.sorted.bam"
    stderr_path = out_dir / "bowtie2.stderr.txt"

    if fq_info["layout"] == "paired":
        align_cmd = [
            "bowtie2",
            "--very-sensitive-local",
            "-p", str(threads),
            "-x", str(index_prefix),
            "-1", ",".join(fq_info["r1"]),
            "-2", ",".join(fq_info["r2"]),
        ]
    else:
        align_cmd = [
            "bowtie2",
            "--very-sensitive-local",
            "-p", str(threads),
            "-x", str(index_prefix),
            "-U", ",".join(fq_info["single"]),
        ]

    p1 = subprocess.Popen(align_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2 = subprocess.run(
        ["samtools", "sort", "-@", str(threads), "-o", str(bam), "-"],
        stdin=p1.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    p1.stdout.close()
    stderr = p1.stderr.read()
    rc = p1.wait()
    stderr_path.write_bytes(stderr or b"")

    if rc != 0:
        raise PipelineError(f"bowtie2 failed: {stderr[:2000]}")
    if p2.returncode != 0:
        raise PipelineError(f"samtools sort failed: {p2.stderr[:2000]}")

    run_cmd(["samtools", "index", str(bam)], log_path=out_dir / "samtools.log")
    bt2_metrics = _parse_bowtie2_stderr(stderr.decode("utf-8", errors="ignore"))
    return bam, bt2_metrics


def align_long_reads_minimap2(
    ref_fna: Path,
    fastqs: List[Path],
    out_dir: Path,
    threads: int = 8,
    preset: str = "map-ont",
) -> Path:
    ensure_dir(out_dir)
    require_tools(["minimap2", "samtools"])
    bam = out_dir / "aligned.sorted.bam"

    mm_args = ["-ax", preset]
    p1 = subprocess.Popen(
        ["minimap2", "-t", str(threads), *mm_args, str(ref_fna), *[str(x) for x in fastqs]],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    p2 = subprocess.run(
        ["samtools", "sort", "-@", str(threads), "-o", str(bam), "-"],
        stdin=p1.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    p1.stdout.close()
    err = p1.stderr.read()
    rc = p1.wait()
    (out_dir / "minimap2.stderr.txt").write_bytes(err or b"")

    if rc != 0:
        raise PipelineError(f"minimap2 failed: {err[:2000]}")
    if p2.returncode != 0:
        raise PipelineError(f"samtools sort failed: {p2.stderr[:2000]}")

    run_cmd(["samtools", "index", str(bam)], log_path=out_dir / "samtools.log")
    return bam


def collect_alignment_metrics(bam: Path, out_dir: Path) -> Dict[str, Any]:
    ensure_dir(out_dir)
    cp1 = run_cmd(["samtools", "flagstat", str(bam)], log_path=out_dir / "samtools.log")
    cp2 = run_cmd(["samtools", "idxstats", str(bam)], log_path=out_dir / "samtools.log")
    cp3 = run_cmd(["samtools", "stats", str(bam)], log_path=out_dir / "samtools.log")
    (out_dir / "flagstat.txt").write_text(cp1.stdout, encoding="utf-8")
    (out_dir / "idxstats.txt").write_text(cp2.stdout, encoding="utf-8")
    (out_dir / "samtools.stats.txt").write_text(cp3.stdout, encoding="utf-8")
    return _parse_flagstat_text(cp1.stdout)


# -----------------------------------------------------------------------------
# Counting and TPM
# -----------------------------------------------------------------------------

def _bam_is_paired(bam: Path) -> bool:
    require_tools(["samtools"])
    total_cp = run_cmd(["samtools", "view", "-c", str(bam)], check=True)
    paired_cp = run_cmd(["samtools", "view", "-c", "-f", "1", str(bam)], check=True)
    try:
        total = int(total_cp.stdout.strip() or "0")
        paired = int(paired_cp.stdout.strip() or "0")
    except Exception:
        return False
    return total > 0 and (paired / total) >= 0.5


def run_featurecounts(
    bam: Path,
    gff: Path,
    out_dir: Path,
    threads: int = 8,
    locus_map_tsv: Path | None = None,
    stranded: int = 0,
) -> Path:
    ensure_dir(out_dir)
    require_tools(["featureCounts"])
    out_txt = out_dir / "featureCounts.txt"
    paired = _bam_is_paired(bam)
    pair_args = ["-p", "-B", "-C"] if paired else []
    strand_args = ["-s", str(stranded)]

    # SAF is often more robust for bacterial annotations than trusting GFF attributes.
    if locus_map_tsv is not None and Path(locus_map_tsv).exists():
        locus = pd.read_csv(locus_map_tsv, sep="\t", dtype=str).fillna("")
        need = {"locus_tag", "seqid", "start", "end"}
        if need.issubset(set(locus.columns)):
            saf = out_dir / "locus.saf"
            if "strand" in locus.columns:
                strand_series = locus["strand"].astype(str)
            else:
                strand_series = pd.Series(["."] * len(locus), index=locus.index)
            saf_df = pd.DataFrame({
                "GeneID": locus["locus_tag"].astype(str),
                "Chr": locus["seqid"].astype(str),
                "Start": locus["start"].astype(int),
                "End": locus["end"].astype(int),
                "Strand": strand_series,
            })
            saf_df["Strand"] = saf_df["Strand"].where(saf_df["Strand"].isin(["+", "-"]), ".")
            saf_df.to_csv(saf, sep="\t", index=False)
            cmd = [
                "featureCounts",
                "-T", str(threads),
                "-F", "SAF",
                "-a", str(saf),
                "-o", str(out_txt),
                *pair_args,
                *strand_args,
                str(bam),
            ]
            run_cmd(cmd, log_path=out_dir / "featureCounts.log")
            return out_txt

    ftype, gid = detect_featurecounts_keys(gff)
    cmd = [
        "featureCounts",
        "-T", str(threads),
        "-F", "GFF",
        "-a", str(gff),
        "-o", str(out_txt),
        "-t", ftype,
        "-g", gid,
        *pair_args,
        *strand_args,
        str(bam),
    ]
    run_cmd(cmd, log_path=out_dir / "featureCounts.log")
    return out_txt


def featurecounts_to_tpm(fc_txt: Path, locus_map_tsv: Path, out_dir: Path) -> Tuple[Path, Path]:
    ensure_dir(out_dir)
    locus = pd.read_csv(locus_map_tsv, sep="\t", dtype=str).dropna(subset=["locus_tag"]).drop_duplicates(subset=["locus_tag"])
    lens = locus.set_index("locus_tag")["length"].astype(float).to_dict()

    lines = fc_txt.read_text(encoding="utf-8", errors="ignore").splitlines()
    data_lines = [ln for ln in lines if ln and not ln.startswith("#")]
    df = pd.read_csv(io.StringIO("\n".join(data_lines)), sep="\t")
    count_col = df.columns[-1]
    df = df.rename(columns={"Geneid": "locus_tag", count_col: "count"})
    df["locus_tag"] = df["locus_tag"].astype(str)

    if "Length" in df.columns:
        df["length"] = df["Length"].astype(float)
    else:
        df["length"] = df["locus_tag"].map(lens).astype(float)

    df = df.dropna(subset=["count", "length"])
    df["count"] = df["count"].astype(float)
    df["length"] = df["length"].astype(float).clip(lower=1.0)
    df["rpk"] = df["count"] / (df["length"] / 1000.0)
    denom = float(df["rpk"].sum())
    df["tpm"] = 0.0 if denom <= 0 else (df["rpk"] / denom * 1e6)

    counts_path = out_dir / "counts.tsv"
    tpm_path = out_dir / "tpm.tsv"
    df[["locus_tag", "count", "length"]].to_csv(counts_path, sep="\t", index=False)
    df[["locus_tag", "tpm"]].to_csv(tpm_path, sep="\t", index=False)
    return counts_path, tpm_path


# -----------------------------------------------------------------------------
# One-run workflow
# -----------------------------------------------------------------------------

def rnaseq_one_run_standard(
    run_id: str,
    biosample: str,
    genome_fna: Path,
    gff: Path,
    locus_map: Path,
    ena_ftp: str,
    ena_md5: str,
    out_root: Path,
    fastq_root: Path | None = None,
    threads: int = 8,
    platform_hint: str = "",
    model_hint: str = "",
    do_trim: bool = True,
    stranded: int = 0,
    force: bool = False,
) -> Dict[str, Any]:
    run_dir = ensure_dir(out_root / "runs" / run_id)
    done = run_dir / ".done"
    metrics_path = run_dir / "run_metrics.json"
    if done.exists() and not force and metrics_path.exists():
        return json.loads(metrics_path.read_text(encoding="utf-8"))

    fastq_base = Path(fastq_root) if fastq_root else (out_root / "fastq")
    fastq_dir = ensure_dir(fastq_base / run_id)

    existing = sorted(
        list(fastq_dir.glob("*.fastq.gz")) +
        list(fastq_dir.glob("*.fastq")) +
        list(fastq_dir.glob("*.fq.gz")) +
        list(fastq_dir.glob("*.fq"))
    )

    if existing:
        fastqs = existing
        dl_method = "LOCAL"
    else:
        fastqs, dl_method = download_fastqs_resilient(
            run_id=run_id,
            out_dir=fastq_dir,
            threads=threads,
            ena_ftp_field=ena_ftp,
            ena_md5_field=ena_md5,
        )

    if not fastqs:
        raise PipelineError(f"No FASTQ files found for {run_id}")

    raw_qc = fastq_sample_stats(fastqs[0], max_reads=50000)
    (run_dir / "raw_fastq_sample_stats.json").write_text(json.dumps(raw_qc, indent=2), encoding="utf-8")

    tech = _detect_run_technology(platform_hint=platform_hint, model_hint=model_hint, read_qc=raw_qc)

    out: Dict[str, Any] = {
        "run_id": run_id,
        "biosample_accession": biosample,
        "download_method": dl_method,
        "technology": tech,
        "fastqs_raw": [str(x) for x in fastqs],
        "raw_qc": raw_qc,
        "status": "running",
    }

    # 1) Raw QC
    qc_dir = ensure_dir(run_dir / "qc")
    if tech == "illumina":
        run_fastqc(fastqs, qc_dir / "raw_fastqc", threads=min(threads, 8))
    else:
        if _tool_exists("NanoPlot"):
            run_nanoplot(fastqs, qc_dir / "raw_nanoplot", threads=min(threads, 8))

    # 2) Trimming / filtering
    trim_dir = ensure_dir(run_dir / "trim")
    trimmed_fastqs: List[Path]
    trim_metrics: Dict[str, Any] = {}
    warnings: List[str] = []

    if tech == "illumina":
        if do_trim:
            trimmed_fastqs, fastp_json, fastp_html = run_fastp(fastqs, trim_dir, threads=threads)
            trim_metrics = parse_fastp_summary(fastp_json)
            out["fastp_json"] = str(fastp_json)
            out["fastp_html"] = str(fastp_html)
        else:
            trimmed_fastqs = fastqs
        run_fastqc(trimmed_fastqs, qc_dir / "trimmed_fastqc", threads=min(threads, 8))
    else:
        if do_trim and _tool_exists("chopper"):
            trimmed_fastqs = run_chopper_ont(fastqs, trim_dir)
        else:
            trimmed_fastqs = fastqs
            if do_trim and not _tool_exists("chopper"):
                warnings.append("ONT_TRIMMING_SKIPPED_CHOPPER_NOT_FOUND")
        if _tool_exists("NanoPlot"):
            run_nanoplot(trimmed_fastqs, qc_dir / "trimmed_nanoplot", threads=min(threads, 8))
        elif tech == "ont":
            warnings.append("ONT_QC_NANOPLOT_NOT_FOUND")

    out["fastqs_trimmed"] = [str(x) for x in trimmed_fastqs]
    out["trim_metrics"] = trim_metrics

    # 3) Alignment
    align_dir = ensure_dir(run_dir / "alignment")
    if tech == "illumina":
        bam, aligner_metrics = align_short_reads_bowtie2(genome_fna, trimmed_fastqs, align_dir, threads=threads)
        out["aligner"] = "bowtie2"
        out["aligner_metrics"] = aligner_metrics
    else:
        bam = align_long_reads_minimap2(genome_fna, trimmed_fastqs, align_dir, threads=threads, preset="map-ont")
        out["aligner"] = "minimap2"
        out["aligner_metrics"] = {}

    mapping = collect_alignment_metrics(bam, align_dir)
    out["bam"] = str(bam)
    out["mapping"] = mapping

    mp = mapping.get("mapped_pct")
    if mp is not None and float(mp) < 20:
        warnings.append("LOW_MAPPING_PCT<20")
    if raw_qc.get("q_mean") is not None and float(raw_qc["q_mean"]) < 7:
        warnings.append("LOW_MEAN_Q<7")
    if raw_qc.get("n_rate_mean") is not None and float(raw_qc["n_rate_mean"]) > 0.05:
        warnings.append("HIGH_N_RATE>5%")

    # 4) Quantification
    counts_dir = ensure_dir(run_dir / "counts")
    fc_txt = run_featurecounts(bam, gff, counts_dir, threads=threads, locus_map_tsv=locus_map, stranded=stranded)
    counts_tsv, tpm_tsv = featurecounts_to_tpm(fc_txt, locus_map, counts_dir)
    out["featurecounts_txt"] = str(fc_txt)
    out["counts_tsv"] = str(counts_tsv)
    out["tpm_tsv"] = str(tpm_tsv)

    out["warnings"] = warnings
    out["status"] = "ok"
    metrics_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    done.write_text("ok\n", encoding="utf-8")
    return out


# -----------------------------------------------------------------------------
# Batch runner
# -----------------------------------------------------------------------------

def run_multiqc(out_dir: Path) -> Path | None:
    if not _tool_exists("multiqc"):
        return None
    qc_root = ensure_dir(out_dir / "multiqc")
    run_cmd(["multiqc", str(out_dir), "-o", str(qc_root), "-f"], log_path=qc_root / "multiqc.log")
    report = qc_root / "multiqc_report.html"
    return report if report.exists() else None


def run_rnaseq_standard(
    manifest_tsv: Path,
    genomes_dir: Path,
    out_dir: Path,
    fastq_root: Path | None = None,
    threads: int = 8,
    max_runs_per_biosample: int = 999,
    do_trim: bool = True,
    stranded: int = 0,
    force: bool = False,
    make_multiqc: bool = True,
) -> Path:
    ensure_dir(out_dir)
    man = pd.read_csv(manifest_tsv, sep="\t", dtype=str).fillna("")
    results: List[Dict[str, Any]] = []

    for _, r in man.iterrows():
        bs = str(r.get("biosample_accession", "")).strip()
        runs = safe_split_semi(r.get("run_list", ""))
        if not bs or bs.lower() == "nan" or not runs:
            continue

        runs = runs[:max_runs_per_biosample]
        bs_dir = Path(genomes_dir) / bs
        genome_fna = bs_dir / "genome.fna"
        gff = bs_dir / "annotation.gff"
        locus_map = bs_dir / "locus_map.tsv"
        if not genome_fna.exists() or not gff.exists() or not locus_map.exists():
            raise PipelineError(f"Missing genome files for {bs}. Run Module B first.")

        ena_ftp = r.get("ena_fastq_ftp", "")
        ena_md5 = r.get("ena_fastq_md5", "")
        platform_hint = r.get("instrument_platform", "") or r.get("platform", "")
        model_hint = r.get("instrument_model", "")

        for run_id in runs:
            try:
                res = rnaseq_one_run_standard(
                    run_id=run_id,
                    biosample=bs,
                    genome_fna=genome_fna,
                    gff=gff,
                    locus_map=locus_map,
                    ena_ftp=ena_ftp,
                    ena_md5=ena_md5,
                    out_root=out_dir,
                    fastq_root=fastq_root,
                    threads=threads,
                    platform_hint=platform_hint,
                    model_hint=model_hint,
                    do_trim=do_trim,
                    stranded=stranded,
                    force=force,
                )
            except Exception as e:
                res = {
                    "run_id": run_id,
                    "biosample_accession": bs,
                    "status": "fail",
                    "error": str(e)[:4000],
                }
            results.append(res)

    df = pd.DataFrame(results)
    out_path = out_dir / "rnaseq_run_index.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    log(f"RNA-seq index: {out_path}")

    if make_multiqc:
        try:
            report = run_multiqc(out_dir)
            if report:
                log(f"MultiQC report: {report}")
        except Exception as e:
            log(f"MultiQC skipped/failed: {e}")

    return out_path


# Backwards-compatible name expected by the original CLI.
run_rnaseq = run_rnaseq_standard
