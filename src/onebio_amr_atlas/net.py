from __future__ import annotations
import os, time, random
import re, shutil
from pathlib import Path
from typing import List, Tuple
import requests
from .utils import ensure_dir, log, PipelineError, md5_file, run_cmd, which

NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def http_get(url: str, out: Path, timeout: int = 120, retries: int = 6) -> None:
    """HTTP GET with retries + NCBI-friendly headers.

    - Writes to a temporary .part file then renames (avoids partial corrupt outputs).
    - On some networks NCBI may return 403 for the ftp host; we retry with the download host.
    - If requests keeps failing and curl is available, fall back to curl.
    """
    ensure_dir(out.parent)

    ua = os.environ.get("USER_AGENT") or os.environ.get("NCBI_TOOL") or "onebio/1.0"
    headers = {"User-Agent": ua, "Accept": "*/*"}

    def _alt_host(u: str) -> str:
        return u.replace("https://ftp.ncbi.nlm.nih.gov/", "https://download.ncbi.nlm.nih.gov/")

    part = out.with_suffix(out.suffix + ".part")

    last_err: Exception | None = None
    cur_url = url

    for i in range(max(1, retries)):
        try:
            with requests.get(cur_url, stream=True, timeout=timeout, headers=headers, allow_redirects=True) as r:
                if r.status_code == 403 and "ftp.ncbi.nlm.nih.gov" in cur_url:
                    alt = _alt_host(cur_url)
                    if alt != cur_url:
                        cur_url = alt
                        raise requests.HTTPError("403 on ftp host; retrying with download host", response=r)
                r.raise_for_status()
                with part.open("wb") as f:
                    for chunk in r.iter_content(chunk_size=2**20):
                        if chunk:
                            f.write(chunk)
            if part.exists():
                part.replace(out)
            return
        except Exception as e:
            last_err = e
            sleep_s = min(30.0, 1.0 * (2**i)) + random.random() * 0.25
            time.sleep(sleep_s)

    if which("curl"):
        log(f"Retrying with curl: {cur_url}")
        run_cmd(["curl", "-L", "-A", ua, "-o", str(out), cur_url], log_path=out.parent / "http_get.log")
        if out.exists() and out.stat().st_size > 0:
            return

    if part.exists():
        try:
            part.unlink()
        except Exception:
            pass

    raise last_err  # type: ignore


def gunzip_to(src_gz: Path, dst: Path) -> None:
    import gzip
    ensure_dir(dst.parent)
    with gzip.open(src_gz, "rb") as f_in, dst.open("wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

def ncbi_assembly_esearch_uid(assembly_acc: str) -> str:
    params = {"db": "assembly", "term": assembly_acc, "retmode": "json"}
    r = requests.get(f"{NCBI_EUTILS_BASE}/esearch.fcgi", params=params, timeout=60)
    r.raise_for_status()
    js = r.json()
    ids = js.get("esearchresult", {}).get("idlist", [])
    if not ids:
        raise PipelineError(f"NCBI esearch found no Assembly UID for {assembly_acc}")
    return ids[0]

def ncbi_assembly_ftp_paths(assembly_acc: str) -> Tuple[str, str]:
    uid = ncbi_assembly_esearch_uid(assembly_acc)
    params = {"db": "assembly", "id": uid, "retmode": "json"}
    r = requests.get(f"{NCBI_EUTILS_BASE}/esummary.fcgi", params=params, timeout=60)
    r.raise_for_status()
    js = r.json()
    doc = js.get("result", {}).get(uid, {})
    gb = doc.get("ftppath_genbank", "") or doc.get("FtpPath_GenBank", "") or ""
    rs = doc.get("ftppath_refseq", "") or doc.get("FtpPath_RefSeq", "") or ""
    return gb, rs

def ncbi_download_assembly_and_gff(assembly_acc: str, out_dir: Path, prefer: str = "genbank") -> Tuple[Path, Path]:
    ensure_dir(out_dir)
    # If we already have the decompressed outputs, reuse them (avoids re-downloading).
    genome_fna = out_dir / "genome.fna"
    annot_gff = out_dir / "annotation.gff"
    if genome_fna.exists() and annot_gff.exists() and genome_fna.stat().st_size > 0 and annot_gff.stat().st_size > 0:
        return genome_fna, annot_gff

    gb, rs = ncbi_assembly_ftp_paths(assembly_acc)
    ftp = gb if (prefer == "genbank" and gb) else (rs if rs else gb)
    if not ftp:
        raise PipelineError(f"Could not resolve FTP path for {assembly_acc}")

    # NCBI gives ftp:// paths; requests can't fetch ftp://. Use https instead.
    if ftp.startswith('ftp://'):
        ftp = 'https://' + ftp[len('ftp://'):]

    base = ftp.rstrip("/").split("/")[-1]
    fna_gz = out_dir / f"{base}_genomic.fna.gz"
    gff_gz = out_dir / f"{base}_genomic.gff.gz"
    md5_txt = out_dir / "md5checksums.txt"

    # md5 list (best effort)
    try:
        http_get(f"{ftp}/md5checksums.txt", md5_txt)
    except Exception:
        md5_txt = None

    http_get(f"{ftp}/{fna_gz.name}", fna_gz)
    http_get(f"{ftp}/{gff_gz.name}", gff_gz)

    if md5_txt and md5_txt.exists():
        md5_map = {}
        for line in md5_txt.read_text(encoding="utf-8", errors="ignore").splitlines():
            parts = line.strip().split()
            if len(parts) >= 2:
                md5_map[parts[1].lstrip("*")] = parts[0]
        for p in [fna_gz, gff_gz]:
            exp = md5_map.get(p.name)
            if exp:
                got = md5_file(p)
                if got != exp:
                    raise PipelineError(f"MD5 mismatch for {p.name}: expected {exp}, got {got}")

    genome_fna = out_dir / "genome.fna"
    annot_gff = out_dir / "annotation.gff"
    gunzip_to(fna_gz, genome_fna)
    gunzip_to(gff_gz, annot_gff)
    return genome_fna, annot_gff

def ena_download_fastqs(ena_ftp_field: str, ena_md5_field: str, out_dir: Path) -> List[Path]:
    ensure_dir(out_dir)
    urls = []
    for part in re.split(r"[;,\s]+", str(ena_ftp_field or "").strip()):
        if not part or str(part).lower() == "nan":
            continue
        if part.startswith("ftp://") or part.startswith("http://") or part.startswith("https://"):
            urls.append(part)
        else:
            urls.append("https://" + part.lstrip("/"))
    md5s = [m for m in re.split(r"[;,\s]+", str(ena_md5_field or "").strip()) if m and m.lower() != "nan"]

    out_paths: List[Path] = []
    for i, url in enumerate(urls):
        name = url.split("/")[-1]
        out = out_dir / name
        if not out.exists():
            log(f"Downloading ENA: {url}")
            http_get(url, out, timeout=240)
        if i < len(md5s):
            exp = md5s[i]
            got = md5_file(out)
            if got != exp:
                raise PipelineError(f"ENA MD5 mismatch for {out.name}: expected {exp}, got {got}")
        out_paths.append(out)
    return out_paths

def sra_download_fastqs(run_id: str, out_dir: Path, threads: int = 8) -> List[Path]:
    ensure_dir(out_dir)
    if which("prefetch") is None or which("fasterq-dump") is None:
        raise PipelineError("SRA Toolkit not found (prefetch/fasterq-dump). Provide ENA fastqs or install SRA Toolkit.")
    run_cmd(["prefetch", run_id, "-O", str(out_dir)], log_path=out_dir/"sra_download.log")
    sra_path = out_dir / run_id / f"{run_id}.sra"
    if not sra_path.exists():
        alt = out_dir / f"{run_id}.sra"
        sra_path = alt if alt.exists() else sra_path
    if not sra_path.exists():
        raise PipelineError(f"prefetch succeeded but .sra not found for {run_id}")
    run_cmd(["fasterq-dump", str(sra_path), "-O", str(out_dir), "-e", str(threads)], log_path=out_dir/"sra_download.log")

    fastqs = sorted(out_dir.glob(f"{run_id}*.fastq"))
    if not fastqs:
        fastqs = sorted(out_dir.glob(f"{run_id}*.fq"))

    gz_paths = []
    if which("pigz"):
        for fq in fastqs:
            run_cmd(["pigz", "-f", "-p", str(threads), str(fq)], log_path=out_dir/"sra_download.log")
            gz_paths.append(Path(str(fq) + ".gz"))
    else:
        gz_paths = fastqs
    return gz_paths

# -----------------------------------------------------------------------------
# Resilient FASTQ resolver used by the browser UI / unified runner
# -----------------------------------------------------------------------------
def ena_lookup_fastqs_by_run(run_id: str, out_dir: Path | None = None) -> tuple[str, str]:
    """Resolve FASTQ FTP paths for a run accession through ENA Portal.

    Returns (fastq_ftp, fastq_md5). Empty strings mean ENA did not return FASTQs.
    This is deliberately small and dependency-light because download failures are
    already enough theatre without adding more moving parts.
    """
    import pandas as pd
    from io import StringIO

    run_id = str(run_id or "").strip().upper()
    if not run_id:
        return "", ""

    cache_file = None
    if out_dir is not None:
        ensure_dir(Path(out_dir))
        cache_file = Path(out_dir) / f"{run_id}.ena_lookup.tsv"
        if cache_file.exists() and cache_file.stat().st_size > 0:
            try:
                df = pd.read_csv(cache_file, sep="\t", dtype=str).fillna("")
                if not df.empty:
                    return str(df.iloc[0].get("fastq_ftp", "")), str(df.iloc[0].get("fastq_md5", ""))
            except Exception:
                pass

    params = {
        "accession": run_id,
        "result": "read_run",
        "fields": "run_accession,fastq_ftp,fastq_md5,library_strategy,secondary_sample_accession,sample_accession",
        "format": "tsv",
        "download": "false",
    }
    try:
        r = requests.get(ENA_PORTAL, params=params, timeout=120)
        r.raise_for_status()
        txt = r.text.strip()
        if not txt:
            return "", ""
        df = pd.read_csv(StringIO(txt), sep="\t", dtype=str).fillna("")
        if cache_file is not None:
            df.to_csv(cache_file, sep="\t", index=False)
        if df.empty:
            return "", ""
        fastq_ftp = ";".join([x for x in df.get("fastq_ftp", pd.Series(dtype=str)).astype(str).tolist() if x and x.lower() != "nan"])
        fastq_md5 = ";".join([x for x in df.get("fastq_md5", pd.Series(dtype=str)).astype(str).tolist() if x and x.lower() != "nan"])
        return fastq_ftp, fastq_md5
    except Exception as e:
        log(f"ENA lookup failed for {run_id}: {e}")
        return "", ""


def download_fastqs_resilient(
    run_id: str,
    out_dir: Path,
    threads: int = 8,
    ena_ftp_field: str = "",
    ena_md5_field: str = "",
) -> tuple[List[Path], str]:
    """Try every practical route before declaring a run dead.

    Order:
      1. ENA URLs already present in the checker/manifest
      2. ENA lookup by run accession, then ENA HTTPS download
      3. SRA Toolkit prefetch + fasterq-dump

    Returns (fastq_paths, method_label).
    """
    ensure_dir(out_dir)

    # Existing local FASTQs win immediately. Useful for demo mode and reruns.
    existing = sorted(
        list(out_dir.glob("*.fastq.gz")) +
        list(out_dir.glob("*.fastq")) +
        list(out_dir.glob("*.fq.gz")) +
        list(out_dir.glob("*.fq"))
    )
    if existing:
        return existing, "LOCAL"

    if ena_ftp_field and str(ena_ftp_field).strip() and str(ena_ftp_field).lower() != "nan":
        try:
            fq = ena_download_fastqs(str(ena_ftp_field), str(ena_md5_field), out_dir)
            if fq:
                return fq, "ENA_MANIFEST"
        except Exception as e:
            log(f"ENA manifest download failed for {run_id}: {e}")

    try:
        ena_ftp, ena_md5 = ena_lookup_fastqs_by_run(run_id, out_dir=out_dir)
        if ena_ftp:
            fq = ena_download_fastqs(ena_ftp, ena_md5, out_dir)
            if fq:
                return fq, "ENA_RUN_LOOKUP"
    except Exception as e:
        log(f"ENA run lookup/download failed for {run_id}: {e}")

    try:
        fq = sra_download_fastqs(run_id, out_dir, threads=threads)
        if fq:
            return fq, "SRA_TOOLKIT"
    except Exception as e:
        log(f"SRA Toolkit download failed for {run_id}: {e}")

    return [], "FAILED"
