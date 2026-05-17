import os
import time
import random
import json
import hashlib
import shutil
import zipfile
import re
import argparse
from io import StringIO
from urllib.parse import quote
from typing import Optional, List, Tuple
from datetime import datetime

import pandas as pd
import requests
from tqdm import tqdm

# ============================================================
# Multi-source Conservative RNA-seq + AMR checker (union mode)
#
# Goal:
# - Discover phenotype candidates from BV-BRC + NCBI AST Browser + EBI AMR Portal
# - Take the union of those phenotype candidates
# - Then validate genome + RNA-seq conservatively
#
# Key behaviour:
# - RNA-seq remains STRONG BioSample-only via SRA + ENA
# - BV-BRC still provides the rich genome metadata backbone where available
# - External-only phenotype candidates are retained as synthetic rows so you can see
#   whether NCBI/EBI add extra RNA-seq-positive hits beyond BV-BRC AMR
# ============================================================

SPECIES_INPUTS_DEFAULT = (
    "Enterococcus faecium, Staphylococcus aureus, Klebsiella pneumoniae, "
    "Acinetobacter baumannii, Pseudomonas aeruginosa, Enterobacter"
)

RESET_CACHE_FOR_EACH_SPECIES = False
ONLY_GENOMES_WITH_AMR = True
PRINT_STEP_STATS = True
PRINT_EVERY_N_TAXON_IDS = 25
SHOW_TQDM_BARS = True
BVBRC_VERBOSE = False

MAX_BVBRC_GENOMES = None
MAX_BVBRC_AMR_RECORDS_PER_TAXON = None
MAX_BIOSAMPLES_PER_SPECIES = None
MAX_SRA_RECORDS = None

BVBRC_PAGE = 25000
SRA_EFETCH_PAGE = 500

NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "").strip()
NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "").strip()
RPS = 2 if not NCBI_API_KEY else 6
MIN_INTERVAL = 1.0 / RPS

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
BVBRC_API = "https://www.bv-brc.org/api/"
ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
EBI_AMR_RELEASES_URL = "https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/"

RUN_RE = re.compile(r"(?:SRR|ERR|DRR)\d+", re.IGNORECASE)
BS_TOKEN_RE = re.compile(r"(?:SAMN|SAMEA|SAMD)\d+", re.IGNORECASE)
ASM_RE = re.compile(r"(?:GCA|GCF)_\d+\.\d+", re.IGNORECASE)
YYYY_MM_RE = re.compile(r"(20\d{2}-\d{2})")

sess = requests.Session()
sess.headers.update({"User-Agent": "rnaseq_amr_checker_multisource_union/8.0"})
_last_call = 0.0

SOURCE_BVBRC = "BV-BRC"
SOURCE_NCBI = "NCBI_AST_BROWSER"
SOURCE_EBI = "EBI_AMR_PORTAL"

LAB_POS = re.compile(
    r"(?:disk|disc)\s*diffusion|kirby\s*[- ]?bauer|broth\s*dilution|agar\s*dilution|"
    r"\bmic\b|minimum\s*inhibitory\s*concentration|e-?\s*test|etest|"
    r"\bvitek\b|\bphoenix\b|sensititre|micronaut|"
    r"\bclsi\b|\beucast\b|"
    r"\bast\b|antimicrobial\s*susceptibility\s*test|phenotypic|phenotype|"
    r"zone\s*diameter|inhibition\s*zone",
    re.IGNORECASE,
)

NONLAB_NEG = re.compile(
    r"predicted|prediction|in\s*silico|computational|genotype|genomic|"
    r"resistance\s*gene|amr\s*gene|gene\s*presence|"
    r"amrfinder|resfinder|card|rgi|abricate|unifire|mettannotator|"
    r"k-?mer|machine\s*learning|\bml\b|model|"
    r"homology|blast|annotation|assembly|variant|mutation|snp|"
    r"rule[- ]based|inferred",
    re.IGNORECASE,
)

PHENO_CANON = {
    "resistant": "R",
    "r": "R",
    "non-susceptible": "R",
    "nonsusceptible": "R",
    "non susceptible": "R",
    "susceptible": "S",
    "s": "S",
    "susceptible-dose dependent": "S",
    "susceptible dose dependent": "S",
    "sdd": "S",
    "intermediate": "I",
    "i": "I",
}


def norm_species_name(name: str) -> str:
    name = " ".join(str(name).strip().split())
    if not name:
        raise ValueError("Empty species name.")
    words = name.split()
    words[0] = words[0].capitalize()
    for i in range(1, len(words)):
        words[i] = words[i].lower()
    return " ".join(words)


def normalize_text(x: str) -> str:
    return " ".join(str(x or "").strip().split())


def normalize_antibiotic(x: str) -> str:
    return normalize_text(str(x or "").lower())


def norm_pheno(x) -> str:
    return PHENO_CANON.get(normalize_text(str(x or "").lower()), "")


def rql_val(s: str) -> str:
    return quote(str(s), safe="")


def _params(extra: dict) -> dict:
    p = dict(extra)
    p["tool"] = "rnaseq_amr_checker_multisource_union"
    if NCBI_EMAIL:
        p["email"] = NCBI_EMAIL
    if NCBI_API_KEY:
        p["api_key"] = NCBI_API_KEY
    return p


def pick_first_token(x: str, token_re: re.Pattern) -> str:
    x = str(x or "")
    m = token_re.search(x)
    return m.group(0).upper() if m else ""


def extract_run_ids(x: str) -> list:
    return sorted(set([m.group(0).upper() for m in RUN_RE.finditer(str(x or ""))]))


def union_run_strings(*vals) -> str:
    out = set()
    for v in vals:
        for rid in extract_run_ids(v):
            out.add(rid)
    return ";".join(sorted(out))


def force_text_ids(df: pd.DataFrame, cols: List[str]) -> None:
    for c in cols:
        if c in df.columns:
            df[c] = df[c].fillna("").astype(str).str.strip()
        else:
            df[c] = ""


def safe_to_csv(df: pd.DataFrame, path: str, index: bool = False) -> str:
    try:
        df.to_csv(path, index=index)
        return path
    except PermissionError:
        base, ext = os.path.splitext(path)
        alt = f"{base}_{datetime.now().strftime('%Y%m%d_%H%M%S')}{ext}"
        df.to_csv(alt, index=index)
        return alt


def eutils_get(endpoint: str, params: dict, timeout: int = 60, max_retries: int = 8) -> requests.Response:
    global _last_call
    url = EUTILS + endpoint
    last_r = None
    last_exc = None
    payload = _params(params)
    term_len = len(str(payload.get("term", "")))
    use_post = (term_len > 1200)
    for attempt in range(max_retries):
        now = time.monotonic()
        wait = (_last_call + MIN_INTERVAL) - now
        if wait > 0:
            time.sleep(wait)
        try:
            if use_post:
                r = sess.post(url, data=payload, timeout=timeout)
            else:
                r = sess.get(url, params=payload, timeout=timeout)
            last_r = r
            _last_call = time.monotonic()

            if r.status_code == 429 or (500 <= r.status_code < 600):
                ra = r.headers.get("Retry-After")
                if ra and ra.isdigit():
                    sleep_s = int(ra)
                else:
                    sleep_s = min(60, 2 ** attempt) + random.uniform(0, 0.5)
                time.sleep(sleep_s)
                continue

            r.raise_for_status()
            return r
        except requests.RequestException as e:
            last_exc = e
            sleep_s = min(60, 2 ** attempt) + random.uniform(0, 0.5)
            time.sleep(sleep_s)
            continue

    if last_r is not None:
        last_r.raise_for_status()
        return last_r
    if last_exc is not None:
        raise last_exc
    raise RuntimeError("NCBI eutils request failed with no response")


def http_get_retry(url: str, timeout: int = 120, max_retries: int = 10, backoff_base: float = 1.6) -> requests.Response:
    last = None
    for attempt in range(max_retries):
        try:
            r = sess.get(url, timeout=timeout)
            if r.status_code in (429, 500, 502, 503, 504):
                ra = r.headers.get("Retry-After")
                if ra and ra.isdigit():
                    sleep_s = int(ra)
                else:
                    sleep_s = min(90, (backoff_base ** attempt)) + random.uniform(0, 0.8)
                time.sleep(sleep_s)
                continue
            r.raise_for_status()
            return r
        except requests.RequestException as e:
            last = e
            sleep_s = min(90, (backoff_base ** attempt)) + random.uniform(0, 0.8)
            time.sleep(sleep_s)
    raise last


def bvbrc_count(endpoint: str, rql_core: str, timeout: int = 120, label: Optional[str] = None) -> int:
    url = f"{BVBRC_API}{endpoint}/?{rql_core}&limit(0,0)&http_accept=application/solr+json"
    r = http_get_retry(url, timeout=timeout)
    j = r.json()
    n = int(j.get("response", {}).get("numFound", 0))
    if BVBRC_VERBOSE:
        tag = f" | {label}" if label else ""
        print(f"[BV-BRC]{tag} {endpoint} expected total: {n}")
    return n


def bvbrc_fetch_all(endpoint: str, rql_core: str, select_fields: list,
                    page: int = BVBRC_PAGE, max_rows: Optional[int] = None,
                    timeout: int = 120, label: Optional[str] = None) -> Tuple[pd.DataFrame, int]:
    expected = bvbrc_count(endpoint, rql_core, timeout=timeout, label=label)
    out = []
    start = 0
    while True:
        if max_rows is not None and start >= max_rows:
            break
        this_limit = page if max_rows is None else min(page, max_rows - start)
        url = (
            f"{BVBRC_API}{endpoint}/?{rql_core}"
            f"&select({','.join(select_fields)})"
            f"&limit({this_limit},{start})"
            f"&http_accept=application/json"
        )
        r = http_get_retry(url, timeout=timeout)
        batch = r.json()
        if not batch:
            break
        out.extend(batch)
        got = len(batch)
        start += got
        if got < this_limit:
            break
    df = pd.DataFrame(out)
    for col in select_fields:
        if col not in df.columns:
            df[col] = pd.Series(dtype="object")
    df = df[select_fields + [c for c in df.columns if c not in select_fields]]
    return df, expected


def ena_filereport(accessions: List[str], fields: List[str], chunk_size: int = 25, timeout: int = 120) -> pd.DataFrame:
    if not accessions:
        return pd.DataFrame()
    frames = []
    for i in range(0, len(accessions), chunk_size):
        chunk = accessions[i:i + chunk_size]
        params = {
            "accession": ",".join(chunk),
            "result": "read_run",
            "fields": ",".join(fields),
            "format": "tsv",
            "download": "false",
        }
        r = sess.get(ENA_PORTAL, params=params, timeout=timeout)
        if r.status_code in (429, 500, 502, 503, 504):
            r = http_get_retry(r.url, timeout=timeout)
        r.raise_for_status()
        txt = r.text.strip()
        if not txt:
            continue
        df = pd.read_csv(StringIO(txt), sep="\t")
        if not df.empty:
            frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def build_sra_runinfo_by_biosamples(biosamples: List[str], chunk_size: int = 20,
                                    max_records: Optional[int] = None) -> pd.DataFrame:
    if not biosamples:
        return pd.DataFrame()
    dfs = []
    fetched_total = 0

    def _fetch_chunk_recursive(chunk: List[str], depth: int = 0) -> List[pd.DataFrame]:
        if not chunk:
            return []
        term = "(" + " OR ".join([f"{bs}[BioSample]" for bs in chunk]) + ") AND rna seq[Strategy]"
        try:
            j0 = eutils_get("esearch.fcgi", {
                "db": "sra", "term": term, "retmode": "json", "retmax": 0, "usehistory": "y"
            }, timeout=120).json()
        except requests.RequestException:
            if len(chunk) == 1:
                return []
            mid = len(chunk) // 2
            return _fetch_chunk_recursive(chunk[:mid], depth + 1) + _fetch_chunk_recursive(chunk[mid:], depth + 1)

        es = j0.get("esearchresult", {})
        total = int(es.get("count", 0))
        webenv = es.get("webenv", "")
        query_key = es.get("querykey", "")
        if total == 0:
            return []

        local_frames = []
        to_fetch = total
        if max_records is not None:
            remaining = max_records - fetched_total
            if remaining <= 0:
                return []
            to_fetch = min(to_fetch, remaining)

        for start in range(0, to_fetch, SRA_EFETCH_PAGE):
            this_max = min(SRA_EFETCH_PAGE, to_fetch - start)
            try:
                txt = eutils_get("efetch.fcgi", {
                    "db": "sra",
                    "query_key": query_key,
                    "WebEnv": webenv,
                    "rettype": "runinfo",
                    "retmode": "text",
                    "retstart": start,
                    "retmax": this_max,
                }, timeout=120).text.strip()
            except requests.RequestException:
                if len(chunk) == 1:
                    return local_frames
                mid = len(chunk) // 2
                return _fetch_chunk_recursive(chunk[:mid], depth + 1) + _fetch_chunk_recursive(chunk[mid:], depth + 1)

            if txt:
                df = pd.read_csv(StringIO(txt))
                local_frames.append(df)
        return local_frames

    itr = range(0, len(biosamples), chunk_size)
    itr = tqdm(itr, desc="SRA BioSample chunks", disable=(not SHOW_TQDM_BARS))
    for i in itr:
        chunk = biosamples[i:i + chunk_size]
        frames = _fetch_chunk_recursive(chunk)
        for df in frames:
            dfs.append(df)
            fetched_total += len(df)
            if max_records is not None and fetched_total >= max_records:
                break
        if max_records is not None and fetched_total >= max_records:
            break

    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def run_map_biosample_only_sra_ena(candidates: pd.DataFrame, cache_dir: str,
                                   max_biosamples: Optional[int] = None,
                                   max_sra_records: Optional[int] = None) -> pd.DataFrame:
    g = candidates.copy()
    g["biosample_accession"] = g["biosample_accession"].apply(lambda x: pick_first_token(x, BS_TOKEN_RE))
    biosamples = sorted([x for x in g["biosample_accession"].unique().tolist() if x])
    if max_biosamples is not None:
        biosamples = biosamples[:max_biosamples]

    runinfo = build_sra_runinfo_by_biosamples(biosamples, chunk_size=60, max_records=max_sra_records)
    sra_bs_map = pd.DataFrame(columns=["biosample_accession", "rnaseq_run_ids_sra"])
    if not runinfo.empty:
        for col in ["Run", "BioSample", "LibraryStrategy"]:
            if col not in runinfo.columns:
                runinfo[col] = ""
        runinfo["Run"] = runinfo["Run"].fillna("").astype(str).str.strip().str.upper()
        runinfo["BioSample"] = runinfo["BioSample"].fillna("").astype(str).str.strip().str.upper()
        runinfo["LibraryStrategy"] = runinfo["LibraryStrategy"].fillna("").astype(str).str.strip().str.lower()
        ri = runinfo[(runinfo["BioSample"] != "") & (runinfo["Run"] != "") & (runinfo["LibraryStrategy"].str.contains("rna"))].copy()
        sra_bs_map = (
            ri.groupby("BioSample")["Run"]
            .apply(lambda s: ";".join(sorted(set([x for x in s if RUN_RE.search(str(x))]))[:500]))
            .reset_index()
            .rename(columns={"BioSample": "biosample_accession", "Run": "rnaseq_run_ids_sra"})
        )

    ena_fields = ["run_accession", "sample_accession", "secondary_sample_accession", "library_strategy", "fastq_ftp", "fastq_md5"]
    ena = ena_filereport(biosamples, fields=ena_fields, chunk_size=25, timeout=120)
    ena_bs_map = pd.DataFrame(columns=["biosample_accession", "rnaseq_run_ids_ena", "ena_fastq_ftp", "ena_fastq_md5"])
    if not ena.empty:
        for c in ena_fields:
            if c not in ena.columns:
                ena[c] = ""
        ena["library_strategy"] = ena["library_strategy"].fillna("").astype(str).str.lower()
        ena = ena[ena["library_strategy"].str.contains("rna")].copy()
        ena["biosample_accession"] = ena["secondary_sample_accession"].apply(lambda x: pick_first_token(x, BS_TOKEN_RE))
        mask_empty = ena["biosample_accession"].eq("")
        ena.loc[mask_empty, "biosample_accession"] = ena.loc[mask_empty, "sample_accession"].apply(lambda x: pick_first_token(x, BS_TOKEN_RE))
        ena["run_accession"] = ena["run_accession"].fillna("").astype(str).str.strip().str.upper()
        ena["fastq_ftp"] = ena["fastq_ftp"].fillna("").astype(str).str.strip()
        ena["fastq_md5"] = ena["fastq_md5"].fillna("").astype(str).str.strip()
        ena = ena[(ena["biosample_accession"] != "") & (ena["run_accession"] != "")]
        ena_bs_map = (
            ena.groupby("biosample_accession")
            .agg(
                rnaseq_run_ids_ena=("run_accession", lambda s: ";".join(sorted(set([x for x in s if RUN_RE.search(str(x))]))[:500])),
                ena_fastq_ftp=("fastq_ftp", lambda s: ";".join(sorted(set([x for x in s if x]))[:50])),
                ena_fastq_md5=("fastq_md5", lambda s: ";".join(sorted(set([x for x in s if x]))[:50])),
            )
            .reset_index()
        )

    bs_map = sra_bs_map.merge(ena_bs_map, on="biosample_accession", how="outer")
    bs_map["rnaseq_run_ids"] = bs_map.apply(lambda r: union_run_strings(r.get("rnaseq_run_ids_sra", ""), r.get("rnaseq_run_ids_ena", "")), axis=1)
    bs_map = bs_map[["biosample_accession", "rnaseq_run_ids", "rnaseq_run_ids_sra", "rnaseq_run_ids_ena", "ena_fastq_ftp", "ena_fastq_md5"]].copy()
    bs_map.to_parquet(os.path.join(cache_dir, "biosample_run_map_SRA_ENA_STRONG.parquet"), index=False)
    return bs_map


def source_lab_mask(df: pd.DataFrame, source_db: str) -> pd.Series:
    n = len(df)
    if n == 0:
        return pd.Series([], dtype=bool)

    def col_text(c):
        if c not in df.columns:
            return pd.Series([""] * n, index=df.index)
        return df[c].fillna("").astype(str).str.strip()

    lab_method = col_text("laboratory_typing_method")
    standard = col_text("testing_standard")
    meas = col_text("measurement_value")
    evidence = col_text("evidence")
    source = col_text("source")
    txt = (lab_method + " " + standard + " " + meas + " " + evidence + " " + source).str.lower()

    has_nonlab = txt.str.contains(NONLAB_NEG)
    has_labpos = txt.str.contains(LAB_POS)

    if source_db == SOURCE_BVBRC:
        strong_lab = (lab_method.ne("") | standard.ne("") | meas.ne(""))
        return (~has_nonlab) & (strong_lab | has_labpos)
    return ~has_nonlab


def make_species_cache_dir(cache_root: str, species: str) -> str:
    species_key = species.lower().replace(" ", "_")
    d = os.path.join(cache_root, species_key)
    os.makedirs(d, exist_ok=True)
    os.makedirs(os.path.join(d, "amr_by_taxon"), exist_ok=True)
    return d


def standardize_amr_columns(df: pd.DataFrame, source_db: str) -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame(columns=[
            "genome_id", "biosample_accession", "assembly_accession", "antibiotic", "antibiotic_canonical",
            "resistant_phenotype", "pheno", "laboratory_typing_method", "testing_standard", "measurement_value",
            "evidence", "source", "source_db", "source_record_uid", "raw_source_payload"
        ])
    out = df.copy()
    for c in [
        "genome_id", "biosample_accession", "assembly_accession", "antibiotic", "resistant_phenotype",
        "laboratory_typing_method", "testing_standard", "measurement_value", "evidence", "source",
        "source_record_uid", "raw_source_payload"
    ]:
        if c not in out.columns:
            out[c] = ""
    force_text_ids(out, [
        "genome_id", "biosample_accession", "assembly_accession", "antibiotic", "resistant_phenotype",
        "laboratory_typing_method", "testing_standard", "measurement_value", "evidence", "source",
        "source_record_uid", "raw_source_payload"
    ])
    out["biosample_accession"] = out["biosample_accession"].apply(lambda x: pick_first_token(x, BS_TOKEN_RE))
    out["assembly_accession"] = out["assembly_accession"].apply(lambda x: pick_first_token(x, ASM_RE))
    out["antibiotic_canonical"] = out["antibiotic"].apply(normalize_antibiotic)
    out["pheno"] = out["resistant_phenotype"].apply(norm_pheno)
    out["source_db"] = source_db
    out = out[source_lab_mask(out, source_db)].copy()
    out = out[(out["antibiotic_canonical"] != "") & (out["pheno"] != "")].copy()
    dedup_cols = [
        "genome_id", "biosample_accession", "assembly_accession", "antibiotic_canonical", "pheno",
        "laboratory_typing_method", "testing_standard", "measurement_value", "source_db"
    ]
    out = out.drop_duplicates(dedup_cols).reset_index(drop=True)
    return out


def detect_latest_ebi_release(timeout: int = 60) -> str:
    r = http_get_retry(EBI_AMR_RELEASES_URL, timeout=timeout)
    hits = sorted(set(YYYY_MM_RE.findall(r.text)))
    if not hits:
        raise RuntimeError("Could not detect any EMBL-EBI AMR portal release directories")
    return hits[-1]


def fetch_ebi_phenotypes(species: str, cache_dir: str, release: str = "auto") -> pd.DataFrame:
    if release == "auto":
        release = detect_latest_ebi_release()

    cache_file = os.path.join(cache_dir, f"ebi_amr_portal_phenotype_{release.replace('-', '_')}_{species.lower().replace(' ', '_')}.parquet")
    if os.path.exists(cache_file):
        return pd.read_parquet(cache_file)

    try:
        import duckdb
    except ImportError as e:
        raise RuntimeError("duckdb is required for EMBL-EBI AMR portal integration. Install with: pip install duckdb") from e

    parquet_url = f"https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/{release}/phenotype.parquet"
    species_norm = norm_species_name(species)
    words = species_norm.split()
    if len(words) >= 2:
        species_clause = "(lower(species) = lower(?) OR lower(organism) = lower(?))"
        params = [species_norm, species_norm]
    else:
        species_clause = "lower(genus) = lower(?)"
        params = [words[0]]

    # Write species-restricted data straight to parquet on disk first.
    # That avoids the worst in-memory blowups from .df() on very large remote tables.
    temp_extract = os.path.join(cache_dir, f"ebi_extract_{release.replace('-', '_')}_{species.lower().replace(' ', '_')}.parquet")
    con = duckdb.connect(database=":memory:")
    con.execute("INSTALL httpfs;")
    con.execute("LOAD httpfs;")
    con.execute("PRAGMA threads=4;")

    sql = f"""
    COPY (
      SELECT
        BioSample_ID AS biosample_accession,
        assembly_ID AS assembly_accession,
        antibiotic_name AS antibiotic,
        resistance_phenotype AS resistant_phenotype,
        laboratory_typing_method,
        ast_standard AS testing_standard,
        measurement AS measurement_value,
        database AS evidence,
        organism,
        species
      FROM read_parquet('{parquet_url}')
      WHERE {species_clause}
        AND (BioSample_ID IS NOT NULL OR assembly_ID IS NOT NULL)
        AND antibiotic_name IS NOT NULL
        AND resistance_phenotype IS NOT NULL
    ) TO '{temp_extract}' (FORMAT PARQUET, COMPRESSION ZSTD)
    """
    con.execute(sql, params)
    con.close()

    df = pd.read_parquet(temp_extract)
    if df.empty:
        return pd.DataFrame()
    df["source"] = "EMBL-EBI AMR Portal phenotype parquet"
    df["raw_source_payload"] = ""
    df["source_record_uid"] = df.apply(lambda r: "|".join([
        SOURCE_EBI,
        str(r.get("biosample_accession", "")),
        str(r.get("assembly_accession", "")),
        str(r.get("antibiotic", "")),
        str(r.get("resistant_phenotype", "")),
        str(r.get("measurement_value", "")),
    ]), axis=1)
    out = standardize_amr_columns(df[[
        "biosample_accession", "assembly_accession", "antibiotic", "resistant_phenotype",
        "laboratory_typing_method", "testing_standard", "measurement_value", "evidence",
        "source", "source_record_uid", "raw_source_payload"
    ]], SOURCE_EBI)
    out.to_parquet(cache_file, index=False)
    return out


def normalize_column_map(cols: List[str]) -> dict:
    return {c: re.sub(r"[^a-z0-9]+", "_", c.strip().lower()).strip("_") for c in cols}


def first_existing(df: pd.DataFrame, candidates: List[str], default: str = "") -> pd.Series:
    for c in candidates:
        if c in df.columns:
            return df[c].fillna("").astype(str)
    return pd.Series([default] * len(df), index=df.index)



def load_ncbi_ast_export(path: str, species: str, keep_raw_payload: bool = False) -> pd.DataFrame:
    """
    Robust NCBI AST Browser export loader.

    Notes:
    - Uses stdlib csv parsing because some pandas builds choke on quoted rows.
    - Species matching is intentionally tolerant because NCBI often stores
      strain-level scientific names like "Escherichia coli O157:H7 ..." rather
      than the bare species string.
    - "Isolate" is NOT treated as an assembly accession fallback.
    """
    import csv

    def _clean(x: str) -> str:
        return " ".join(str(x or "").strip().split()).lower()

    def _species_match(row: dict, species_norm: str) -> bool:
        species_norm = _clean(species_norm)

        organism_group = _clean(row.get("organism_group", ""))
        scientific_name = _clean(row.get("scientific_name", ""))
        organism = _clean(row.get("organism", ""))
        taxgroup_name = _clean(row.get("taxgroup_name", ""))
        species_col = _clean(row.get("species", ""))

        candidates = [organism_group, scientific_name, organism, taxgroup_name, species_col]
        candidates = [c for c in candidates if c]

        aliases = {species_norm}

        if species_norm == "escherichia coli":
            aliases.update({
                "e. coli",
                "e.coli",
                "escherichia coli and shigella",
                "e.coli and shigella",
                "escherichia/shigella",
                "escherichia coli/shigella",
            })

        for c in candidates:
            if c in aliases:
                return True
            if c.startswith(species_norm + " "):
                return True
            if species_norm == "escherichia coli":
                if "escherichia coli" in c:
                    return True
                if "escherichia" in c and "shigella" in c:
                    return True

        return False

    with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
        sample = f.read(10000)
    sep = "\t" if sample.count("\t") > sample.count(",") else ","

    species_norm = norm_species_name(species)
    rows = []
    bad_rows = 0

    with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter=sep)
        if not reader.fieldnames:
            return pd.DataFrame()

        # normalize_column_map already returns {original_header: normalized_header}
        # Do NOT zip it again, otherwise each header maps to itself and the loader
        # silently fails to find keys like organism_group / scientific_name.
        field_map = normalize_column_map(reader.fieldnames)

        for raw in reader:
            try:
                row = {field_map.get(k, k): ("" if v is None else str(v)) for k, v in raw.items()}

                if not _species_match(row, species_norm):
                    continue

                biosample_raw = (
                    row.get("biosample")
                    or row.get("biosample_acc")
                    or row.get("biosample_accession")
                    or row.get("_biosample")
                    or row.get("#biosample")
                    or ""
                )

                assembly_raw = (
                    row.get("asm_acc")
                    or row.get("assembly_accession")
                    or row.get("assembly_acc")
                    or ""
                )

                antibiotic = row.get("antibiotic") or row.get("antibiotic_name") or ""
                phenotype = (
                    row.get("resistance_phenotype")
                    or row.get("phenotype")
                    or row.get("result")
                    or row.get("susceptibility_phenotype")
                    or ""
                )
                lab_method = (
                    row.get("laboratory_typing_method")
                    or row.get("laboratory_typing_platform")
                    or row.get("typing_method")
                    or row.get("ast_method")
                    or row.get("method")
                    or ""
                )
                standard = (
                    row.get("testing_standard")
                    or row.get("ast_standard")
                    or row.get("breakpoint_standard")
                    or row.get("standard")
                    or ""
                )
                measurement_sign = row.get("measurement_sign") or row.get("sign") or ""
                mic = row.get("mic_mg_l") or row.get("mic") or row.get("measurement_value") or row.get("measurement") or ""
                disk = row.get("disk_diffusion_mm") or row.get("disk_diffusion") or ""
                measurement = mic if str(mic).strip() else disk

                biosample = pick_first_token(biosample_raw, BS_TOKEN_RE)
                assembly = pick_first_token(assembly_raw, ASM_RE)

                out_row = {
                    "biosample_accession": biosample,
                    "assembly_accession": assembly,
                    "antibiotic": antibiotic,
                    "resistant_phenotype": phenotype,
                    "laboratory_typing_method": lab_method,
                    "testing_standard": standard,
                    "measurement_value": f"{measurement_sign}{measurement}".strip(),
                    "evidence": "NCBI AST Browser export",
                    "source": "NCBI AST Browser export",
                    "raw_source_payload": json.dumps(row, ensure_ascii=False) if keep_raw_payload else "",
                }
                out_row["source_record_uid"] = "|".join([
                    SOURCE_NCBI,
                    str(out_row.get("biosample_accession", "")),
                    str(out_row.get("assembly_accession", "")),
                    str(out_row.get("antibiotic", "")),
                    str(out_row.get("resistant_phenotype", "")),
                    str(out_row.get("measurement_value", "")),
                ])
                rows.append(out_row)

            except Exception:
                bad_rows += 1
                continue

    if not rows:
        try:
            log(f"Loaded NCBI AST Browser export rows | rows=0 | path={path} | bad_rows_skipped={bad_rows}")
        except Exception:
            pass
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    df = standardize_amr_columns(df, SOURCE_NCBI).drop_duplicates().reset_index(drop=True)

    try:
        log(f"Loaded NCBI AST Browser export rows | rows={len(df)} | path={path} | bad_rows_skipped={bad_rows}")
    except Exception:
        pass

    return df



def try_bigquery_import():
    try:
        from google.cloud import bigquery  # type: ignore
        return bigquery
    except Exception:
        return None


def fetch_ncbi_ast_bigquery(species: str, project_id: str, cache_dir: str) -> pd.DataFrame:
    cache_file = os.path.join(cache_dir, f"ncbi_ast_bigquery_{species.lower().replace(' ', '_')}.parquet")
    if os.path.exists(cache_file):
        return pd.read_parquet(cache_file)
    bigquery = try_bigquery_import()
    if bigquery is None:
        raise RuntimeError("google-cloud-bigquery is required for NCBI AST BigQuery integration. Install with: pip install google-cloud-bigquery")
    client = bigquery.Client(project=project_id)
    species_norm = norm_species_name(species)
    words = species_norm.split()
    where_clause = "LOWER(i.taxgroup_name) = LOWER(@species)" if len(words) >= 2 else "STARTS_WITH(LOWER(i.taxgroup_name), LOWER(@species))"
    sql = f"""
    SELECT
      i.biosample_acc AS biosample_accession,
      i.asm_acc AS assembly_accession,
      TO_JSON_STRING(ast_item) AS ast_json
    FROM `ncbi-pathogen-detect.pdbrowser.isolates` AS i,
    UNNEST(i.AST_phenotypes) AS ast_item
    WHERE {where_clause}
    """
    job_config = bigquery.QueryJobConfig(query_parameters=[bigquery.ScalarQueryParameter("species", "STRING", species_norm)])
    df = client.query(sql, job_config=job_config).result().to_dataframe()
    if df.empty:
        return pd.DataFrame()
    rows = []
    for _, r in df.iterrows():
        try:
            ast = json.loads(r.get("ast_json", "{}") or "{}")
        except Exception:
            ast = {}
        rows.append({
            "biosample_accession": str(r.get("biosample_accession", "")),
            "assembly_accession": str(r.get("assembly_accession", "")),
            "antibiotic": str(ast.get("antibiotic", ast.get("antibiotic_name", ""))),
            "resistant_phenotype": str(ast.get("phenotype", ast.get("resistance_phenotype", ast.get("result", "")))),
            "laboratory_typing_method": str(ast.get("laboratory_typing_method", ast.get("method", ""))),
            "testing_standard": str(ast.get("testing_standard", ast.get("ast_standard", ast.get("standard", "")))),
            "measurement_value": " ".join([str(ast.get("measurement_sign", "")).strip(), str(ast.get("measurement", ast.get("mic", ast.get("disk_diffusion", "")))).strip()]).strip(),
            "evidence": "NCBI Pathogen Detection BigQuery isolates.AST_phenotypes",
            "source": "NCBI AST Browser / Pathogen Detection BigQuery",
            "raw_source_payload": "",
        })
    out = pd.DataFrame(rows)
    out["source_record_uid"] = out.apply(lambda r: "|".join([
        SOURCE_NCBI,
        str(r.get("biosample_accession", "")),
        str(r.get("assembly_accession", "")),
        str(r.get("antibiotic", "")),
        str(r.get("resistant_phenotype", "")),
        str(r.get("measurement_value", "")),
    ]), axis=1)
    out = standardize_amr_columns(out, SOURCE_NCBI)
    out.to_parquet(cache_file, index=False)
    return out


def build_bvbrc_raw_amr(genomes: pd.DataFrame, cache_dir: str, stat, warn) -> Tuple[pd.DataFrame, pd.DataFrame]:
    taxon_ids = sorted([int(x) for x in pd.to_numeric(genomes["taxon_id"], errors="coerce").dropna().unique().tolist()])
    stat("BV-BRC AMR stage start", True, taxon_ids=len(taxon_ids), print_every_n=PRINT_EVERY_N_TAXON_IDS)
    amr_fields = [
        "genome_id", "antibiotic", "resistant_phenotype",
        "laboratory_typing_method", "testing_standard", "measurement_value",
        "evidence", "source"
    ]
    taxon_summary_rows = []
    raw_parts = []
    total_taxa = len(taxon_ids)
    for idx, tid in enumerate(taxon_ids, start=1):
        do_print_taxon = idx == 1 or idx == total_taxa or (PRINT_EVERY_N_TAXON_IDS and idx % PRINT_EVERY_N_TAXON_IDS == 0)
        pq = os.path.join(cache_dir, "amr_by_taxon", f"taxon_{tid}_raw.parquet")
        meta = os.path.join(cache_dir, "amr_by_taxon", f"taxon_{tid}_raw.meta.json")
        if os.path.exists(pq) and os.path.exists(meta):
            with open(meta, "r", encoding="utf-8") as f:
                m = json.load(f)
            if m.get("complete") is True:
                df = pd.read_parquet(pq)
                raw_parts.append(df)
                taxon_summary_rows.append({
                    "taxon_id": tid,
                    "taxon_idx": idx,
                    "expected_raw": m.get("expected_raw"),
                    "fetched_raw": m.get("fetched_raw"),
                    "lab_only_kept": m.get("lab_only_kept"),
                    "from_cache": True,
                })
                stat("BV-BRC taxon AMR (cache)", force_print=do_print_taxon, obey_global=False, taxon_id=tid, taxon_idx=idx, taxon_total=total_taxa, lab_kept=int(m.get("lab_only_kept", 0)))
                continue

        stat("BV-BRC taxon AMR (fetch)", force_print=do_print_taxon, obey_global=False, taxon_id=tid, taxon_idx=idx, taxon_total=total_taxa)
        rql_amr = f"eq(taxon_id,{tid})"
        amr, expected = bvbrc_fetch_all("genome_amr", rql_amr, amr_fields, max_rows=MAX_BVBRC_AMR_RECORDS_PER_TAXON, label=f"taxon_id={tid}")
        fetched = len(amr)
        force_text_ids(amr, amr_fields)
        amr["source_record_uid"] = amr.apply(lambda r: "|".join([
            SOURCE_BVBRC,
            str(r.get("genome_id", "")),
            str(r.get("antibiotic", "")),
            str(r.get("resistant_phenotype", "")),
            str(r.get("laboratory_typing_method", "")),
            str(r.get("testing_standard", "")),
            str(r.get("measurement_value", "")),
        ]), axis=1)
        amr["raw_source_payload"] = ""
        amr = standardize_amr_columns(amr, SOURCE_BVBRC)
        raw_parts.append(amr)
        amr.to_parquet(pq, index=False)
        with open(meta, "w", encoding="utf-8") as f:
            json.dump({
                "taxon_id": tid,
                "expected_raw": expected,
                "fetched_raw": fetched,
                "lab_only_kept": int(len(amr)),
                "complete": True,
            }, f)
        taxon_summary_rows.append({
            "taxon_id": tid,
            "taxon_idx": idx,
            "expected_raw": expected,
            "fetched_raw": fetched,
            "lab_only_kept": int(len(amr)),
            "from_cache": False,
        })
        stat("BV-BRC taxon AMR done", force_print=do_print_taxon, obey_global=False, taxon_id=tid, taxon_idx=idx, expected=expected, fetched=fetched, lab_kept=len(amr))
    raw = pd.concat(raw_parts, ignore_index=True) if raw_parts else pd.DataFrame()
    tax = pd.DataFrame(taxon_summary_rows)
    return raw, tax


def summarise_amr_per_genome(ab_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gid, g in ab_df.groupby("genome_id"):
        mp = dict(zip(g["antibiotic_canonical"], g["pheno"]))
        conflicts = sorted([a for a, p in mp.items() if len(p) > 1])
        resistant_any = sorted([a for a, p in mp.items() if "R" in p])
        susceptible_any = sorted([a for a, p in mp.items() if "S" in p])
        intermediate_any = sorted([a for a, p in mp.items() if "I" in p])
        source_list = sorted(set([x for x in g["source_db"].astype(str).tolist() if x]))
        rows.append({
            "genome_id": str(gid),
            "amr_antibiotics_count": len(mp),
            "amr_resistant_count": len(resistant_any),
            "amr_susceptible_count": len(susceptible_any),
            "amr_intermediate_count": len(intermediate_any),
            "amr_conflict_count": len(conflicts),
            "amr_conflict_antibiotics": ";".join(conflicts),
            "antibiotics_resistant": ";".join(resistant_any),
            "antibiotics_susceptible": ";".join(susceptible_any),
            "antibiotics_intermediate": ";".join(intermediate_any),
            "amr_records_total": int(len(g)),
            "amr_sources_count": len(source_list),
            "amr_sources": ";".join(source_list),
            "amr_records_bvbrc": int((g["source_db"] == SOURCE_BVBRC).sum()),
            "amr_records_ncbi_ast": int((g["source_db"] == SOURCE_NCBI).sum()),
            "amr_records_ebi_amr_portal": int((g["source_db"] == SOURCE_EBI).sum()),
        })
    return pd.DataFrame(rows)



def apply_bvbrc_genome_quality_filters(genomes: pd.DataFrame, stat, warn, write_cache_path: Optional[str] = None) -> pd.DataFrame:
    out = genomes.copy()
    if "genome_status" in out.columns:
        gs_before = len(out)
        out["genome_status_norm"] = out["genome_status"].fillna("").astype(str).str.strip().str.lower()
        out = out[out["genome_status_norm"].isin({"complete", "wgs"})].copy()
        stat("Filter genome_status in {Complete,WGS}", True, before=gs_before, after=len(out), removed=gs_before - len(out))
    else:
        warn("genome_status missing -> skipping Complete/WGS filter")

    if "genome_quality" in out.columns:
        gq_before = len(out)
        out["genome_quality_norm"] = out["genome_quality"].fillna("").astype(str).str.strip().str.lower()
        out = out[out["genome_quality_norm"].eq("good")].copy()
        stat("Filter genome_quality == Good", True, before=gq_before, after=len(out), removed=gq_before - len(out))
    elif "genome_quality_flags" in out.columns:
        gqf_before = len(out)
        out["genome_quality_flags_norm"] = out["genome_quality_flags"].fillna("").astype(str).str.strip().str.lower()
        out = out[(out["genome_quality_flags_norm"].eq("")) | (out["genome_quality_flags_norm"].str.contains("good"))].copy()
        stat("Filter genome_quality_flags ~ good/empty", True, before=gqf_before, after=len(out), removed=gqf_before - len(out))
    else:
        warn("genome_quality + flags missing -> skipping Good filter")

    if write_cache_path:
        out.to_parquet(write_cache_path, index=False)
        stat("Saved filtered genomes cache", True, rows=len(out), path=write_cache_path)
    return out


def fetch_species_genomes(species: str, cache_dir: str, stat, warn) -> pd.DataFrame:
    genomes_file = os.path.join(cache_dir, "genomes.parquet")
    filtered_genomes_file = os.path.join(cache_dir, "genomes.filtered.good_complete_wgs.parquet")
    if os.path.exists(filtered_genomes_file):
        genomes = pd.read_parquet(filtered_genomes_file)
        stat("Loaded filtered genomes from cache", True, rows=len(genomes))
        return genomes
    if os.path.exists(genomes_file):
        genomes = pd.read_parquet(genomes_file)
        stat("Loaded genomes from cache", True, rows=len(genomes))
        return apply_bvbrc_genome_quality_filters(genomes, stat, warn, write_cache_path=filtered_genomes_file)

    stat("Fetching BV-BRC genomes by species string (full genome universe)", True, species=species)
    rql_species = f"eq(species,{rql_val(species)})"
    required_fields = ["genome_id", "genome_name", "strain", "taxon_id", "biosample_accession", "bioproject_accession"]
    accession_fields = ["assembly_accession", "genbank_accessions", "refseq_accessions", "sra_accession"]
    optional_fields = ["genome_status", "genome_quality", "genome_quality_flags"]
    try:
        genomes, _ = bvbrc_fetch_all("genome", rql_species, required_fields + optional_fields + accession_fields, max_rows=MAX_BVBRC_GENOMES, label=f"species={species}")
    except Exception as e:
        warn("Genome fetch with accession fields failed; retrying without them", error=str(e)[:200])
        try:
            genomes, _ = bvbrc_fetch_all("genome", rql_species, required_fields + optional_fields, max_rows=MAX_BVBRC_GENOMES, label=f"species={species}")
        except Exception as e2:
            warn("Optional genome fields failed; retrying required fields only", error=str(e2)[:200])
            genomes, _ = bvbrc_fetch_all("genome", rql_species, required_fields, max_rows=MAX_BVBRC_GENOMES, label=f"species={species}")
    if genomes.empty:
        raise ValueError(f"No genomes returned for species='{species}'.")
    force_text_ids(genomes, required_fields + accession_fields)
    genomes["taxon_id"] = pd.to_numeric(genomes["taxon_id"], errors="coerce")
    genomes["biosample_accession"] = genomes["biosample_accession"].apply(lambda x: pick_first_token(x, BS_TOKEN_RE))
    genomes["genbank_assembly_accession"] = genomes["genbank_accessions"].apply(lambda x: pick_first_token(x, ASM_RE))
    genomes["refseq_assembly_accession"] = genomes["refseq_accessions"].apply(lambda x: pick_first_token(x, ASM_RE))
    genomes["assembly_accession_best"] = genomes["assembly_accession"].fillna("").astype(str).str.strip()
    mask_empty = genomes["assembly_accession_best"].eq("")
    genomes.loc[mask_empty, "assembly_accession_best"] = genomes.loc[mask_empty, "refseq_assembly_accession"]
    mask_empty = genomes["assembly_accession_best"].eq("")
    genomes.loc[mask_empty, "assembly_accession_best"] = genomes.loc[mask_empty, "genbank_assembly_accession"]
    genomes["candidate_origin"] = "BV-BRC_genome_universe"
    genomes["candidate_is_external_stub"] = False
    genomes.to_parquet(genomes_file, index=False)
    stat("Saved raw genomes cache", True, rows=len(genomes))
    genomes = apply_bvbrc_genome_quality_filters(genomes, stat, warn, write_cache_path=filtered_genomes_file)
    return genomes


def source_union_label(vals: List[str]) -> str:
    vals = sorted(set([str(v) for v in vals if str(v).strip()]))
    return ";".join(vals)


def build_external_stub_candidates(genomes: pd.DataFrame, external_union: pd.DataFrame, species: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if external_union.empty:
        return pd.DataFrame(), pd.DataFrame()

    g_bs = set([x for x in genomes["biosample_accession"].astype(str).tolist() if x])
    g_asm = set([x for x in genomes["assembly_accession_best"].astype(str).tolist() if x])

    ext = external_union.copy()
    ext["biosample_accession"] = ext["biosample_accession"].apply(lambda x: pick_first_token(x, BS_TOKEN_RE))
    ext["assembly_accession"] = ext["assembly_accession"].apply(lambda x: pick_first_token(x, ASM_RE))
    ext["matches_existing_genome"] = ext["biosample_accession"].isin(g_bs) | ext["assembly_accession"].isin(g_asm)

    unmatched = ext[~ext["matches_existing_genome"]].copy()
    if unmatched.empty:
        return pd.DataFrame(), ext

    def key_for_row(r):
        bs = str(r.get("biosample_accession", "")).strip()
        asm = str(r.get("assembly_accession", "")).strip()
        if bs:
            return f"BS::{bs}"
        if asm:
            return f"ASM::{asm}"
        return ""

    unmatched["stub_key"] = unmatched.apply(key_for_row, axis=1)
    unmatched = unmatched[unmatched["stub_key"] != ""].copy()
    if unmatched.empty:
        return pd.DataFrame(), ext

    rows = []
    for stub_key, g in unmatched.groupby("stub_key"):
        bs_vals = [x for x in g["biosample_accession"].astype(str).tolist() if x]
        asm_vals = [x for x in g["assembly_accession"].astype(str).tolist() if x]
        biosample = sorted(set(bs_vals))[0] if bs_vals else ""
        assembly = sorted(set(asm_vals))[0] if asm_vals else ""
        srcs = sorted(set(g["source_db"].astype(str).tolist()))
        if biosample:
            synthetic_id = f"EXT::{biosample}"
            display_name = f"External candidate {biosample}"
        elif assembly:
            synthetic_id = f"EXT::{assembly}"
            display_name = f"External candidate {assembly}"
        else:
            digest = hashlib.md5(stub_key.encode("utf-8")).hexdigest()[:12]
            synthetic_id = f"EXT::{digest}"
            display_name = f"External candidate {digest}"
        rows.append({
            "genome_id": synthetic_id,
            "genome_name": display_name,
            "strain": "",
            "taxon_id": "",
            "biosample_accession": biosample,
            "bioproject_accession": "",
            "assembly_accession": assembly,
            "genbank_accessions": "",
            "refseq_accessions": "",
            "sra_accession": "",
            "genbank_assembly_accession": "",
            "refseq_assembly_accession": "",
            "assembly_accession_best": assembly,
            "candidate_origin": "external_only_" + "+".join([s.lower() for s in srcs]),
            "candidate_is_external_stub": True,
            "candidate_source_db_union": ";".join(srcs),
            "species": species,
        })
    stub_df = pd.DataFrame(rows)
    return stub_df, ext


def map_external_amr_to_candidates(candidates: pd.DataFrame, ext_df: pd.DataFrame) -> pd.DataFrame:
    if ext_df.empty:
        return ext_df.copy()
    c_bs = candidates[["genome_id", "biosample_accession"]].drop_duplicates().copy()
    c_asm = candidates[["genome_id", "assembly_accession_best"]].drop_duplicates().rename(columns={"assembly_accession_best": "assembly_accession"}).copy()
    ext = ext_df.drop(columns=["genome_id"], errors="ignore").copy()
    ext["biosample_accession"] = ext["biosample_accession"].apply(lambda x: pick_first_token(x, BS_TOKEN_RE))
    ext["assembly_accession"] = ext["assembly_accession"].apply(lambda x: pick_first_token(x, ASM_RE))
    m1 = ext.merge(c_bs, on="biosample_accession", how="left")
    matched_bs = m1[m1["genome_id"].fillna("") != ""].copy()
    unmatched = m1[m1["genome_id"].fillna("") == ""].drop(columns=["genome_id"], errors="ignore").copy()
    m2 = unmatched.merge(c_asm, on="assembly_accession", how="left")
    matched_asm = m2[m2["genome_id"].fillna("") != ""].copy()
    out = pd.concat([matched_bs, matched_asm], ignore_index=True)
    out = out.drop_duplicates([
        "genome_id", "biosample_accession", "assembly_accession", "antibiotic_canonical", "pheno",
        "laboratory_typing_method", "testing_standard", "measurement_value", "source_db"
    ]).reset_index(drop=True)
    return out


def build_multisource_candidates(species: str, genomes: pd.DataFrame, bvbrc_raw: pd.DataFrame,
                                 enable_ebi: bool, enable_ncbi: bool,
                                 ebi_release: str, ncbi_ast_export: str, ncbi_gcp_project: str,
                                 cache_dir: str, stat, warn):
    source_summary_rows = []
    raw_parts = []

    # BV-BRC raw is already genome-mapped.
    raw_parts.append(bvbrc_raw)
    source_summary_rows.append({
        "source_db": SOURCE_BVBRC,
        "raw_rows": int(len(bvbrc_raw)),
        "mapped_rows": int(len(bvbrc_raw)),
        "discovery_role": "seed+union",
    })

    external_raw_parts = []

    if enable_ebi:
        try:
            ebi_raw = fetch_ebi_phenotypes(species, cache_dir=cache_dir, release=ebi_release)
            stat("Loaded EBI phenotype rows", True, rows=len(ebi_raw))
            external_raw_parts.append(ebi_raw)
        except Exception as e:
            warn("EMBL-EBI AMR Portal integration failed", error=str(e)[:300])

    if enable_ncbi:
        ncbi_ok = False
        if ncbi_ast_export:
            try:
                ncbi_raw = load_ncbi_ast_export(ncbi_ast_export, species)
                stat("Loaded NCBI AST Browser export rows", True, rows=len(ncbi_raw), path=ncbi_ast_export)
                external_raw_parts.append(ncbi_raw)
                ncbi_ok = True
            except Exception as e:
                warn("NCBI AST Browser export integration failed", error=str(e)[:300])
        if (not ncbi_ok) and ncbi_gcp_project:
            try:
                ncbi_raw = fetch_ncbi_ast_bigquery(species, ncbi_gcp_project, cache_dir=cache_dir)
                stat("Loaded NCBI AST BigQuery rows", True, rows=len(ncbi_raw), project=ncbi_gcp_project)
                external_raw_parts.append(ncbi_raw)
                ncbi_ok = True
            except Exception as e:
                warn("NCBI AST Browser BigQuery integration failed", error=str(e)[:300])
        if not ncbi_ok and not ncbi_ast_export and not ncbi_gcp_project:
            warn("NCBI AST integration enabled but no source succeeded", hint="Provide --ncbi-ast-export <AST Browser csv/tsv> OR --ncbi-gcp-project <your_gcp_project>")

    external_union = pd.concat(external_raw_parts, ignore_index=True) if external_raw_parts else pd.DataFrame()
    if not external_union.empty:
        external_union = external_union.drop_duplicates([
            "biosample_accession", "assembly_accession", "antibiotic_canonical", "pheno",
            "laboratory_typing_method", "testing_standard", "measurement_value", "source_db"
        ]).reset_index(drop=True)

    stub_df, external_annot = build_external_stub_candidates(genomes, external_union, species)
    candidate_table = pd.concat([genomes.copy(), stub_df], ignore_index=True, sort=False) if not stub_df.empty else genomes.copy()
    candidate_table["species"] = species
    if "candidate_source_db_union" not in candidate_table.columns:
        candidate_table["candidate_source_db_union"] = ""

    mapped_external = map_external_amr_to_candidates(candidate_table, external_union)

    if not external_union.empty:
        for src, sdf in external_union.groupby("source_db"):
            mapped_count = int((mapped_external["source_db"] == src).sum()) if not mapped_external.empty else 0
            unmatched_src = sdf.copy()
            if not mapped_external.empty:
                unmatched_src = sdf.merge(
                    mapped_external[["source_record_uid"]].drop_duplicates(),
                    on="source_record_uid", how="left", indicator=True
                )
                unmatched_src = unmatched_src[unmatched_src["_merge"] == "left_only"].drop(columns=["_merge"])
            source_summary_rows.append({
                "source_db": src,
                "raw_rows": int(len(sdf)),
                "mapped_rows": mapped_count,
                "external_only_rows_unmapped": int(len(unmatched_src)),
                "external_stub_candidates": int(len(stub_df[stub_df["candidate_source_db_union"].fillna("").str.contains(src, regex=False)])) if not stub_df.empty else 0,
                "discovery_role": "union_only",
            })

    raw_merged = pd.concat([bvbrc_raw, mapped_external], ignore_index=True) if not mapped_external.empty else bvbrc_raw.copy()
    if not raw_merged.empty:
        raw_merged = raw_merged.drop_duplicates([
            "genome_id", "biosample_accession", "assembly_accession", "antibiotic_canonical", "pheno",
            "laboratory_typing_method", "testing_standard", "measurement_value", "source_db"
        ]).reset_index(drop=True)

    return candidate_table, raw_merged, pd.DataFrame(source_summary_rows), stub_df



def source_slug_v2(source_db: str) -> str:
    return {
        SOURCE_BVBRC: "BVBRC",
        SOURCE_NCBI: "NCBI_AST",
        SOURCE_EBI: "EBI_AMR",
    }.get(source_db, re.sub(r"[^A-Za-z0-9]+", "_", str(source_db)).strip("_") or "SOURCE")


def make_isolate_key_v2(biosample_accession: str = "", assembly_accession: str = "", genome_id: str = "") -> str:
    bs = pick_first_token(biosample_accession, BS_TOKEN_RE)
    asm = pick_first_token(assembly_accession, ASM_RE)
    gid = str(genome_id or "").strip()
    if bs:
        return f"BS::{bs}"
    if asm:
        return f"ASM::{asm}"
    if gid:
        return f"GENOME::{gid}"
    return ""


def assembly_version_rank_v2(acc: str) -> int:
    acc = str(acc or "").strip()
    m = re.search(r"\.(\d+)$", acc)
    return int(m.group(1)) if m else -1


def build_canonical_genome_universe_v2(genomes: pd.DataFrame, species: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    g = genomes.copy()
    required = [
        "genome_id", "genome_name", "strain", "taxon_id", "biosample_accession", "bioproject_accession",
        "assembly_accession", "genbank_accessions", "refseq_accessions", "sra_accession",
        "genbank_assembly_accession", "refseq_assembly_accession", "assembly_accession_best",
    ]
    for c in required:
        if c not in g.columns:
            g[c] = ""
    force_text_ids(g, required)
    g["biosample_accession"] = g["biosample_accession"].apply(lambda x: pick_first_token(x, BS_TOKEN_RE))
    for c in ["assembly_accession", "genbank_assembly_accession", "refseq_assembly_accession", "assembly_accession_best"]:
        g[c] = g[c].apply(lambda x: pick_first_token(x, ASM_RE))
    mask_empty = g["assembly_accession_best"].eq("")
    g.loc[mask_empty, "assembly_accession_best"] = g.loc[mask_empty, "assembly_accession"]
    mask_empty = g["assembly_accession_best"].eq("")
    g.loc[mask_empty, "assembly_accession_best"] = g.loc[mask_empty, "refseq_assembly_accession"]
    mask_empty = g["assembly_accession_best"].eq("")
    g.loc[mask_empty, "assembly_accession_best"] = g.loc[mask_empty, "genbank_assembly_accession"]

    g["isolate_key"] = g.apply(
        lambda r: make_isolate_key_v2(
            biosample_accession=r.get("biosample_accession", ""),
            assembly_accession=r.get("assembly_accession_best", ""),
            genome_id=r.get("genome_id", ""),
        ),
        axis=1,
    )
    g["has_assembly_v2"] = g["assembly_accession_best"].ne("")
    g["refseq_pref_v2"] = g["assembly_accession_best"].str.startswith("GCF_", na=False)
    g["assembly_version_rank_v2"] = g["assembly_accession_best"].apply(assembly_version_rank_v2)

    grp = g.groupby("isolate_key", dropna=False).agg(
        same_isolate_genome_rows=("genome_id", "size"),
        same_isolate_genome_ids=("genome_id", lambda s: ";".join(sorted(set([str(x) for x in s if str(x).strip()])))),
        same_isolate_assemblies=("assembly_accession_best", lambda s: ";".join(sorted(set([str(x) for x in s if str(x).strip()])))),
    ).reset_index()

    g_sorted = g.sort_values(
        ["isolate_key", "has_assembly_v2", "refseq_pref_v2", "assembly_version_rank_v2", "genome_id"],
        ascending=[True, False, False, False, True],
    ).copy()

    canonical = g_sorted.drop_duplicates("isolate_key", keep="first").copy()
    canonical = canonical.merge(grp, on="isolate_key", how="left")
    canonical["candidate_origin"] = "BV-BRC_canonical_genome"
    canonical["candidate_is_external_stub"] = False
    canonical["candidate_source_db_union"] = ""
    canonical["species"] = species
    canonical["canonical_pick_reason"] = canonical.apply(
        lambda r: "best assembly accession/version within isolate"
        if int(r.get("same_isolate_genome_rows", 1) or 1) > 1 else "single genome row for isolate",
        axis=1,
    )

    return g_sorted, canonical


def build_source_lookup_v2(all_genomes_with_keys: pd.DataFrame, canonical_genomes: pd.DataFrame) -> dict:
    genome_to_key = dict(zip(all_genomes_with_keys["genome_id"].astype(str), all_genomes_with_keys["isolate_key"].astype(str)))
    bs_map = (
        all_genomes_with_keys[["biosample_accession", "isolate_key"]]
        .copy()
    )
    bs_map = bs_map[bs_map["biosample_accession"].astype(str) != ""].drop_duplicates("biosample_accession")
    biosample_to_key = dict(zip(bs_map["biosample_accession"].astype(str), bs_map["isolate_key"].astype(str)))

    asm_map = (
        all_genomes_with_keys[["assembly_accession_best", "isolate_key"]]
        .copy()
        .rename(columns={"assembly_accession_best": "assembly_accession"})
    )
    asm_map = asm_map[asm_map["assembly_accession"].astype(str) != ""].drop_duplicates("assembly_accession")
    assembly_to_key = dict(zip(asm_map["assembly_accession"].astype(str), asm_map["isolate_key"].astype(str)))

    key_to_candidate_id = dict(zip(canonical_genomes["isolate_key"].astype(str), canonical_genomes["genome_id"].astype(str)))
    return {
        "genome_to_key": genome_to_key,
        "biosample_to_key": biosample_to_key,
        "assembly_to_key": assembly_to_key,
        "existing_keys": set([str(x) for x in canonical_genomes["isolate_key"].astype(str).tolist() if str(x).strip()]),
        "key_to_candidate_id": key_to_candidate_id,
    }


def resolve_external_isolate_key_v2(biosample_accession: str, assembly_accession: str, lookup: dict) -> str:
    bs = pick_first_token(biosample_accession, BS_TOKEN_RE)
    asm = pick_first_token(assembly_accession, ASM_RE)
    if bs and bs in lookup["biosample_to_key"]:
        return lookup["biosample_to_key"][bs]
    if asm and asm in lookup["assembly_to_key"]:
        return lookup["assembly_to_key"][asm]
    if bs:
        return f"BS::{bs}"
    if asm:
        return f"ASM::{asm}"
    return ""


def annotate_source_raw_with_isolates_v2(raw_df: pd.DataFrame, source_db: str, lookup: dict) -> pd.DataFrame:
    cols = [
        "genome_id", "biosample_accession", "assembly_accession", "antibiotic", "antibiotic_canonical",
        "resistant_phenotype", "pheno", "laboratory_typing_method", "testing_standard", "measurement_value",
        "evidence", "source", "source_db", "source_record_uid", "raw_source_payload",
        "original_source_genome_id", "isolate_key", "mapped_to_existing_bvbrc_genome",
    ]
    if raw_df is None or raw_df.empty:
        return pd.DataFrame(columns=cols)

    out = raw_df.copy()
    for c in cols:
        if c not in out.columns:
            out[c] = ""
    force_text_ids(out, [
        "genome_id", "biosample_accession", "assembly_accession", "source_record_uid",
        "original_source_genome_id",
    ])
    out["biosample_accession"] = out["biosample_accession"].apply(lambda x: pick_first_token(x, BS_TOKEN_RE))
    out["assembly_accession"] = out["assembly_accession"].apply(lambda x: pick_first_token(x, ASM_RE))
    if source_db == SOURCE_BVBRC:
        out["original_source_genome_id"] = out["genome_id"].astype(str)
        out["isolate_key"] = out["genome_id"].astype(str).map(lookup["genome_to_key"]).fillna("")
        mask = out["isolate_key"].eq("")
        if mask.any():
            out.loc[mask, "isolate_key"] = out.loc[mask].apply(
                lambda r: make_isolate_key_v2(
                    biosample_accession=r.get("biosample_accession", ""),
                    assembly_accession=r.get("assembly_accession", ""),
                    genome_id=r.get("genome_id", ""),
                ),
                axis=1,
            )
    else:
        out["original_source_genome_id"] = ""
        out["isolate_key"] = out.apply(
            lambda r: resolve_external_isolate_key_v2(
                biosample_accession=r.get("biosample_accession", ""),
                assembly_accession=r.get("assembly_accession", ""),
                lookup=lookup,
            ),
            axis=1,
        )

    out = out[out["isolate_key"].astype(str) != ""].copy()
    out["mapped_to_existing_bvbrc_genome"] = out["isolate_key"].isin(lookup["existing_keys"])
    dedup_cols = [
        "isolate_key", "biosample_accession", "assembly_accession", "antibiotic_canonical", "pheno",
        "laboratory_typing_method", "testing_standard", "measurement_value", "source_db",
    ]
    out = out.drop_duplicates(dedup_cols).reset_index(drop=True)
    return out[cols].copy()


def build_external_stub_candidates_v2(external_annotated_frames: List[pd.DataFrame], existing_keys: set, species: str) -> pd.DataFrame:
    if not external_annotated_frames:
        return pd.DataFrame()
    ext = pd.concat([df for df in external_annotated_frames if df is not None and not df.empty], ignore_index=True)
    if ext.empty:
        return pd.DataFrame()

    ext = ext[~ext["isolate_key"].isin(existing_keys)].copy()
    if ext.empty:
        return pd.DataFrame()

    rows = []
    for isolate_key, g in ext.groupby("isolate_key"):
        biosamples = sorted(set([x for x in g["biosample_accession"].astype(str).tolist() if x]))
        assemblies = sorted(set([x for x in g["assembly_accession"].astype(str).tolist() if x]))
        srcs = sorted(set([x for x in g["source_db"].astype(str).tolist() if x]))
        biosample = biosamples[0] if biosamples else ""
        assembly = assemblies[0] if assemblies else ""

        if not biosample and str(isolate_key).startswith("BS::"):
            biosample = str(isolate_key).split("BS::", 1)[1]
        if not assembly and str(isolate_key).startswith("ASM::"):
            assembly = str(isolate_key).split("ASM::", 1)[1]

        if biosample:
            synthetic_id = f"EXT::{biosample}"
            display_name = f"External isolate {biosample}"
        elif assembly:
            synthetic_id = f"EXT::{assembly}"
            display_name = f"External isolate {assembly}"
        else:
            digest = hashlib.md5(str(isolate_key).encode("utf-8")).hexdigest()[:12]
            synthetic_id = f"EXT::{digest}"
            display_name = f"External isolate {digest}"

        rows.append({
            "genome_id": synthetic_id,
            "genome_name": display_name,
            "strain": "",
            "taxon_id": "",
            "biosample_accession": biosample,
            "bioproject_accession": "",
            "assembly_accession": assembly,
            "genbank_accessions": "",
            "refseq_accessions": "",
            "sra_accession": "",
            "genbank_assembly_accession": assembly,
            "refseq_assembly_accession": assembly if str(assembly).startswith("GCF_") else "",
            "assembly_accession_best": assembly,
            "isolate_key": isolate_key,
            "same_isolate_genome_rows": 0,
            "same_isolate_genome_ids": "",
            "same_isolate_assemblies": assembly,
            "candidate_origin": "External_only_stub",
            "candidate_is_external_stub": True,
            "candidate_source_db_union": ";".join(srcs),
            "species": species,
            "canonical_pick_reason": "external isolate stub (no BV-BRC genome row for isolate)",
        })
    return pd.DataFrame(rows)


def assign_candidate_ids_v2(raw_iso_df: pd.DataFrame, candidate_universe: pd.DataFrame) -> pd.DataFrame:
    if raw_iso_df is None or raw_iso_df.empty:
        out = raw_iso_df.copy() if raw_iso_df is not None else pd.DataFrame()
        if "genome_id" not in out.columns:
            out["genome_id"] = ""
        return out
    key_to_candidate_id = dict(zip(candidate_universe["isolate_key"].astype(str), candidate_universe["genome_id"].astype(str)))
    cand_meta = candidate_universe[["genome_id", "candidate_origin", "candidate_is_external_stub"]].drop_duplicates().copy()
    out = raw_iso_df.copy()
    out["genome_id"] = out["isolate_key"].astype(str).map(key_to_candidate_id).fillna("")
    out = out[out["genome_id"].astype(str) != ""].copy()
    out = out.drop_duplicates([
        "genome_id", "isolate_key",
        "antibiotic_canonical", "pheno", "laboratory_typing_method", "testing_standard",
        "measurement_value", "source_db"
    ]).reset_index(drop=True)
    out = out.merge(cand_meta, on="genome_id", how="left")
    return out


def build_presence_flags_v2(source_raw_map: dict) -> pd.DataFrame:
    ids = set()
    for df in source_raw_map.values():
        if df is not None and not df.empty:
            ids.update([str(x) for x in df["genome_id"].astype(str).tolist() if str(x).strip()])
    if not ids:
        return pd.DataFrame(columns=["genome_id", "present_in_bvbrc_amr", "present_in_ncbi_ast", "present_in_ebi_amr", "candidate_source_db_union"])

    base = pd.DataFrame({"genome_id": sorted(ids)})
    for src, col in [
        (SOURCE_BVBRC, "present_in_bvbrc_amr"),
        (SOURCE_NCBI, "present_in_ncbi_ast"),
        (SOURCE_EBI, "present_in_ebi_amr"),
    ]:
        src_ids = set()
        df = source_raw_map.get(src)
        if df is not None and not df.empty:
            src_ids = set([str(x) for x in df["genome_id"].astype(str).tolist() if str(x).strip()])
        base[col] = base["genome_id"].astype(str).isin(src_ids)

    def mk_union(r):
        vals = []
        if bool(r.get("present_in_bvbrc_amr", False)):
            vals.append(SOURCE_BVBRC)
        if bool(r.get("present_in_ncbi_ast", False)):
            vals.append(SOURCE_NCBI)
        if bool(r.get("present_in_ebi_amr", False)):
            vals.append(SOURCE_EBI)
        return ";".join(vals)

    base["candidate_source_db_union"] = base.apply(mk_union, axis=1)
    return base


OUTPUT_COLS_V2 = [
    "genome_id", "genome_name", "strain", "taxon_id", "biosample_accession", "bioproject_accession",
    "assembly_accession", "genbank_accessions", "refseq_accessions", "sra_accession",
    "genbank_assembly_accession", "refseq_assembly_accession", "assembly_accession_best",
    "isolate_key", "same_isolate_genome_rows", "same_isolate_genome_ids", "same_isolate_assemblies",
    "candidate_origin", "candidate_is_external_stub", "candidate_source_db_union",
    "present_in_bvbrc_amr", "present_in_ncbi_ast", "present_in_ebi_amr",
    "genome_available", "rnaseq_available", "final_conservative_hit", "is_extra_over_bvbrc_amr_baseline",
    "evidence_tier", "rnaseq_run_count", "rnaseq_run_ids", "rnaseq_run_ids_sra", "rnaseq_run_ids_ena",
    "ena_fastq_ftp", "ena_fastq_md5", "access_method",
    "amr_antibiotics_count", "amr_resistant_count", "amr_susceptible_count", "amr_intermediate_count",
    "amr_conflict_count", "amr_conflict_antibiotics", "antibiotics_resistant", "antibiotics_susceptible",
    "antibiotics_intermediate", "amr_records_total", "amr_sources_count", "amr_sources",
    "amr_records_bvbrc", "amr_records_ncbi_ast", "amr_records_ebi_amr_portal"
]


def build_candidate_output_table_v2(candidate_subset: pd.DataFrame,
                                    mapped_raw: pd.DataFrame,
                                    bs_map: pd.DataFrame,
                                    baseline_isolate_keys: set,
                                    presence_flags: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    out = candidate_subset.copy()
    if out.empty:
        out = pd.DataFrame(columns=OUTPUT_COLS_V2)
        return out

    out = out.merge(bs_map, on="biosample_accession", how="left")
    out["rnaseq_available"] = out["rnaseq_run_ids"].fillna("").astype(str).str.contains(RUN_RE)
    out["evidence_tier"] = out["rnaseq_available"].map(lambda v: "Strong (same BioSample: SRA+ENA)" if v else "None")
    out["access_method"] = out["rnaseq_available"].map(lambda v: "Direct BioSample match → SRA Toolkit (prefetch/fasterq-dump) or ENA FTP FASTQ" if v else "")
    out["rnaseq_run_count"] = out["rnaseq_run_ids"].fillna("").apply(lambda s: len(extract_run_ids(s)))
    out["genome_available"] = (~out["candidate_is_external_stub"].fillna(False)) | out["assembly_accession_best"].fillna("").astype(str).ne("")

    if mapped_raw is None or mapped_raw.empty:
        amr_sum = pd.DataFrame(columns=[
            "genome_id", "amr_antibiotics_count", "amr_resistant_count", "amr_susceptible_count",
            "amr_intermediate_count", "amr_conflict_count", "amr_conflict_antibiotics",
            "antibiotics_resistant", "antibiotics_susceptible", "antibiotics_intermediate",
            "amr_records_total", "amr_sources_count", "amr_sources",
            "amr_records_bvbrc", "amr_records_ncbi_ast", "amr_records_ebi_amr_portal"
        ])
    else:
        amr_sum = summarise_amr_per_genome(
            mapped_raw[["genome_id", "antibiotic_canonical", "pheno", "source_db"]].drop_duplicates().copy()
        )
    out = out.merge(amr_sum, on="genome_id", how="left")

    if presence_flags is None:
        presence_flags = build_presence_flags_v2({})
    out = out.merge(presence_flags, on="genome_id", how="left", suffixes=("", "_pf"))
    for c in ["present_in_bvbrc_amr", "present_in_ncbi_ast", "present_in_ebi_amr"]:
        if c not in out.columns:
            out[c] = False
        out[c] = out[c].fillna(False)

    if "candidate_source_db_union_pf" in out.columns:
        out["candidate_source_db_union"] = out["candidate_source_db_union_pf"].fillna(out["candidate_source_db_union"])
        out = out.drop(columns=["candidate_source_db_union_pf"])

    out["final_conservative_hit"] = out["genome_available"].fillna(False) & out["rnaseq_available"].fillna(False) & out["amr_records_total"].fillna(0).astype(float).gt(0)
    out["is_extra_over_bvbrc_amr_baseline"] = ~out["isolate_key"].astype(str).isin(set([str(x) for x in baseline_isolate_keys if str(x).strip()]))

    for c in OUTPUT_COLS_V2:
        if c not in out.columns:
            out[c] = ""
    out = out[OUTPUT_COLS_V2].copy()
    force_text_ids(out, [
        "genome_id", "biosample_accession", "bioproject_accession", "assembly_accession_best", "assembly_accession",
        "genbank_assembly_accession", "refseq_assembly_accession", "isolate_key", "same_isolate_genome_ids", "same_isolate_assemblies"
    ])
    return out


def source_discovery_summary_rows_v2(source_raw_map: dict,
                                     source_tables: dict,
                                     union_table: pd.DataFrame,
                                     baseline_isolate_keys: set,
                                     candidate_universe: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for src in [SOURCE_BVBRC, SOURCE_NCBI, SOURCE_EBI]:
        raw_df = source_raw_map.get(src, pd.DataFrame())
        tab = source_tables.get(src, pd.DataFrame())
        isolate_keys = set([str(x) for x in raw_df.get("isolate_key", pd.Series(dtype=str)).astype(str).tolist() if str(x).strip()]) if raw_df is not None and not raw_df.empty else set()
        biosamples = set([str(x) for x in raw_df.get("biosample_accession", pd.Series(dtype=str)).astype(str).tolist() if str(x).strip()]) if raw_df is not None and not raw_df.empty else set()
        rows.append({
            "source_db": src,
            "raw_phenotype_rows": int(len(raw_df)) if raw_df is not None else 0,
            "unique_isolates_with_ast": int(len(isolate_keys)),
            "unique_biosamples_with_ast": int(len(biosamples)),
            "source_validated_rows": int(len(tab)) if tab is not None else 0,
            "source_rnaseq_true": int(tab["rnaseq_available"].sum()) if tab is not None and not tab.empty else 0,
            "source_final_conservative_hits": int(tab["final_conservative_hit"].sum()) if tab is not None and not tab.empty else 0,
            "new_isolates_over_bvbrc_baseline": int(len([k for k in isolate_keys if k not in baseline_isolate_keys])),
        })
    union_keys = set([str(x) for x in union_table.get("isolate_key", pd.Series(dtype=str)).astype(str).tolist() if str(x).strip()]) if union_table is not None and not union_table.empty else set()
    rows.append({
        "source_db": "FINAL_UNION",
        "raw_phenotype_rows": int(sum([(len(df) if df is not None else 0) for df in source_raw_map.values()])),
        "unique_isolates_with_ast": int(len(union_keys)),
        "unique_biosamples_with_ast": int(union_table["biosample_accession"].astype(str).replace("", pd.NA).dropna().nunique()) if union_table is not None and not union_table.empty else 0,
        "source_validated_rows": int(len(union_table)) if union_table is not None else 0,
        "source_rnaseq_true": int(union_table["rnaseq_available"].sum()) if union_table is not None and not union_table.empty else 0,
        "source_final_conservative_hits": int(union_table["final_conservative_hit"].sum()) if union_table is not None and not union_table.empty else 0,
        "new_isolates_over_bvbrc_baseline": int(len([k for k in union_keys if k not in baseline_isolate_keys])),
    })
    return pd.DataFrame(rows)


def build_multisource_source_tables_v2(species: str,
                                       genomes: pd.DataFrame,
                                       bvbrc_raw: pd.DataFrame,
                                       enable_ebi: bool,
                                       enable_ncbi: bool,
                                       ebi_release: str,
                                       ncbi_ast_export: str,
                                       ncbi_gcp_project: str,
                                       cache_dir: str,
                                       stat,
                                       warn) -> tuple:
    all_genomes_with_keys, canonical_genomes = build_canonical_genome_universe_v2(genomes, species)
    lookup = build_source_lookup_v2(all_genomes_with_keys, canonical_genomes)
    baseline_isolate_keys = set()

    source_raw_original = {SOURCE_BVBRC: bvbrc_raw.copy() if bvbrc_raw is not None else pd.DataFrame()}
    external_source_frames = []

    if enable_ebi:
        try:
            ebi_raw = fetch_ebi_phenotypes(species, cache_dir=cache_dir, release=ebi_release)
            stat("Loaded EBI phenotype rows", True, rows=len(ebi_raw))
            source_raw_original[SOURCE_EBI] = ebi_raw
            external_source_frames.append((SOURCE_EBI, ebi_raw))
        except Exception as e:
            warn("EMBL-EBI AMR Portal integration failed", error=str(e)[:300])
            source_raw_original[SOURCE_EBI] = pd.DataFrame()

    if enable_ncbi:
        ncbi_ok = False
        if ncbi_ast_export:
            try:
                ncbi_raw = load_ncbi_ast_export(ncbi_ast_export, species)
                stat("Loaded NCBI AST Browser export rows", True, rows=len(ncbi_raw), path=ncbi_ast_export)
                source_raw_original[SOURCE_NCBI] = ncbi_raw
                external_source_frames.append((SOURCE_NCBI, ncbi_raw))
                ncbi_ok = True
            except Exception as e:
                warn("NCBI AST Browser export integration failed", error=str(e)[:300])
        if (not ncbi_ok) and ncbi_gcp_project:
            try:
                ncbi_raw = fetch_ncbi_ast_bigquery(species, ncbi_gcp_project, cache_dir=cache_dir)
                stat("Loaded NCBI AST BigQuery rows", True, rows=len(ncbi_raw), project=ncbi_gcp_project)
                source_raw_original[SOURCE_NCBI] = ncbi_raw
                external_source_frames.append((SOURCE_NCBI, ncbi_raw))
                ncbi_ok = True
            except Exception as e:
                warn("NCBI AST Browser BigQuery integration failed", error=str(e)[:300])
        if not ncbi_ok:
            if not ncbi_ast_export and not ncbi_gcp_project:
                warn("NCBI AST integration enabled but no source succeeded", hint="Provide --ncbi-ast-export <AST Browser csv/tsv> OR --ncbi-gcp-project <your_gcp_project>")
            source_raw_original[SOURCE_NCBI] = source_raw_original.get(SOURCE_NCBI, pd.DataFrame())
    else:
        source_raw_original[SOURCE_NCBI] = pd.DataFrame()

    annotated_source_raw = {}
    annotated_source_raw[SOURCE_BVBRC] = annotate_source_raw_with_isolates_v2(source_raw_original.get(SOURCE_BVBRC, pd.DataFrame()), SOURCE_BVBRC, lookup)
    baseline_isolate_keys = set([str(x) for x in annotated_source_raw[SOURCE_BVBRC].get("isolate_key", pd.Series(dtype=str)).astype(str).tolist() if str(x).strip()])
    external_annotated_frames = []
    for src, raw_df in external_source_frames:
        ann = annotate_source_raw_with_isolates_v2(raw_df, src, lookup)
        annotated_source_raw[src] = ann
        external_annotated_frames.append(ann)
    for src in [SOURCE_NCBI, SOURCE_EBI]:
        if src not in annotated_source_raw:
            annotated_source_raw[src] = pd.DataFrame()

    stub_df = build_external_stub_candidates_v2(external_annotated_frames, lookup["existing_keys"], species)
    candidate_universe = canonical_genomes.copy()
    if not stub_df.empty:
        candidate_universe = pd.concat([candidate_universe, stub_df], ignore_index=True, sort=False)
    candidate_universe["species"] = species
    candidate_universe = candidate_universe.drop_duplicates("isolate_key", keep="first").copy()

    key_to_candidate_id = dict(zip(candidate_universe["isolate_key"].astype(str), candidate_universe["genome_id"].astype(str)))
    lookup["key_to_candidate_id"] = key_to_candidate_id

    mapped_source_raw = {}
    for src in [SOURCE_BVBRC, SOURCE_NCBI, SOURCE_EBI]:
        mapped_source_raw[src] = assign_candidate_ids_v2(annotated_source_raw.get(src, pd.DataFrame()), candidate_universe)

    raw_union = pd.concat([df for df in mapped_source_raw.values() if df is not None and not df.empty], ignore_index=True) if any(df is not None and not df.empty for df in mapped_source_raw.values()) else pd.DataFrame()
    if not raw_union.empty:
        raw_union = raw_union.drop_duplicates([
            "genome_id", "isolate_key", "biosample_accession", "assembly_accession",
            "antibiotic_canonical", "pheno", "laboratory_typing_method", "testing_standard",
            "measurement_value", "source_db"
        ]).reset_index(drop=True)

    return candidate_universe, mapped_source_raw, raw_union, stub_df, baseline_isolate_keys


def run_one_species(species_input: str, cache_root: str, out_root: str,
                    enable_ebi: bool, enable_ncbi: bool,
                    ebi_release: str, ncbi_ast_export: str, ncbi_gcp_project: str):
    species = norm_species_name(species_input)
    stats = []

    def stat(msg, force_print=False, obey_global=True, **kwargs):
        row = {"message": msg}
        row.update(kwargs)
        stats.append(row)
        if (PRINT_STEP_STATS and obey_global) or force_print:
            tail = " | ".join([f"{k}={v}" for k, v in kwargs.items()])
            print(f"[{len(stats):04d}] {msg}" + (f" | {tail}" if tail else ""))

    def warn(msg, **kwargs):
        tail = " | ".join([f"{k}={v}" for k, v in kwargs.items()])
        print(f"[!!] {msg}" + (f" | {tail}" if tail else ""))

    cache_dir = make_species_cache_dir(cache_root, species)
    if RESET_CACHE_FOR_EACH_SPECIES and os.path.exists(cache_dir):
        shutil.rmtree(cache_dir)
        cache_dir = make_species_cache_dir(cache_root, species)
    print(f"\n[Species] {species} | cache={cache_dir}")

    genomes = fetch_species_genomes(species, cache_dir, stat, warn)
    stat("Genome universe ready", True, rows=len(genomes), biosamples=int((genomes["biosample_accession"].astype(str) != "").sum()))

    bvbrc_raw, taxon_summary = build_bvbrc_raw_amr(genomes, cache_dir, stat, warn)
    stat("BV-BRC raw AMR ready", True, rows=len(bvbrc_raw), baseline_bvbrc_amr_genomes=int(bvbrc_raw["genome_id"].astype(str).replace("", pd.NA).dropna().nunique()) if not bvbrc_raw.empty else 0)

    candidate_universe, mapped_source_raw, raw_union, stub_df, baseline_isolate_keys = build_multisource_source_tables_v2(
        species=species,
        genomes=genomes,
        bvbrc_raw=bvbrc_raw,
        enable_ebi=enable_ebi,
        enable_ncbi=enable_ncbi,
        ebi_release=ebi_release,
        ncbi_ast_export=ncbi_ast_export,
        ncbi_gcp_project=ncbi_gcp_project,
        cache_dir=cache_dir,
        stat=stat,
        warn=warn,
    )
    stat(
        "Canonical isolate universe ready",
        True,
        canonical_candidates=len(candidate_universe),
        external_stub_candidates=len(stub_df),
        baseline_bvbrc_isolates=len(baseline_isolate_keys),
        raw_union_rows=len(raw_union),
    )

    union_candidate_ids = sorted(set([str(x) for x in raw_union.get("genome_id", pd.Series(dtype=str)).astype(str).tolist() if str(x).strip()]))
    union_candidates = candidate_universe[candidate_universe["genome_id"].astype(str).isin(set(union_candidate_ids))].copy() if union_candidate_ids else candidate_universe.iloc[0:0].copy()

    bs_map_file = os.path.join(cache_dir, "biosample_run_map_SRA_ENA_STRONG_canonical_union.parquet")
    if os.path.exists(bs_map_file):
        bs_map = pd.read_parquet(bs_map_file)
        stat("Loaded STRONG SRA+ENA BioSample map from cache", True, rows=len(bs_map))
    else:
        stat("Building STRONG SRA+ENA BioSample map (canonical union candidates)", True)
        bs_map = run_map_biosample_only_sra_ena(union_candidates, cache_dir=cache_dir, max_biosamples=MAX_BIOSAMPLES_PER_SPECIES, max_sra_records=MAX_SRA_RECORDS)
        bs_map.to_parquet(bs_map_file, index=False)
        stat("Saved STRONG SRA+ENA BioSample map", True, rows=len(bs_map))
    force_text_ids(bs_map, ["biosample_accession", "rnaseq_run_ids", "rnaseq_run_ids_sra", "rnaseq_run_ids_ena", "ena_fastq_ftp", "ena_fastq_md5"])

    presence_flags = build_presence_flags_v2(mapped_source_raw)
    source_tables = {}
    for src in [SOURCE_BVBRC, SOURCE_NCBI, SOURCE_EBI]:
        raw_df = mapped_source_raw.get(src, pd.DataFrame())
        src_ids = sorted(set([str(x) for x in raw_df.get("genome_id", pd.Series(dtype=str)).astype(str).tolist() if str(x).strip()])) if raw_df is not None and not raw_df.empty else []
        src_candidates = candidate_universe[candidate_universe["genome_id"].astype(str).isin(set(src_ids))].copy() if src_ids else candidate_universe.iloc[0:0].copy()
        source_tables[src] = build_candidate_output_table_v2(
            candidate_subset=src_candidates,
            mapped_raw=raw_df,
            bs_map=bs_map,
            baseline_isolate_keys=baseline_isolate_keys,
            presence_flags=presence_flags,
        )
        stat(
            f"Built source table {src}",
            True,
            rows=len(source_tables[src]),
            rnaseq_true=int(source_tables[src]["rnaseq_available"].sum()) if not source_tables[src].empty else 0,
            final_hits=int(source_tables[src]["final_conservative_hit"].sum()) if not source_tables[src].empty else 0,
        )

    final_union = build_candidate_output_table_v2(
        candidate_subset=union_candidates,
        mapped_raw=raw_union,
        bs_map=bs_map,
        baseline_isolate_keys=baseline_isolate_keys,
        presence_flags=presence_flags,
    )
    stat(
        "Built final union table",
        True,
        rows=len(final_union),
        unique_biosamples=int(final_union["biosample_accession"].astype(str).replace("", pd.NA).dropna().nunique()) if not final_union.empty else 0,
        rnaseq_true=int(final_union["rnaseq_available"].sum()) if not final_union.empty else 0,
        final_hits=int(final_union["final_conservative_hit"].sum()) if not final_union.empty else 0,
        true_extra_isolates=int((~final_union["isolate_key"].astype(str).isin(set(baseline_isolate_keys))).sum()) if not final_union.empty else 0,
    )

    raw_debug = raw_union.copy()
    if not raw_debug.empty:
        raw_debug = raw_debug.merge(
            candidate_universe[[
                "genome_id", "isolate_key", "candidate_origin", "candidate_is_external_stub",
                "same_isolate_genome_rows", "same_isolate_genome_ids", "same_isolate_assemblies"
            ]].drop_duplicates(["genome_id", "isolate_key"]),
            on=["genome_id", "isolate_key"],
            how="left",
        )

    source_summary = source_discovery_summary_rows_v2(
        source_raw_map=mapped_source_raw,
        source_tables=source_tables,
        union_table=final_union,
        baseline_isolate_keys=baseline_isolate_keys,
        candidate_universe=candidate_universe,
    )

    discovery_summary = pd.DataFrame([{
        "species": species,
        "bvbrc_genome_universe_rows": int(len(genomes)),
        "canonical_candidate_universe": int(len(candidate_universe)),
        "baseline_bvbrc_isolates_with_ast": int(len(baseline_isolate_keys)),
        "external_only_stub_isolates": int(len(stub_df)),
        "final_union_isolates_with_ast": int(len(final_union)),
        "final_union_unique_biosamples": int(final_union["biosample_accession"].astype(str).replace("", pd.NA).dropna().nunique()) if not final_union.empty else 0,
        "final_union_rnaseq_true": int(final_union["rnaseq_available"].sum()) if not final_union.empty else 0,
        "final_union_final_hits": int(final_union["final_conservative_hit"].sum()) if not final_union.empty else 0,
        "true_extra_isolates_over_bvbrc_baseline": int((~final_union["isolate_key"].astype(str).isin(set(baseline_isolate_keys))).sum()) if not final_union.empty else 0,
        "true_extra_rnaseq_isolates_over_bvbrc_baseline": int(final_union[(~final_union["isolate_key"].astype(str).isin(set(baseline_isolate_keys))) & (final_union["rnaseq_available"])]["isolate_key"].nunique()) if not final_union.empty else 0,
        "true_extra_final_hits_over_bvbrc_baseline": int(final_union[(~final_union["isolate_key"].astype(str).isin(set(baseline_isolate_keys))) & (final_union["final_conservative_hit"])]["isolate_key"].nunique()) if not final_union.empty else 0,
    }])

    final_hits_only = final_union[final_union["final_conservative_hit"]].copy()

    os.makedirs(out_root, exist_ok=True)
    species_key = species.replace(" ", "_")

    paths = []
    final_union_csv = safe_to_csv(final_union, os.path.join(out_root, f"{species_key}_FINAL_union_complete_dedup.csv"), index=False)
    paths.append(final_union_csv)
    final_hits_csv = safe_to_csv(final_hits_only, os.path.join(out_root, f"{species_key}_FINAL_hits_only.csv"), index=False)
    paths.append(final_hits_csv)

    for src in [SOURCE_BVBRC, SOURCE_NCBI, SOURCE_EBI]:
        slug = source_slug_v2(src)
        p = safe_to_csv(source_tables.get(src, pd.DataFrame(columns=OUTPUT_COLS_V2)), os.path.join(out_root, f"{species_key}_{slug}_validated.csv"), index=False)
        paths.append(p)

    raw_debug_csv = safe_to_csv(raw_debug, os.path.join(out_root, f"{species_key}_raw_amr_records_canonicalised.csv"), index=False)
    source_summary_csv = safe_to_csv(source_summary, os.path.join(out_root, f"{species_key}_source_summary.csv"), index=False)
    discovery_csv = safe_to_csv(discovery_summary, os.path.join(out_root, f"{species_key}_discovery_summary.csv"), index=False)
    dbg_csv = safe_to_csv(pd.DataFrame(stats), os.path.join(out_root, f"{species_key}_debug_stats.csv"), index=False)
    tax_csv = safe_to_csv(taxon_summary, os.path.join(out_root, f"{species_key}_amr_taxon_summary.csv"), index=False)
    paths.extend([raw_debug_csv, source_summary_csv, discovery_csv, dbg_csv, tax_csv])

    stat(
        "Saved outputs",
        True,
        final_union_csv=final_union_csv,
        final_hits_csv=final_hits_csv,
        source_summary_csv=source_summary_csv,
        discovery_csv=discovery_csv,
        raw_debug_csv=raw_debug_csv,
    )
    return {
        "primary_csv": final_union_csv,
        "all_outputs": paths,
    }


def parse_species_list(species_arg: str, species_file: Optional[str]) -> List[str]:
    if species_file:
        with open(species_file, "r", encoding="utf-8") as f:
            lines = [ln.strip() for ln in f.readlines()]
        return [ln for ln in lines if ln and not ln.startswith("#")]
    if species_arg:
        return [s.strip() for s in species_arg.split(",") if s.strip()]
    return [s.strip() for s in SPECIES_INPUTS_DEFAULT.split(",") if s.strip()]


def zip_outputs(out_root: str, paths: List[str], zip_path: str) -> None:
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for f in paths:
            if os.path.exists(f):
                z.write(f, arcname=os.path.relpath(f, out_root))


def _read_csv_as_text(path: str) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False, na_filter=False)


def build_combined_rnaseq_true_csv(species_csvs: List[Tuple[str, str]], out_root: str) -> str:
    keep_frames = []
    template_cols = None
    for sp, path in species_csvs:
        df = _read_csv_as_text(path)
        if template_cols is None:
            template_cols = ["species"] + df.columns.tolist()
        rn = df["rnaseq_available"].astype(str).str.strip().str.lower().isin({"true", "1", "yes"})
        df_true = df[rn].copy()
        if not df_true.empty:
            df_true.insert(0, "species", sp)
            keep_frames.append(df_true)
        print(f"[Combined] {sp} | total={len(df)} | rnaseq_true={int(rn.sum())}")
    combined = pd.concat(keep_frames, ignore_index=True) if keep_frames else pd.DataFrame(columns=template_cols if template_cols else ["species"])
    combined_true_csv = os.path.join(out_root, f"rnaseq_true_combined_{int(time.time())}.csv")
    combined.to_csv(combined_true_csv, index=False)
    print("Saved combined rnaseq_true CSV:", combined_true_csv, "| rows:", len(combined))
    return combined_true_csv


def _safe_delete(paths: List[str]) -> None:
    for p in paths:
        try:
            if p and os.path.exists(p):
                os.remove(p)
        except Exception:
            pass



def main() -> None:
    ap = argparse.ArgumentParser(description="Conservative RNA-seq + AMR checker (multisource complete isolate mode)")
    ap.add_argument("--species", type=str, default="", help="Comma-separated species list")
    ap.add_argument("--species-file", type=str, default="", help="Text file with one species per line")
    ap.add_argument("--cache-root", type=str, default="./cache", help="Cache directory")
    ap.add_argument("--out-root", type=str, default="./outputs", help="Output directory")
    ap.add_argument("--disable-ebi", action="store_true", help="Disable EMBL-EBI AMR portal integration")
    ap.add_argument("--disable-ncbi", action="store_true", help="Disable NCBI AST integration")
    ap.add_argument("--ebi-release", type=str, default=os.environ.get("EBI_AMR_RELEASE", "auto"), help="EMBL-EBI AMR release like 2025-12 or 'auto'")
    ap.add_argument("--ncbi-ast-export", type=str, default=os.environ.get("NCBI_AST_EXPORT", ""), help="Path to AST Browser CSV/TSV export")
    ap.add_argument("--ncbi-gcp-project", type=str, default=os.environ.get("NCBI_GCP_PROJECT", ""), help="GCP project for querying public NCBI BigQuery tables")
    args = ap.parse_args()

    species_file = args.species_file.strip() or None
    species_list = parse_species_list(args.species, species_file)
    print("Species list:", species_list)
    os.makedirs(args.cache_root, exist_ok=True)
    os.makedirs(args.out_root, exist_ok=True)

    all_outputs: List[str] = []
    species_primary_csvs: List[Tuple[str, str]] = []
    for s in species_list:
        print("\n" + "=" * 90)
        print("RUNNING SPECIES:", s)
        print("=" * 90)
        try:
            result = run_one_species(
                s, cache_root=args.cache_root, out_root=args.out_root,
                enable_ebi=not args.disable_ebi, enable_ncbi=not args.disable_ncbi,
                ebi_release=args.ebi_release, ncbi_ast_export=args.ncbi_ast_export.strip(), ncbi_gcp_project=args.ncbi_gcp_project.strip(),
            )
        except Exception as e:
            print(f"[!!] FAILED species='{s}' | error={type(e).__name__}: {str(e)[:400]}")
            continue
        all_outputs.extend(result["all_outputs"])
        species_primary_csvs.append((norm_species_name(s), result["primary_csv"]))

    if len(species_primary_csvs) > 1:
        print("\n" + "=" * 90)
        print("BUILDING COMBINED rnaseq_available==True CSV")
        print("=" * 90)
        combined_true_csv = build_combined_rnaseq_true_csv(species_primary_csvs, args.out_root)
        all_outputs.append(combined_true_csv)
        ts = time.strftime("%Y%m%d_%H%M%S")
        zip_path = os.path.join(args.out_root, f"rnaseq_amr_outputs_{ts}.zip")
        zip_outputs(args.out_root, all_outputs, zip_path)
        print("Zipped outputs:", zip_path)
        _safe_delete(all_outputs)
        print("Cleaned up individual CSVs (multi-species run) → outputs folder now has the ZIP only.")
    else:
        print("\nSingle species run → CSV outputs only (no zip, no combined file).")

    print("\nDONE.")
    print("Outputs folder:", os.path.abspath(args.out_root))
    print("Cache folder:", os.path.abspath(args.cache_root))


if __name__ == "__main__":
    main()
