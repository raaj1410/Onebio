from __future__ import annotations
from pathlib import Path
from typing import Tuple
import os
import re
import time
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError
import xml.etree.ElementTree as ET
import json
import pandas as pd
from .utils import ensure_dir, log, PipelineError, safe_split_semi

REQUIRED_COLS = ["biosample_accession", "assembly_accession_best"]


def _http_get(url: str, params: dict, timeout: int = 30, retries: int = 3) -> bytes:
    """Stdlib-only HTTP GET helper with retries."""
    full = f"{url}?{urlencode({k: v for k, v in params.items() if v is not None})}"
    ua = os.environ.get("USER_AGENT", "biosample-meta-fetch/1.0")
    last_err = None
    for i in range(retries):
        try:
            req = Request(full, headers={"User-Agent": ua})
            with urlopen(req, timeout=timeout) as r:
                return r.read()
        except (HTTPError, URLError, TimeoutError) as e:
            last_err = e
            time.sleep(1.0 + i * 1.5)
    raise last_err  # type: ignore


def _norm_key(k: str) -> str:
    k = (k or "").strip().lower()
    k = re.sub(r"[\s\-\/]+", "_", k)
    k = re.sub(r"[^a-z0-9_]+", "", k)
    return k


def _pick_first(d: dict, keys: list[str]) -> str:
    for k in keys:
        v = d.get(k)
        if v:
            return str(v)
    return ""


def fetch_biosample_metadata_ncbi(
    biosample_acc: str,
    *,
    email: str = "",
    tool: str = "amr_atlas",
    api_key: str = "",
    sleep_s: float = 0.34,
) -> dict:
    """Fetch BioSample attributes from NCBI (E-utilities) using a SAMN accession.

    Returns a dict with selected normalised fields plus raw attributes JSON.
    """
    biosample_acc = str(biosample_acc).strip()
    if not biosample_acc or biosample_acc.lower() == "nan":
        return {}

    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    # 1) resolve SAMN accession -> numeric BioSample ID
    esearch_params = {
        "db": "biosample",
        "term": biosample_acc,
        "retmode": "json",
        "tool": tool,
        "email": email or None,
        "api_key": api_key or None,
    }
    esearch = json.loads(_http_get(f"{base}/esearch.fcgi", esearch_params).decode("utf-8"))
    idlist = esearch.get("esearchresult", {}).get("idlist", []) or []
    if not idlist:
        return {"biosample_accession": biosample_acc, "fetch_status": "not_found"}
    bs_id = idlist[0]

    # 2) fetch BioSample XML
    efetch_params = {
        "db": "biosample",
        "id": bs_id,
        "retmode": "xml",
        "tool": tool,
        "email": email or None,
        "api_key": api_key or None,
    }
    xml_bytes = _http_get(f"{base}/efetch.fcgi", efetch_params)
    time.sleep(max(0.0, sleep_s))

    root = ET.fromstring(xml_bytes)
    bs = root.find(".//BioSample")
    if bs is None:
        return {"biosample_accession": biosample_acc, "fetch_status": "xml_parse_failed"}

    # Organism + taxonomy
    org = bs.find(".//Description/Organism")
    organism_name = org.text.strip() if (org is not None and org.text) else ""
    taxid = org.attrib.get("taxonomy_id", "") if org is not None else ""

    title_el = bs.find(".//Description/Title")
    title = title_el.text.strip() if (title_el is not None and title_el.text) else ""

    # Collect all attributes
    raw_attrs: dict[str, str] = {}
    for a in bs.findall(".//Attributes/Attribute"):
        key = a.attrib.get("attribute_name") or a.attrib.get("harmonized_name") or ""
        key_n = _norm_key(key)
        val = (a.text or "").strip()
        if not key_n or not val:
            continue
        if key_n in raw_attrs and raw_attrs[key_n] != val:
            raw_attrs[key_n] = f"{raw_attrs[key_n]};{val}"
        else:
            raw_attrs[key_n] = val

    host = _pick_first(raw_attrs, ["host", "host_scientific_name", "host_common_name"])
    isolation_source = _pick_first(raw_attrs, ["isolation_source"])
    body_site = _pick_first(raw_attrs, ["body_site", "anatomical_site", "host_tissue_sampled"])
    geo = _pick_first(raw_attrs, ["geo_loc_name", "geographic_location", "country"])
    collection_date = _pick_first(raw_attrs, ["collection_date", "date_of_isolation", "isolation_date"])
    strain = _pick_first(raw_attrs, ["strain", "isolate", "isolation_id"])

    collection_year = ""
    if collection_date:
        m = re.search(r"(19\d{2}|20\d{2})", collection_date)
        if m:
            collection_year = m.group(1)

    host_l = host.lower()
    host_is_human = int(("homo sapiens" in host_l) or (host_l == "human") or (raw_attrs.get("host_taxid") == "9606"))
    iso_l = f"{isolation_source} {body_site}".lower()
    isolation_is_blood = int(any(x in iso_l for x in ["blood", "bloodstream", "bacteremia", "septicaemia", "sepsis"]))

    return {
        "biosample_accession": biosample_acc,
        "organism_name": organism_name,
        "taxon_id": taxid,
        "title": title,
        "host": host,
        "isolation_source": isolation_source,
        "body_site": body_site,
        "geo_loc_name": geo,
        "collection_date": collection_date,
        "collection_year": collection_year,
        "strain": strain,
        "isolate": strain,
        "strain_or_isolate": strain,
        "host_is_human": host_is_human,
        "isolation_is_blood": isolation_is_blood,
        "raw_attributes_json": json.dumps(raw_attrs, ensure_ascii=False),
        "fetch_status": "ok",
    }

def build_manifest(xlsx: Path, out_dir: Path) -> Tuple[Path, Path, Path]:
    ensure_dir(out_dir)
    df = pd.read_excel(xlsx)

    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise PipelineError(f"XLSX missing required columns: {missing}")

    df["biosample_accession"] = df["biosample_accession"].astype(str).str.strip()
    df["assembly_accession_best"] = df["assembly_accession_best"].astype(str).str.strip()

    if "rnaseq_run_ids" not in df.columns:
        df["rnaseq_run_ids"] = ""

    df["run_list"] = df["rnaseq_run_ids"].apply(lambda s: ";".join(safe_split_semi(str(s))))

    for c in ["ena_fastq_ftp", "ena_fastq_md5"]:
        if c not in df.columns:
            df[c] = ""

    keep_cols = [
        "biosample_accession",
        "bioproject_accession",
        "taxon_id",
        "strain",
        "genome_id",
        "assembly_accession_best",
        "evidence_tier",
        "rnaseq_run_count",
        "run_list",
        "ena_fastq_ftp",
        "ena_fastq_md5",
        "access_method",
    ]
    for c in keep_cols:
        if c not in df.columns:
            df[c] = ""

    manifest = df[keep_cols].copy()

    def split_ab(s):
        return [x.strip() for x in str(s or "").split(";") if x.strip() and str(x).lower() != "nan"]

    ast_rows = []
    for _, r in df.iterrows():
        bs = str(r["biosample_accession"]).strip()
        for ab in split_ab(r.get("antibiotics_resistant", "")):
            ast_rows.append({"biosample_accession": bs, "antibiotic": ab, "phenotype": "R"})
        for ab in split_ab(r.get("antibiotics_susceptible", "")):
            ast_rows.append({"biosample_accession": bs, "antibiotic": ab, "phenotype": "S"})
        for ab in split_ab(r.get("antibiotics_intermediate", "")):
            ast_rows.append({"biosample_accession": bs, "antibiotic": ab, "phenotype": "I"})
    ast_long = pd.DataFrame(ast_rows).drop_duplicates()

    # --- BioSample metadata table (1 row per BioSample) ---
    # XLSX may not contain provenance fields like host/isolation source/collection date.
    # We enrich by fetching BioSample attributes from NCBI using the SAMN accession.
    biosample_meta = manifest.copy()

    # add derived flags from sheet
    biosample_meta["rnaseq_run_n"] = biosample_meta["run_list"].apply(
        lambda s: len([x for x in safe_split_semi(str(s)) if x.strip() and str(x).lower() != "nan"])
    )
    biosample_meta["has_rnaseq"] = (biosample_meta["rnaseq_run_n"] > 0).astype(int)
    biosample_meta = biosample_meta.drop_duplicates(subset=["biosample_accession"], keep="first")

    # Fetch + cache
    cache_dir = out_dir / "_cache"
    ensure_dir(cache_dir)
    cache_path = cache_dir / "biosample_meta_cache.json"
    cache: dict = {}
    if cache_path.exists():
        try:
            cache = json.loads(cache_path.read_text(encoding="utf-8"))
        except Exception:
            cache = {}

    email = os.environ.get("NCBI_EMAIL", "")
    api_key = os.environ.get("NCBI_API_KEY", "")
    tool = os.environ.get("NCBI_TOOL", "amr_atlas")
    sleep_s = float(os.environ.get("NCBI_SLEEP", "0.34"))

    fetched_rows = []
    for acc in biosample_meta["biosample_accession"].tolist():
        acc = str(acc).strip()
        if not acc or acc.lower() == "nan":
            continue
        if acc in cache and isinstance(cache.get(acc), dict) and cache[acc].get("fetch_status") == "ok":
            meta = cache[acc]
        else:
            try:
                meta = fetch_biosample_metadata_ncbi(acc, email=email, api_key=api_key, tool=tool, sleep_s=sleep_s)
            except Exception as e:
                meta = {"biosample_accession": acc, "fetch_status": f"error:{type(e).__name__}"}
            cache[acc] = meta
        fetched_rows.append(meta)

    try:
        cache_path.write_text(json.dumps(cache, indent=2, ensure_ascii=False), encoding="utf-8")
    except Exception:
        pass

    if fetched_rows:
        fetched_df = pd.DataFrame(fetched_rows).drop_duplicates(subset=["biosample_accession"], keep="first")
        biosample_meta = biosample_meta.merge(fetched_df, on="biosample_accession", how="left", suffixes=("", "_ncbi"))

        # Prefer NCBI values when sheet has blanks
        for col in ["taxon_id", "strain", "isolate"]:
            nc = f"{col}_ncbi"
            if col in biosample_meta.columns and nc in biosample_meta.columns:
                base = biosample_meta[col].astype(str)
                base = base.where(~base.str.lower().isin(["", "nan", "none"]), "")
                biosample_meta[col] = base.where(base != "", biosample_meta[nc])
                biosample_meta = biosample_meta.drop(columns=[nc])

        # keep the rest as new columns (organism_name/host/isolation_source/etc.)
        drop_ncbi = [c for c in biosample_meta.columns if c.endswith("_ncbi")]
        if drop_ncbi:
            biosample_meta = biosample_meta.drop(columns=drop_ncbi)

    prov = {
        "source_xlsx": str(xlsx),
        "rows": int(len(df)),
        "columns_in_xlsx": list(df.columns),
        "biosample_meta_tsv": str(out_dir / "biosample_meta.tsv"),
        "note": "BioSample anchored manifest and AST long built from xlsx.",
    }

    manifest_path = out_dir / "manifest.tsv"
    ast_path = out_dir / "ast_long.tsv"
    biosample_meta_path = out_dir / "biosample_meta.tsv"
    prov_path = out_dir / "provenance.json"

    manifest.to_csv(manifest_path, sep="\t", index=False)
    ast_long.to_csv(ast_path, sep="\t", index=False)
    biosample_meta.to_csv(biosample_meta_path, sep="\t", index=False)
    prov_path.write_text(json.dumps(prov, indent=2), encoding="utf-8")

    log(f"Manifest: {manifest_path}")
    log(f"AST long: {ast_path}")
    log(f"BioSample meta: {biosample_meta_path}")
    log(f"Provenance: {prov_path}")
    return manifest_path, ast_path, prov_path
