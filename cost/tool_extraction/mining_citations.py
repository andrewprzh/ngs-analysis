"""Enrich the bio.tools catalog with citation counts.

For each tool, extract DOI/PMID from its 'publication' field, then look up
citation counts. OpenAlex (free, no key) is the default source; Scopus is used
as an optional enrichment when SCOPUS_API_KEY is set in the environment.

Lookups are cached to disk so the script is resumable across runs.

Usage:
    python mining_citations.py
    SCOPUS_API_KEY=... python mining_citations.py   # also query Scopus
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

import pandas as pd
import requests

INPUT_FILE = Path("bio_tools/bio_tools_multiomics.csv")
OUTPUT_FILE = Path("bio_tools/bio_tools_multiomics_citations.csv")
CACHE_FILE = Path("bio_tools/.citation_cache.json")

OPENALEX_URL = "https://api.openalex.org/works"
SCOPUS_URL = "https://api.elsevier.com/content/search/scopus"

# Polite OpenAlex usage: an email puts you in the faster "polite pool".
CONTACT_EMAIL = os.getenv("OPENALEX_EMAIL", "")
SCOPUS_API_KEY = os.getenv("SCOPUS_API_KEY", "")

REQUEST_TIMEOUT = 60
MAX_RETRIES = 4
RETRY_SLEEP = 2.0
QUERY_SLEEP = 0.2

# Manual DOI/PMID fixes for tools whose publication metadata is missing/wrong.
MANUAL_FIXES = {
    "mixomics": {"doi": "10.1371/journal.pcbi.1005752", "pmid": "29099853"},
}


# --------------------------------------------------------------------------- #
# Cache
# --------------------------------------------------------------------------- #
def load_cache() -> dict:
    if CACHE_FILE.exists():
        return json.loads(CACHE_FILE.read_text(encoding="utf-8"))
    return {}


def save_cache(cache: dict) -> None:
    CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)
    CACHE_FILE.write_text(json.dumps(cache), encoding="utf-8")


# --------------------------------------------------------------------------- #
# ID extraction
# --------------------------------------------------------------------------- #
def parse_json_field(value) -> list | dict | None:
    """Parse a CSV-stored JSON field, tolerating doubled quotes and NaN."""
    if value is None or (isinstance(value, float)) or not str(value).strip():
        return None
    text = str(value).replace('""', '"')
    try:
        return json.loads(text)
    except (json.JSONDecodeError, ValueError):
        return None


def extract_ids(publication_field) -> tuple[str | None, str | None]:
    """Return (doi, pmid) from a tool's publication field; first of each."""
    pubs = parse_json_field(publication_field)
    if not pubs:
        return None, None
    if isinstance(pubs, dict):
        pubs = [pubs]

    doi = pmid = None
    for pub in pubs:
        if not isinstance(pub, dict):
            continue
        if doi is None and pub.get("doi"):
            doi = str(pub["doi"]).strip()
        if pmid is None and pub.get("pmid"):
            pmid = str(pub["pmid"]).strip()
    return doi, pmid


# --------------------------------------------------------------------------- #
# Citation lookups
# --------------------------------------------------------------------------- #
def get_json(url: str, params: dict, headers: dict | None = None) -> dict | None:
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            resp = requests.get(
                url, params=params, headers=headers, timeout=REQUEST_TIMEOUT
            )
            if resp.status_code == 404:
                return None
            resp.raise_for_status()
            return resp.json()
        except requests.RequestException as exc:
            print(f"    request failed (attempt {attempt}/{MAX_RETRIES}): {exc}")
            time.sleep(RETRY_SLEEP * attempt)
    return None


def get_openalex_citations(doi: str | None, pmid: str | None) -> int | None:
    """Look up cited_by_count from OpenAlex by DOI (preferred) or PMID."""
    if doi:
        ident = f"https://doi.org/{doi}"
    elif pmid:
        ident = f"pmid:{pmid}"
    else:
        return None

    params = {}
    if CONTACT_EMAIL:
        params["mailto"] = CONTACT_EMAIL
    data = get_json(f"{OPENALEX_URL}/{ident}", params)
    if data is None:
        return None
    return data.get("cited_by_count")


def get_scopus_citations(doi: str | None, pmid: str | None) -> int | None:
    """Look up citedby-count from Scopus; requires SCOPUS_API_KEY."""
    if not SCOPUS_API_KEY:
        return None
    if doi:
        query = f"DOI({doi})"
    elif pmid:
        query = f"PMID({pmid})"
    else:
        return None

    data = get_json(
        SCOPUS_URL,
        params={"query": query, "count": 1},
        headers={"X-ELS-APIKey": SCOPUS_API_KEY, "Accept": "application/json"},
    )
    if data is None:
        return None
    entries = data.get("search-results", {}).get("entry", [])
    if not entries:
        return None
    count = entries[0].get("citedby-count")
    return int(count) if count is not None else None


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main() -> None:
    if not INPUT_FILE.exists():
        raise SystemExit(f"Input file not found: {INPUT_FILE}")
    if "publication" not in pd.read_csv(INPUT_FILE, nrows=0).columns:
        raise SystemExit("No 'publication' column found in input file.")

    df = pd.read_csv(INPUT_FILE)
    use_scopus = bool(SCOPUS_API_KEY)
    print(f"Loaded {len(df)} tools. Scopus enrichment: {'on' if use_scopus else 'off'}.\n")

    # Extract DOI/PMID
    ids = df["publication"].apply(extract_ids)
    df["doi"] = [i[0] for i in ids]
    df["pmid"] = [i[1] for i in ids]

    # Apply manual fixes by (case-insensitive substring) tool name
    if "name" in df.columns:
        name_lower = df["name"].fillna("").str.lower()
        for key, fix in MANUAL_FIXES.items():
            mask = name_lower.str.contains(key, regex=False)
            df.loc[mask, "doi"] = fix["doi"]
            df.loc[mask, "pmid"] = fix["pmid"]

    cache = load_cache()
    openalex_counts: list[int | None] = []
    scopus_counts: list[int | None] = []

    for i, row in df.iterrows():
        doi, pmid = row["doi"], row["pmid"]
        tool_id = row.get("biotoolsID", i)
        key = f"{doi}|{pmid}"
        print(f"Checking {i + 1} of {len(df)}: {tool_id}")

        entry = cache.get(key, {})

        if "openalex" not in entry:
            entry["openalex"] = get_openalex_citations(doi, pmid)
            time.sleep(QUERY_SLEEP)
        if use_scopus and "scopus" not in entry:
            entry["scopus"] = get_scopus_citations(doi, pmid)
            time.sleep(QUERY_SLEEP)

        cache[key] = entry
        save_cache(cache)  # incremental -> resumable

        openalex_counts.append(entry.get("openalex"))
        scopus_counts.append(entry.get("scopus"))

    df["openalex_citations"] = openalex_counts
    if use_scopus:
        df["scopus_citations"] = scopus_counts

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_FILE, index=False)
    print(f"\nDone. Saved to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
