"""Harvest multi-omics tools from the bio.tools registry.

Reads search terms from search_terms.txt, queries the bio.tools API, paginates
through every result page, deduplicates by biotoolsID, flattens nested fields to
JSON strings, and writes a flat CSV.

Usage:
    python mining_biotools.py
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd
import requests

API_URL = "https://bio.tools/api/t/"
SEARCH_TERMS_FILE = Path("search_terms.txt")
OUTPUT_FILE = Path("bio_tools/bio_tools_multiomics.csv")

REQUEST_TIMEOUT = 60
MAX_RETRIES = 4
RETRY_SLEEP = 2.0
PAGE_SLEEP = 0.3  # polite pause between pages


def load_search_terms(path: Path) -> list[str]:
    """Read terms from file, dropping blank lines and '#' comments, deduped."""
    terms: list[str] = []
    seen: set[str] = set()
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].strip()
        if line and line.lower() not in seen:
            seen.add(line.lower())
            terms.append(line)
    return terms


def get_json(url: str, params: dict) -> dict | None:
    """GET a URL with retries; return parsed JSON or None on failure."""
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            resp = requests.get(url, params=params, timeout=REQUEST_TIMEOUT)
            resp.raise_for_status()
            return resp.json()
        except requests.RequestException as exc:
            print(f"  request failed (attempt {attempt}/{MAX_RETRIES}): {exc}")
            time.sleep(RETRY_SLEEP * attempt)
    return None


def fetch_query(term: str) -> list[dict]:
    """Fetch all pages of results for a single search term."""
    rows: list[dict] = []
    page = 1
    while True:
        print(f"Query: {term} | Page: {page}")
        data = get_json(API_URL, {"q": term, "page": page, "format": "json"})
        if data is None:
            break

        hits = data.get("list") or []
        if not hits:
            break

        for hit in hits:
            hit["_matched_query"] = term  # keep track of which term found it
        rows.extend(hits)

        # bio.tools returns "next" as e.g. "?page=2" while more pages remain.
        if not data.get("next"):
            break
        page += 1
        time.sleep(PAGE_SLEEP)
    return rows


def flatten_value(value):
    """Render nested (list/dict) values as JSON strings; pass scalars through."""
    if isinstance(value, (list, dict)):
        return json.dumps(value, ensure_ascii=False) if value else None
    return value


def main() -> None:
    if not SEARCH_TERMS_FILE.exists():
        raise SystemExit(f"Search terms file not found: {SEARCH_TERMS_FILE}")

    terms = load_search_terms(SEARCH_TERMS_FILE)
    print(f"Loaded {len(terms)} unique search terms.\n")

    all_rows: list[dict] = []
    for term in terms:
        all_rows.extend(fetch_query(term))

    if not all_rows:
        raise SystemExit("No results returned from bio.tools API.")

    df = pd.DataFrame(all_rows)

    if "biotoolsID" in df.columns:
        df = df.drop_duplicates(subset="biotoolsID", keep="first")
    else:
        print("Warning: 'biotoolsID' column not found; using full-row unique.")
        df = df.drop_duplicates()

    df = df.applymap(flatten_value)

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_FILE, index=False)
    print(f"\nDone. Saved {len(df)} unique tools to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
