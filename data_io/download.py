"""
data_io/download.py
====================
Fetch and cache raw meteorological data for Niort-Souche (WMO 07330).

Priority:
  1. Meteo-France 'Donnees climatologiques horaires' on data.gouv.fr
     - files are grouped by multi-year PERIOD (not one-per-month).
       We resolve the actual file list from the dataset API before downloading.
  2. Infoclimat HTML archive, one page per calendar day (fallback).

Cache contract: every downloaded payload is written to data/raw/ once;
subsequent calls return CACHED without network access.
"""

import gzip
import json
import logging
import time
from datetime import datetime, timedelta
from pathlib import Path

import requests

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import RAW, STATION_WMO

PROJECT_ROOT = Path(RAW).resolve().parent.parent


def _relpath(path) -> str:
    """Return a path relative to the project root, never absolute.

    Keeps the cached-download manifest free of the machine's directory layout.
    """
    try:
        return str(Path(path).resolve().relative_to(PROJECT_ROOT))
    except ValueError:
        return Path(path).name


log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (compatible; SainteSolineEvap/1.0; research use)"
    )
}
TIMEOUT = 30
MAX_RETRIES = 3


def _get(url: str, stream: bool = False, **kwargs):
    """GET with retries and exponential backoff."""
    for attempt in range(MAX_RETRIES):
        try:
            r = requests.get(url, headers=HEADERS, timeout=TIMEOUT,
                             stream=stream, **kwargs)
            if r.status_code == 200:
                return r
            if r.status_code == 404:
                log.warning("404 %s", url)
                return None
            log.warning("HTTP %s for %s (attempt %d)", r.status_code, url, attempt + 1)
        except requests.RequestException as e:
            log.warning("Request error %s: %s (attempt %d)", url, e, attempt + 1)
        time.sleep(2 ** attempt)
    return None


def month_range(start: str, end: str) -> list:
    """Return list of (year, month) tuples covering start...end (inclusive)."""
    t0 = datetime.strptime(start[:10], "%Y-%m-%d")
    t1 = datetime.strptime(end[:10], "%Y-%m-%d")
    result = []
    cur = t0.replace(day=1)
    while cur <= t1:
        result.append((cur.year, cur.month))
        if cur.month == 12:
            cur = cur.replace(year=cur.year + 1, month=1)
        else:
            cur = cur.replace(month=cur.month + 1)
    return result


# ---------------------------------------------------------------------------
# Source 1 -- data.gouv.fr 'Donnees climatologiques horaires'
# ---------------------------------------------------------------------------

DATAGOUV_DATASET_ID = "6569b4473bedf2e7abad3b72"
DATAGOUV_API = (
    f"https://www.data.gouv.fr/api/1/datasets/{DATAGOUV_DATASET_ID}/"
)


def resolve_datagouv_resources(dept: str = "79", cache_dir: Path = RAW):
    """
    Query the data.gouv API once to get the list of resource URLs for the
    requested department. Cache the manifest so repeat calls are offline.
    Returns list of {title, url} dicts, or [] on failure.
    """
    cache_file = cache_dir / f"datagouv_resources_H_{dept}.json"
    if cache_file.exists() and cache_file.stat().st_size > 0:
        log.info("CACHED datagouv resource list %s", cache_file)
        with open(cache_file) as f:
            return json.load(f)

    log.info("Resolving datagouv resources for dept %s ...", dept)
    r = _get(DATAGOUV_API)
    if r is None:
        log.error("Cannot reach data.gouv API -- datagouv source unavailable")
        return []

    try:
        data = r.json()
    except Exception as e:
        log.error("JSON parse error from data.gouv API: %s", e)
        return []

    resources = []
    for res in data.get("resources", []):
        title = res.get("title", "")
        url   = res.get("url", "")
        if f"departement_{dept}" in title and url.endswith(".csv.gz"):
            resources.append({"title": title, "url": url})

    log.info("Found %d resource(s) matching H_%s", len(resources), dept)
    cache_dir.mkdir(parents=True, exist_ok=True)
    with open(cache_file, "w") as f:
        json.dump(resources, f, indent=2)
    return resources


def fetch_datagouv_file(resource: dict, cache_dir: Path = RAW):
    """Download one csv.gz resource and cache it. Return local path or None."""
    title = resource["title"]
    url   = resource["url"]
    safe  = "".join(c if c.isalnum() or c in "-_." else "_" for c in title)
    dest  = cache_dir / f"{safe}.csv.gz"
    if dest.exists() and dest.stat().st_size > 0:
        log.info("CACHED %s", dest)
        return dest

    log.info("Downloading datagouv %s -> %s", url, dest)
    cache_dir.mkdir(parents=True, exist_ok=True)
    r = _get(url, stream=True)
    if r is None:
        return None
    with open(dest, "wb") as f:
        for chunk in r.iter_content(65536):
            f.write(chunk)
    log.info("Saved %d bytes -> %s", dest.stat().st_size, dest)
    return dest


def fetch_datagouv_all(dept: str = "79", cache_dir: Path = RAW,
                       start_year: int = 2025) -> list:
    """
    Fetch period files for the department that overlap with start_year or later.
    Only downloads relevant files (avoids gigabytes of historical data).
    Returns list of local paths.
    """
    resources = resolve_datagouv_resources(dept, cache_dir)
    paths = []
    for res in resources:
        title = res.get("title", "")
        # Extract period from title (e.g. '2025-2026')
        # Download only if the period end year >= start_year
        import re as _re
        m = _re.search(r"(\d{4})-(\d{4})", title)
        if m:
            end_yr = int(m.group(2))
            if end_yr < start_year:
                log.info("Skipping historical file %s (before %d)", title, start_year)
                continue
        p = fetch_datagouv_file(res, cache_dir)
        if p:
            paths.append(p)
    return paths


# ---------------------------------------------------------------------------
# Source 2 -- Infoclimat HTML (fallback, per calendar day)
# ---------------------------------------------------------------------------

INFOCLIMAT_BASE = (
    "https://www.infoclimat.fr/observations-meteo/archives"
    "/{day:02d}/{month_fr}/{year}/niort-souche/07330.html"
)
MONTHS_FR = {
    1: "janvier", 2: "fevrier", 3: "mars", 4: "avril",
    5: "mai", 6: "juin", 7: "juillet", 8: "aout",
    9: "septembre", 10: "octobre", 11: "novembre", 12: "decembre",
}


def fetch_infoclimat_day(date_str: str, cache_dir: Path = RAW):
    """
    Fetch one day's Infoclimat HTML archive page and cache it.
    date_str: 'YYYY-MM-DD'
    Returns local .html path or None.
    """
    d = datetime.strptime(date_str, "%Y-%m-%d")
    month_fr = MONTHS_FR[d.month]
    url  = INFOCLIMAT_BASE.format(day=d.day, month_fr=month_fr, year=d.year)
    dest = cache_dir / f"infoclimat_{date_str}.html"

    if dest.exists() and dest.stat().st_size > 200:
        log.info("CACHED %s", dest)
        return dest

    log.info("Downloading Infoclimat %s ...", url)
    cache_dir.mkdir(parents=True, exist_ok=True)
    r = _get(url)
    if r is None:
        return None

    dest.write_bytes(r.content)

    text = dest.read_text(errors="replace")
    if "aucune" in text.lower() and len(text) < 5000:
        log.warning("Infoclimat: no data for %s", date_str)
        dest.unlink()
        return None

    log.info("Saved %d bytes -> %s", dest.stat().st_size, dest)
    return dest


# ---------------------------------------------------------------------------
# Orchestrator: ensure_raw
# ---------------------------------------------------------------------------

def _all_dates(start: str, end: str) -> list:
    t0 = datetime.strptime(start[:10], "%Y-%m-%d")
    t1 = datetime.strptime(end[:10], "%Y-%m-%d")
    dates = []
    cur = t0
    while cur <= t1:
        dates.append(cur.strftime("%Y-%m-%d"))
        cur += timedelta(days=1)
    return dates


def ensure_raw(start: str, end: str) -> dict:
    """
    Guarantee raw files cover [start, end].

    Strategy:
      1. Try to download all data.gouv.fr period files for dept 79.
         Scan them to determine which dates for our station are present.
      2. For any date still missing, attempt infoclimat.
      3. Write manifest and return it.
    """
    RAW.mkdir(parents=True, exist_ok=True)
    manifest_path = RAW / "manifest.json"

    dates_needed = set(_all_dates(start, end))

    # --- Try datagouv --------------------------------------------------
    start_year = int(start[:4])
    datagouv_paths = fetch_datagouv_all(dept="79", cache_dir=RAW, start_year=start_year)
    datagouv_covered = set()

    if datagouv_paths:
        import csv
        from config import STATION_NATIONAL as _STA_NAT
        for p in datagouv_paths:
            try:
                with gzip.open(p, "rt", encoding="latin-1") as f:
                    reader = csv.DictReader(f, delimiter=";")
                    for row in reader:
                        sta = row.get("NUM_POSTE", row.get("POSTE", ""))
                        # Match on either WMO id substring or national id
                        if STATION_WMO not in sta and _STA_NAT not in sta:
                            continue
                        date_raw = row.get("AAAAMMJJHH", row.get("DATE", ""))
                        if len(date_raw) >= 8:
                            y, m, d = date_raw[:4], date_raw[4:6], date_raw[6:8]
                            day_str = f"{y}-{m}-{d}"
                            if day_str in dates_needed:
                                datagouv_covered.add(day_str)
            except Exception as e:
                log.warning("Could not scan %s for dates: %s", p, e)

    # --- Infoclimat fallback for gaps -----------------------------------
    dates_missing_after_datagouv = dates_needed - datagouv_covered
    infoclimat_paths = []
    infoclimat_covered = set()
    still_missing = set()

    for date_str in sorted(dates_missing_after_datagouv):
        p = fetch_infoclimat_day(date_str, cache_dir=RAW)
        if p:
            infoclimat_paths.append(str(p))
            infoclimat_covered.add(date_str)
        else:
            still_missing.add(date_str)

    covered = datagouv_covered | infoclimat_covered

    manifest = {
        "start": start,
        "end": end,
        "source_used": (
            "datagouv+infoclimat" if (datagouv_paths and infoclimat_covered) else
            "datagouv" if datagouv_paths else
            "infoclimat" if infoclimat_covered else
            "none"
        ),
        "datagouv_files": [_relpath(p) for p in datagouv_paths],
        "infoclimat_files": [_relpath(p) for p in infoclimat_paths],
        "covered_dates": sorted(covered),
        "missing": sorted(still_missing),
    }

    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    log.info("Manifest written to %s", manifest_path)
    log.info("Covered %d/%d dates; missing %d",
             len(covered), len(dates_needed), len(still_missing))
    return manifest
