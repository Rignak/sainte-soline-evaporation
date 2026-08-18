"""
data_io/parse.py
================
Parse raw cached files (datagouv csv.gz and Infoclimat HTML) into a
canonical tidy hourly DataFrame, then assemble and quality-control it.

Canonical schema (one row per UTC hour, DatetimeIndex named 'time'):
  T_C      -- air temperature, degrees Celsius
  RH_pct   -- relative humidity, percent
  u10_ms   -- wind speed at 10 m, m/s
  Rs_Wm2   -- global (shortwave) radiation, W/m2 (may be NaN at night)
  N_frac   -- cloud cover fraction 0-1 (if available, else NaN)
  source   -- 'datagouv' or 'infoclimat'

Gap handling:
  * Reindex to complete hourly DatetimeIndex.
  * Linearly interpolate gaps <= 3 h for T_C, RH_pct, u10_ms.
  * For Rs_Wm2: interpolate gaps <= 2 h ONLY during daylight; force 0 when
    extraterrestrial radiation Ra = 0 (night check per latitude/doy/hour).
  * Gaps > those limits are left as NaN and REPORTED, not silently filled.
  * gap_fraction per field = (interpolated + NaN) / total rows, returned in
    gap_report dict and written to results/gap_report.json.
"""

import gzip
import json
import logging
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))
from config import RAW, RESULTS, STATION_WMO, STATION_NATIONAL, LAT_DEG

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

# Canonical column names
COLS = ["T_C", "RH_pct", "u10_ms", "Rs_Wm2", "N_frac"]


# ---------------------------------------------------------------------------
# parse_datagouv
# ---------------------------------------------------------------------------

def parse_datagouv(path: Path) -> pd.DataFrame:
    """
    Parse one data.gouv.fr H_79_*.csv.gz file.
    Filters to STATION_WMO.
    Returns tidy DataFrame with 'time' as DatetimeIndex (UTC).

    Key columns in source (semicolon-delimited, Latin-1):
      NUM_POSTE or POSTE -- station id (7-digit, matching 0733000x)
      AAAAMMJJHH         -- UTC timestamp YYYYMMDDHH
      T                  -- temperature (deg C, or K in some versions)
      U                  -- relative humidity (%)
      FF                 -- wind speed at 10 m (m/s, factor 1/10 in some files)
      GLO                -- global radiation (J/cm2 per hour; 1 J/cm2 = 0.01 MJ/m2)
      N                  -- cloud cover (oktas 0-8)
    """
    log.info("Parsing datagouv %s", path)
    try:
        with gzip.open(path, "rt", encoding="latin-1") as f:
            df_raw = pd.read_csv(f, sep=";", dtype=str, low_memory=False)
    except Exception as e:
        log.error("Failed to read %s: %s", path, e)
        return pd.DataFrame()

    log.info("  Columns: %s", list(df_raw.columns)[:20])

    # Identify station column
    sta_col = None
    for c in ["NUM_POSTE", "POSTE", "numer_sta"]:
        if c in df_raw.columns:
            sta_col = c
            break
    if sta_col is None:
        log.error("No station column found in %s", path)
        return pd.DataFrame()

    # Filter station: try both WMO id (07330) and national id (79191005)
    mask = (df_raw[sta_col].str.contains(STATION_WMO, na=False) |
            df_raw[sta_col].str.contains(STATION_NATIONAL, na=False))
    df = df_raw[mask].copy()
    if df.empty:
        log.warning("No rows for station %s in %s", STATION_WMO, path)
        return pd.DataFrame()
    log.info("  Rows for station: %d", len(df))

    # Date column
    date_col = None
    for c in ["AAAAMMJJHH", "DATE", "date"]:
        if c in df.columns:
            date_col = c
            break
    if date_col is None:
        log.error("No date column in %s", path)
        return pd.DataFrame()

    # Parse timestamp: YYYYMMDDHH -> UTC
    def parse_ts(s):
        s = str(s).strip()
        if len(s) >= 10:
            try:
                return pd.Timestamp(s[:4] + "-" + s[4:6] + "-" + s[6:8] + " " + s[8:10] + ":00")
            except Exception:
                pass
        return pd.NaT

    df["time"] = df[date_col].map(parse_ts)
    df = df.dropna(subset=["time"]).set_index("time")

    out = pd.DataFrame(index=df.index)

    # Temperature: T in degC (some files store K, check range)
    if "T" in df.columns:
        t_vals = pd.to_numeric(df["T"], errors="coerce")
        # If mean > 100 assume Kelvin -> convert
        if t_vals.dropna().mean() > 100:
            t_vals = t_vals - 273.15
        out["T_C"] = t_vals
    else:
        out["T_C"] = np.nan

    # Humidity
    if "U" in df.columns:
        out["RH_pct"] = pd.to_numeric(df["U"], errors="coerce")
    else:
        out["RH_pct"] = np.nan

    # Wind at 10 m in m/s
    # FF may be stored as tenths (divide by 10) or as m/s directly
    if "FF" in df.columns:
        ff = pd.to_numeric(df["FF"], errors="coerce")
        # Heuristic: if median > 30, likely tenths
        if ff.dropna().median() > 30:
            ff = ff / 10.0
        out["u10_ms"] = ff
    else:
        out["u10_ms"] = np.nan

    # Global radiation: GLO in J/cm2/h -> W/m2 = (J/cm2 * 10000 J/m2 per J/cm2) / 3600 s
    glo_col = None
    for c in ["GLO", "GLO2", "glo"]:
        if c in df.columns:
            glo_col = c
            break
    if glo_col:
        glo = pd.to_numeric(df[glo_col], errors="coerce")
        out["Rs_Wm2"] = glo * 10000.0 / 3600.0
    else:
        out["Rs_Wm2"] = np.nan

    # Cloud cover: N in oktas 0-8 -> N_frac 0-1
    if "N" in df.columns:
        n_raw = pd.to_numeric(df["N"], errors="coerce")
        # oktas 0-8
        out["N_frac"] = n_raw / 8.0
    else:
        out["N_frac"] = np.nan

    out["source"] = "datagouv"
    return out.sort_index()


# ---------------------------------------------------------------------------
# parse_infoclimat
# ---------------------------------------------------------------------------

def parse_infoclimat(path: Path) -> pd.DataFrame:
    """
    Parse one Infoclimat HTML file (one calendar day).
    Uses regex to extract the data table rows.
    Returns tidy DataFrame with 'time' as DatetimeIndex (UTC).

    Expected Infoclimat columns (may vary):
      Heure  -- time as 'HH:MM' (UTC in the archive pages)
      Temperature -- degC
      Vent   -- 'NN km/h' (possibly with gust info)
      Humidite -- %
      Radiations -- W/m2
    """
    log.info("Parsing infoclimat %s", path)
    text = path.read_text(errors="replace")

    # Extract date from filename
    m = re.search(r"infoclimat_(\d{4}-\d{2}-\d{2})\.html", path.name)
    if not m:
        log.error("Cannot extract date from filename %s", path.name)
        return pd.DataFrame()
    date_str = m.group(1)

    # Try to parse HTML table with pandas
    try:
        import io
        tables = pd.read_html(io.StringIO(text), decimal=",", thousands=" ")
    except Exception as e:
        log.error("pd.read_html failed on %s: %s", path, e)
        return pd.DataFrame()

    if not tables:
        log.warning("No tables found in %s", path)
        return pd.DataFrame()

    # Find the table with the most numeric content (observation table)
    best = max(tables, key=lambda t: t.shape[0] * t.shape[1])
    log.info("  Using table shape %s", best.shape)
    log.info("  Columns: %s", list(best.columns))

    out_rows = []
    for _, row in best.iterrows():
        # Find time column
        time_val = None
        for col in best.columns:
            val = str(row[col])
            if re.match(r"^\d{1,2}[h:]\d{2}$", val.strip()):
                time_val = val.strip()
                break
            # sometimes just hour integer
            if re.match(r"^\d{1,2}$", val.strip()) and int(val.strip()) < 24:
                time_val = val.strip() + ":00"
                break

        if time_val is None:
            continue

        # Normalise time string
        time_val = time_val.replace("h", ":").replace("H", ":")
        parts = time_val.split(":")
        if len(parts) == 2:
            hh, mm = parts
        else:
            continue
        try:
            ts = pd.Timestamp(f"{date_str} {int(hh):02d}:{mm}")
        except Exception:
            continue

        rec = {"time": ts}

        # Try to extract T, RH, wind, radiation by scanning columns
        cols_lower = {str(c).lower(): c for c in best.columns}

        # Temperature
        for k in ["temp", "tempe", "temperature"]:
            for ck, cname in cols_lower.items():
                if k in ck:
                    v = pd.to_numeric(str(row[cname]).replace(",", "."), errors="coerce")
                    if pd.notna(v) and -50 < v < 60:
                        rec["T_C"] = v
                    break
            if "T_C" in rec:
                break

        # Humidity
        for ck, cname in cols_lower.items():
            if "humid" in ck or "hum" == ck:
                v = pd.to_numeric(str(row[cname]).replace(",", ".").replace("%", ""), errors="coerce")
                if pd.notna(v) and 0 <= v <= 100:
                    rec["RH_pct"] = v
                break

        # Wind: column may contain strings like '14 km/h', extract number
        for ck, cname in cols_lower.items():
            if "vent" in ck or "wind" in ck or "ff" == ck:
                raw_wind = str(row[cname])
                nums = re.findall(r"[\d]+(?:[.,]\d+)?", raw_wind)
                if nums:
                    v = float(nums[0].replace(",", "."))
                    # Convert km/h -> m/s
                    rec["u10_ms"] = v * 0.27778
                break

        # Radiation
        for ck, cname in cols_lower.items():
            if "radia" in ck or "glo" in ck or "solar" in ck or "solaire" in ck:
                v = pd.to_numeric(str(row[cname]).replace(",", "."), errors="coerce")
                if pd.notna(v) and v >= 0:
                    rec["Rs_Wm2"] = v
                break

        rec["N_frac"] = np.nan  # not typically in infoclimat HTML
        rec["source"] = "infoclimat"
        out_rows.append(rec)

    if not out_rows:
        log.warning("No parseable rows in %s", path)
        return pd.DataFrame()

    df = pd.DataFrame(out_rows).set_index("time").sort_index()
    for c in COLS:
        if c not in df.columns:
            df[c] = np.nan
    log.info("  Parsed %d rows from %s", len(df), path)
    return df


# ---------------------------------------------------------------------------
# load_hourly
# ---------------------------------------------------------------------------

def _night_mask(index: pd.DatetimeIndex, lat_rad: float) -> np.ndarray:
    """Return boolean mask True where extraterrestrial radiation Ra <= 0 (night)."""
    from physics.radiation import Ra_hourly
    import math
    lon_deg = -0.40142  # from config
    ra = np.array([
        Ra_hourly(
            ts.timetuple().tm_yday,
            lat_rad,
            ts.hour + 0.5,  # mid-hour UTC
            lon_deg,
        )
        for ts in index
    ])
    return ra <= 0.0


def load_hourly(start: str, end: str):
    """
    Assemble one UTC-hourly DataFrame for [start, end], from cached raw files.
    Returns (df, gap_report) where gap_report is a dict.

    Steps:
      1. Read all datagouv files from RAW; concatenate; filter to [start, end].
      2. Read all infoclimat HTML files from RAW for dates in [start, end].
      3. Merge: datagouv preferred; infoclimat fills gaps.
      4. Reindex to full hourly grid; interpolate short gaps; report.
    """
    import math

    t_start = pd.Timestamp(start)
    t_end   = pd.Timestamp(end)
    lat_rad = math.radians(LAT_DEG)

    # --- datagouv ---
    dg_frames = []
    for p in sorted(RAW.glob("*.csv.gz")):
        df = parse_datagouv(p)
        if df.empty:
            continue
        df = df[(df.index >= t_start) & (df.index <= t_end)]
        if not df.empty:
            dg_frames.append(df)

    # --- infoclimat ---
    ic_frames = []
    for p in sorted(RAW.glob("infoclimat_*.html")):
        m = re.search(r"infoclimat_(\d{4}-\d{2}-\d{2})\.html", p.name)
        if not m:
            continue
        day_ts = pd.Timestamp(m.group(1))
        if not (t_start.date() <= day_ts.date() <= t_end.date()):
            continue
        df = parse_infoclimat(p)
        if df.empty:
            continue
        df = df[(df.index >= t_start) & (df.index <= t_end)]
        if not df.empty:
            ic_frames.append(df)

    if not dg_frames and not ic_frames:
        raise RuntimeError(
            "No data found in data/raw/ for the requested period. "
            "Run ensure_raw() first."
        )

    # Combine datagouv frames
    dg = pd.concat(dg_frames).sort_index() if dg_frames else pd.DataFrame(columns=COLS + ["source"])
    # Remove duplicate timestamps (keep first)
    dg = dg[~dg.index.duplicated(keep="first")]

    ic = pd.concat(ic_frames).sort_index() if ic_frames else pd.DataFrame(columns=COLS + ["source"])
    ic = ic[~ic.index.duplicated(keep="first")]

    # Build full hourly index
    full_idx = pd.date_range(start=t_start, end=t_end, freq="h")

    # Merge: start from infoclimat, override with datagouv where available
    merged = ic.reindex(full_idx) if not ic.empty else pd.DataFrame(index=full_idx, columns=COLS + ["source"])
    if not dg.empty:
        for col in COLS:
            if col in dg.columns:
                merged.loc[dg.index.intersection(full_idx), col] = dg.loc[dg.index.intersection(full_idx), col]
        merged.loc[dg.index.intersection(full_idx), "source"] = "datagouv"

    # Record which hours came from each source before filling
    original_valid = {}
    for c in COLS:
        original_valid[c] = merged[c].notna()

    # ---- Gap filling ----
    # For T_C, RH_pct, u10_ms: linear interp up to 3 h
    for col in ["T_C", "RH_pct", "u10_ms"]:
        if col in merged.columns:
            merged[col] = merged[col].interpolate(method="time", limit=3, limit_direction="both")

    # For Rs_Wm2: force 0 at night, interpolate <= 2 h during day
    night = _night_mask(full_idx, lat_rad)
    merged["Rs_Wm2"] = pd.to_numeric(merged["Rs_Wm2"], errors="coerce")
    merged.loc[night, "Rs_Wm2"] = 0.0
    merged["Rs_Wm2"] = merged["Rs_Wm2"].interpolate(method="time", limit=2, limit_direction="both")
    merged.loc[night, "Rs_Wm2"] = 0.0  # re-enforce night zeros

    # N_frac: fill with NaN (not critical; only used as fallback)
    # Source column: fill forward
    merged["source"] = merged["source"].ffill().bfill().fillna("none")

    # ---- Gap report ----
    total = len(merged)
    gap_report = {}
    for c in COLS:
        if c not in merged.columns:
            gap_report[c] = {"gap_fraction": 1.0, "n_nan": total, "n_interpolated": 0}
            continue
        n_orig_valid = int(original_valid[c].sum())
        n_now_valid  = int(merged[c].notna().sum())
        n_interp     = max(0, n_now_valid - n_orig_valid)
        n_nan_final  = total - n_now_valid
        gap_report[c] = {
            "gap_fraction": round((total - n_orig_valid) / total, 4),
            "n_nan_remaining": n_nan_final,
            "n_interpolated": n_interp,
        }
        if n_nan_final > 0:
            log.warning("Field %s: %d hours remain NaN after interpolation (reported gap)", c, n_nan_final)

    gap_report["total_hours"] = total
    gap_report["sources_used"] = (
        ("datagouv" if not dg.empty else "") + ("+" if not dg.empty and not ic.empty else "") +
        ("infoclimat" if not ic.empty else "")
    )

    RESULTS.mkdir(parents=True, exist_ok=True)
    with open(RESULTS / "gap_report.json", "w") as f:
        json.dump(gap_report, f, indent=2)
    log.info("Gap report: %s", gap_report)

    # Ensure numeric
    for c in COLS:
        if c in merged.columns:
            merged[c] = pd.to_numeric(merged[c], errors="coerce")

    return merged, gap_report
