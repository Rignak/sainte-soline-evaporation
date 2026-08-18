"""
analysis.py
===========
Main analysis pipeline for Sainte-Soline open-water evaporation.

Two timeframes:
  T1: 2026-06-15 00:00 UTC -> 2026-07-20 23:00 UTC
  T2: 2026-06-15 00:00 UTC -> 2026-08-15 23:00 UTC

Computation:
  Penman is applied to each hourly record independently (hourly evaporation),
  then hourly values are summed to daily and cumulative totals.

Note on condensation:
  Negative hourly E (condensation / dew) is returned as-is from penman_E.
  Daily and cumulative totals include negatives.
  A second 'clipped' cumulative (E>=0 only) is reported alongside, clearly labelled.

Note on incomplete days (reviewer fix #2):
  UTC days with fewer than 24 valid input hours are flagged in the daily CSV
  with 'n_valid_hours' < 24. Their daily E is marked as partial and excluded from
  the mean daily E computation (to avoid downward bias from treating zero-padding
  as real data). They remain in the daily CSV for transparency.
"""

import json
import logging
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import pandas as pd

ROOT = Path(__file__).parent.resolve()
sys.path.insert(0, str(ROOT))

from config import (AREA_M2, T1_START, T1_END, T2_START, T2_END,
                    RESULTS, ALT_M, LAT_DEG, E_MIN_SANITY, E_MAX_SANITY)
from data_io.download import ensure_raw
from data_io.parse import load_hourly
from physics.penman import penman_E, wind_2m
from physics.radiation import Ra_hourly, Rn_hourly_series, Ra_daily, Rso_fn, es_kpa, ea_kpa, Rnl_hourly

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

LAT_RAD = math.radians(LAT_DEG)
RESULTS.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Plot styling helpers
# ---------------------------------------------------------------------------

def _style_axes(ax, minor=False):
    """Apply consistent grid styling to an axes: grid behind data, major
    (+ optional minor) gridlines."""
    ax.set_axisbelow(True)
    ax.grid(True, which="major", alpha=0.35, linewidth=0.6)
    if minor:
        ax.minorticks_on()
        ax.grid(True, which="minor", alpha=0.15, linewidth=0.4)


def _style_date_axis(ax, fig):
    """Readable, non-overlapping date ticks."""
    locator = mdates.AutoDateLocator()
    formatter = mdates.ConciseDateFormatter(locator)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)
    fig.autofmt_xdate()


# ---------------------------------------------------------------------------
# Hourly pipeline
# ---------------------------------------------------------------------------

def compute_hourly(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute hourly Penman evaporation and radiation components.
    Adds columns: u2, Rs_MJ, Ra_MJ, Rso_MJ, Rn_MJ, E_hourly_mm.
    Rows with any NaN in required inputs produce NaN E_hourly_mm.
    """
    out = df.copy()
    n = len(out)
    idx = out.index

    T   = out["T_C"].values.astype(float)
    RH  = out["RH_pct"].values.astype(float)
    u10 = out["u10_ms"].values.astype(float)
    Rs  = out["Rs_Wm2"].values.astype(float)

    # Wind at 2 m
    out["u2"] = wind_2m(u10)

    # Radiation
    doy_arr  = np.array([ts.timetuple().tm_yday for ts in idx])
    hour_arr = np.array([ts.hour for ts in idx])

    Ra_arr  = np.array([Ra_hourly(doy_arr[i], LAT_RAD, hour_arr[i] + 0.5) for i in range(n)])
    Rso_arr = Rso_fn(Ra_arr)
    Rs_MJ   = Rs * 3600.0 / 1e6

    out["Rs_MJ"]  = Rs_MJ
    out["Ra_MJ"]  = Ra_arr
    out["Rso_MJ"] = Rso_arr

    # Net radiation (vectorised, with carry-forward for night cloud-correction)
    Rn_arr = Rn_hourly_series(T, RH, Rs, doy_arr, hour_arr, LAT_RAD)
    out["Rn_MJ"] = Rn_arr

    # Penman per hour: period_hours=1
    E_arr = np.full(n, np.nan)
    for i in range(n):
        if any(math.isnan(v) for v in [T[i], RH[i], u10[i]]):
            continue
        E_arr[i] = penman_E(
            Rn_MJ=float(Rn_arr[i]),
            T_C=float(T[i]),
            RH_pct=float(RH[i]),
            u2=float(out["u2"].iloc[i]),
            period_hours=1,
            z=ALT_M,
        )
    out["E_hourly_mm"] = E_arr

    return out


# ---------------------------------------------------------------------------
# Daily aggregation
# ---------------------------------------------------------------------------

def aggregate_daily(hourly: pd.DataFrame) -> pd.DataFrame:
    """
    Sum hourly E to daily (UTC day).
    Reviewer fix #2: flag incomplete days (n_valid_hours < 24).
    Returns DataFrame indexed by date with E_mm, Rn_day_MJ, n_valid_hours, partial_day.
    """
    df = hourly.copy()
    df["date"] = df.index.date

    def agg_day(g):
        n_valid = int(g["E_hourly_mm"].notna().sum())
        E_sum   = g["E_hourly_mm"].sum(min_count=1)  # NaN if all NaN
        Rn_sum  = g["Rn_MJ"].sum(min_count=1)
        return pd.Series({
            "E_mm": E_sum,
            "Rn_day_MJ": Rn_sum,
            "n_valid_hours": n_valid,
            "partial_day": n_valid < 24,
        })

    daily = df.groupby("date").apply(agg_day)
    daily.index = pd.to_datetime(daily.index)
    return daily


# ---------------------------------------------------------------------------
# Per-timeframe aggregation
# ---------------------------------------------------------------------------

def per_timeframe(hourly: pd.DataFrame, start: str, end: str, label: str) -> dict:
    """
    Slice hourly data to [start, end], compute daily/cumulative evaporation,
    write CSVs.

    Reviewer fix #9: report actual T2 coverage (first/last valid hour, missing days).
    """
    t_start = pd.Timestamp(start)
    t_end   = pd.Timestamp(end)
    sub = hourly[(hourly.index >= t_start) & (hourly.index <= t_end)].copy()

    # --- Actual coverage ---
    valid_e = sub["E_hourly_mm"].dropna()
    first_valid = str(valid_e.index[0]) if len(valid_e) else "none"
    last_valid  = str(valid_e.index[-1]) if len(valid_e) else "none"

    # Missing days: UTC days with 0 valid E
    sub["date"] = sub.index.date
    valid_by_day = sub.groupby("date")["E_hourly_mm"].apply(lambda x: x.notna().sum())
    missing_days = sorted([str(d) for d in valid_by_day[valid_by_day == 0].index])

    # --- Daily aggregation ---
    daily = aggregate_daily(sub)

    daily["E_m3"] = daily["E_mm"] / 1000.0 * AREA_M2

    # Cumulative (including negatives)
    daily["cum_mm"] = daily["E_mm"].cumsum()
    daily["cum_m3"] = daily["E_m3"].cumsum()

    # Clipped cumulative (E >= 0 only, clearly labelled)
    daily["E_mm_clipped"] = daily["E_mm"].clip(lower=0)
    daily["cum_mm_clipped"] = daily["E_mm_clipped"].cumsum()

    # Write daily CSV
    daily.to_csv(RESULTS / f"daily_{label}.csv", float_format="%.4f")
    log.info("Written daily_%s.csv (%d rows)", label, len(daily))

    # Write cumulative CSV
    cum_cols = [c for c in daily.columns if c.startswith("cum")]
    daily[cum_cols].to_csv(RESULTS / f"cumulative_{label}.csv", float_format="%.4f")

    # --- Statistics on complete days only ---
    complete = daily[~daily["partial_day"]].copy() if "partial_day" in daily.columns else daily.copy()
    n_partial = int(daily["partial_day"].sum()) if "partial_day" in daily.columns else 0

    mean_e = float(complete["E_mm"].mean(skipna=True))
    total_mm = float(daily["E_mm"].sum(skipna=True))
    total_m3 = total_mm / 1000.0 * AREA_M2

    # Condensation stats (from hourly)
    neg_mask = sub["E_hourly_mm"] < 0
    condensation_hours = int(neg_mask.sum())
    condensation_mm    = float(sub.loc[neg_mask, "E_hourly_mm"].sum())

    summary = {
        "label": label,
        "start": start,
        "end": end,
        "actual_first_valid_hour": first_valid,
        "actual_last_valid_hour":  last_valid,
        "missing_days": missing_days,
        "n_days_total": len(daily),
        "n_days_partial": n_partial,
        "mean_daily_mm": round(mean_e, 3),
        "total_mm": round(total_mm, 2),
        "total_m3": round(total_m3, 1),
        "condensation": {
            "hours": condensation_hours,
            "total_mm": round(condensation_mm, 3),
        },
    }

    # Sanity check
    m = summary["mean_daily_mm"]
    if not math.isnan(m):
        if not (E_MIN_SANITY <= m <= E_MAX_SANITY):
            log.warning(
                "SANITY BAND EXCEEDED: %s mean daily E = %.2f mm/day "
                "(expected %.1f-%.1f mm/day). "
                "Check for data gaps or physics errors.",
                label, m, E_MIN_SANITY, E_MAX_SANITY
            )

    return summary, daily


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def make_plots(hourly: pd.DataFrame, daily_T1: pd.DataFrame, daily_T2: pd.DataFrame,
               summary_T1: dict, summary_T2: dict):
    """Produce all required figures."""

    # 1. Hourly time series (T2 window)
    fig, ax = plt.subplots(figsize=(14, 4.5), dpi=140)
    e = hourly["E_hourly_mm"]
    ax.plot(hourly.index, e, lw=0.5, color="steelblue", label="Hourly evaporation")
    # Mark condensation (negative) hours in plain words
    neg = hourly[hourly["E_hourly_mm"] < 0]
    if not neg.empty:
        ax.scatter(neg.index, neg["E_hourly_mm"], s=1, color="purple",
                   label=f"Condensation hours (n={len(neg)})")
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xlabel("Date (UTC)")
    ax.set_ylabel("Evaporation (mm/h)")
    ax.set_title("Hourly open-water evaporation — Niort-Souche (WMO 07330), Sainte-Soline")
    ax.legend(fontsize=8, loc="upper right")
    _style_axes(ax, minor=True)
    _style_date_axis(ax, fig)
    fig.tight_layout()
    fig.savefig(RESULTS / "fig_hourly_timeseries.png", dpi=140)
    plt.close(fig)

    # 2 & 3. Daily evaporation + cumulative for each timeframe
    for label, daily, summary in [("T1", daily_T1, summary_T1), ("T2", daily_T2, summary_T2)]:
        dates = daily.index
        # Titles name the actual span; "T1"/"T2" alone is meaningless to a reader.
        span = f"{dates[0]:%d %b %Y} \u2013 {dates[-1]:%d %b %Y}"

        # Daily evaporation (single series, no legend needed)
        fig, ax = plt.subplots(figsize=(12, 4.5), dpi=140)
        ax.bar(dates, daily["E_mm"], width=0.8, color="steelblue", alpha=0.85)
        ax.axhline(0, color="k", lw=0.7)
        ax.set_ylabel("Evaporation (mm/day)")
        ax.set_xlabel("Date (UTC)")
        ax.set_title(f"Daily open-water evaporation — Sainte-Soline, {span}")
        _style_axes(ax, minor=True)
        _style_date_axis(ax, fig)
        fig.tight_layout()
        fig.savefig(RESULTS / f"fig_daily_{label}.png", dpi=140)
        plt.close(fig)

        # Cumulative evaporation, mm (left) and m3 (right, consistent scale)
        fig, ax = plt.subplots(figsize=(12, 4.5), dpi=140)
        ax.plot(dates, daily["cum_mm"], color="steelblue", lw=1.5)
        ax.set_ylabel("Cumulative evaporation (mm)")
        ax.set_xlabel("Date (UTC)")
        ax.set_title(f"Cumulative open-water evaporation — Sainte-Soline, {span}")
        _style_axes(ax, minor=True)

        ax_vol = ax.twinx()
        mn, mx = ax.get_ylim()
        ax_vol.set_ylim(mn / 1000.0 * AREA_M2, mx / 1000.0 * AREA_M2)
        ax_vol.set_ylabel("Cumulative volume (m3)")
        ax_vol.grid(False)

        _style_date_axis(ax, fig)
        fig.tight_layout()
        fig.savefig(RESULTS / f"fig_cumulative_{label}.png", dpi=140)
        plt.close(fig)

    log.info("All figures saved to %s", RESULTS)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    log.info("=== Sainte-Soline evaporation analysis ===")

    # Step 1: ensure data downloaded
    log.info("Ensuring raw data for T2 window (superset of T1)...")
    manifest = ensure_raw(T2_START, T2_END)
    log.info("Manifest source: %s, covered: %d dates, missing: %d",
             manifest["source_used"], len(manifest["covered_dates"]), len(manifest["missing"]))
    if manifest["missing"]:
        log.warning("Missing dates: %s", manifest["missing"])

    # Step 2: parse hourly data
    log.info("Loading hourly data for T2...")
    hourly_raw, gap_report = load_hourly(T2_START, T2_END)

    # Step 3: compute hourly evaporation
    log.info("Computing hourly evaporation...")
    hourly = compute_hourly(hourly_raw)

    # Save full hourly CSV
    hourly.to_csv(RESULTS / "hourly_evaporation.csv", float_format="%.6f")
    log.info("Saved hourly_evaporation.csv (%d rows)", len(hourly))

    # Condensation global report
    neg_mask = hourly["E_hourly_mm"] < 0
    n_neg = int(neg_mask.sum())
    sum_neg = float(hourly.loc[neg_mask, "E_hourly_mm"].sum())
    log.info("Condensation: %d hours, total = %.3f mm (negative E, NOT clipped)", n_neg, sum_neg)

    # Step 4: per-timeframe
    log.info("Computing T1...")
    summary_T1, daily_T1 = per_timeframe(hourly, T1_START, T1_END, "T1")
    log.info("Computing T2...")
    summary_T2, daily_T2 = per_timeframe(hourly, T2_START, T2_END, "T2")

    # Consistency check: T1 head of T2
    t1_dates = set(str(d.date()) for d in daily_T1.index)
    t2_dates = set(str(d.date()) for d in daily_T2.index)
    assert t1_dates.issubset(t2_dates), "T1 dates not subset of T2"
    cum_t1_end = float(daily_T1["cum_mm"].iloc[-1])
    t1_end_in_t2 = daily_T2.loc[daily_T2.index == daily_T1.index[-1], "cum_mm"]
    if not t1_end_in_t2.empty:
        cum_t2_at_t1end = float(t1_end_in_t2.iloc[0])
        diff = abs(cum_t1_end - cum_t2_at_t1end)
        log.info("Consistency check: T1 final cum=%.3f mm, T2 at T1-end=%.3f mm, diff=%.4f",
                 cum_t1_end, cum_t2_at_t1end, diff)
        assert diff < 0.01, f"T1/T2 cumulative inconsistency: {diff}"

    assert summary_T2["total_mm"] >= summary_T1["total_mm"], \
        "T2 cumulative should be >= T1"

    # Step 5: plots
    log.info("Generating plots...")
    make_plots(hourly, daily_T1, daily_T2, summary_T1, summary_T2)

    # Step 6: write summary.json
    summary = {
        "gap_report": gap_report,
        "data_manifest": {
            "source_used": manifest["source_used"],
            "n_covered_dates": len(manifest["covered_dates"]),
            "n_missing_dates": len(manifest["missing"]),
            "missing_dates": manifest["missing"],
        },
        "T1": summary_T1,
        "T2": summary_T2,
        "condensation_T2": {
            "hours": n_neg,
            "total_mm": round(sum_neg, 3),
            "note": "Negative E (condensation/dew) retained in totals; not clipped.",
        },
    }

    with open(RESULTS / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    log.info("summary.json written")

    # Print headline numbers
    print("\n=== HEADLINE RESULTS ===")
    for label, s in [("T1", summary_T1), ("T2", summary_T2)]:
        print(f"\n  [{label}] {s['start']} -> {s['end']}")
        print(f"    Actual coverage: {s['actual_first_valid_hour']} -> {s['actual_last_valid_hour']}")
        if s["missing_days"]:
            print(f"    Missing days: {s['missing_days']}")
        print(f"    Mean daily E = {s['mean_daily_mm']:.3f} mm/day")
        print(f"    Total E      = {s['total_mm']:.2f} mm  = {s['total_m3']:.0f} m3")
        print(f"    Condensation: {s['condensation']['hours']} h, {s['condensation']['total_mm']:.3f} mm")

    print(f"\n  Gap report: {gap_report}")
    print(f"\n  Result files: {RESULTS}")

    return summary


if __name__ == "__main__":
    main()
