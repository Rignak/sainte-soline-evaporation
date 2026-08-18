"""
config.py — centralised constants for Sainte-Soline open-water evaporation.
All timestamps are UTC (tz-naive, documented as UTC throughout).
"""

from pathlib import Path

ROOT = Path(__file__).parent.resolve()
RAW  = ROOT / "data" / "raw"
RESULTS = ROOT / "results"

# --- Site ----------------------------------------------------------------
AREA_M2   = 160894.0          # Sainte-Soline reservoir surface, m²
STATION_WMO = "07330"         # Niort-Souche, Météo-France WMO id
LAT_DEG   = 46.3114           # ° N  (source: Météo-France station metadata)
LON_DEG   = -0.40142          # ° E  (negative = west)
ALT_M     = 62.0              # m above MSL

# --- Physics -------------------------------------------------------------
ALBEDO    = 0.08              # open-water albedo (FAO-56 p.45)
SIGMA     = 4.903e-9          # MJ K⁻⁴ m⁻² day⁻¹ (Stefan-Boltzmann, FAO-56 eq.39)
GSC       = 0.0820            # MJ m⁻² min⁻¹ solar constant (FAO-56 eq.21)
A_S, B_S  = 0.25, 0.50       # Angstrom coefficients (FAO-56 default)
RSO_A, RSO_B = 0.75, 2e-5    # clear-sky Rs0 = (RSO_A + RSO_B*z)*Ra  (FAO-56 eq.37)

# --- Timeframes ----------------------------------------------------------
T1_START = "2026-06-15 00:00"; T1_END = "2026-07-20 23:00"   # UTC
T2_START = "2026-06-15 00:00"; T2_END = "2026-08-15 23:00"   # UTC

# --- Acceptance sanity band (mm/day) ------------------------------------
E_MIN_SANITY = 2.0
E_MAX_SANITY = 8.0

# Meteo-France national station code (used in data.gouv.fr files)
# WMO 07330 maps to national code 79191005 (NIORT, dept 79)
STATION_NATIONAL = "79191005"
