"""
physics/radiation.py
====================
Pure functions for computing solar radiation components and net radiation.
All energy in MJ m-2 per period (hourly or daily).

References:
  Allen, R.G. et al. (1998). FAO Irrigation and Drainage Paper 56 (FAO-56).
  Equations numbered per FAO-56 unless otherwise noted.

Longitude sign convention (IMPORTANT -- reviewer fix #5):
  FAO-56 uses Lz = longitude of the centre of the time zone (0 for UTC),
  and Lm = longitude of the measurement site in DEGREES WEST for west-of-UTC.
  Niort-Souche is at 0.40142 deg W -> Lm = +0.40142 (positive west in FAO convention).
  In this code, LON_DEG (from config) is stored as -0.40142 (negative = west, geographic convention).
  Inside Ra_hourly we therefore use: b = -LON_DEG (i.e. +0.40142) for the FAO Lm term.
  A unit test (test_radiation.py) asserts that sum(Ra_hourly, 24h) == Ra_daily within 3%.
"""

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
from config import GSC, SIGMA, ALBEDO, RSO_A, RSO_B, A_S, B_S, ALT_M, LON_DEG

LAT_RAD = math.radians(46.3114)


def inv_rel_distance(doy):
    """FAO-56 eq.23: inverse relative distance Earth-Sun."""
    return 1.0 + 0.033 * np.cos(2 * np.pi / 365 * doy)


def solar_declination(doy):
    """FAO-56 eq.24: solar declination (radians)."""
    return 0.409 * np.sin(2 * np.pi / 365 * doy - 1.39)


def sunset_hour_angle(lat_rad, delta):
    """FAO-56 eq.25: sunset hour angle (radians)."""
    arg = -np.tan(lat_rad) * np.tan(delta)
    arg = np.clip(arg, -1.0, 1.0)
    return np.arccos(arg)


def Ra_daily(doy, lat_rad=LAT_RAD):
    """
    FAO-56 eq.21: daily extraterrestrial radiation [MJ m-2 day-1].
    doy: day of year (1-365)
    """
    dr = inv_rel_distance(doy)
    delta = solar_declination(doy)
    ws = sunset_hour_angle(lat_rad, delta)
    ra = (
        24.0 * 60.0 / np.pi * GSC * dr
        * (ws * np.sin(lat_rad) * np.sin(delta)
           + np.cos(lat_rad) * np.cos(delta) * np.sin(ws))
    )
    return np.maximum(ra, 0.0)


def Ra_hourly(doy, lat_rad, t_mid_hour_utc, lon_deg=LON_DEG):
    """
    FAO-56 eqs 28-33: hourly extraterrestrial radiation [MJ m-2 h-1].

    doy:             day of year (1-365)
    lat_rad:         latitude in radians
    t_mid_hour_utc:  mid-point of the hour in UTC decimal hours (e.g. 13.5 for 13:00-14:00)
    lon_deg:         longitude in geographic convention (negative = west)

    FAO longitude convention (positive west):
      Lm = -lon_deg  (converts geographic negative-west to FAO positive-west)
      Lz = 0         (UTC time zone)
    """
    # FAO-56 eq.32: seasonal correction Sc (hours)
    b_fao = 2 * np.pi * (doy - 81) / 364.0
    Sc = 0.1645 * np.sin(2 * b_fao) - 0.1255 * np.cos(b_fao) - 0.025 * np.sin(b_fao)

    # FAO-56 eq.31: solar time (hours) from standard time
    # Lm is longitude WEST of Greenwich (positive west); Lz=0 for UTC
    Lm = -lon_deg  # geographic deg E -> FAO deg W
    Lz = 0.0
    t_solar = t_mid_hour_utc + (Lm - Lz) / 15.0 + Sc

    # FAO-56 eq.29: hour angle at mid-point of the hour (radians)
    omega = np.pi / 12.0 * (t_solar - 12.0)

    # Hour angles at start and end of hour interval (FAO-56 eq.30)
    omega1 = omega - np.pi / 24.0
    omega2 = omega + np.pi / 24.0

    dr = inv_rel_distance(doy)
    delta = solar_declination(doy)

    # Clip to daylight interval (FAO-56 eq.28)
    ws = sunset_hour_angle(lat_rad, delta)
    omega1 = np.clip(omega1, -ws, ws)
    omega2 = np.clip(omega2, -ws, ws)
    omega1 = np.minimum(omega1, omega2)  # ensure omega1 <= omega2

    ra = (
        12.0 * 60.0 / np.pi * GSC * dr
        * ((omega2 - omega1) * np.sin(lat_rad) * np.sin(delta)
           + np.cos(lat_rad) * np.cos(delta) * (np.sin(omega2) - np.sin(omega1)))
    )
    return np.maximum(ra, 0.0)


def Rso_fn(Ra, z=ALT_M):
    """FAO-56 eq.37: clear-sky solar radiation [same units as Ra]."""
    return (RSO_A + RSO_B * z) * Ra


def es_kpa(T_C):
    """FAO-56 eq.11: saturation vapour pressure [kPa]."""
    return 0.6108 * np.exp(17.27 * T_C / (T_C + 237.3))


def ea_kpa(T_C, RH_pct):
    """FAO-56 eq.17: actual vapour pressure [kPa]."""
    return es_kpa(T_C) * RH_pct / 100.0


def Rns(Rs_MJ):
    """FAO-56 eq.38: net shortwave radiation [MJ m-2 per period]."""
    return (1.0 - ALBEDO) * Rs_MJ


def Rnl_hourly(T_C, ea, Rs_MJ_h, Rso_MJ_h, _carry_ratio=None):
    """
    FAO-56 eq.39 adapted for hourly: net longwave radiation [MJ m-2 h-1].
    Positive value means upward (loss) -- Rn = Rns - Rnl.

    Night treatment (Rso_MJ_h ~ 0): the (Rs/Rso) cloud-correction term
    cannot be computed directly; per FAO-56 recommendation we carry forward
    the cloud-correction from the last daylight hour. If _carry_ratio is
    provided (scalar or array) it is used for those night hours.
    When vectorised, the caller should pass Rs/Rso for each hour (with the
    night values substituted by carry-forward) directly as _carry_ratio.
    """
    sigma_h = SIGMA / 24.0  # MJ m-2 h-1 K-4
    T_K = T_C + 273.16

    ea_arr = np.asarray(ea, dtype=float)
    f_ea   = 0.34 - 0.14 * np.sqrt(np.maximum(ea_arr, 0.0))

    if _carry_ratio is not None:
        f_sun = 1.35 * np.asarray(_carry_ratio, dtype=float) - 0.35
    else:
        # direct computation (daytime only; night will give odd values)
        rso_safe = np.where(Rso_MJ_h > 0, Rso_MJ_h, 1e-6)
        ratio = np.minimum(Rs_MJ_h / rso_safe, 1.0)
        f_sun = 1.35 * ratio - 0.35

    f_sun = np.clip(f_sun, 0.05, 1.0)  # physical bounds
    rnl = sigma_h * T_K**4 * f_ea * f_sun
    return np.maximum(rnl, 0.0)


def Rn_hourly_series(T_C_arr, RH_arr, Rs_Wm2_arr, doy_arr, hour_arr, lat_rad=LAT_RAD):
    """
    Vectorised net radiation over an array of hourly observations.

    Inputs (1-D numpy arrays, same length):
      T_C_arr    -- temperature (deg C)
      RH_arr     -- relative humidity (%)
      Rs_Wm2_arr -- global radiation (W/m2)
      doy_arr    -- day of year (int)
      hour_arr   -- hour of day UTC (int, 0-23); mid-hour = hour + 0.5 used inside

    Returns Rn_MJ_arr [MJ m-2 h-1], may be negative at night.
    """
    n = len(T_C_arr)
    Rs_MJ   = Rs_Wm2_arr * 3600.0 / 1e6  # W/m2 -> MJ/m2/h
    Ra_arr  = np.array([Ra_hourly(doy_arr[i], lat_rad, hour_arr[i] + 0.5) for i in range(n)])
    Rso_arr = Rso_fn(Ra_arr)
    ea_arr  = ea_kpa(T_C_arr, RH_arr)

    # Build carry-forward cloud-correction ratio
    rso_safe = np.where(Rso_arr > 0, Rso_arr, np.nan)
    ratio_raw = np.where(Rso_arr > 0, np.minimum(Rs_MJ / rso_safe, 1.0), np.nan)

    # Forward-fill NaN (night) with last daylight value; if no prior day, use 0.5
    carry_ratio = np.full(n, 0.5)
    last_valid = 0.5
    for i in range(n):
        if not np.isnan(ratio_raw[i]):
            last_valid = ratio_raw[i]
        carry_ratio[i] = last_valid

    rnl_arr = Rnl_hourly(T_C_arr, ea_arr, Rs_MJ, Rso_arr, _carry_ratio=carry_ratio)
    rns_arr = Rns(Rs_MJ)
    rn_arr  = rns_arr - rnl_arr  # may be negative at night
    return rn_arr
