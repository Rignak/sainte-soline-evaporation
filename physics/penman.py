"""
physics/penman.py
=================
Pure functions implementing the Penman (1948) open-water combination equation.

Wind function used:
  f(u2) = 0.26 * (1 + 0.54 * u2)  [mm day-1 kPa-1]

Source / citation:
  Penman, H.L. (1948). Natural evaporation from open water, bare soil and grass.
  Proceedings of the Royal Society of London. Series A, 193(1032), 120-145.
  Wind function as formulated in:
  FAO Irrigation and Drainage Paper 24 (Doorenbos & Pruitt, 1977), and
  Valiantzas, J.D. (2006). Simplified versions for the Penman evaporation
  equation using routine weather data. J. Hydrol., 331, 690-702.

10 m -> 2 m wind conversion:
  u2 = u10 * 4.87 / ln(67.8 * 10 - 5.42)   [FAO-56 eq.47]
  = u10 * 4.87 / ln(672.58)
  = u10 * 4.87 / 6.5093
  ~ u10 * 0.7480

Assumptions documented (reviewer fix #6):
  - Psychrometric constant gamma is derived from altitude-modelled pressure
    (FAO-56 eq.7-8), NOT from hourly measured pressure. This is standard FAO-56
    practice. Station pressure (if available) could be substituted for each hour
    for marginal improvement; that path is not taken here.
  - Water heat storage in the reservoir body is NEGLECTED (stated assumption).
    This is appropriate for shallow/well-mixed reservoirs over short periods but
    introduces a phase shift at daily timescales.
  - Negative evaporation (condensation / dew) is RETURNED AS-IS (not clipped).
    Callers are responsible for reporting condensation statistics.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
from config import ALT_M
from physics.radiation import es_kpa, ea_kpa


def delta_svp(T_C):
    """
    FAO-56 eq.13: slope of saturation vapour pressure curve [kPa deg C-1].
    """
    return 4098.0 * es_kpa(T_C) / (T_C + 237.3) ** 2


def psychrometric(z=ALT_M):
    """
    FAO-56 eqs.7-8: psychrometric constant [kPa deg C-1].
    Derived from altitude-modelled atmospheric pressure.
    See module docstring for notes on using measured vs modelled pressure.
    """
    P = 101.3 * ((293.0 - 0.0065 * z) / 293.0) ** 5.26  # kPa
    return 0.000665 * P


def latent_heat(T_C):
    """Latent heat of vaporisation [MJ kg-1]. Approximation: lambda = 2.501 - 0.002361*T."""
    return 2.501 - 0.002361 * T_C


def wind_2m(u10):
    """
    FAO-56 eq.47: convert wind speed from 10 m to 2 m height.
    u2 = u10 * 4.87 / ln(67.8 * 10 - 5.42)
    """
    return u10 * 4.87 / np.log(67.8 * 10.0 - 5.42)


def f_wind(u2):
    """
    Penman (1948) wind function [mm day-1 kPa-1]:
      f(u2) = 0.26 * (1 + 0.54 * u2)
    Cited: Penman 1948 / FAO Paper 24 / Valiantzas 2006.
    """
    return 0.26 * (1.0 + 0.54 * u2)


def penman_E(Rn_MJ, T_C, RH_pct, u2, period_hours, z=ALT_M):
    """
    Penman (1948) open-water evaporation over `period_hours` [mm].

    Parameters
    ----------
    Rn_MJ       : net radiation for the period [MJ m-2]
    T_C         : air temperature [deg C]
    RH_pct      : relative humidity [%]
    u2          : wind speed at 2 m [m/s]
    period_hours: duration [hours] (1 for hourly, 24 for daily)
    z           : elevation [m], default from config

    Returns
    -------
    E [mm] -- NOT clipped; negative values = condensation (dew deposition).

    Physics notes:
      - Rn_MJ / lambda gives mm because 1 kg m-2 = 1 mm and 1 MJ = 1e6 J.
      - f(u2) is in mm day-1 kPa-1; scaled by (period_hours / 24) for sub-daily.
      - The aerodynamic term (Ea component) uses Penman's wind function, not PM.
      - Water heat storage G is assumed zero (stated assumption in module docstring).
    """
    Delta = delta_svp(T_C)
    gamma = psychrometric(z)
    lam   = latent_heat(T_C)
    _es   = es_kpa(T_C)
    _ea   = ea_kpa(T_C, RH_pct)
    vpd   = _es - _ea

    # Radiation term: [mm per period]
    Erad  = (Delta / (Delta + gamma)) * (Rn_MJ / lam)

    # Aerodynamic term: f(u2) in mm/day/kPa, scale to period
    Ea_per_period = f_wind(u2) * vpd * (period_hours / 24.0)

    # Aerodynamic contribution to E
    Eaero = (gamma / (Delta + gamma)) * Ea_per_period

    E = Erad + Eaero  # may be negative
    return E
