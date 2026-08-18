"""
tests/test_radiation.py
=======================
Unit tests for physics/radiation.py.
Validates Ra_daily, Ra_hourly (sum vs daily), and Rn sign at night/day.
"""

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from physics.radiation import (
    Ra_daily, Ra_hourly, Rn_hourly_series, solar_declination, inv_rel_distance
)

LAT_RAD = math.radians(46.3114)
LON_DEG = -0.40142
DOY_JUNE15 = 166  # approx day of year for 2026-06-15


def test_Ra_daily_midsummer():
    """Ra at lat 46.3 N, doy 166 should be ~41-42 MJ/m2/day."""
    Ra = float(Ra_daily(DOY_JUNE15, LAT_RAD))
    print(f"Ra_daily(doy={DOY_JUNE15}, lat=46.3N) = {Ra:.2f} MJ/m2/day")
    assert 38 < Ra < 44, f"Ra_daily out of expected range: {Ra}"


def test_Ra_hourly_sum_vs_daily():
    """Sum of 24 Ra_hourly values should match Ra_daily within 3%."""
    Ra_d = float(Ra_daily(DOY_JUNE15, LAT_RAD))
    Ra_h_sum = sum(
        Ra_hourly(DOY_JUNE15, LAT_RAD, h + 0.5, LON_DEG)
        for h in range(24)
    )
    print(f"Ra_daily={Ra_d:.4f}  Ra_hourly_sum={Ra_h_sum:.4f}")
    rel_diff = abs(Ra_h_sum - Ra_d) / Ra_d
    assert rel_diff < 0.03, f"Ra_hourly sum vs daily: relative diff {rel_diff:.4f} > 3%"


def test_Ra_hourly_night_zero():
    """Ra_hourly at midnight should be 0 (sun below horizon)."""
    Ra_night = Ra_hourly(DOY_JUNE15, LAT_RAD, 0.5, LON_DEG)  # 00:00-01:00 UTC
    print(f"Ra_hourly at 00:30 UTC = {Ra_night:.6f} (expect 0)")
    assert Ra_night == 0.0, f"Ra at midnight should be 0, got {Ra_night}"


def test_Rn_sign_night_positive_day():
    """
    Night: Rn should be negative (net longwave loss, Rs=0).
    Peak day: Rn should be positive (>1.5 MJ/m2/h with Rs~681 W/m2).
    """
    import numpy as np

    # Single night hour
    T_night = np.array([18.0])
    RH_night = np.array([70.0])
    Rs_night = np.array([0.0])
    doy = np.array([DOY_JUNE15])
    hour_night = np.array([1])  # 01:00 UTC

    Rn_night = Rn_hourly_series(T_night, RH_night, Rs_night, doy, hour_night, LAT_RAD)
    print(f"Rn at night (01:00 UTC): {Rn_night[0]:.4f} MJ/m2/h (expect < 0)")
    assert Rn_night[0] < 0, f"Expected negative Rn at night, got {Rn_night[0]}"

    # Peak daytime hour (~14:00 UTC) with high radiation
    T_day = np.array([33.7])
    RH_day = np.array([15.0])
    Rs_day = np.array([681.0])  # W/m2 from sample data
    hour_day = np.array([14])   # 14:00-15:00 UTC

    Rn_day = Rn_hourly_series(T_day, RH_day, Rs_day, doy, hour_day, LAT_RAD)
    print(f"Rn at peak day (14:00 UTC, Rs=681 W/m2): {Rn_day[0]:.4f} MJ/m2/h (expect >1.5)")
    assert Rn_day[0] > 1.5, f"Expected Rn > 1.5 MJ/m2/h at peak day, got {Rn_day[0]}"


def test_seasonal_correction_sign():
    """
    Verify solar time correction is sensible:
    Niort is 0.40 deg W, so solar noon is ~1.6 min after UTC noon -- tiny effect.
    Test that Ra_hourly at hour 12 > Ra_hourly at hour 6 on mid-summer day.
    """
    Ra_noon = Ra_hourly(DOY_JUNE15, LAT_RAD, 12.5, LON_DEG)
    Ra_morning = Ra_hourly(DOY_JUNE15, LAT_RAD, 6.5, LON_DEG)
    print(f"Ra at 12:30 UTC = {Ra_noon:.4f}, at 06:30 UTC = {Ra_morning:.4f}")
    assert Ra_noon > Ra_morning, "Ra should be higher at solar noon than morning"


if __name__ == "__main__":
    test_Ra_daily_midsummer()
    test_Ra_hourly_sum_vs_daily()
    test_Ra_hourly_night_zero()
    test_Rn_sign_night_positive_day()
    test_seasonal_correction_sign()
    print("All radiation tests passed.")
