"""
tests/test_penman.py
====================
Unit tests for physics/penman.py.
Hand-computed FAO-56 style worked example with stated tolerance.

Reference calculation (documented here, not from external source):
  Inputs: T=25 degC, ea=1.5 kPa -> RH = ea/es(25)*100, u2=2 m/s,
          Rn=10 MJ/m2/day, z=0 m.
  es(25) = 0.6108 * exp(17.27*25/(25+237.3)) = 3.1674 kPa
  RH = 1.5/3.1674*100 = 47.36 %
  Delta = 4098*3.1674/(25+237.3)^2 = 0.18882 kPa/degC
  gamma = 0.000665 * 101.3 * ((293-0)/293)^5.26 = 0.0674 kPa/degC
  lambda = 2.501 - 0.002361*25 = 2.4419 MJ/kg
  Erad = (0.18882/(0.18882+0.0674)) * (10/2.4419) = 0.7369 * 4.095 = 3.018 mm/day
  Ea = 0.26*(1+0.54*2)*(3.1674-1.5) = 0.541*1.6*1.6674 = ... let us compute step by step
  f_wind(2) = 0.26*(1+1.08) = 0.26*2.08 = 0.5408 mm/day/kPa
  vpd = 3.1674 - 1.5 = 1.6674 kPa
  Ea_full = 0.5408*1.6674 = 0.9015 mm/day
  Eaero = (0.0674/0.25622) * 0.9015 = 0.2630 * 0.9015 = 0.2372 mm/day
  E_total = 3.018 + 0.2372 = 3.255 mm/day

Tolerance: +/-3% of 3.255 = +/-0.098 mm/day, so abs(result - 3.255) < 0.10.
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from physics.penman import penman_E, wind_2m, f_wind, delta_svp, psychrometric
from physics.radiation import es_kpa, ea_kpa


def test_hand_computed_example():
    """
    Hand-computed reference: T=25 C, ea=1.5 kPa, u2=2 m/s, Rn=10 MJ/day, z=0.
    Tolerance: abs(E - 3.255) < 0.10 mm/day (3%).
    """
    T = 25.0
    ea_target = 1.5  # kPa
    es = float(es_kpa(T))
    RH = ea_target / es * 100.0
    u2 = 2.0
    Rn_MJ = 10.0
    z = 0.0

    E = penman_E(Rn_MJ=Rn_MJ, T_C=T, RH_pct=RH, u2=u2, period_hours=24, z=z)
    expected = 3.255
    tolerance = 0.10

    print(f"E computed: {E:.4f} mm/day  (expected ~{expected})")
    assert abs(E - expected) < tolerance, (
        f"Hand-computed Penman: E={E:.4f}, expected {expected} +/-{tolerance}"
    )


def test_wind_2m_conversion():
    """u2 = u10 * 4.87 / ln(67.8*10 - 5.42). For u10=10 -> ~7.480; u10=1 -> ~0.748."""
    factor = 4.87 / math.log(67.8 * 10 - 5.42)
    assert abs(wind_2m(10.0) - 10.0 * factor) < 1e-9
    assert abs(wind_2m(1.0) - 1.0 * factor) < 1e-9
    print(f"wind_2m(10)={wind_2m(10.0):.4f}, wind_2m(1)={wind_2m(1.0):.4f}, factor={factor:.4f}")


def test_night_condensation_negative():
    """Night case: negative Rn + high RH -> condensation returned as negative E."""
    E = penman_E(Rn_MJ=-0.3, T_C=18.0, RH_pct=90.0, u2=1.0, period_hours=1, z=62.0)
    print(f"Night E (should be < 0): {E:.6f} mm")
    assert E < 0, f"Expected negative E (condensation) at night, got {E}"


def test_positive_daytime():
    """Sunny day with strong radiation -> positive evaporation."""
    E = penman_E(Rn_MJ=2.0, T_C=30.0, RH_pct=30.0, u2=3.0, period_hours=1)
    print(f"Daytime E (should be > 0): {E:.6f} mm/h")
    assert E > 0, f"Expected positive E during daytime, got {E}"


def test_intermediate_values():
    """Check intermediate physics values."""
    T = 25.0
    Delta = float(delta_svp(T))
    gamma = psychrometric(0)
    es = float(es_kpa(T))
    print(f"es(25)={es:.4f} kPa (expect ~3.167)")
    print(f"Delta(25)={Delta:.5f} kPa/C (expect ~0.1888)")
    print(f"gamma(z=0)={gamma:.5f} kPa/C (expect ~0.0674)")
    assert abs(es - 3.1674) < 0.01
    assert abs(Delta - 0.1888) < 0.002
    assert abs(gamma - 0.0674) < 0.001


if __name__ == "__main__":
    test_hand_computed_example()
    test_wind_2m_conversion()
    test_night_condensation_negative()
    test_positive_daytime()
    test_intermediate_values()
    print("All penman tests passed.")
