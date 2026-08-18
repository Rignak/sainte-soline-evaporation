# Sainte-Soline Open-Water Evaporation — Penman 1948

Computation of open-water evaporation over the Sainte-Soline reservoir
(Deux-Sevres, France) using hourly meteorological data and the Penman (1948)
combination equation.

---

## Station and Site

| Parameter | Value |
|-----------|-------|
| Station | Niort-Souche (NIORT, Meteo-France) |
| WMO id | 07330 |
| National id | 79191005 |
| Latitude | 46.3114 deg N |
| Longitude | 0.40142 deg W |
| Altitude | 62 m a.s.l. |
| Reservoir area | 160 894 m2 |

---

## Timeframes

| Label | Start (UTC) | End (UTC) |
|-------|-------------|-----------|
| T1 | 2026-06-15 00:00 | 2026-07-20 23:00 |
| T2 | 2026-06-15 00:00 | 2026-08-15 23:00 |

---

## Data Source

**Primary source (used exclusively):** Meteo-France 'Donnees climatologiques
de base - horaires' via data.gouv.fr, dataset id
`6569b4473bedf2e7abad3b72`, file `HOR_departement_79_periode_2025-2026`
(H_79_latest-2025-2026.csv.gz). Hourly UTC observations.

The data.gouv.fr API (dataset 6569b27d… old id) returned HTTP 404; the
correct current dataset id is 6569b4473bedf2e7abad3b72. An Infoclimat HTML
fallback was implemented but not needed: data.gouv covered all 62 dates.

**Fields used:**

| Field | Unit in source | Conversion |
|-------|---------------|------------|
| T | deg C | direct |
| U | % | direct |
| FF | m/s (10 m) | -> 2 m via FAO-56 eq.47 |
| GLO | J/cm2/h | * 10000/3600 -> W/m2 |

**Gap report (T2, 1488 hours):**

| Field | Gap fraction | Interpolated hours | NaN remaining |
|-------|-------------|-------------------|---------------|
| T_C | 0.00 % | 0 | 0 |
| RH_pct | 0.00 % | 0 | 0 |
| u10_ms | 0.00 % | 0 | 0 |
| Rs_Wm2 | 0.00 % | 0 | 0 |
| N_frac | 100.00 % | 0 | 1488 |

N (cloud cover) is entirely absent in the data.gouv hourly file for this
station. It is only needed as an Angstrom fallback for Rs when GLO is missing;
since GLO is complete, N_frac = NaN causes no impact.

---

## Method

### Penman (1948) Combination Equation

```
E = [Delta * Rn/lambda + gamma * Ea] / (Delta + gamma)   [mm per period]
```

where:
- Delta = slope of saturation vapour pressure curve (FAO-56 eq.13)
- gamma = psychrometric constant (FAO-56 eqs.7-8, altitude-derived)
- lambda = latent heat of vaporisation (2.501 - 0.002361*T, MJ/kg)
- Rn = net radiation (MJ/m2/period), see below
- Ea = f(u2) * (es - ea) = aerodynamic evaporative demand

### Wind Function (cited)

```
f(u2) = 0.26 * (1 + 0.54 * u2)   [mm/day/kPa]
```

Source: Penman (1948) *Proc. R. Soc. A* 193:120, as expressed in FAO
Irrigation and Drainage Paper 24 (Doorenbos and Pruitt 1977) and
Valiantzas (2006, J. Hydrol. 331:690).

### Wind Height Conversion

10 m to 2 m: u2 = u10 * 4.87 / ln(67.8 * 10 - 5.42) = u10 * 0.7480
(FAO-56 eq.47).

### Net Radiation

Rn = Rns - Rnl (may be negative at night).

Rns = (1 - albedo) * Rs, albedo = 0.08 (open water, FAO-56 p.45).

Rnl computed from FAO-56 eq.39 (hourly form with sigma/24 and
carry-forward of cloud-correction for night hours).

Extraterrestrial radiation Ra: hourly formulation, FAO-56 eqs.28-33.
Longitude sign convention: Lm = 0.40142 deg W (positive west, FAO convention),
Lz = 0 (UTC).

Clear-sky radiation: Rso = (0.75 + 2e-5 * z) * Ra (FAO-56 eq.37).

Rs is taken directly from GLO (J/cm2/h -> W/m2) at every hour. Night hours
(Ra = 0) have Rs = 0 enforced; no Angstrom fallback was needed.

### Assumptions

1. **Water heat storage neglected.** The thermal storage term G = 0 is
   standard for short-period (daily/hourly) open-water Penman and is
   explicitly assumed here. It introduces a phase lag at daily scale.

2. **Night longwave carry-forward.** The cloud-correction factor (Rs/Rso) in
   the Rnl formula cannot be evaluated when Rso = 0 (night). The last valid
   daytime ratio is carried forward; this affects the sign/magnitude of nightly
   Rnl but Rn at night is dominated by the sigma*T^4 upward emission anyway.

3. **Condensation retained.** Negative hourly E (condensation / dewfall) is
   returned as-is and included in daily and cumulative totals. A second
   "clipped" cumulative (E >= 0 only) is reported alongside in the CSVs.

4. **Psychrometric constant.** Gamma is derived from altitude-modelled
   pressure (FAO-56 eq.8), not from hourly measured station pressure (pres).
   This is standard FAO-56 practice; the difference is < 1%.

5. **Interpolation limits.** T, RH, u10: linear interpolation up to 3-h gaps.
   Rs: up to 2-h gaps during daylight only. Gaps exceeding those limits remain
   NaN and are reported. For this dataset no such gaps existed.

---

## Computation

Penman is applied to each of the 1488 hourly records independently. Daily
evaporation is the sum of 24 hourly values, and cumulative evaporation is the
running sum of daily values.

---

## Headline Results

| | Mean daily | Total | Volume |
|--|-----------|-------|--------|
| **T1** | 6.660 mm/day | 239.77 mm | 38 578 m3 |
| **T2** | 6.135 mm/day | 380.35 mm | 61 196 m3 |

Condensation (negative hourly E):
- T1: 292 hours, -18.075 mm total
- T2: 532 hours, -33.006 mm total
These occur almost exclusively at night (Rs=0, negative Rn) and are included
in all totals. A clipped cumulative curve (zeroing condensation hours) is
provided separately in the CSVs.

---

## Files

| File | Description |
|------|-------------|
| `config.py` | All constants (area, station, physics) |
| `data_io/download.py` | Fetch + cache raw data; data.gouv API + Infoclimat fallback |
| `data_io/parse.py` | Parse CSV.gz and HTML -> tidy hourly DataFrame |
| `physics/radiation.py` | Ra, Rso, Rns, Rnl, Rn (pure functions) |
| `physics/penman.py` | Penman combination, wind function, unit conversions |
| `analysis.py` | Full pipeline: hourly evaporation, aggregation, plots |
| `tests/test_penman.py` | Hand-computed unit test (tolerance 3%) |
| `tests/test_radiation.py` | Ra daily/hourly, Rn sign tests |
| `results/hourly_evaporation.csv` | 1488 hourly rows with all intermediate columns |
| `results/daily_T1.csv` | 36 daily rows of evaporation totals |
| `results/daily_T2.csv` | 62 daily rows of evaporation totals |
| `results/cumulative_T1.csv` | Cumulative mm and m3 |
| `results/cumulative_T2.csv` | Cumulative mm and m3 |
| `results/summary.json` | Headline numbers, gap report, condensation stats |
| `results/fig_hourly_timeseries.png` | Hourly E time series with condensation marked |
| `results/fig_daily_T1.png` / `fig_daily_T2.png` | Daily evaporation (mm/day) |
| `results/fig_cumulative_T1/T2.png` | Cumulative curves (mm left, m3 right) |

---

## References

- Allen, R.G., Pereira, L.S., Raes, D., Smith, M. (1998). FAO Irrigation and
  Drainage Paper 56: Crop evapotranspiration. FAO, Rome.
- Penman, H.L. (1948). Natural evaporation from open water, bare soil and
  grass. Proc. R. Soc. A, 193(1032), 120-145.
- Valiantzas, J.D. (2006). Simplified versions for the Penman evaporation
  equation using routine weather data. J. Hydrology, 331, 690-702.
