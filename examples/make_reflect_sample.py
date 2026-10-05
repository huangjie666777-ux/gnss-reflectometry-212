"""Generate a synthetic GNSS-IR reflectometry sample with a known water level.

Station: lat 30 N, lon 114 E, h 50 m. Antenna phase-center height above the
gauge zero Z = 12.0 m; true reflector height H = 8.0 m, so the known water
level is Z - H = 4.0 m.

Satellite geometry is built directly from azimuth/elevation sweeps at a fixed
slant range, so elevations are exact by construction. S1C follows the static
planar-reflector model |1 + 0.4*exp(i*4*pi*H*sin(elev)/lam1)| on top of a
smooth quadratic background in sin(elevation), plus 0.15 dB noise.

  G01 rising arc, G02 setting arc, G03 rising with an S1C gap (arc break),
  G04 rising-then-setting (reversal split), G05 S1C always blank.

Writes examples/obs_reflect.rnx (S1C only) and examples/eph_reflect.sp3.
"""

import datetime as dt
import math
import os

import numpy as np

GPS_EPOCH = dt.datetime(1980, 1, 6, tzinfo=dt.timezone.utc)
C = 299792458.0
F1 = 1575.42e6
LAM1 = C / F1

LAT, LON, HGT = math.radians(30.0), math.radians(114.0), 50.0
A = 6378137.0
F = 1 / 298.257223563
E2 = F * (2 - F)
N = A / math.sqrt(1 - E2 * math.sin(LAT) ** 2)
RX = np.array([(N + HGT) * math.cos(LAT) * math.cos(LON),
               (N + HGT) * math.cos(LAT) * math.sin(LON),
               (N * (1 - E2) + HGT) * math.sin(LAT)])

Z_ANTENNA_M = 12.0       # antenna phase center above gauge zero
H_TRUE_M = 8.0           # reflector height -> water level 4.0 m
WATER_LEVEL_M = Z_ANTENNA_M - H_TRUE_M
RANGE_M = 20200e3
REFL_COEF = 0.4

N_EPOCHS = 90
STEP_S = 20
T0 = dt.datetime(2024, 6, 1, 1, 0, 0, tzinfo=dt.timezone.utc)


def sat_ecef(elev_deg: float, az_deg: float) -> np.ndarray:
    e = math.radians(elev_deg)
    a = math.radians(az_deg)
    enu = np.array([math.cos(e) * math.sin(a),
                    math.cos(e) * math.cos(a),
                    math.sin(e)])
    slat, clat = math.sin(LAT), math.cos(LAT)
    slon, clon = math.sin(LON), math.cos(LON)
    rot = np.array([[-slon, -slat * clon, clat * clon],
                    [clon, -slat * slon, clat * slon],
                    [0.0, clat, slat]])
    return RX + RANGE_M * (rot @ enu)


def elevation_deg(sat_idx: int, epoch: float) -> float:
    """Elevation profile per satellite (degrees)."""
    f = epoch / (N_EPOCHS - 1)
    if sat_idx in (0, 2):
        return 4.0 + 23.0 * f            # G01/G03 rising 4 -> 27
    if sat_idx == 1:
        return 27.0 - 23.0 * f           # G02 setting 27 -> 4
    if sat_idx == 3:
        # G04: smooth rise-and-set 8 -> 14 -> 8 deg (reversal at mid-window)
        return 8.0 + 6.0 * math.sin(math.pi * epoch / (N_EPOCHS - 1))
    return 15.0                          # G05 (S1C always blank anyway)


AZIMUTHS = [70.0, 250.0, 120.0, 300.0, 10.0]


def s1c_dbhz(elev_deg: float, rng) -> float:
    x = math.sin(math.radians(elev_deg))
    base = 25.0 + 30.0 * (x - 0.1) - 40.0 * (x - 0.1) ** 2
    phase = 4.0 * math.pi * H_TRUE_M * x / LAM1
    amp = base + 6.0 * math.cos(phase)  # smooth background + multipath oscillation
    return 20.0 * math.log10(amp) + rng.normal(0.0, 0.15)


def main():
    outdir = os.path.dirname(os.path.abspath(__file__))
    rng = np.random.default_rng(7)

    # --- SP3-c: 300 s spacing, generous margin for 8-node windows ---
    sp3_start = T0 - dt.timedelta(seconds=1200)
    sp3_n = 15
    lines = [
        f"#cP{sp3_start.year:4d} {sp3_start.month:2d} {sp3_start.day:2d} "
        f"{sp3_start.hour:2d} {sp3_start.minute:2d} {sp3_start.second:11.8f} "
        f"{sp3_n:7d}    u+U  IGS14 FIT  AIUB",
        "## 0000      0.00000000    300.00000000   00000   0.0000000000000",
        "+    5   G01G02G03G04G05  0  0  0  0  0  0  0  0  0  0  0",
        "++          0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0",
        "%c M  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc",
        "%f  0.0000000  0.000000000  0.00000000000  0.000000000000000",
        "%i    0    0    0    0         0         0         0         0",
        "/* SYNTHETIC REFLECTOMETRY EPHEMERIS - NOT REAL NAVIGATION DATA",
    ]
    for k in range(sp3_n):
        t = sp3_start + dt.timedelta(seconds=300 * k)
        epoch_frac = (t - T0).total_seconds() / STEP_S
        lines.append(f"*  {t.year:4d} {t.month:2d} {t.day:2d} {t.hour:2d} "
                     f"{t.minute:2d} {t.second:11.8f}")
        for s in range(5):
            p = sat_ecef(elevation_deg(s, epoch_frac), AZIMUTHS[s]) / 1000.0
            lines.append(f"PG{s + 1:02d}{p[0]:14.6f}{p[1]:14.6f}{p[2]:14.6f}"
                         f"{0.0:14.6f}")
    lines.append("EOF")
    with open(os.path.join(outdir, "eph_reflect.sp3"), "w") as fh:
        fh.write("\n".join(lines) + "\n")

    # --- RINEX 3.04 obs, S1C only ---
    lines = [
        f"{'3.04':<9}{'':11}{'O':<20}{'G':<20}RINEX VERSION / TYPE",
        f"{'SYNTH-REFL':<20}{'TEST':<20}{'20240601 000000 UTC':<20}PGM / RUN BY / DATE",
        f"{RX[0]:14.4f}{RX[1]:14.4f}{RX[2]:14.4f}{'':18}APPROX POSITION XYZ",
        f"{0:6d}{'':54}RCV CLOCK OFFS APPL",
        f"{'G':<1}{1:5d} {'S1C':<4}{'':49}SYS / # / OBS TYPES",
        f"{T0.year:6d}{T0.month:6d}{T0.day:6d}{T0.hour:6d}{T0.minute:6d}"
        f"{T0.second:12.7f}GPS{'':9}TIME OF FIRST OBS",
        f"{'':60}END OF HEADER",
    ]
    for e in range(N_EPOCHS):
        t = T0 + dt.timedelta(seconds=STEP_S * e)
        lines.append(f"> {t.year:4d} {t.month:02d} {t.day:02d} {t.hour:02d} "
                     f"{t.minute:02d}{t.second:11.7f}  0{5:3d}{'':6}")
        for s in range(5):
            prn = f"G{s + 1:02d}"
            blank = " " * 16
            if s == 4:
                lines.append(prn + blank)          # G05: all-empty record
                continue
            if s == 2 and 30 <= e <= 34:
                lines.append(prn + blank)          # G03: S1C gap -> arc break
                continue
            lines.append(prn + f"{s1c_dbhz(elevation_deg(s, e), rng):14.3f}  ")
    with open(os.path.join(outdir, "obs_reflect.rnx"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"wrote examples/obs_reflect.rnx, examples/eph_reflect.sp3 "
          f"(Z={Z_ANTENNA_M} m, H={H_TRUE_M} m, water level={WATER_LEVEL_M} m)")


if __name__ == "__main__":
    main()
