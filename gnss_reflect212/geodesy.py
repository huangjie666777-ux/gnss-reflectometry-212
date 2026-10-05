"""WGS84 coordinate conversions and GPS time helpers."""

from __future__ import annotations

import datetime as dt
import math

import numpy as np

WGS84_A = 6378137.0
WGS84_F = 1.0 / 298.257223563
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)
OMEGA_E = 7.2921151467e-5  # Earth rotation rate, rad/s
C_LIGHT = 299792458.0  # m/s

GPS_EPOCH = dt.datetime(1980, 1, 6, tzinfo=dt.timezone.utc)


def ecef_to_geodetic(xyz: np.ndarray) -> tuple[float, float, float]:
    """ECEF (m) -> (lat deg, lon deg, ellipsoidal height m), WGS84."""
    x, y, z = float(xyz[0]), float(xyz[1]), float(xyz[2])
    lon = math.atan2(y, x)
    p = math.hypot(x, y)
    lat = math.atan2(z, p * (1.0 - WGS84_E2))
    for _ in range(50):
        n = WGS84_A / math.sqrt(1.0 - WGS84_E2 * math.sin(lat) ** 2)
        h = p / math.cos(lat) - n
        new_lat = math.atan2(z, p * (1.0 - WGS84_E2 * n / (n + h)))
        if abs(new_lat - lat) < 1e-14:
            lat = new_lat
            break
        lat = new_lat
    n = WGS84_A / math.sqrt(1.0 - WGS84_E2 * math.sin(lat) ** 2)
    h = p / math.cos(lat) - n
    return math.degrees(lat), math.degrees(lon), h


def geodetic_to_ecef(lat_deg: float, lon_deg: float, h: float) -> np.ndarray:
    """(lat deg, lon deg, ellipsoidal height m) -> ECEF (m), WGS84."""
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    n = WGS84_A / math.sqrt(1.0 - WGS84_E2 * math.sin(lat) ** 2)
    return np.array([
        (n + h) * math.cos(lat) * math.cos(lon),
        (n + h) * math.cos(lat) * math.sin(lon),
        (n * (1.0 - WGS84_E2) + h) * math.sin(lat),
    ])


def epoch_to_seconds(t: dt.datetime) -> float:
    """Seconds since GPS epoch (continuous, week-rollover safe)."""
    return (t - GPS_EPOCH).total_seconds()
