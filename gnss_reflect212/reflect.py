"""Single-frequency GNSS-IR reflectometry from S1C signal-to-noise ratios.

Geometry model: the water surface is a static horizontal plane; atmospheric
refraction is neglected. The SNR series of one rising/setting arc oscillates
in sin(elevation) with phase 4*pi*H*sin(elev)/lam1, where H is the vertical
reflector height (antenna phase center above the water surface). Water level
relative to the gauge zero is Z - H, with Z the antenna phase-center height
above the gauge zero.

Because sin(elevation) is not sampled uniformly, no FFT is used: heights are
scanned on a regular grid and, for each trial H, cos/sin of the phase plus a
constant are fitted by least squares; the height minimizing the residual sum
of squares wins (ties resolve to the smaller H).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .geodesy import C_LIGHT, elevation_azimuth, epoch_to_seconds
from .interp import satellite_position
from .rinex import RinexData
from .sp3 import Sp3Data

F1 = 1575.42e6  # GPS L1 frequency, Hz
LAM1 = C_LIGHT / F1

OBS_TYPE = "S1C"
MIN_ELEV_DEG = 5.0
MAX_ELEV_DEG = 25.0
MAX_GAP_S = 120.0
MIN_ARC_POINTS = 12
MIN_ARC_SPAN_DEG = 5.0
MAX_GRID_POINTS = 5001


@dataclass
class _Sample:
    time_iso: str
    t_sec: float
    s1c: float | None = None
    elev_deg: float | None = None
    problem: str | None = None  # filter/break reason


def _collect_samples(obs: RinexData, sp3: Sp3Data, prn: str,
                     station_ecef: np.ndarray,
                     lat_deg: float, lon_deg: float) -> list[_Sample]:
    samples: list[_Sample] = []
    rec = sp3.sats.get(prn)
    for ep in obs.epochs:
        t_sec = epoch_to_seconds(ep.time)
        iso = ep.time.isoformat()
        entry = ep.obs.get(prn)
        if entry is None:
            samples.append(_Sample(iso, t_sec,
                                   problem="satellite not observed at epoch"))
            continue
        if OBS_TYPE not in entry:
            samples.append(_Sample(iso, t_sec,
                                   problem="missing S1C signal-to-noise observation"))
            continue
        s1c = entry[OBS_TYPE][0]
        if rec is None:
            samples.append(_Sample(iso, t_sec, s1c=s1c,
                                   problem="satellite not present in SP3 file"))
            continue
        pos, reason = satellite_position(rec, t_sec)
        if pos is None:
            samples.append(_Sample(iso, t_sec, s1c=s1c, problem=reason))
            continue
        elev, _az = elevation_azimuth(station_ecef, lat_deg, lon_deg, pos)
        if not (MIN_ELEV_DEG <= elev <= MAX_ELEV_DEG):
            samples.append(_Sample(
                iso, t_sec, s1c=s1c, elev_deg=elev,
                problem=f"elevation {elev:.2f} deg outside "
                        f"[{MIN_ELEV_DEG:.0f}, {MAX_ELEV_DEG:.0f}] deg filter"))
            continue
        samples.append(_Sample(iso, t_sec, s1c=s1c, elev_deg=elev))
    return samples


def _split_arcs(samples: list[_Sample]) -> list[list[_Sample]]:
    """Split valid samples into monotonic rising/setting elevation arcs.

    Breaks on filtered samples (missing S1C, no ephemeris, elevation outside
    the 5-25 deg window), time gaps > MAX_GAP_S, and rising/setting reversals.
    """
    arcs: list[list[_Sample]] = []
    cur: list[_Sample] = []
    direction = 0  # +1 rising, -1 setting, 0 undetermined
    for s in samples:
        if s.problem is not None:
            if cur:
                arcs.append(cur)
                cur = []
            direction = 0
            continue
        if cur:
            if s.t_sec - cur[-1].t_sec > MAX_GAP_S:
                arcs.append(cur)
                cur = []
                direction = 0
            else:
                d = (s.elev_deg > cur[-1].elev_deg) - (s.elev_deg < cur[-1].elev_deg)
                if d != 0:
                    if direction == 0:
                        direction = d
                    elif d != direction:
                        arcs.append(cur)
                        cur = []
                        direction = d
        cur.append(s)
    if cur:
        arcs.append(cur)
    return arcs


def invert_arc(elev_deg: np.ndarray, s1c: np.ndarray,
               grid_m: np.ndarray) -> dict:
    """Estimate reflector height from one arc; returns a result dict."""
    x = np.sin(np.radians(elev_deg))
    amp = np.power(10.0, s1c / 20.0)  # dB -> linear amplitude

    # remove quadratic trend in sin(elevation)
    trend_design = np.column_stack([np.ones_like(x), x, x * x])
    coef, *_ = np.linalg.lstsq(trend_design, amp, rcond=None)
    if np.linalg.matrix_rank(trend_design) < 3:
        return {"status": "failed",
                "reason": "rank-deficient quadratic trend fit "
                          "(sin(elevation) span too small)"}
    y = amp - trend_design @ coef
    if float(np.std(y)) <= 1e-9 * max(1.0, float(np.std(amp))):
        return {"status": "failed",
                "reason": "no residual oscillation after detrending"}

    rss = np.empty(len(grid_m))
    coefs = np.empty((len(grid_m), 3))
    for i, h in enumerate(grid_m):
        phase = 4.0 * np.pi * h * x / LAM1
        design = np.column_stack([np.cos(phase), np.sin(phase), np.ones_like(x)])
        if np.linalg.matrix_rank(design) < 3:
            return {"status": "failed",
                    "reason": f"rank-deficient phase fit at H={h:.4f} m"}
        c, *_ = np.linalg.lstsq(design, y, rcond=None)
        r = y - design @ c
        rss[i] = float(r @ r)
        coefs[i] = c

    best = int(np.argmin(rss))  # first minimum -> smaller H on ties
    c = coefs[best]
    warnings: list[str] = []
    if best == 0 or best == len(grid_m) - 1:
        warnings.append("optimum at reflector-height grid boundary; "
                        "consider widening the search interval")
    return {
        "status": "ok",
        "reflector_height_m": float(grid_m[best]),
        "amplitude": float(math.hypot(c[0], c[1])),
        "residual_rms": float(math.sqrt(rss[best] / len(x))),
        "warnings": warnings,
        "height_scan": {"h_m": [float(v) for v in grid_m],
                        "rss": [float(v) for v in rss]},
    }


def compute_reflectometry(obs: RinexData, sp3: Sp3Data, station_ecef: np.ndarray,
                          lat_deg: float, lon_deg: float, z_m: float,
                          h_min: float, h_max: float, h_step: float) -> dict:
    """Per-satellite arc splitting and reflector-height inversion."""
    n_grid = int(round((h_max - h_min) / h_step)) + 1
    grid_m = h_min + h_step * np.arange(n_grid)

    all_prns = sorted({prn for ep in obs.epochs for prn in ep.obs})
    out_sats: dict[str, dict] = {}
    for prn in all_prns:
        samples = _collect_samples(obs, sp3, prn, station_ecef, lat_deg, lon_deg)
        arcs = _split_arcs(samples)
        arc_of: dict[int, int] = {}  # id(sample) -> arc number
        arc_results: list[dict] = []
        for arc_no, arc in enumerate(arcs, start=1):
            for s in arc:
                arc_of[id(s)] = arc_no
            direction = ("rising" if arc[-1].elev_deg > arc[0].elev_deg
                         else "setting")
            result: dict = {
                "arc_id": f"{prn}-A{arc_no}",
                "start": arc[0].time_iso,
                "end": arc[-1].time_iso,
                "direction": direction,
                "n_samples": len(arc),
                "obs_type": OBS_TYPE,
                "elevation_span_deg": float(arc[-1].elev_deg - arc[0].elev_deg)
                if direction == "rising"
                else float(arc[0].elev_deg - arc[-1].elev_deg),
            }
            span = max(s.elev_deg for s in arc) - min(s.elev_deg for s in arc)
            if len(arc) < MIN_ARC_POINTS:
                result.update(status="failed",
                              reason=f"arc has {len(arc)} samples "
                                     f"(< {MIN_ARC_POINTS})")
            elif span < MIN_ARC_SPAN_DEG:
                result.update(status="failed",
                              reason=f"elevation span {span:.2f} deg "
                                     f"(< {MIN_ARC_SPAN_DEG:.0f} deg)")
            else:
                inv = invert_arc(np.array([s.elev_deg for s in arc]),
                                 np.array([s.s1c for s in arc]), grid_m)
                result.update(inv)
                if inv["status"] == "ok":
                    result["water_level_m"] = z_m - inv["reflector_height_m"]
            arc_results.append(result)

        records = []
        for s in samples:
            record = {
                "time": s.time_iso,
                "status": "ok" if s.problem is None else "excluded",
                "reason": s.problem,
                "arc_id": (f"{prn}-A{arc_of[id(s)]}"
                           if s.problem is None and id(s) in arc_of else None),
                "s1c_dbhz": s.s1c,
                "elevation_deg": s.elev_deg,
            }
            records.append(record)
        out_sats[prn] = {"n_arcs": len(arcs), "records": records,
                         "arcs": arc_results}

    return {
        "antenna_height_z_m": z_m,
        "reflector_grid": {"min_m": h_min, "max_m": h_max,
                           "step_m": h_step, "n_points": n_grid},
        "elevation_filter_deg": [MIN_ELEV_DEG, MAX_ELEV_DEG],
        "satellites": out_sats,
    }
