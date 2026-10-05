import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from fastapi.testclient import TestClient

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from gnss_reflect212.geodesy import geodetic_to_ecef
from gnss_reflect212.interp import satellite_position
from gnss_reflect212.reflect import compute_reflectometry
from gnss_reflect212.rinex import parse_rinex
from gnss_reflect212.sp3 import parse_sp3
from main import app


@pytest.fixture(scope="module", autouse=True)
def synthetic():
    subprocess.run([sys.executable, str(ROOT / "examples" / "make_reflect_sample.py")],
                   check=True, cwd=ROOT)


@pytest.fixture(scope="module")
def refl_result():
    obs = parse_rinex((ROOT / "examples" / "obs_reflect.rnx").read_text())
    eph = parse_sp3((ROOT / "examples" / "eph_reflect.sp3").read_text())
    station = geodetic_to_ecef(30.0, 114.0, 50.0)
    return compute_reflectometry(obs, eph, station, 30.0, 114.0,
                                 12.0, 1.0, 15.0, 0.01)


def test_s1c_only_header_accepted():
    obs = parse_rinex((ROOT / "examples" / "obs_reflect.rnx").read_text())
    assert obs.obs_types_gps == ["S1C"]
    assert len(obs.epochs) == 90


def test_all_empty_satellite_record_kept():
    obs = parse_rinex((ROOT / "examples" / "obs_reflect.rnx").read_text())
    # G05 has blank S1C in every epoch but must not vanish
    assert all("G05" in ep.obs for ep in obs.epochs)
    assert all(obs_epoch.obs["G05"] == {} for obs_epoch in obs.epochs)


def test_sp3_no_interpolation_across_missing_nodes():
    text = (ROOT / "examples" / "eph_reflect.sp3").read_text()
    lines = text.splitlines()
    # drop all G01 records at one epoch -> a gap in the time series
    out, skipped = [], 0
    for i, ln in enumerate(lines):
        if ln.startswith("PG01") and skipped < 1 and i > 20:
            skipped += 1
            continue
        out.append(ln)
    eph = parse_sp3("\n".join(out) + "\n")
    rec = eph.sats["G01"]
    t_gap = rec.times[1] + 1.0  # window here straddles the removed node
    pos, reason = satellite_position(rec, t_gap)
    assert pos is None and "missing SP3 nodes" in reason
    # away from the gap, interpolation still works
    pos, reason = satellite_position(rec, rec.times[8] + 1.0)
    assert pos is not None


def test_recovers_known_water_level(refl_result):
    for prn in ("G01", "G02"):
        arc = refl_result["satellites"][prn]["arcs"][0]
        assert arc["status"] == "ok"
        assert abs(arc["reflector_height_m"] - 8.0) <= 0.05
        assert abs(arc["water_level_m"] - 4.0) <= 0.05
        assert arc["obs_type"] == "S1C"
        assert arc["n_samples"] >= 12
        assert arc["amplitude"] > 0.0
        assert arc["residual_rms"] > 0.0
        assert len(arc["height_scan"]["h_m"]) == 1401
    assert refl_result["satellites"]["G01"]["arcs"][0]["direction"] == "rising"
    assert refl_result["satellites"]["G02"]["arcs"][0]["direction"] == "setting"


def test_arc_breaks_on_s1c_gap_and_reversal(refl_result):
    g03 = refl_result["satellites"]["G03"]
    assert g03["n_arcs"] == 2  # S1C gap at epochs 30-34
    assert all(a["status"] == "ok" for a in g03["arcs"])
    gap_records = [r for r in g03["records"] if r["status"] == "excluded"
                   and "S1C" in (r["reason"] or "")]
    assert len(gap_records) == 5
    g04 = refl_result["satellites"]["G04"]
    assert g04["n_arcs"] == 2  # rise/set reversal splits the arc
    assert [a["direction"] for a in g04["arcs"]] == ["rising", "setting"]
    assert all(abs(a["reflector_height_m"] - 8.0) <= 0.05
               for a in g04["arcs"] if a["status"] == "ok")


def test_all_empty_satellite_reported(refl_result):
    g05 = refl_result["satellites"]["G05"]
    assert g05["n_arcs"] == 0
    assert len(g05["records"]) == 90
    assert all(r["status"] == "excluded" and "S1C" in r["reason"]
               for r in g05["records"])


def _post_reflect(client, **overrides):
    f1 = open(ROOT / "examples" / "obs_reflect.rnx", "rb")
    f2 = open(ROOT / "examples" / "eph_reflect.sp3", "rb")
    data = {"station_lat_deg": "30.0", "station_lon_deg": "114.0",
            "station_height_m": "50.0", "antenna_height_z_m": "12.0",
            "reflector_min_m": "1.0", "reflector_max_m": "15.0",
            "reflector_step_m": "0.01"}
    data.update(overrides)
    r = client.post("/reflect",
                    files={"rinex": ("obs.rnx", f1), "sp3": ("eph.sp3", f2)},
                    data=data)
    f1.close()
    f2.close()
    return r


def test_reflect_endpoint():
    client = TestClient(app)
    r = _post_reflect(client)
    assert r.status_code == 200
    body = r.json()
    arc = body["satellites"]["G01"]["arcs"][0]
    assert arc["status"] == "ok"
    assert abs(arc["water_level_m"] - 4.0) <= 0.05
    assert body["reflector_grid"]["n_points"] == 1401


def test_reflect_rejects_bad_grid():
    client = TestClient(app)
    assert _post_reflect(client, reflector_min_m="-1.0").status_code == 422
    assert _post_reflect(client, reflector_max_m="0.5").status_code == 422
    assert _post_reflect(client, reflector_step_m="0.0").status_code == 422
    # 1..15 m at 0.002 m -> 7001 points (> 5001) rejected
    r = _post_reflect(client, reflector_step_m="0.002")
    assert r.status_code == 422 and "5001" in r.json()["detail"]


def test_position_rejects_s1c_only_file():
    client = TestClient(app)
    f1 = open(ROOT / "examples" / "obs_reflect.rnx", "rb")
    f2 = open(ROOT / "examples" / "eph_reflect.sp3", "rb")
    r = client.post("/position",
                    files={"rinex": ("obs.rnx", f1), "sp3": ("eph.sp3", f2)})
    f1.close()
    f2.close()
    assert r.status_code == 422
    assert "C1C" in r.json()["detail"]
