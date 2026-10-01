"""Shared helpers for the Daymet basin-aggregation tools (EO_data_processing/daymet/).

Plan: docs/plans/2026-09-29-daymet-reprocessing-action-plan.md. The tools compute a
coverage-weighted daily mean per watershed polygon straight from the Daymet V4 R1
North America annual mosaics (one NetCDF-4 file per variable-year), without
resampling the raster: the polygons are reprojected to the Daymet grid instead.

Grid facts are READ FROM EACH FILE and asserted against the constants below. The checks
cover:
  * the shape, the first and last cell centres and the 1 km spacing;
  * the LCC parameters and ellipsoid;
  * the fill value, and the absence of packing;
  * the time units, and the Daymet calendar: 365 strictly daily steps from Jan 1, so
    Dec 31 is absent exactly in leap years.
A file with a different grid, CRS or calendar therefore stops the run instead of
silently misaligning. The CRS, last-cell and calendar asserts were added on 2026-10-01.
Before that they were recorded in the probe JSONs but not enforced; all 276 files of the
2026-09-29 run carry the identical CRS.
"""
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
import warnings
from dataclasses import dataclass

import h5py
import numpy as np
from pyproj import CRS

DAYMET_VARS = ["prcp", "tmin", "tmax", "swe", "vp", "srad"]
FILL = -9999.0
NX, NY = 7814, 8075
CELL = 1000.0
X0_CENTRE, Y0_CENTRE = -4560250.0, 4984000.0   # first column / first (top) row centres
X1_CENTRE, Y1_CENTRE = X0_CENTRE + (NX - 1) * CELL, Y0_CENTRE - (NY - 1) * CELL   # last column / bottom row
# The Daymet V4 R1 NA projection as PROJ parameters (CRS.from_cf of every file of the 2026-09-29
# run). Compared parameter by parameter: CRS.equals() is False against the files' own CRS
# because their datum and ellipsoid are named "undefined".
DAYMET_LCC = {"proj": "lcc", "lat_0": 42.5, "lon_0": -100.0, "lat_1": 25.0, "lat_2": 60.0,
              "x_0": 0.0, "y_0": 0.0, "ellps": "WGS84", "units": "m"}
HERE = os.path.dirname(os.path.abspath(__file__))

# Plausible physical ranges per variable (QA only; a value outside is reported, not altered).
RANGES = {
    "prcp": (0.0, 1000.0),    # mm/day
    "swe": (0.0, 50000.0),    # kg/m2
    "tmin": (-80.0, 50.0),    # degC
    "tmax": (-70.0, 60.0),    # degC
    "vp": (0.0, 10000.0),     # Pa
    "srad": (0.0, 1500.0),    # W/m2 (daylight average)
}


def canon(x):
    """Canonical gage id: digits -> leading zeros stripped; otherwise upper-cased.

    Same rule as the geometry rebuild and the Julia runner's join id. Join ids across
    products on this form, never on re-padded strings (some staged ids are un-padded).
    """
    s = str(x).strip()
    return str(int(s)) if s.isdigit() else s.upper()


def _attr(v):
    if isinstance(v, bytes):
        return v.decode()
    if isinstance(v, np.ndarray):
        return v.tolist() if v.size > 1 else v.item()
    if isinstance(v, np.generic):
        return v.item()
    return v


@dataclass
class DaymetGrid:
    x: np.ndarray            # cell-centre x (m), ascending
    y: np.ndarray            # cell-centre y (m), descending (row 0 = north)
    crs: CRS
    var: str
    shape: tuple             # (nt, ny, nx)
    chunks: tuple            # (ct, cy, cx)
    compression: str
    shuffle: bool
    hdf5_fillvalue: float
    time: np.ndarray         # days since 1950-01-01 (noon values)
    dates: np.ndarray        # datetime64[D]

    @property
    def bounds(self):
        """(xmin, ymin, xmax, ymax) of the cell edges."""
        return (self.x[0] - CELL / 2, self.y[-1] - CELL / 2,
                self.x[-1] + CELL / 2, self.y[0] + CELL / 2)


def read_grid(path, var=None):
    """Open a Daymet NA file and return its grid, asserting the expected geometry."""
    with h5py.File(path, "r") as f:
        if var is None:
            cands = [v for v in DAYMET_VARS + ["dayl"] if v in f]
            if len(cands) != 1:
                raise ValueError(f"{path}: cannot infer the variable (found {cands})")
            var = cands[0]
        d = f[var]
        x = f["x"][:].astype(np.float64)
        y = f["y"][:].astype(np.float64)
        t = f["time"][:].astype(np.float64)
        tunits = _attr(f["time"].attrs.get("units", b""))
        cf = {k: _attr(v) for k, v in f["lambert_conformal_conic"].attrs.items()}
        for bad in ("scale_factor", "add_offset"):
            if bad in d.attrs:
                raise ValueError(f"{path}: {var} carries {bad}; packed data is not handled")
        fv = _attr(d.attrs.get("_FillValue", np.nan))
        if float(fv) != FILL:
            raise ValueError(f"{path}: _FillValue {fv} != {FILL}")
        grid = DaymetGrid(
            x=x, y=y, crs=CRS.from_cf(cf), var=var, shape=d.shape, chunks=d.chunks,
            compression=str(d.compression), shuffle=bool(d.shuffle),
            hdf5_fillvalue=float(d.fillvalue), time=t,
            dates=(np.datetime64("1950-01-01") + np.floor(t).astype("timedelta64[D]")),
        )
    if not tunits.startswith("days since 1950-01-01"):
        raise ValueError(f"{path}: unexpected time units {tunits!r}")
    bad = crs_mismatch(grid.crs)
    if bad:
        raise ValueError(f"{path}: CRS differs from the Daymet NA LCC: {bad}")
    if grid.shape[1:] != (NY, NX):
        raise ValueError(f"{path}: grid {grid.shape[1:]} != {(NY, NX)}")
    if not (np.allclose(np.diff(x), CELL, rtol=0, atol=1e-3) and np.allclose(np.diff(y), -CELL, rtol=0, atol=1e-3)):
        raise ValueError(f"{path}: x/y spacing is not +/-{CELL} m")
    if (x[0], y[0]) != (X0_CENTRE, Y0_CENTRE):
        raise ValueError(f"{path}: grid origin {(x[0], y[0])} != {(X0_CENTRE, Y0_CENTRE)}")
    if (x[-1], y[-1]) != (X1_CENTRE, Y1_CENTRE):
        raise ValueError(f"{path}: last cell centre {(x[-1], y[-1])} != {(X1_CENTRE, Y1_CENTRE)}")
    if grid.shape[0] != 365:
        raise ValueError(f"{path}: {grid.shape[0]} time steps (Daymet years have 365)")
    first = str(grid.dates[0])
    if first[5:] != "01-01" or not (np.diff(grid.dates.astype("int64")) == 1).all():
        # with 365 steps this is the Daymet calendar: Dec 31 absent exactly in leap years
        raise ValueError(f"{path}: time axis is not 365 consecutive days from Jan 1 (starts {first}, "
                         f"ends {grid.dates[-1]})")
    return grid


def crs_mismatch(crs):
    """The PROJ parameters of `crs` that differ from DAYMET_LCC (empty dict = same projection)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")      # to_dict() warns that a PROJ string is lossy
        got = crs.to_dict()
    out = {}
    for k, v in DAYMET_LCC.items():
        g = got.get(k)
        same = (isinstance(v, float) and isinstance(g, (int, float)) and abs(g - v) < 1e-9) or g == v
        if not same:
            out[k] = (g, v)
    extra = set(got) - set(DAYMET_LCC) - {"no_defs", "type"}
    for k in sorted(extra):
        out[k] = (got[k], None)
    return out


def cell_true_area_km2(cell_id, crs, block=4_000_000):
    """True (ellipsoidal) area in km2 of the 1 km Daymet cells with these flat indices.

    The grid is conformal (LCC), so a projected 1 km2 cell covers 1 / areal_scale km2 on
    the ellipsoid (areal_scale = k^2, 0.91-1.21 over the grid's populated area). Weighting
    by coverage x true area gives the area-weighted basin mean an equal-area weight
    generator (gdptools, the co-authors' tool) produces; weighting by coverage alone
    over-weights cells where k < 1. Evaluated at the cell centre (k varies by < 1e-4
    across one cell).
    """
    from pyproj import Proj, Transformer
    cell_id = np.asarray(cell_id, dtype=np.int64)
    out = np.empty(cell_id.shape, dtype=np.float64)
    to_ll = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
    proj = Proj(crs)
    for s in range(0, cell_id.size, block):
        c = cell_id[s:s + block]
        x = X0_CENTRE + (c % NX) * CELL
        y = Y0_CENTRE - (c // NX) * CELL
        lon, lat = to_ll.transform(x, y)
        out[s:s + block] = 1.0 / np.asarray(proj.get_factors(lon, lat).areal_scale)
    return out


def md5sum(path, block=1 << 24):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(block), b""):
            h.update(b)
    return h.hexdigest()


def sha256sum(path, block=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(block), b""):
            h.update(b)
    return h.hexdigest()


def software_versions():
    import exactextract, geopandas, pandas, pyarrow, pyogrio, pyproj, scipy, shapely
    return {
        "python": sys.version.split()[0], "platform": platform.platform(),
        "machine": platform.machine(), "numpy": np.__version__, "h5py": h5py.__version__,
        "hdf5": h5py.version.hdf5_version, "scipy": scipy.__version__,
        "pandas": pandas.__version__, "pyarrow": pyarrow.__version__,
        "geopandas": geopandas.__version__, "shapely": shapely.__version__,
        "geos": shapely.geos_version_string, "pyogrio": pyogrio.__version__,
        "gdal": getattr(pyogrio, "__gdal_version_string__", "unknown"),
        "pyproj": pyproj.__version__, "proj": pyproj.proj_version_str,
        "exactextract": getattr(exactextract, "__version__", "unknown"),
    }


def git_state(path=HERE):
    """Commit of the checkout holding `path`, and whether `path` has uncommitted changes."""
    def run(*args):
        r = subprocess.run(["git", "-C", path] + list(args), capture_output=True, text=True)
        return r.stdout.strip() if r.returncode == 0 else None
    return {"commit": run("rev-parse", "HEAD"),
            "daymet_tools_dirty": bool(run("status", "--porcelain", "--", path))}


def utc_now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def write_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=2, default=str)
    os.replace(tmp, path)
