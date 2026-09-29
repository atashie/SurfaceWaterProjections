"""Shared helpers for the Daymet basin-aggregation tools (EO_data_processing/daymet/).

Plan: docs/plans/2026-09-29-daymet-reprocessing-action-plan.md. The tools compute a
coverage-weighted daily mean per watershed polygon straight from the Daymet V4 R1
North America annual mosaics (one NetCDF-4 file per variable-year), without
resampling the raster: the polygons are reprojected to the Daymet grid instead.

Grid facts are READ FROM EACH FILE and asserted against the constants below, so a
file with a different grid or CRS stops the run instead of silently misaligning.
"""
import hashlib
import json
import os
import platform
import sys
from dataclasses import dataclass

import h5py
import numpy as np
from pyproj import CRS

DAYMET_VARS = ["prcp", "tmin", "tmax", "swe", "vp", "srad"]
FILL = -9999.0
NX, NY = 7814, 8075
CELL = 1000.0
X0_CENTRE, Y0_CENTRE = -4560250.0, 4984000.0   # first column / first (top) row centres

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
    if grid.shape[1:] != (NY, NX):
        raise ValueError(f"{path}: grid {grid.shape[1:]} != {(NY, NX)}")
    if not (np.allclose(np.diff(x), CELL) and np.allclose(np.diff(y), -CELL)):
        raise ValueError(f"{path}: x/y spacing is not +/-{CELL} m")
    if (x[0], y[0]) != (X0_CENTRE, Y0_CENTRE):
        raise ValueError(f"{path}: grid origin {(x[0], y[0])} != {(X0_CENTRE, Y0_CENTRE)}")
    if grid.shape[0] != 365:
        raise ValueError(f"{path}: {grid.shape[0]} time steps (Daymet years have 365)")
    return grid


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
    import exactextract, geopandas, pandas, pyarrow, pyproj, scipy, shapely
    return {
        "python": sys.version.split()[0], "platform": platform.platform(),
        "machine": platform.machine(), "numpy": np.__version__, "h5py": h5py.__version__,
        "hdf5": h5py.version.hdf5_version, "scipy": scipy.__version__,
        "pandas": pandas.__version__, "pyarrow": pyarrow.__version__,
        "geopandas": geopandas.__version__, "shapely": shapely.__version__,
        "pyproj": pyproj.__version__, "proj": pyproj.proj_version_str,
        "exactextract": getattr(exactextract, "__version__", "unknown"),
    }


def write_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=2, default=str)
    os.replace(tmp, path)
