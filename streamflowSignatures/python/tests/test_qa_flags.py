"""flagged_for_high_na denominator — regression for the 2026-09-04 finding.

Mirrors julia/test/test_qa_high_na.jl. The denominator is defined by NAME (the 16
statistic suffixes + the per-gage scalar registry, config qa_qc.high_na_denominator);
metadata, unregistered diagnostics and flag columns never enter it. Before the fix
Python counted every non-flag column (incl. string metadata), Julia's runner path
counted metadata only, and rpkg counted mean/median keys only.
"""
import numpy as np
import pandas as pd

from streamflow_signatures.qa_qc import (
    compute_qa_flags,
    get_flag_columns,
    high_na_denominator_columns,
)
from streamflow_signatures.config import QAQC_HIGH_NA_SUFFIXES, QAQC_HIGH_NA_SCALARS
from streamflow_signatures.io import (
    ALL_STAT_SUFFIXES,
    PER_GAGE_SCALARS,
    expected_output_columns,
)


def test_manifest_matches_schema_registries():
    assert list(QAQC_HIGH_NA_SUFFIXES) == list(ALL_STAT_SUFFIXES)
    assert list(QAQC_HIGH_NA_SCALARS) == list(PER_GAGE_SCALARS)
    assert len(QAQC_HIGH_NA_SUFFIXES) == 16 and len(QAQC_HIGH_NA_SCALARS) == 21


def test_full_product_header_selects_1621_columns():
    header = expected_output_columns()
    assert len(header) == 1653
    assert len(high_na_denominator_columns(header)) == 1621
    assert not any(c.startswith("flagged_") for c in high_na_denominator_columns(header))


def _frame():
    return pd.DataFrame({
        "gage_id": ["g1", "g2", "g3", "g4"],
        "Qann_mean": [1.0, 1.0, 1.0, np.nan],
        "Qann_median": [1.0, np.nan, 1.0, np.nan],
        "Q5_pettitt_cp_year": [2000.0, np.nan, np.nan, np.nan],
        "elasticity_static": [1.2, np.nan, np.nan, np.nan],   # registered scalar
        "swe_max_mean": [np.nan, np.nan, 3.0, np.nan],        # NA-filled family counts
        # metadata: excluded from the denominator whatever their NA state
        "latitude": [40.0, 41.0, 42.0, 43.0],
        "NDAMS_2009": [np.nan, np.nan, np.nan, 3.0],
        "MAJ_DDENS_2009": [np.nan, np.nan, np.nan, 1.0],
        "num_water_years": [20, 21, 22, 23],
        "area_normalized": [True, True, True, False],
        "gage_type": ["Canada", "Canada", "USGS", "USGS"],
        "CLASS": [None, None, "Ref", None],
    })


def test_denominator_is_signature_columns_only():
    df = _frame()
    assert set(high_na_denominator_columns(df.columns)) == {
        "Qann_mean", "Qann_median", "Q5_pettitt_cp_year", "elasticity_static", "swe_max_mean"}
    out = compute_qa_flags(df)
    # NA fractions over the 5 signature columns: g1 1/5, g2 4/5, g3 3/5, g4 5/5
    assert out["flagged_for_high_na"].tolist() == [False, True, True, True]
    # a second pass must not count the flag columns themselves
    again = compute_qa_flags(out.drop(columns=get_flag_columns()))
    assert again["flagged_for_high_na"].tolist() == out["flagged_for_high_na"].tolist()


def test_metadata_only_frame_never_flags():
    df = pd.DataFrame({"gage_id": ["a", "b"], "latitude": [1.0, 2.0],
                       "NDAMS_2009": [np.nan, np.nan], "CLASS": [None, None]})
    assert high_na_denominator_columns(df.columns) == []
    assert compute_qa_flags(df)["flagged_for_high_na"].tolist() == [False, False]
