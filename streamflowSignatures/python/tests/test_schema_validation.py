"""
Tests for the io.py schema model and validate_schema().

Pins the August 2026 product arithmetic:

    1,653 = 20 metadata + 100 bases x 16 suffixes + 21 scalars + 12 QA flags

and exercises both validation modes: strict (complete production CSV — every
expected column required, nothing unexpected) and non-strict (partial outputs
such as a bare calculate_all_signatures dict must not fail).
"""

import numpy as np
import pandas as pd
import pytest

from streamflow_signatures.io import (
    ALL_STAT_SUFFIXES,
    EXPECTED_SIGNATURE_BASES,
    METADATA_COLUMNS,
    PER_GAGE_SCALARS,
    PETTITT_SUFFIXES,
    QA_FLAG_COLUMNS,
    STAT_SUFFIXES,
    expected_output_columns,
    validate_schema,
)


def _full_header_df():
    """A one-row DataFrame with exactly the full expected column set."""
    cols = expected_output_columns()
    return pd.DataFrame([[np.nan] * len(cols)], columns=cols)


class TestSchemaArithmetic:
    def test_component_counts(self):
        assert len(METADATA_COLUMNS) == 20
        assert len(EXPECTED_SIGNATURE_BASES) == 100
        assert len(STAT_SUFFIXES) == 8
        assert len(PETTITT_SUFFIXES) == 8
        assert len(ALL_STAT_SUFFIXES) == 16
        assert len(PER_GAGE_SCALARS) == 21
        assert len(QA_FLAG_COLUMNS) == 12

    def test_no_duplicates_within_components(self):
        assert len(set(METADATA_COLUMNS)) == 20
        assert len(set(EXPECTED_SIGNATURE_BASES)) == 100
        assert len(set(ALL_STAT_SUFFIXES)) == 16
        assert len(set(PER_GAGE_SCALARS)) == 21
        assert len(set(QA_FLAG_COLUMNS)) == 12

    def test_total_is_1653(self):
        expected = expected_output_columns()
        assert len(expected) == 1653
        assert len(set(expected)) == 1653  # no collisions across components
        assert 20 + 100 * 16 + 21 + 12 == 1653

    def test_expected_contains_known_columns(self):
        expected = set(expected_output_columns())
        # Metadata (production names, not the stale pre-Aug-2026 ones)
        assert "basin_area" in expected
        assert "area_normalized" in expected
        assert "start_water_year" in expected
        assert "basin_area_km2" not in expected
        assert "processing_status" not in expected
        # Ordinary stat + Pettitt field of the same base
        assert "Qann_mean" in expected
        assert "Qann_pettitt_cp_year" in expected
        # Scalar, flag, drought/snow families
        assert "recession_alpha_point_cloud_linear_reservoir" in expected
        assert "drought_threshold_fixed_p2" in expected
        assert "flagged_for_high_na" in expected
        assert "swe_max_median" in expected
        assert "drought_deficit_fixed_p30_pettitt_post_mk_pval" in expected


class TestValidateSchemaStrict:
    def test_full_header_passes(self):
        result = validate_schema(_full_header_df(), strict=True)
        assert result["valid"] is True
        assert result["missing_metadata"] == []
        assert result["missing_signatures"] == []
        assert result["missing_columns"] == []
        assert result["extra_columns"] == []
        assert result["n_expected"] == 1653

    @pytest.mark.parametrize(
        "dropped",
        ["gage_id", "Qann_pettitt_pval", "elasticity_static", "flagged_for_high_na"],
    )
    def test_dropping_one_column_fails(self, dropped):
        df = _full_header_df().drop(columns=[dropped])
        result = validate_schema(df, strict=True)
        assert result["valid"] is False
        assert result["missing_columns"] == [dropped]
        assert result["extra_columns"] == []

    def test_adding_one_bogus_column_fails(self):
        df = _full_header_df()
        df["bogus_column"] = 1.0
        result = validate_schema(df, strict=True)
        assert result["valid"] is False
        assert result["missing_columns"] == []
        assert result["extra_columns"] == ["bogus_column"]


class TestValidateSchemaNonStrict:
    def test_partial_output_passes(self):
        """A bare signatures dict (8 stats only, no metadata / Pettitt /
        scalars / flags) must validate in the default non-strict mode."""
        cols = [f"{base}{suffix}"
                for base in EXPECTED_SIGNATURE_BASES for suffix in STAT_SUFFIXES]
        df = pd.DataFrame([[np.nan] * len(cols)], columns=cols)
        result = validate_schema(df)
        assert result["valid"] is True
        assert result["missing_signatures"] == []
        # The gaps are still reported, just not fatal
        assert result["missing_metadata"] == METADATA_COLUMNS
        assert len(result["missing_columns"]) == 1653 - len(cols)
        assert result["extra_columns"] == []

    def test_entirely_absent_base_fails_non_strict(self):
        cols = [f"{base}{suffix}"
                for base in EXPECTED_SIGNATURE_BASES for suffix in STAT_SUFFIXES
                if base != "Qann"]
        df = pd.DataFrame([[np.nan] * len(cols)], columns=cols)
        result = validate_schema(df)
        assert result["valid"] is False
        assert result["missing_signatures"] == ["Qann"]

    def test_full_header_passes_non_strict_too(self):
        assert validate_schema(_full_header_df())["valid"] is True
