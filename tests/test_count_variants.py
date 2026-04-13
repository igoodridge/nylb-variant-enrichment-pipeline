"""Tests for scripts/count_variants.py — count_reads, calculate_frequencies,
calculate_enrichment, and write_report."""

import math
from unittest.mock import MagicMock, patch

import pytest

from count_variants import (
    calculate_enrichment,
    calculate_frequencies,
    count_reads,
    write_report,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _mock_read(is_unmapped: bool, reference_name: str | None) -> MagicMock:
    read = MagicMock()
    read.is_unmapped = is_unmapped
    read.reference_name = reference_name
    return read


def _mock_bam(reads: list) -> MagicMock:
    bam = MagicMock()
    bam.__iter__ = MagicMock(return_value=iter(reads))
    bam.__enter__ = MagicMock(return_value=bam)
    bam.__exit__ = MagicMock(return_value=False)
    return bam


# ---------------------------------------------------------------------------
# count_reads
# ---------------------------------------------------------------------------


class TestCountReads:
    def test_counts_mapped_reads_per_reference(self):
        reads = [
            _mock_read(False, "var1"),
            _mock_read(False, "var1"),
            _mock_read(False, "var2"),
        ]
        with patch("count_variants.pysam.AlignmentFile", return_value=_mock_bam(reads)):
            assert count_reads("dummy.bam") == {"var1": 2, "var2": 1}

    def test_unmapped_reads_are_excluded(self):
        reads = [
            _mock_read(True, None),
            _mock_read(False, "var1"),
        ]
        with patch("count_variants.pysam.AlignmentFile", return_value=_mock_bam(reads)):
            assert count_reads("dummy.bam") == {"var1": 1}

    def test_empty_bam_returns_empty_dict(self):
        with patch("count_variants.pysam.AlignmentFile", return_value=_mock_bam([])):
            assert count_reads("dummy.bam") == {}

    def test_all_unmapped_returns_empty_dict(self):
        reads = [_mock_read(True, None), _mock_read(True, None)]
        with patch("count_variants.pysam.AlignmentFile", return_value=_mock_bam(reads)):
            assert count_reads("dummy.bam") == {}


# ---------------------------------------------------------------------------
# calculate_frequencies
# ---------------------------------------------------------------------------


class TestCalculateFrequencies:
    def test_frequencies_sum_to_one(self):
        freqs = calculate_frequencies({"a": 3, "b": 7})
        assert math.isclose(sum(freqs.values()), 1.0)

    def test_single_reference_has_frequency_one(self):
        assert calculate_frequencies({"only": 100}) == {"only": 1.0}

    def test_equal_counts_give_equal_frequencies(self):
        freqs = calculate_frequencies({"x": 50, "y": 50})
        assert math.isclose(freqs["x"], 0.5)
        assert math.isclose(freqs["y"], 0.5)

    def test_proportional_values_are_correct(self):
        freqs = calculate_frequencies({"a": 1, "b": 3})
        assert math.isclose(freqs["a"], 0.25)
        assert math.isclose(freqs["b"], 0.75)

    def test_all_keys_are_preserved(self):
        counts = {"a": 1, "b": 2, "c": 3}
        assert set(calculate_frequencies(counts).keys()) == {"a", "b", "c"}


# ---------------------------------------------------------------------------
# calculate_enrichment
# ---------------------------------------------------------------------------


class TestCalculateEnrichment:
    @pytest.fixture()
    def basic_counts(self):
        pre = {"var1": 100, "var2": 200, "var3": 700}
        post = {"var1": 500, "var2": 200, "var3": 300}
        return pre, post

    def test_one_row_per_unique_reference(self, basic_counts):
        pre, post = basic_counts
        assert len(calculate_enrichment(pre, post)) == 3

    def test_rows_contain_required_fields(self, basic_counts):
        pre, post = basic_counts
        required = {"variant", "pre_count", "post_count", "pre_freq", "post_freq", "enrichment_ratio"}
        for row in calculate_enrichment(pre, post):
            assert set(row.keys()) == required

    def test_results_sorted_by_variant_name(self, basic_counts):
        pre, post = basic_counts
        names = [r["variant"] for r in calculate_enrichment(pre, post)]
        assert names == sorted(names)

    def test_enriched_variant_has_ratio_above_one(self):
        # a goes from 10 % → 90 %, so should be enriched
        pre = {"a": 1, "b": 9}
        post = {"a": 9, "b": 1}
        rows = {r["variant"]: r for r in calculate_enrichment(pre, post)}
        assert rows["a"]["enrichment_ratio"] > 1
        assert rows["b"]["enrichment_ratio"] < 1

    def test_absent_from_pre_gives_infinite_enrichment(self):
        pre = {"a": 100}
        post = {"a": 50, "new": 50}
        rows = {r["variant"]: r for r in calculate_enrichment(pre, post)}
        assert rows["new"]["enrichment_ratio"] == float("inf")
        assert rows["new"]["pre_count"] == 0

    def test_absent_from_post_gives_zero_enrichment(self):
        pre = {"a": 100, "b": 100}
        post = {"a": 200}
        rows = {r["variant"]: r for r in calculate_enrichment(pre, post)}
        assert rows["b"]["post_count"] == 0
        assert rows["b"]["enrichment_ratio"] == 0.0

    def test_original_counts_are_preserved(self):
        pre = {"x": 10}
        post = {"x": 20}
        row = calculate_enrichment(pre, post)[0]
        assert row["pre_count"] == 10
        assert row["post_count"] == 20


# ---------------------------------------------------------------------------
# write_report
# ---------------------------------------------------------------------------


class TestWriteReport:
    def test_header_row_is_correct(self, tmp_path, sample_enrichment_rows):
        out = tmp_path / "report.tsv"
        write_report(sample_enrichment_rows, str(out))
        header = out.read_text().splitlines()[0]
        assert header == "variant\tpre_count\tpost_count\tpre_freq\tpost_freq\tenrichment_ratio"

    def test_correct_number_of_data_rows(self, tmp_path, sample_enrichment_rows):
        out = tmp_path / "report.tsv"
        write_report(sample_enrichment_rows, str(out))
        lines = out.read_text().splitlines()
        assert len(lines) == len(sample_enrichment_rows) + 1  # +1 for header

    def test_data_values_are_tab_separated(self, tmp_path, sample_enrichment_rows):
        out = tmp_path / "report.tsv"
        write_report(sample_enrichment_rows, str(out))
        first_data = out.read_text().splitlines()[1].split("\t")
        assert first_data[0] == "var1"
        assert first_data[1] == "100"
        assert first_data[4] == "0.2"
        assert first_data[5] == "2.0"

    def test_creates_missing_parent_directories(self, tmp_path, sample_enrichment_rows):
        nested = tmp_path / "a" / "b" / "report.tsv"
        write_report(sample_enrichment_rows, str(nested))
        assert nested.exists()

    def test_empty_rows_writes_header_only(self, tmp_path):
        out = tmp_path / "empty.tsv"
        write_report([], str(out))
        lines = out.read_text().splitlines()
        assert len(lines) == 1
        assert lines[0].startswith("variant")
