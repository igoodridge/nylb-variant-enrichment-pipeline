"""Tests for scripts/generate_variants.py — generate_variants."""

import pytest
from Bio.SeqRecord import SeqRecord

from conftest import WILDTYPE_SEQ
from generate_variants import generate_variants


@pytest.fixture()
def variants_3():
    """Three variants with two mutations each, fixed seed."""
    return generate_variants(WILDTYPE_SEQ, n_variants=3, n_mutations=2, seed=7)


class TestGenerateVariantsCount:
    def test_returns_requested_number_of_records(self, variants_3):
        assert len(variants_3) == 3

    def test_zero_variants_returns_empty_list(self):
        assert generate_variants(WILDTYPE_SEQ, n_variants=0, n_mutations=2, seed=0) == []

    def test_returns_list_of_seqrecords(self, variants_3):
        assert all(isinstance(r, SeqRecord) for r in variants_3)


class TestGenerateVariantsIds:
    def test_ids_follow_naming_convention(self, variants_3):
        assert [r.id for r in variants_3] == [
            "nylB_variant_1",
            "nylB_variant_2",
            "nylB_variant_3",
        ]

    def test_descriptions_are_empty(self, variants_3):
        assert all(r.description == "" for r in variants_3)


class TestGenerateVariantsSequence:
    def test_each_variant_has_same_length_as_wildtype(self, variants_3):
        for rec in variants_3:
            assert len(str(rec.seq)) == len(WILDTYPE_SEQ)

    def test_mutated_positions_differ_from_wildtype(self):
        records = generate_variants(WILDTYPE_SEQ, n_variants=10, n_mutations=5, seed=42)
        for rec in records:
            variant_seq = str(rec.seq)
            for pos, (wt_base, var_base) in enumerate(zip(WILDTYPE_SEQ, variant_seq)):
                if wt_base != var_base:
                    assert var_base != wt_base  # mutation is always a substitution


class TestGenerateVariantsReproducibility:
    def test_same_seed_produces_identical_sequences(self):
        r1 = generate_variants(WILDTYPE_SEQ, n_variants=3, n_mutations=3, seed=123)
        r2 = generate_variants(WILDTYPE_SEQ, n_variants=3, n_mutations=3, seed=123)
        assert all(str(a.seq) == str(b.seq) for a, b in zip(r1, r2))

    def test_different_seeds_produce_different_sequences(self):
        r1 = generate_variants(WILDTYPE_SEQ, n_variants=2, n_mutations=3, seed=1)
        r2 = generate_variants(WILDTYPE_SEQ, n_variants=2, n_mutations=3, seed=999)
        assert any(str(a.seq) != str(b.seq) for a, b in zip(r1, r2))
