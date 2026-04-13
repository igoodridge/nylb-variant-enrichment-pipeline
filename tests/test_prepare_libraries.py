"""Tests for scripts/prepare_libraries.py — write_individual_fastas."""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from prepare_libraries import write_individual_fastas


class TestWriteIndividualFastas:
    def test_creates_one_file_per_record(self, tmp_path, seqrecords):
        write_individual_fastas(seqrecords, tmp_path)
        assert len(list(tmp_path.glob("*.fasta"))) == len(seqrecords)

    def test_files_are_named_by_record_id(self, tmp_path, seqrecords):
        write_individual_fastas(seqrecords, tmp_path)
        for rec in seqrecords:
            assert (tmp_path / f"{rec.id}.fasta").exists()

    def test_creates_output_directory_when_missing(self, tmp_path, seqrecords):
        nested = tmp_path / "deep" / "nested"
        write_individual_fastas(seqrecords, nested)
        assert nested.is_dir()

    def test_fasta_content_contains_id_and_sequence(self, tmp_path):
        rec = SeqRecord(Seq("AATTCCGG"), id="myseq", description="")
        write_individual_fastas([rec], tmp_path)
        content = (tmp_path / "myseq.fasta").read_text()
        assert ">myseq" in content
        assert "AATTCCGG" in content

    def test_empty_records_list_writes_no_files(self, tmp_path):
        write_individual_fastas([], tmp_path)
        assert list(tmp_path.glob("*.fasta")) == []
