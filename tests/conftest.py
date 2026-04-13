"""Shared fixtures for the test suite."""

import sys
import textwrap
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Make scripts/ importable from all test modules
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))


WILDTYPE_SEQ = "ATGCATGCATGCATGCATGCATGCATGCATGC"  # 32-bp synthetic sequence

MINIMAL_CONFIG_YAML = textwrap.dedent("""\
    paths:
      wildtype: "ref/wt.fasta"
      variants: "ref/variants.fasta"
      pre_selection: "data/pre.fastq"
      post_selection: "data/post.fastq"
      individual_fastas: "data/individual"

    simulation:
      n_variants: 5
      n_mutations: 3
      seed: 99
      base_quantity: "100x"
      winner_quantity: "300x"
      winners:
        - 1
        - 3

    samples:
      - pre_selection
      - post_selection
""")


@pytest.fixture()
def config_file(tmp_path):
    """Write a minimal config.yaml and return its path."""
    cfg = tmp_path / "config.yaml"
    cfg.write_text(MINIMAL_CONFIG_YAML)
    return cfg


@pytest.fixture()
def seqrecords():
    """A small list of SeqRecord objects for use in prepare_libraries tests."""
    return [
        SeqRecord(Seq("ATGC"), id="seq_a", description=""),
        SeqRecord(Seq("TTTT"), id="seq_b", description=""),
        SeqRecord(Seq("CCCC"), id="seq_c", description=""),
    ]


@pytest.fixture()
def sample_enrichment_rows():
    """Pre-built enrichment result rows for write_report tests."""
    return [
        {
            "variant": "var1",
            "pre_count": 100,
            "post_count": 200,
            "pre_freq": 0.1,
            "post_freq": 0.2,
            "enrichment_ratio": 2.0,
        },
        {
            "variant": "var2",
            "pre_count": 900,
            "post_count": 800,
            "pre_freq": 0.9,
            "post_freq": 0.8,
            "enrichment_ratio": 0.889,
        },
    ]
