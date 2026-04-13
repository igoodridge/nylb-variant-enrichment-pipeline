"""Tests for scripts/config.py — load_config."""

import textwrap
from pathlib import Path

import pytest

from config import Config, PathConfig, SimulationConfig, load_config


class TestLoadConfig:
    def test_returns_config_object(self, config_file):
        assert isinstance(load_config(str(config_file)), Config)

    def test_path_fields_are_path_objects(self, config_file):
        paths = load_config(str(config_file)).paths
        assert isinstance(paths.wildtype, Path)
        assert isinstance(paths.variants, Path)
        assert isinstance(paths.pre_selection, Path)
        assert isinstance(paths.post_selection, Path)
        assert isinstance(paths.individual_fastas, Path)

    def test_path_values(self, config_file):
        paths = load_config(str(config_file)).paths
        assert paths.wildtype == Path("ref/wt.fasta")
        assert paths.variants == Path("ref/variants.fasta")
        assert paths.pre_selection == Path("data/pre.fastq")
        assert paths.post_selection == Path("data/post.fastq")
        assert paths.individual_fastas == Path("data/individual")

    def test_simulation_values(self, config_file):
        sim = load_config(str(config_file)).simulation
        assert sim.n_variants == 5
        assert sim.n_mutations == 3
        assert sim.seed == 99
        assert sim.base_quantity == "100x"
        assert sim.winner_quantity == "300x"
        assert sim.winners == [1, 3]

    def test_samples_list(self, config_file):
        assert load_config(str(config_file)).samples == ["pre_selection", "post_selection"]

    def test_missing_samples_key_defaults_to_empty_list(self, tmp_path):
        cfg = tmp_path / "config.yaml"
        cfg.write_text(textwrap.dedent("""\
            paths:
              wildtype: "wt.fasta"
              variants: "v.fasta"
              pre_selection: "pre.fastq"
              post_selection: "post.fastq"
              individual_fastas: "ind"
            simulation:
              n_variants: 1
              n_mutations: 1
              seed: 0
              base_quantity: "10x"
              winner_quantity: "10x"
              winners: []
        """))
        assert load_config(str(cfg)).samples == []
