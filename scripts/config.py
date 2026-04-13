from dataclasses import dataclass
from pathlib import Path

import yaml


@dataclass
class SimulationConfig:
    """Configuration for variant generation and library simulation.

    Attributes:
        n_variants: Number of variant sequences to generate.
        n_mutations: Number of point mutations per variant.
        seed: Random seed for reproducibility.
        base_depth: Sequencing depth for non-winner variants.
        winner_depth: Sequencing depth for winner variants.
        winners: List of variant indices designated as winners.
    """
    n_variants: int
    n_mutations: int
    seed: int
    base_depth: int
    winner_depth: int
    winners: list[int]


@dataclass
class PathConfig:
    """Configuration for input and output file paths.

    Attributes:
        wildtype: Path to wildtype FASTA file.
        variants: Path to generated variants FASTA file.
        pre_selection: Path to pre-selection library FASTA file.
        post_selection: Path to post-selection library FASTA file.
    """
    wildtype: Path
    variants: Path
    pre_selection: Path
    post_selection: Path


@dataclass
class Config:
    """Top-level configuration object.

    Attributes:
        paths: File path configuration.
        simulation: Simulation parameter configuration.
    """
    paths: PathConfig
    simulation: SimulationConfig


def load_config(config_path: str = "config.yaml") -> Config:
    """Load and parse configuration from a YAML file.

    Args:
        config_path: Path to the YAML configuration file.

    Returns:
        Parsed Config object.
    """
    with open(config_path, "r") as f:
        raw = yaml.safe_load(f)

    paths = PathConfig(
        wildtype=Path(raw["paths"]["wildtype"]),
        variants=Path(raw["paths"]["variants"]),
        pre_selection=Path(raw["paths"]["pre_selection"]),
        post_selection=Path(raw["paths"]["post_selection"]),
    )

    simulation = SimulationConfig(
        n_variants=raw["simulation"]["n_variants"],
        n_mutations=raw["simulation"]["n_mutations"],
        seed=raw["simulation"]["seed"],
        base_depth=raw["simulation"]["base_depth"],
        winner_depth=raw["simulation"]["winner_depth"],
        winners=raw["simulation"]["winners"],
    )

    return Config(paths=paths, simulation=simulation)