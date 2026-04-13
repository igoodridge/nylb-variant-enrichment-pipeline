from dataclasses import dataclass, field
from pathlib import Path

import yaml


@dataclass
class SimulationConfig:
    """Configuration for variant generation and library simulation.

    Attributes:
        n_variants: Number of variant sequences to generate.
        n_mutations: Number of point mutations per variant.
        seed: Random seed for reproducibility.
        base_quantity: Badread quantity for non-winner variants.
        winner_quantity: Badread quantity for winner variants.
        winners: List of variant indices designated as winners.
    """
    n_variants: int
    n_mutations: int
    seed: int
    base_quantity: str
    winner_quantity: str
    winners: list[int]


@dataclass
class PathConfig:
    """Configuration for input and output file paths.

    Attributes:
        wildtype: Path to wildtype FASTA file.
        variants: Path to generated variants FASTA file.
        pre_selection: Path to pre-selection FASTQ file.
        post_selection: Path to post-selection FASTQ file.
        individual_fastas: Directory for individual per-variant FASTA files.
    """
    wildtype: Path
    variants: Path
    pre_selection: Path
    post_selection: Path
    individual_fastas: Path


@dataclass
class Config:
    """Top-level configuration object.

    Attributes:
        paths: File path configuration.
        simulation: Simulation parameter configuration.
        samples: List of sample names to process.
    """
    paths: PathConfig
    simulation: SimulationConfig
    samples: list[str] = field(default_factory=list)


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
        individual_fastas=Path(raw["paths"]["individual_fastas"]),
    )

    simulation = SimulationConfig(
        n_variants=raw["simulation"]["n_variants"],
        n_mutations=raw["simulation"]["n_mutations"],
        seed=raw["simulation"]["seed"],
        base_quantity=raw["simulation"]["base_quantity"],
        winner_quantity=raw["simulation"]["winner_quantity"],
        winners=raw["simulation"]["winners"],
    )

    return Config(
        paths=paths,
        simulation=simulation,
        samples=raw.get("samples", []),
    )