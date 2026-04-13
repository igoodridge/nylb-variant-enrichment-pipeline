import logging
import random
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parent))

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from config import load_config

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def generate_variants(sequence: str, n_variants: int, n_mutations: int, seed: int) -> list[SeqRecord]:
    """Generate a list of mutant sequences from a wildtype sequence.

    Each variant is created by introducing a fixed number of random single
    nucleotide substitutions at randomly selected positions. Substitutions
    are guaranteed to differ from the wildtype base at each position.

    Args:
        sequence: Wildtype nucleotide sequence as a string.
        n_variants: Number of variant sequences to generate.
        n_mutations: Number of point mutations to introduce per variant.
        seed: Random seed for reproducibility.

    Returns:
        List of SeqRecord objects representing the mutant sequences.
    """
    random.seed(seed)
    len_seq = len(sequence)
    records = []

    for i in range(n_variants):
        mutant = list(sequence)
        for _ in range(n_mutations):
            pos = random.randint(0, len_seq - 1)
            base = random.choice(["A", "T", "G", "C"])
            while base == mutant[pos]:
                base = random.choice(["A", "T", "G", "C"])
            logger.info(f"Variant {i + 1}: position {pos} changed from {mutant[pos]} to {base}")
            mutant[pos] = base
        name = f"nylB_variant_{i + 1}"
        records.append(SeqRecord(Seq("".join(mutant)), id=name, description=""))
        logger.info(f"Generated variant {name}")

    return records


def main():
    """Main entry point for variant generation pipeline.

    Reads a wildtype FASTA sequence, generates synthetic variants according
    to parameters defined in config.yaml, and writes all variants to a
    single output FASTA file.
    """
    config = load_config()

    logger.info(f"Reading wildtype sequence from {config.paths.wildtype}")
    record = SeqIO.read(config.paths.wildtype, "fasta")
    sequence = str(record.seq)
    logger.info(f"Wildtype sequence length: {len(sequence)}bp")

    records = generate_variants(
        sequence,
        config.simulation.n_variants,
        config.simulation.n_mutations,
        config.simulation.seed,
    )

    SeqIO.write(records, config.paths.variants, "fasta")
    logger.info(f"Written {len(records)} variants to {config.paths.variants}")


if __name__ == "__main__":
    main()