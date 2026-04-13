import logging
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parent))

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from config import load_config

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def write_library(records: list[SeqRecord], depths: dict, output_path: Path) -> None:
    """Write a FASTA library with Badread-compatible depth annotations.

    Args:
        records: List of SeqRecord objects to write.
        depths: Dictionary mapping sequence IDs to depth values.
        output_path: Path to write the output FASTA file.
    """
    output_records = []
    for rec in records:
        depth = depths.get(rec.id, 10)
        new_rec = SeqRecord(rec.seq, id=rec.id, description=f"depth={depth}")
        output_records.append(new_rec)
        logger.info(f"Added {rec.id} with depth={depth}")

    SeqIO.write(output_records, output_path, "fasta")
    logger.info(f"Written {len(output_records)} sequences to {output_path}")


def main():
    """Generate pre- and post-selection FASTA libraries for Badread simulation."""
    config = load_config()

    wildtype = SeqIO.read(config.paths.wildtype, "fasta")
    variants = list(SeqIO.parse(config.paths.variants, "fasta"))
    all_records = [wildtype] + variants
    logger.info(f"Loaded {len(all_records)} sequences")

    # Pre-selection: all sequences at equal depth
    pre_depths = {rec.id: config.simulation.base_depth for rec in all_records}
    write_library(all_records, pre_depths, config.paths.pre_selection)

    # Post-selection: winner variants at higher depth
    post_depths = {rec.id: config.simulation.base_depth for rec in all_records}
    for winner in config.simulation.winners:
        name = f"nylB_variant_{winner}"
        post_depths[name] = config.simulation.winner_depth
        logger.info(f"Set winner {name} to depth={config.simulation.winner_depth}")
    write_library(all_records, post_depths, config.paths.post_selection)


if __name__ == "__main__":
    main()