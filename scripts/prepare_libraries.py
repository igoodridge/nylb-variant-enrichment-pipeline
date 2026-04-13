import logging
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

sys.path.append(str(Path(__file__).parent))
from config import load_config

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def write_individual_fastas(
    records: list[SeqRecord],
    output_dir: Path,
    ) -> None:
    """Write individual FASTA files for each sequence.

    Args:
        records: List of SeqRecord objects to write.
        output_dir: Directory to write individual FASTA files.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    for rec in records:
        out_path = output_dir / f"{rec.id}.fasta"
        SeqIO.write(rec, out_path, "fasta")
        logger.info(f"Written {rec.id} to {out_path}")


def main() -> None:
    """Generate individual FASTA files for pre- and post-selection libraries."""
    config = load_config()

    wildtype = SeqIO.read(config.paths.wildtype, "fasta")
    wildtype.id = "nylB_wildtype"
    wildtype.description = ""
    variants = list(SeqIO.parse(config.paths.variants, "fasta"))
    all_records = [wildtype] + variants
    logger.info(f"Loaded {len(all_records)} sequences")

    pre_dir = Path(config.paths.individual_fastas) / "pre_selection"
    post_dir = Path(config.paths.individual_fastas) / "post_selection"

    write_individual_fastas(all_records, pre_dir)
    write_individual_fastas(all_records, post_dir)
    logger.info("Written individual FASTAs for pre- and post-selection libraries")


if __name__ == "__main__":
    main()