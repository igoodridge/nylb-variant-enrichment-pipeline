import argparse
import logging
from pathlib import Path
import sys

import pysam

sys.path.append(str(Path(__file__).parent))
from config import load_config

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def count_reads(bam_path: str) -> dict[str, int]:
    """Count reads mapping to each reference sequence in a BAM file.

    Args:
        bam_path: Path to the BAM file.

    Returns:
        Dictionary mapping reference sequence names to read counts.
    """
    counts = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.is_unmapped:
                continue
            ref = read.reference_name
            counts[ref] = counts.get(ref, 0) + 1
    logger.info(f"Counted reads in {bam_path}: {sum(counts.values())} total mapped reads")
    return counts


def calculate_frequencies(counts: dict[str, int]) -> dict[str, float]:
    """Calculate relative frequency of each reference sequence.

    Args:
        counts: Dictionary mapping reference names to read counts.

    Returns:
        Dictionary mapping reference names to frequencies.
    """
    total = sum(counts.values())
    return {ref: count / total for ref, count in counts.items()}


def calculate_enrichment(
    pre_counts: dict[str, int],
    post_counts: dict[str, int]
) -> list[dict]:
    """Calculate enrichment ratios between post- and pre-selection libraries.

    Args:
        pre_counts: Read counts from pre-selection BAM.
        post_counts: Read counts from post-selection BAM.

    Returns:
        List of dicts with per-variant counts, frequencies, and enrichment ratios.
    """
    pre_freqs = calculate_frequencies(pre_counts)
    post_freqs = calculate_frequencies(post_counts)

    all_refs = sorted(set(pre_counts.keys()) | set(post_counts.keys()))
    rows = []
    for ref in all_refs:
        pre_count = pre_counts.get(ref, 0)
        post_count = post_counts.get(ref, 0)
        pre_freq = pre_freqs.get(ref, 0)
        post_freq = post_freqs.get(ref, 0)
        if pre_freq > 0:
            enrichment = post_freq / pre_freq
        else:
            enrichment = 9999.0   # sentinel: variant absent pre-selection
        rows.append({
            "variant": ref,
            "pre_count": pre_count,
            "post_count": post_count,
            "pre_freq": round(pre_freq, 4),
            "post_freq": round(post_freq, 4),
            "enrichment_ratio": round(enrichment, 3),
        })
        logger.info(f"{ref}: enrichment ratio = {enrichment:.3f}")
    return rows


def write_report(rows: list[dict], output_path: str) -> None:
    """Write enrichment report to a TSV file.

    Args:
        rows: List of per-variant result dicts.
        output_path: Path to write the TSV output.
    """
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    headers = ["variant", "pre_count", "post_count", "pre_freq", "post_freq", "enrichment_ratio"]
    with open(output_path, "w") as f:
        f.write("\t".join(headers) + "\n")
        for row in rows:
            f.write("\t".join(str(row[h]) for h in headers) + "\n")
    logger.info(f"Written enrichment report to {output_path}")


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    Returns:
        Parsed arguments namespace.
    """
    parser = argparse.ArgumentParser(description="Count variant read frequencies and calculate enrichment.")
    parser.add_argument("--pre", required=True, help="Path to pre-selection BAM file")
    parser.add_argument("--post", required=True, help="Path to post-selection BAM file")
    parser.add_argument("--output", required=True, help="Path to output TSV file")
    return parser.parse_args()


def main() -> None:
    """Main entry point for variant enrichment analysis."""
    args = parse_args()

    logger.info("Counting reads in pre-selection BAM")
    pre_counts = count_reads(args.pre)

    logger.info("Counting reads in post-selection BAM")
    post_counts = count_reads(args.post)

    logger.info("Calculating enrichment ratios")
    rows = calculate_enrichment(pre_counts, post_counts)

    write_report(rows, args.output)


if __name__ == "__main__":
    main()