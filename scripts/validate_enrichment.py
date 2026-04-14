import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

sys.path.append(str(Path(__file__).parent))
from config import load_config

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def load_enrichment_report(report_path: str) -> pd.DataFrame:
    """Load enrichment report from a TSV file.

    Args:
        report_path: Path to the enrichment report TSV.

    Returns:
        DataFrame containing enrichment results.
    """
    df = pd.read_csv(report_path, sep="\t")
    df = df.sort_values("enrichment_ratio", ascending=False).reset_index(drop=True)
    return df


def validate(df: pd.DataFrame, winners: list[str], expected_ratio: float) -> dict:
    """Validate enrichment results against known ground truth.

    Args:
        df: Enrichment report dataframe sorted by enrichment ratio.
        winners: List of variant IDs expected to be enriched.
        expected_ratio: Expected enrichment ratio (winner_quantity / base_quantity).

    Returns:
        Dictionary containing validation results.
    """
    n_winners = len(winners)
    top_variants = df.head(n_winners)["variant"].tolist()

    correctly_identified = all(w in top_variants for w in winners)
    n_correct = sum(1 for w in winners if w in top_variants)

    winner_rows = df[df["variant"].isin(winners)]
    observed_ratios = dict(zip(winner_rows["variant"], winner_rows["enrichment_ratio"]))
    mean_observed = winner_rows["enrichment_ratio"].mean()

    return {
        "passed": correctly_identified,
        "n_winners": n_winners,
        "n_correctly_identified": n_correct,
        "expected_ratio": expected_ratio,
        "mean_observed_ratio": round(mean_observed, 3),
        "observed_ratios": observed_ratios,
        "top_variants": top_variants,
    }


def write_validation_report(results: dict, output_path: str) -> None:
    """Write validation report to a text file.

    Args:
        results: Dictionary of validation results.
        output_path: Path to write the validation report.
    """
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    status = "PASSED" if results["passed"] else "FAILED"

    lines = [
        "=" * 50,
        f"ENRICHMENT VALIDATION REPORT — {status}",
        "=" * 50,
        f"Winners correctly identified: {results['n_correctly_identified']}/{results['n_winners']}",
        f"Expected enrichment ratio:    {results['expected_ratio']:.1f}x",
        f"Mean observed ratio:          {results['mean_observed_ratio']:.3f}x",
        "",
        "Per-winner observed ratios:",
    ]

    for variant, ratio in sorted(results["observed_ratios"].items()):
        lines.append(f"  {variant}: {ratio:.3f}x")

    lines += [
        "",
        f"Top {results['n_winners']} variants by enrichment ratio:",
    ]
    for i, variant in enumerate(results["top_variants"], 1):
        lines.append(f"  {i}. {variant}")

    with open(output_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    logger.info(f"Written validation report to {output_path}")
    logger.info(f"Validation {status}")


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    Returns:
        Parsed arguments namespace.
    """
    parser = argparse.ArgumentParser(description="Validate enrichment pipeline against known ground truth.")
    parser.add_argument("--report", required=True, help="Path to enrichment report TSV")
    parser.add_argument("--output", required=True, help="Path to output validation report")
    return parser.parse_args()


def main() -> None:
    """Main entry point for enrichment validation."""
    args = parse_args()
    config = load_config()

    logger.info(f"Loading enrichment report from {args.report}")
    df = load_enrichment_report(args.report)

    winners = [f"nylB_variant_{w}" for w in config.simulation.winners]
    expected_ratio = float(config.simulation.winner_quantity.replace("x", "")) / \
                     float(config.simulation.base_quantity.replace("x", ""))
    logger.info(f"Expected enrichment ratio: {expected_ratio:.1f}x")

    results = validate(df, winners, expected_ratio)
    write_validation_report(results, args.output)


if __name__ == "__main__":
    main()