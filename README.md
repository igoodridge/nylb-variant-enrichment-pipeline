# NylB Variant Enrichment Pipeline

A Nanopore sequencing analysis pipeline for tracking enzyme variant enrichment in directed evolution experiments. Built as a demonstration of bioinformatics pipeline development using NylB (nylon-degrading enzyme) as a model system.

## Background

Directed evolution is a technique used to engineer proteins with improved properties. A gene of interest is mutagenised to create a library of variants, which are then subjected to selection pressure. Variants with improved activity are enriched in the surviving population. By sequencing the population before and after selection with Oxford Nanopore Technology, we can identify which variants were enriched — and therefore which mutations improved enzyme function.

This pipeline simulates that process using the NylB nylonase gene from *Pseudarthrobacter oxydans* as a reference, generating synthetic variants and simulating realistic Nanopore reads to demonstrate the full analysis workflow.

## Directory Structure

nylb-variant-enrichment-pipeline/
```
nylb-variant-enrichment-pipeline/
├── reference/                  # Wildtype and variant FASTA sequences
├── data/
│   └── simulated/
│       ├── individual/         # Per-variant FASTA and FASTQ files
│       │   ├── pre_selection/
│       │   └── post_selection/
│       ├── pre_selection.fastq # Merged pre-selection reads
│       └── post_selection.fastq# Merged post-selection reads
├── results/
│   ├── qc/                     # NanoStat quality reports
│   ├── trimmed/                # Porechop-trimmed reads
│   ├── aligned/                # Minimap2 BAM files
│   └── report/                 # Enrichment report
├── logs/                       # Per-rule logs
├── scripts/
│   ├── config.py               # Configuration dataclasses
│   ├── generate_variants.py    # Synthetic variant generation
│   ├── prepare_libraries.py    # Per-variant FASTA preparation
│   └── count_variants.py       # Enrichment analysis
├── Snakefile                   # Snakemake workflow
├── config.yaml                 # Pipeline configuration
├── environment.yml             # Conda environment
└── README.md
```

## Installation

Clone the repository and create the conda environment:

```bash
git clone https://github.com/isabellegoodridge/nylb-variant-enrichment-pipeline.git
cd nylb-variant-enrichment-pipeline
conda env create -f environment.yml
conda activate nylb-pipeline
pip install git+https://github.com/rrwick/Badread.git
```

## Usage

Run each step in order from the project root directory.

### 1. Generate synthetic variants

Creates 10 synthetic NylB variants with random point mutations from the wildtype reference sequence:

```bash
python3 scripts/generate_variants.py
```

Output: `reference/nylB_variants.fasta`

### 2. Prepare selection libraries

Creates individual per-variant FASTA files for pre- and post-selection libraries:

```bash
python3 scripts/prepare_libraries.py
```

Output: `data/simulated/individual/pre_selection/` and `data/simulated/individual/post_selection/`

### 3. Run the Snakemake workflow

Runs all remaining steps automatically — Nanopore read simulation, QC, trimming, alignment, and enrichment analysis:

```bash
snakemake --cores 4 --scheduler greedy
```

The workflow will:
- Simulate Nanopore reads per variant using Badread (winners at 100x, others at 10x)
- Merge per-variant FASTQs into pre- and post-selection libraries
- Run QC with NanoStat
- Trim adapters with Porechop
- Align reads to reference variants with Minimap2
- Generate enrichment report

Output: `results/report/enrichment_report.tsv`

## Configuration

All parameters are defined in `config.yaml`:

```yaml
paths:
  wildtype: "reference/nylB_wildtype.fasta"
  variants: "reference/nylB_variants.fasta"
  pre_selection: "data/simulated/pre_selection.fastq"
  post_selection: "data/simulated/post_selection.fastq"
  individual_fastas: "data/simulated/individual"

simulation:
  n_variants: 10
  n_mutations: 15
  seed: 42
  base_quantity: "10x"
  winner_quantity: "100x"
  winners:
    - 2
    - 5
    - 8

samples:
  - pre_selection
  - post_selection
```

Note: `n_mutations` is set to 15 to ensure variants are sufficiently distinct for unambiguous read alignment. With fewer mutations, sequences are too similar (~99.8% identical at 2 mutations in a 1179bp gene) and Minimap2 distributes reads ambiguously across variants, flattening the enrichment signal.

## Expected output

The enrichment report (`results/report/enrichment_report.tsv`) should show variants 2, 5, and 8 with enrichment ratios of approximately 10, reflecting their 10x higher sequencing depth in the post-selection library. All other variants should have ratios close to 1.

## Dependencies

- Python 3.12
- Biopython
- PyYAML
- pysam
- Snakemake
- NanoStat
- Porechop
- Minimap2
- Samtools
- Badread (installed via pip)

See `environment.yml` for full details.