# NylB Variant Enrichment Pipeline

A Nanopore sequencing analysis pipeline for tracking enzyme variant enrichment in directed evolution experiments. Built as a demonstration of bioinformatics pipeline development using NylB (nylon-degrading enzyme) as a model system.

## Background

Directed evolution is a technique used to engineer proteins with improved properties. A gene of interest is mutagenised to create a library of variants, which are then subjected to selection pressure. Variants with improved activity are enriched in the surviving population. By sequencing the population before and after selection with Oxford Nanopore Technology, we can identify which variants were enriched — and therefore which mutations improved enzyme function.

This pipeline simulates that process using the NylB nylonase gene from *Pseudarthrobacter oxydans* as a reference, generating synthetic variants and simulating realistic Nanopore reads to demonstrate the full analysis workflow.

## Directory Structure

nylb-variant-enrichment-pipeline/
├── reference/              # Wildtype and variant FASTA sequences
├── data/
│   ├── simulated/          # Badread-simulated Nanopore reads
│   └── raw/                # Raw sequencing data (if applicable)
├── scripts/                # Analysis scripts
│   ├── config.py           # Configuration dataclasses
│   ├── generate_variants.py
│   └── prepare_libraries.py
├── results/                # Pipeline outputs
├── config.yaml             # Pipeline configuration
├── environment.yml         # Conda environment
└── README.md

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

Creates pre- and post-selection FASTA libraries with Badread-compatible depth annotations. Variants 2, 5, and 8 are designated as enriched winners in the post-selection library:

```bash
python3 scripts/prepare_libraries.py
```

Output: `data/simulated/pre_selection.fasta`, `data/simulated/post_selection.fasta`

### 3. Simulate Nanopore reads

Simulate realistic Nanopore reads from each library using Badread:

```bash
badread simulate --reference data/simulated/pre_selection.fasta --quantity 50x > data/simulated/pre_selection.fastq
badread simulate --reference data/simulated/post_selection.fasta --quantity 50x > data/simulated/post_selection.fastq
```

Output: `data/simulated/pre_selection.fastq`, `data/simulated/post_selection.fastq`

## Configuration

All parameters are defined in `config.yaml`:

```yaml
paths:
  wildtype: "reference/nylB_wildtype.fasta"
  variants: "reference/nylB_variants.fasta"
  pre_selection: "data/simulated/pre_selection.fasta"
  post_selection: "data/simulated/post_selection.fasta"

simulation:
  n_variants: 10
  n_mutations: 2
  seed: 42
  base_depth: 10
  winner_depth: 100
  winners:
    - 2
    - 5
    - 8
```

## Dependencies

- Python 3.12
- Biopython
- PyYAML
- Badread

See `environment.yml` for full details.