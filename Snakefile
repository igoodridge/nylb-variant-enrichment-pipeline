import yaml
from pathlib import Path

with open("config.yaml") as f:
    cfg = yaml.safe_load(f)

SAMPLES = cfg["samples"]
WINNERS = [f"nylB_variant_{w}" for w in cfg["simulation"]["winners"]]
N_VARIANTS = cfg["simulation"]["n_variants"]
ALL_VARIANTS = ["nylB_wildtype"] + [f"nylB_variant_{i}" for i in range(1, N_VARIANTS + 1)]
BASE_QUANTITY = cfg["simulation"]["base_quantity"]
WINNER_QUANTITY = cfg["simulation"]["winner_quantity"]

wildcard_constraints:
    library="pre_selection|post_selection",
    variant="|".join(ALL_VARIANTS)


rule all:
    input:
        expand("results/qc/{sample}_stats.txt", sample=SAMPLES),
        expand("results/trimmed/{sample}.fastq", sample=SAMPLES),
        expand("results/aligned/{sample}.bam", sample=SAMPLES),
        "results/report/enrichment_report.tsv"


rule simulate_reads:
    input:
        "data/simulated/individual/{library}/{variant}.fasta"
    output:
        "data/simulated/individual/{library}/{variant}.fastq"
    log:
        "logs/simulate/{library}/{variant}.log"
    wildcard_constraints:
        library="[^/]+",
        variant="[^/]+"
    params:
        quantity=lambda wildcards: WINNER_QUANTITY if (wildcards.variant in WINNERS and wildcards.library == "post_selection") else BASE_QUANTITY
    shell:
        "badread simulate --reference {input} --quantity {params.quantity} > {output} 2> {log}"


rule merge_fastqs:
    input:
        expand("data/simulated/individual/{{library}}/{variant}.fastq", variant=ALL_VARIANTS)
    output:
        "data/simulated/{library}.fastq"
    log:
        "logs/merge/{library}.log"
    shell:
        "cat {input} > {output} 2> {log}"


rule qc:
    input:
        "data/simulated/{sample}.fastq"
    output:
        "results/qc/{sample}_stats.txt"
    log:
        "logs/qc/{sample}.log"
    shell:
        "NanoStat --fastq {input} > {output} 2> {log}"


rule trim:
    input:
        "data/simulated/{sample}.fastq"
    output:
        "results/trimmed/{sample}.fastq"
    log:
        "logs/trim/{sample}.log"
    shell:
        "porechop -i {input} -o {output} 2> {log}"


rule align:
    input:
        reads="results/trimmed/{sample}.fastq",
        reference="reference/nylB_variants.fasta"
    output:
        "results/aligned/{sample}.bam"
    log:
        "logs/align/{sample}.log"
    shell:
        "minimap2 -a {input.reference} {input.reads} 2> {log} | samtools sort > {output} && samtools index {output}"


rule count_variants:
    input:
        expand("results/aligned/{sample}.bam", sample=SAMPLES)
    output:
        "results/report/enrichment_report.tsv"
    log:
        "logs/count_variants.log"
    shell:
        "python3 scripts/count_variants.py --pre results/aligned/pre_selection.bam --post results/aligned/post_selection.bam --output {output} 2> {log}"