FROM mambaorg/micromamba:latest

WORKDIR /pipeline

COPY environment.yml requirements.txt ./

USER root
RUN apt-get update && apt-get install -y git && rm -rf /var/lib/apt/lists/*
USER mambauser

RUN micromamba env create -f environment.yml -y

RUN micromamba run -n nylb-pipeline pip install -r requirements.txt

COPY --chown=mambauser:mambauser . .

RUN mkdir -p reference data/simulated results logs && \
    chmod -R 777 reference data/simulated results logs

CMD ["micromamba", "run", "-n", "nylb-pipeline", "bash", "-c", \
     "python3 scripts/generate_variants.py && \
      python3 scripts/prepare_libraries.py && \
      snakemake --cores 4 --scheduler greedy"]