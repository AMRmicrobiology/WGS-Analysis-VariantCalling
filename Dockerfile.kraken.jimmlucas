FROM mambaorg/micromamba:1.5.5

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/envs/base/bin:/opt/conda/bin:$PATH

# Paquetes necesarios para el entrypoint y Kraken2
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    kraken2=2.1.3 aria2 pigz coreutils grep sed gawk curl \
    && micromamba clean -a -y

ENV KRAKEN2_DEFAULT_DB=/kraken2-db
USER root
RUN apt-get update && apt-get install -y procps
USER mambauser

COPY --chmod=0755 kraken2-entrypoint.sh /usr/local/bin/kraken2-entrypoint.sh

ENTRYPOINT []