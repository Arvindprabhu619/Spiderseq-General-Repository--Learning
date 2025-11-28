# Dockerfile — unified bacterial WGS image
FROM continuumio/miniconda3:latest

ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      build-essential ca-certificates curl wget git unzip bzip2 \
      libbz2-dev liblzma-dev libzstd1 pigz procps && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

COPY environment.yml /tmp/environment.yml

RUN conda install -y -n base -c conda-forge mamba && \
    mamba env create -n wgs -f /tmp/environment.yml && \
    conda clean -afy && \
    rm -f /tmp/environment.yml

SHELL ["conda", "run", "-n", "wgs", "/bin/bash", "-lc"]
ENV PATH=/opt/conda/envs/wgs/bin:$PATH
WORKDIR /work
ENTRYPOINT [ "/bin/bash" ]
