# wgs_bacteria_project — unified bacterial WGS pipeline (Docker)

**Repository purpose**

This repository contains everything needed to build a reproducible Docker image `wgs_bacteria:1.0` that bundles a full bacterial whole-genome sequencing (WGS) toolchain (download, QC, trim, assembly, typing, AMR detection, phylogeny). The README below explains every step — installing Docker on Ubuntu, creating the project files (Dockerfile + `environment.yml`), building and running the image, running the pipeline for a single isolate, expected outputs, troubleshooting, and how to extend or reproduce the work.

> Target audience: bioinformaticians and lab scientists who want a single, reproducible container with common bacterial WGS tools (FastQC, fastp, SPAdes, Prokka/Bakta, mlst, abricate, kraken2, bwa, samtools, freebayes, snippy, iqtree, etc.).

---

## Table of contents

1. Overview
2. Requirements (host)
3. Installing Docker on Ubuntu — quick summary (official guidance)
4. Project structure
5. Files to create (Dockerfile, environment.yml) — full contents
6. Build the Docker image
7. Run & test the container
8. Example pipeline (one isolate): download → QC → trim → assembly → polish → annotate → typing → AMR → phylogeny
9. Directory layout & recommended conventions
10. Troubleshooting & common fixes
11. Reproducibility notes & pinning
12. Extending the image or building a smaller image
13. Contributing, license and citation

---

## 1. Overview

This project produces a single Docker image that contains a Conda environment named `wgs` with all tools required in a bacterial WGS analysis pipeline. The image is intended to be:

* **Reproducible**: pinned package versions where possible in `environment.yml`.
* **Portable**: runs on any machine with Docker installed (Linux, macOS via Docker Desktop, Windows via Docker Desktop).
* **Isolated**: tools run inside the container; your host OS remains unchanged.

## 2. Requirements (host)

* Ubuntu (tested on Ubuntu LTS releases). Docker supports recent Ubuntu versions — see official Docker docs for the supported list.
* Docker Engine / Docker Desktop installed and configured.
* Sufficient disk space (the built image can be several GB depending on packages). Allow at least 10–30 GB for building and working data.
* Enough RAM/CPU for assembly steps (SPAdes, IQ-TREE). For a single bacterial genome a minimum 8–16 GB RAM is recommended; more is better.

## 3. Installing Docker on Ubuntu — quick summary (official guidance)

**Recommended approach** (production-friendly): install Docker Engine from Docker's official APT repository.

**High-level steps** (run on the Ubuntu host):

```bash
sudo apt update
sudo apt install -y ca-certificates curl gnupg lsb-release

# Add Docker's official GPG key
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg

# Add the Docker apt repository (substitute $(lsb_release -cs) if needed)
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] \
  https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" \
  | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

sudo apt update
sudo apt install -y docker-ce docker-ce-cli containerd.io

# Post-installation (optional): run docker without sudo
sudo groupadd docker || true
sudo usermod -aG docker $USER
# log out and log back in or run: newgrp docker
```

**Alternative**: Docker Desktop for Linux (if you prefer a desktop GUI and integrated Kubernetes). For servers or production the Engine installation above is typical.

**References (official docs)**:

* Docker docs — Install Docker Engine on Ubuntu.

---

## 4. Project structure (expected)

Create your project root (example):

```
~/wgs_bacteria_project/
├── Dockerfile
├── environment.yml
├── README.md          # this file
├── LICENSE
├── data/              # place raw FASTQ or .sra here (NOT required to build image)
└── scripts/           # optional helper scripts (run-pipeline.sh, check-env.sh)
```

Tip: keep raw data outside the image and mount it into `/work` when running the container.

---

## 5. Files to create

Below are the exact contents you should place into your project files.

### `Dockerfile`

```Dockerfile
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
```

**Notes**:

* Starts from `continuumio/miniconda3` to keep package management reproducible via Conda.
* System packages installed via `apt-get` are minimal but necessary for compiling some packages and for common utilities.
* `mamba` is used (faster, robust) to create the `wgs` environment from `environment.yml`.
* `SHELL` is switched to run inside the Conda environment so tools are available in next layers.

### `environment.yml`

```yaml
name: wgs
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - mamba
  - sra-tools
  - fastp=0.23.2
  - fastqc=0.12.1
  - multiqc=1.15
  - spades
  - unicycler=0.5.0
  - quast=5.2.0
  - prokka=1.14.6
  - mlst=2.19.0
  - abricate=1.0.1
  - kraken2=2.1.2
  - bwa=0.7.17
  - samtools=1.17
  - bcftools=1.17
  - freebayes=1.3.5
  - mash=2.3
  - fastani=1.33
  - snippy=4.6.0
  - iqtree=2.2.2
  - python=3.11
  - pip
  - pip:
      - pyfastx
```

**Pin versions where possible**. For packages that are known to be sensitive to the toolchain (e.g. SPAdes), you may choose to pin a version or leave Cona to choose a compatible build.

---

## 6. Build the Docker image

From your project root (e.g. `~/wgs_bacteria_project`):

```bash
# build image (runs Dockerfile, installs system deps and conda env)
docker build -t wgs_bacteria:1.0 .
```

**What happens during build**:

* Docker fetches the base image (`continuumio/miniconda3:latest`).
* Installs system packages with `apt-get`.
* Copies `environment.yml` into the image and creates the `wgs` conda environment using `mamba`.
* Cleans cache to reduce size.

**Build tips**:

* If the build fails because of a mis-resolved conda package, edit `environment.yml` to adjust/remove the problematic package or pin alternate versions.
* Building behind a proxy or in an air-gapped environment requires extra steps (not covered here).

---

## 7. Run & test the container

A minimal smoke test — run transient container and confirm tools are present:

```bash
# run interactive container, auto-remove after exit
docker run -it --rm wgs_bacteria:1.0
# you will be dropped into /work with base conda active
# inside container, list conda envs and check versions
conda env list
conda activate wgs
fastqc --version
spades.py --version
prokka --version
```

**Mounting your project/data**

Always mount host directories with `-v` when you want container to access data and persist results:

```bash
docker run -it --rm -v /home/$USER/wgs_bacteria_project:/work wgs_bacteria:1.0
# inside container work dir: /work
```

**Post-installation (optional, host)**
Follow Docker's post-installation steps to run Docker without `sudo` (add your user to the `docker` group). See Docker docs for details.

---

## 8. Example pipeline (single isolate)

Below is a step-by-step example you can follow inside the running container (assumes you mounted your project folder into `/work` and that FASTQ or SRA files are in `/work`):

### 8.0 Prepare (start container)

```bash
# host
docker run -it --rm -v /home/$USER/wgs_bacteria_project:/work wgs_bacteria:1.0
# inside container
cd /work
conda activate wgs
```

### 8.1 Data download (optional — SRA)

If you have SRA accessions (example `SRR36198494`):

```bash
# ensure sra-tools work; if sra-tools inside conda fails due to SSL/CA issues, see Troubleshooting
prefetch SRR36198494
fasterq-dump SRR36198494 --split-files --gzip
# results: SRR36198494_1.fastq.gz SRR36198494_2.fastq.gz
```

### 8.2 Raw QC — FastQC

```bash
mkdir -p fastqc_raw
fastqc SRR36198494_1.fastq.gz SRR36198494_2.fastq.gz -o fastqc_raw
```

### 8.3 Trimming — fastp

```bash
fastp \
  -i SRR36198494_1.fastq.gz \
  -I SRR36198494_2.fastq.gz \
  -o trimmed_1.fastq.gz \
  -O trimmed_2.fastq.gz \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --length_required 50 \
  --thread 4 \
  --html fastp_report.html \
  --json fastp_report.json

# QC trimmed reads
mkdir -p fastqc_trimmed
fastqc trimmed_1.fastq.gz trimmed_2.fastq.gz -o fastqc_trimmed
```

### 8.4 MultiQC summary

```bash
multiqc .
# outputs multiqc_report.html summarizing all FastQC / fastp reports
```

### 8.5 Assembly — SPAdes

```bash
spades.py \
  -1 trimmed_1.fastq.gz \
  -2 trimmed_2.fastq.gz \
  -o spades_output \
  --careful \
  -t 4 \
  -m 16

# main assembly fasta: spades_output/contigs.fasta (and scaffolds.fasta)
```

### 8.6 Assembly QC — QUAST

```bash
quast.py spades_output/contigs.fasta -o quast_results -t 4
```

### 8.7 Polishing — Pilon (external jar) or other

Pilon is Java-based and not always in conda. You can download pilon.jar and run it.

```bash
# index assembly
bwa index spades_output/scaffolds.fasta

# map reads back
bwa mem -t 4 spades_output/scaffolds.fasta trimmed_1.fastq.gz trimmed_2.fastq.gz | samtools sort -o reads.bam
samtools index reads.bam

# run Pilon (example jar location)
java -Xmx16G -jar pilon-1.24.jar --genome spades_output/scaffolds.fasta --frags reads.bam --output polished_scaffolds --threads 4
```

Repeat mapping + Pilon for 1–2 extra rounds if desired.

### 8.8 Annotation — Bakta or Prokka

Example using Prokka (included in environment):

```bash
prokka --outdir prokka_out --prefix SAMPLE polished_scaffolds.fasta --cpus 4
```

Or Bakta (recommended for high-quality bacterial annotation). Bakta requires downloading a database: `bakta_db download --type light --output bakta_db`.

### 8.9 AMR — Abricate

```bash
abricate --setupdb   # run once to fetch databases
abricate --db ncbi polished_scaffolds.fasta > abricate_ncbi.tsv
```

### 8.10 MLST

```bash
mlst polished_scaffolds.fasta > mlst_result.txt
```

### 8.11 Species identity — fastANI

```bash
# fastANI vs reference fasta (if available)
fastani -q polished_scaffolds.fasta -r /work/ref_genomes/closest_ref.fna -o fastani_out.txt
```

### 8.12 SNP calling & phylogeny (snippy + iqtree)

Map multiple genomes to same reference with `snippy`, run `snippy-core` to obtain core alignment, then build ML tree with IQ-TREE.

```bash
# run snippy for each genome (example)
snippy --cpus 4 --ref /work/ref_genomes/WJL_ref.fna --ctgs polished_scaffolds.fasta --outdir snippy_sample

# combine folders into core alignment
after running snippy for all samples:

snippy-core --ref /work/ref_genomes/WJL_ref.fna --prefix core_all snippy_sample snippy_other1 snippy_other2

# build tree with IQ-TREE
iqtree -s core_all.full.aln -m GTR+G -bb 1000 -nt AUTO
```

---

## 9. Directory layout & recommended conventions

Keep results organized with a consistent layout inside the mounted project folder:

```
wgs_bacteria_project/
├── raw_data/            # raw FASTQ or .sra
├── fastqc_raw/
├── trimmed_reads/
├── spades_output/
├── quast_results/
├── polished_scaffolds.fasta
├── prokka_out/
├── abricate_results.tsv
├── ref_genomes/
├── snippy_out/
├── phylogeny/
└── metadata.csv
```

Use symbolic sample IDs and keep a `metadata.tsv` to track sample collection date, source, location, and any lab metadata.

---

## 10. Troubleshooting & common fixes

### 10.1 Conda/mamba environment fails to solve during image build

* Edit `environment.yml` to remove or change the offending package or pin a different version.
* If a Bioconda binary is incompatible, try specifying channel priority (`conda config --set channel_priority strict`) locally or add `--strict-channel-priority` when resolving locally.

### 10.2 sra-tools/NCBI SSL or CA problems inside container

* Update CA certs inside container: `apt-get update && apt-get install -y ca-certificates && update-ca-certificates`.
* Or use `wget` to download SRA assembly summary and fetch FASTA files via HTTP as described in the examples.

### 10.3 FASTQ downloader fails (prefetch/fasterq-dump)

* Try `prefetch` first, then `fasterq-dump`. If `fasterq-dump` errors, use `fastq-dump` or the newer `fasterq-dump` with `--split-files` and `--gzip`.
* If sra-tools from conda keeps failing, remove the package from the `wgs` env and install the official SRA Toolkit binary manually inside the container (instructions provided in container). Example steps are included in the pipeline notes.

### 10.4 SPAdes throws memory error

* Increase `-m` value for SPAdes (MB). Example: `-m 32` for 32 GB RAM.
* Use `--only-assembler` or `--careful` depending on needs.

### 10.5 Docker build slow or fails due to network

* Retry build; mamba/conda solves packages over network and can occasionally fail due to mirrors.
* Run build on a machine with stable network, or create a custom cache of conda packages.

---

## 11. Reproducibility notes & pinning

* Pin versions when stability is required (we pinned many tools in `environment.yml`).
* For absolute reproducibility across time, consider using conda-lock or `mamba env export --file` from a known-good environment to permanently fix package hashes.
* Version containers using tags (e.g. `wgs_bacteria:1.0`) and store the `Dockerfile` and `environment.yml` in your VCS.

---

## 12. Extending the image or making it smaller

**To add packages**: edit `environment.yml` and rebuild the image.

**To reduce image size**:

* Use `micromamba` as base or a smaller conda base image.
* Remove documentation and caches, strip unused files, or move heavy tools (like large databases) out of the image and mount them at runtime.

---

## 13. Contributing, license & citation

* Add a `CONTRIBUTING.md` for contribution guidelines.
* Include a license (e.g., MIT) in `LICENSE`.
* Cite the tools you used in your publications (Prokka, SPAdes, FastQC, QUAST, IQ-TREE, etc.).

---

### Example `LICENSE` snippet (MIT)

```
MIT License

Copyright (c) YEAR AUTHOR

Permission is hereby granted, free of charge, to any person obtaining a copy
...
```

---

### Final notes

This README aims to be a complete, highly actionable reference for building and using a reproducible Docker image for bacterial WGS analysis. If you want, I can also:

* Generate a `run_pipeline.sh` bash script that automates the entire example pipeline with sensible defaults and checks.
* Produce a Snakemake or Nextflow pipeline that executes every step reproducibly for multiple samples.

If you want either, tell me which (bash script, Snakemake, or Nextflow) and I will create it next.
