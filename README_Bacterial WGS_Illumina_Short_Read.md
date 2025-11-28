Unified Bacterial WGS Pipeline using Docker

This repository contains a fully reproducible workflow for bacterial whole-genome sequencing (WGS) analysis using Docker. The pipeline covers:

Data download from NCBI SRA

Quality control (FastQC)

Read trimming (fastp)

Genome assembly (SPAdes)

Assembly polishing (Pilon)

Genome annotation (Bakta / Prokka)

AMR, virulence, plasmid detection

Species identification, MLST, ANI

Core genome SNP phylogeny (Snippy + IQ-TREE)

All tools are installed inside a Docker container, ensuring your host OS remains clean.

Table of Contents

Requirements

Setup: Docker & Project

Create & Build Docker Image

Running the Container

Download & Preprocess SRA Data

Quality Control & Trimming

Genome Assembly & Polishing

Genome Annotation

AMR, Virulence & Plasmid Detection

Species Identification & MLST

Core Genome Alignment & Phylogeny

References

Requirements

Ubuntu 20.04+

Docker (latest stable version)

At least 16–32 GB RAM for assembly/polishing

CPU: 4+ cores recommended

Setup Docker & Project

Install Docker on Ubuntu:

sudo apt update
sudo apt install -y \
    apt-transport-https \
    ca-certificates \
    curl \
    software-properties-common

curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -

sudo add-apt-repository \
   "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"

sudo apt update
sudo apt install -y docker-ce
sudo systemctl enable docker
sudo systemctl start docker
sudo usermod -aG docker $USER


Why: Docker allows us to isolate all WGS tools inside a reproducible container.

Create project directory:

mkdir -p ~/wgs_bacteria_project
cd ~/wgs_bacteria_project

Create & Build Docker Image

Step 1 — Create Dockerfile

cat > Dockerfile <<'EOF'
# Unified bacterial WGS Docker image
FROM continuumio/miniconda3:latest

ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Install system packages
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
EOF


Step 2 — Create environment.yml

cat > environment.yml <<'EOF'
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
EOF


Why: Pinning versions ensures reproducibility. Conda/mamba installs all bioinformatics tools in an isolated environment.

Step 3 — Build Docker image

docker build -t wgs_bacteria:1.0 .


Result: A self-contained Docker image with all tools for bacterial WGS.

Running the Container
docker run -it --rm -v ~/wgs_bacteria_project:/work wgs_bacteria:1.0


Inside container:

conda activate wgs
cd /work


Why: Mounting the host folder allows access to input/output data. The WGS environment is activated.

Download & Preprocess SRA Data
# Download SRA
prefetch SRR36198494
# Convert to gzipped FASTQ
fasterq-dump SRR36198494 --split-files --gzip


Why: Converts SRA format to FASTQ for downstream QC/assembly.

Quality Control & Trimming

Step 1 — FastQC on raw reads

mkdir -p fastqc_output
fastqc SRR36198494_1.fastq.gz SRR36198494_2.fastq.gz -o fastqc_output


Step 2 — Trim reads using fastp

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


Step 3 — FastQC again on trimmed reads

mkdir -p fastqc_trimmed
fastqc trimmed_1.fastq.gz trimmed_2.fastq.gz -o fastqc_trimmed
multiqc .


Why: Ensures reads are high-quality and adapters are removed.

Genome Assembly & Polishing

SPAdes assembly:

spades.py \
  -1 trimmed_1.fastq.gz \
  -2 trimmed_2.fastq.gz \
  -o spades_output \
  --careful \
  -t 4 \
  -m 16


Assembly quality (QUAST):

quast.py spades_output/contigs.fasta -o quast_results/ -t 4


Polishing (Pilon):

bwa index spades_output/scaffolds.fasta
bwa mem -t 4 spades_output/scaffolds.fasta trimmed_1.fastq.gz trimmed_2.fastq.gz | samtools sort -o reads.bam
samtools index reads.bam
pilon --genome spades_output/scaffolds.fasta --frags reads.bam --output polished_scaffolds --threads 4


Why: Corrects small errors in the draft assembly using mapped reads.

Genome Annotation

Bakta (recommended):

conda install -c conda-forge -c bioconda bakta
mkdir -p bakta_db
bakta_db download --type light --output bakta_db
bakta --db /work/bakta_db/db-light --threads 4 --output bakta_output polished_scaffolds.fasta


Optional Prokka:

prokka --outdir prokka_output --prefix polished polished_scaffolds.fasta

AMR, Virulence & Plasmid Detection
conda install -c bioconda abricate
abricate --setupdb
abricate polished_scaffolds.fasta > abricate_results.tsv


Why: Detects antimicrobial resistance genes, virulence factors, and plasmids.

Species Identification & MLST
conda install -y -c bioconda fastani mlst
mlst polished_scaffolds.fasta > mlst_result.txt


Optional ANI calculation:

fastani -q polished_scaffolds.fasta -r ref_genomes/*.fna -o ani_results.txt

Core Genome Alignment & Phylogeny

Snippy SNP calling:

snippy --cpus 4 --ref ref_genomes/WJL_ref.fna --ctgs polished_scaffolds.fasta --outdir snippy_WJL
snippy-core --ref ref_genomes/WJL_ref.fna --prefix core_alignment snippy_WJL snippy_other_refs...


Phylogenetic tree (IQ-TREE):

iqtree -s core_alignment.full.aln -m GTR+G -bb 1000 -nt AUTO


Why: Generates core genome alignment and maximum-likelihood tree with bootstrap support.

References

Docker Documentation

Conda & Mamba

SPAdes

FastQC

fastp

QUAST

Pilon

Bakta

Snippy

IQ-TREE

✅ This README can be saved as README.md in your GitHub repo, along with Dockerfile and environment.yml. It provides full reproducible instructions, commands, and explanations for bacterial WGS using Docker.
