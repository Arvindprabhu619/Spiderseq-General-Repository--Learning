Bacterial WGS Analysis Pipeline Overview
Complete Workflow: Raw Data → Biological Insights
Raw FASTQ Files
       ↓
📊 Quality Control (FastQC, MultiQC)
       ↓
🧹 Preprocessing (Trimming, Filtering)
       ↓
🧩 Genome Assembly (SPAdes, Unicycler, Flye)
       ↓
✅ Assembly Quality Assessment (QUAST, BUSCO)
       ↓
🏷️ Genome Annotation (Prokka, RAST, NCBI PGAP)
       ↓
🔍 Functional Analysis (COG, KEGG, GO)
       ↓
🛡️ Resistance & Virulence Screening (ResFinder, VirulenceFinder)
       ↓
🌳 Comparative Genomics & Phylogeny
       ↓
📝 Results Interpretation & Reporting
Key Analysis Types
1. Single Isolate Analysis
Genome assembly and annotation
Species identification
Resistance/virulence profiling
Functional characterization
2. Comparative Analysis
Multi-isolate comparison
Outbreak investigation
Phylogenetic reconstruction
Pan-genome analysis
3. Population Genomics
Large-scale strain comparison
Evolution and adaptation studies
Geographic distribution analysis
Tools We'll Use
Quality Control
FastQC: Read quality assessment
MultiQC: Aggregate quality reports
Trimmomatic/fastp: Read trimming and filtering
Assembly
SPAdes: General-purpose assembler
Unicycler: Hybrid assembly (short + long reads)
QUAST: Assembly quality metrics
BUSCO: Genome completeness
Annotation
Prokka: Fast prokaryotic annotation
eggNOG-mapper: Functional annotation
ResFinder: Antimicrobial resistance genes
VirulenceFinder: Virulence factor detection
Analysis
Roary: Pan-genome analysis
FastANI: Average nucleotide identity
IQ-TREE: Phylogenetic inference
Snippy: Variant calling and core genome SNPs
Computing Requirements
Minimum
8 GB RAM
4 CPU cores
50 GB free disk space
Recommended
16+ GB RAM
8+ CPU cores
100+ GB free disk space
SSD storage for better I/O performance
Data Management
File Organization
project/
├── 00_raw_data/          # Original FASTQ files
├── 01_qc/               # Quality control reports
├── 02_trimmed/          # Processed reads
├── 03_assembly/         # Assembled genomes
├── 04_annotation/       # Annotated genomes
├── 05_analysis/         # Comparative analysis
├── 06_results/          # Final outputs
└── scripts/             # Analysis scripts
Naming Conventions
Use consistent, descriptive sample names

For bacterial whole-genome sequencing (WGS) analysis, Docker is strongly recommended, Docker for package-level control. Reproductability is high in docker compared to conda
Include metadata (strain, date, source)
Example: Ecoli_K12_MG1655_2024_lab.fastq.gz
