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


Practical WGS Analysis Commands
🔧 Tool Installation and Setup
Conda Environment Setup (Recommended)
# Create dedicated environment
conda create -n wgs_analysis python=3.9
conda activate wgs_analysis

# Install core tools
conda install -c bioconda fastqc multiqc trimmomatic spades quast prokka
conda install -c conda-forge biopython pandas numpy matplotlib seaborn

# Additional useful tools
conda install -c bioconda unicycler snippy roary iqtree fasttree
Docker Alternative (Containerized)
# Pull popular bioinformatics containers
docker pull staphb/fastqc:latest
docker pull staphb/spades:latest
docker pull staphb/prokka:latest
📁 Data Organization Best Practices
# Create project structure
mkdir bacterial_wgs_project
cd bacterial_wgs_project
mkdir -p {00_raw_data,01_qc,02_trimmed,03_assembly,04_annotation,05_analysis,06_results,scripts}

# Example directory structure
bacterial_wgs_project/
├── 00_raw_data/          # Original FASTQ files
│   ├── sample1_R1.fastq.gz
│   ├── sample1_R2.fastq.gz
│   └── samples.txt       # Sample metadata
├── 01_qc/               # Quality control reports
├── 02_trimmed/          # Processed reads
├── 03_assembly/         # Assembled genomes
├── 04_annotation/       # Annotated genomes
├── 05_analysis/         # Comparative analysis
├── 06_results/          # Final outputs
└── scripts/             # Analysis scripts
🔍 Step 1: Quality Control
FastQC - Per-sample quality assessment
# Single sample
fastqc 00_raw_data/sample1_R1.fastq.gz -o 01_qc/

# Multiple samples
fastqc 00_raw_data/*.fastq.gz -o 01_qc/ -t 4

# Batch processing script
for file in 00_raw_data/*_R1.fastq.gz; do
    base=$(basename $file _R1.fastq.gz)
    fastqc 00_raw_data/${base}_R1.fastq.gz 00_raw_data/${base}_R2.fastq.gz -o 01_qc/
done
MultiQC - Aggregate quality report
# Generate combined report
multiqc 01_qc/ -o 01_qc/ -n combined_qc_report
🧹 Step 2: Read Trimming and Filtering
Trimmomatic - Comprehensive read processing
# Basic trimming (remove adapters and low quality)
trimmomatic PE -threads 4 \
    00_raw_data/sample1_R1.fastq.gz 00_raw_data/sample1_R2.fastq.gz \
    02_trimmed/sample1_R1_trimmed.fastq.gz 02_trimmed/sample1_R1_unpaired.fastq.gz \
    02_trimmed/sample1_R2_trimmed.fastq.gz 02_trimmed/sample1_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Advanced trimming with adapter detection
trimmomatic PE -threads 4 \
    00_raw_data/sample1_R1.fastq.gz 00_raw_data/sample1_R2.fastq.gz \
    02_trimmed/sample1_R1_trimmed.fastq.gz 02_trimmed/sample1_R1_unpaired.fastq.gz \
    02_trimmed/sample1_R2_trimmed.fastq.gz 02_trimmed/sample1_R2_unpaired.fastq.gz \
    ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:2:keepBothReads \
    LEADING:3 TRAILING:3 MAXINFO:40:0.8 MINLEN:50
fastp - Modern alternative to Trimmomatic
# Single command trimming with built-in QC
fastp -i 00_raw_data/sample1_R1.fastq.gz -I 00_raw_data/sample1_R2.fastq.gz \
      -o 02_trimmed/sample1_R1_trimmed.fastq.gz -O 02_trimmed/sample1_R2_trimmed.fastq.gz \
      -h 02_trimmed/sample1_fastp.html -j 02_trimmed/sample1_fastp.json \
      --detect_adapter_for_pe --correction --cut_front --cut_tail \
      --thread 4 --qualified_quality_phred 20
🧩 Step 3: Genome Assembly
SPAdes - General purpose assembler
# Standard bacterial assembly
spades.py --pe1-1 02_trimmed/sample1_R1_trimmed.fastq.gz \
          --pe1-2 02_trimmed/sample1_R2_trimmed.fastq.gz \
          -o 03_assembly/sample1_spades \
          --threads 8 --memory 16

# Careful mode (slower but more accurate)
spades.py --pe1-1 02_trimmed/sample1_R1_trimmed.fastq.gz \
          --pe1-2 02_trimmed/sample1_R2_trimmed.fastq.gz \
          -o 03_assembly/sample1_spades_careful \
          --careful --threads 8 --memory 16

# plasmidSPAdes for plasmid assembly
spades.py --pe1-1 02_trimmed/sample1_R1_trimmed.fastq.gz \
          --pe1-2 02_trimmed/sample1_R2_trimmed.fastq.gz \
          -o 03_assembly/sample1_plasmids \
          --plasmid --threads 8 --memory 16
Unicycler - Hybrid assembler (if you have long reads)
# Short reads only
unicycler -1 02_trimmed/sample1_R1_trimmed.fastq.gz \
          -2 02_trimmed/sample1_R2_trimmed.fastq.gz \
          -o 03_assembly/sample1_unicycler \
          --threads 8

# With long reads (Nanopore/PacBio)
unicycler -1 02_trimmed/sample1_R1_trimmed.fastq.gz \
          -2 02_trimmed/sample1_R2_trimmed.fastq.gz \
          -l long_reads/sample1_long.fastq.gz \
          -o 03_assembly/sample1_unicycler_hybrid \
          --threads 8
✅ Step 4: Assembly Quality Assessment
QUAST - Assembly statistics
# Single assembly
quast.py 03_assembly/sample1_spades/contigs.fasta -o 03_assembly/sample1_quast

# Multiple assemblies comparison
quast.py 03_assembly/*/contigs.fasta -o 03_assembly/quast_comparison

# With reference genome
quast.py 03_assembly/sample1_spades/contigs.fasta \
         -r reference_genomes/reference.fasta \
         -g reference_genomes/reference.gff \
         -o 03_assembly/sample1_quast_ref
BUSCO - Genome completeness
# Check for universal bacterial genes
busco -i 03_assembly/sample1_spades/contigs.fasta \
      -l bacteria_odb10 \
      -o sample1_busco \
      -m genome \
      --out_path 03_assembly/
🏷️ Step 5: Genome Annotation
Prokka - Fast prokaryotic annotation
# Basic annotation
prokka --outdir 04_annotation/sample1_prokka \
       --prefix sample1 \
       --genus Escherichia \
       --species coli \
       --strain sample1 \
       --gram neg \
       --cpus 8 \
       03_assembly/sample1_spades/contigs.fasta

# Advanced annotation with custom databases
prokka --outdir 04_annotation/sample1_prokka_advanced \
       --prefix sample1 \
       --genus Escherichia \
       --addgenes --addmrna \
       --rfam --rnammer \
       --cpus 8 \
       03_assembly/sample1_spades/contigs.fasta
🔬 Step 6: Functional and Resistance Analysis
AMRFinderPlus - Antimicrobial resistance detection
# Download database
amrfinder -u

# Run AMR detection
amrfinder --nucleotide 03_assembly/sample1_spades/contigs.fasta \
          --organism Escherichia \
          --threads 8 \
          > 05_analysis/sample1_amr_results.txt
abricate - Multiple resistance database screening
# Screen against multiple databases
abricate --db resfinder 03_assembly/sample1_spades/contigs.fasta > 05_analysis/sample1_resfinder.txt
abricate --db card 03_assembly/sample1_spades/contigs.fasta > 05_analysis/sample1_card.txt
abricate --db vfdb 03_assembly/sample1_spades/contigs.fasta > 05_analysis/sample1_virulence.txt
🌳 Step 7: Comparative Analysis
Snippy - Variant calling and core genome
# Call variants against reference
snippy --cpus 8 --outdir 05_analysis/sample1_snippy \
       --ref reference_genomes/reference.fasta \
       --pe1 02_trimmed/sample1_R1_trimmed.fastq.gz \
       --pe2 02_trimmed/sample1_R2_trimmed.fastq.gz

# Core genome alignment (multiple samples)
snippy-core --ref reference_genomes/reference.fasta \
           --prefix core_alignment \
           05_analysis/*/snippy_results/ \
           > 05_analysis/core_genome.aln
FastANI - Average Nucleotide Identity
# Compare two genomes
fastANI -q 03_assembly/sample1_spades/contigs.fasta \
        -r reference_genomes/reference.fasta \
        -o 05_analysis/sample1_ani.txt

# All-vs-all comparison
fastANI --ql genome_list.txt --rl genome_list.txt -o 05_analysis/all_ani.txt
Roary - Pan-genome analysis
# Create pan-genome from GFF files
roary -p 8 -e -n -v 04_annotation/*/sample*.gff

# Generate core gene alignment
roary -p 8 -e -n -v -a -cd 99 04_annotation/*/sample*.gff
📊 Step 8: Phylogenetic Analysis
IQ-TREE - Maximum likelihood phylogeny
# Build phylogenetic tree from core genome alignment
iqtree -s 05_analysis/core_genome.aln -nt AUTO -bb 1000 -m MFP

# With partition model
iqtree -s 05_analysis/core_genome.aln -spp partitions.txt -bb 1000 -nt AUTO
FastTree - Quick phylogenetic analysis
# Fast approximate ML tree
FastTree -nt -gtr 05_analysis/core_genome.aln > 05_analysis/tree.newick

# With gamma model
FastTree -nt -gtr -gamma 05_analysis/core_genome.aln > 05_analysis/tree_gamma.newick
🔄 Automation Scripts
Batch Processing Script Example
#!/bin/bash
# batch_wgs_analysis.sh

# Sample list
SAMPLES="sample1 sample2 sample3"

for SAMPLE in $SAMPLES; do
    echo "Processing $SAMPLE..."
    
    # Quality control
    fastqc 00_raw_data/${SAMPLE}_R*.fastq.gz -o 01_qc/
    
    # Trimming
    trimmomatic PE -threads 4 \
        00_raw_data/${SAMPLE}_R1.fastq.gz 00_raw_data/${SAMPLE}_R2.fastq.gz \
        02_trimmed/${SAMPLE}_R1_trimmed.fastq.gz 02_trimmed/${SAMPLE}_R1_unpaired.fastq.gz \
        02_trimmed/${SAMPLE}_R2_trimmed.fastq.gz 02_trimmed/${SAMPLE}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    
    # Assembly
    spades.py --pe1-1 02_trimmed/${SAMPLE}_R1_trimmed.fastq.gz \
              --pe1-2 02_trimmed/${SAMPLE}_R2_trimmed.fastq.gz \
              -o 03_assembly/${SAMPLE}_spades \
              --careful --threads 8 --memory 16
    
    # Annotation
    prokka --outdir 04_annotation/${SAMPLE}_prokka \
           --prefix $SAMPLE \
           --cpus 8 \
           03_assembly/${SAMPLE}_spades/contigs.fasta
    
    echo "Completed $SAMPLE"
done

# Generate summary reports
multiqc 01_qc/ -o 06_results/
quast.py 03_assembly/*/contigs.fasta -o 06_results/assembly_comparison
⚡ Performance Tips
Memory and CPU Optimization
SPAdes: Use --memory parameter to limit RAM usage
Trimmomatic: Use -threads parameter for parallel processing
FastQC: Use -t parameter for multiple threads
Prokka: Use --cpus parameter for parallel annotation
Storage Optimization
Compress intermediate files: gzip *.fastq
Remove unnecessary files after each step
Use symbolic links for reference genomes
Archive completed analyses
Troubleshooting Common Issues
Low memory: Reduce SPAdes k-mer sizes or use --careful mode
Long runtimes: Increase CPU cores and check I/O bottlenecks
Poor assemblies: Check read quality and coverage depth
Failed annotations: Verify genome format and gene models
Include metadata (strain, date, source)
Example: Ecoli_K12_MG1655_2024_lab.fastq.gz
