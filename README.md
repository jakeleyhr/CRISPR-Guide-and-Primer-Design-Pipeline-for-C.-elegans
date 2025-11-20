# CRISPR Guide and Primer Design Pipeline for C. elegans

A comprehensive Python-based pipeline for designing sgRNAs and PCR primers for CRISPR-mediated N- and C-terminal protein tagging in *C. elegans*.

## Overview

This pipeline of 7 scripts automates the complete workflow from identification of post-translational protein features to ready-to-order guide and primers sequences:

1. **Identify mature protein sequence** - Identify signal peptides, lipidation sites, etc to identify permissive N- and C-terminal tagging locations for mature proteins
2. **Identify insertion sites** - Convert mature protein termini into genomic coordinates
3. **Design guide RNAs** - Design high-specificity sgRNA sequences with off-target analysis using FlashFry
4. **Extract genomic context** - Extract genomic sequences around insertion sites
5. **Design inner primers** - Design inner homology arm primers with guide/PAM disrupting mutations
6. **Design outer primers** - Design matching outer primers for amplification of homology arms
7. **Design genotyping primers** - Design primers to verify successful insertions

The scripts vary in speed, but all support incremental processing with checkpoint/resume functionality, meaning thousands of genes can easily be gradually processed on personal laptops.

## Requirements

### Core Dependencies
- pandas
- biopython
- pyfaidx
- primer3

### Recommended (for guide off-target detection)
- [FlashFry](https://github.com/mckennalab/FlashFry) - Fast, PAM-aware off-target detection

## Quick Start

### 1. Install dependencies
```bash
pip install pandas biopython pyfaidx primer3
```

### 2. Download genome file
```bash
# Create genome directory
mkdir -p ~/Desktop/wbcel235
cd ~/Desktop/wbcel235

# C. elegans genome (WS235/ce11)
curl -o ce11.fa.gz http://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz
gunzip ce11.fa.gz

cd ..
```

### 3. Set up FlashFry for off-target analysis
```bash
# Download FlashFry jar file
mkdir -p ~/Desktop/flashfry_cel_db
cd ~/Desktop/flashfry_cel_db
curl -L -o FlashFry-assembly-1.15.jar https://github.com/mckennalab/FlashFry/releases/download/1.15/FlashFry-assembly-1.15.jar

java -Xmx8g -jar FlashFry-assembly-1.15.jar \
    index \
    --tmpLocation /tmp \
    --database ce11_spcas9ngg_db \
    --reference ../wbcel235/ce11.fa \
    --enzyme spcas9ngg # change enzyme if needed

cd ..
```

### 4. Run the pipeline
```bash
# Step 1: Identify protein features
python 1_getproteinfeatures.py \
  --input genes_list.csv \
  --output 1.protein_features.csv \
  --flush-every 100

# Step 2: Identify insertion sites
python 2_getinsertionsites.py \
  --input 1.protein_features.csv \
  --output 2.insertion_sites.csv \
  --flush-every 100

# Step 3: Select optimal guides with off-target filtering
python 3_getguidesequences.py \
  --input 2.insertion_sites.csv \
  --output 3.allguideseqs.csv \
  --genome wbcel235/ce11.fa \
  --max-sgrnas 1 \
  --offtarget-mode flashfry \
  --flashfry-db flashfry_cel_db/ce11_spcas9ngg_db \
  --flush-every 100

# Step 4: Extract genomic sequences for primer design
python 4_getgenomicsequencearoundinsertion.py \
  --input 3.allguideseqs.csv \
  --output 4.allguidesandgenomicseqs.csv \
  --genome wbcel235/ce11.fa \
  --flush-every 10000

# Step 5: Design inner homology arm primers with PAM disruption
python 5_designinnerprimers.py \
  --input 4.allguidesandgenomicseqs.csv \
  --output 5.allwormguidesinternalprimers.csv \
  --flush-every 100

# Step 6: Design outer homology arm primers


# Step 7: Design genotyping/inital amplication primers

---

**Note**: This pipeline was developed for *C. elegans* but the core logic should work for other organisms with appropriate genome files. The main species-specific component is the signal peptide prediction in script 1.
