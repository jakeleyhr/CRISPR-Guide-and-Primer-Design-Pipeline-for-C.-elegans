# CRISPR Guide and Primer Design Pipeline for *C. elegans*

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
- BLAT - Fast primer off-target detection

## Quick Start (MacOS)

### 1. Install dependencies
```bash
pip install pandas biopython pyfaidx primer3
```

### 2. Download genome file
```bash
# Create genome directory
mkdir -p ~/Desktop/wbcel235
cd ~/Desktop/wbcel235

# Download C. elegans genome (WS235/ce11)
curl -o ce11.fa.gz http://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz
gunzip ce11.fa.gz

cd ..
```

### 3. Set up FlashFry for guide off-target analysis
```bash
# Download FlashFry jar file
mkdir -p ~/Desktop/flashfry_cel_db
cd ~/Desktop/flashfry_cel_db
curl -L -o FlashFry-assembly-1.15.jar https://github.com/mckennalab/FlashFry/releases/download/1.15/FlashFry-assembly-1.15.jar

Build database (takes a couple of minutes)
java -Xmx8g -jar FlashFry-assembly-1.15.jar \
    index \
    --tmpLocation /tmp \
    --database ce11_spcas9ngg_db \
    --reference ../wbcel235/ce11.fa \
    --enzyme spcas9ngg # change enzyme if needed

cd ..

```
### 4. Download BLAT for primer off-target analysis
```bash
# Download BLAT
mkdir -p ~/bin
cd ~/bin
curl -o blat http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/blat/blat

# Make it executable
chmod +x blat

# Add to PATH permanently (for future terminal windows)
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.zshrc

# Test it works
blat
```

### 5. Run the pipeline
```bash
# Step 1: Identify protein features (~40 genes/min)
python 1_getproteinfeatures.py \
  --input genes_list.csv \
  --output 1.protein_features.csv \
  --flush-every 100

# Step 2: Identify insertion sites (~60 genes/min)
python 2_getinsertionsites.py \
  --input 1.protein_features.csv \
  --output 2.insertion_sites.csv \
  --flush-every 100

# Step 3: Select optimal guides with off-target filtering (~15 sites/min)
python 3_getguidesequences.py \
  --input 2.insertion_sites.csv \
  --output 3.allguideseqs.csv \
  --genome ~/Desktop/wbcel235/ce11.fa \
  --max-sgrnas 1 \
  --offtarget-mode flashfry \
  --flashfry-db ~/Desktop/flashfry_cel_db/ce11_spcas9ngg_db \
  --flashfry-jar ~/Desktop/FlashFry-assembly-1.15.jar \
  --flush-every 100

# Step 4: Extract genomic sequences for primer design (~500 sites/sec)
python 4_getgenomicsequencearoundinsertion.py \
  --input 3.allguideseqs.csv \
  --output 4.allguidesandgenomicseqs.csv \
  --genome ~/Desktop/wbcel235/ce11.fa \
  --flush-every 10000

# Step 5: Design inner homology arm primers with PAM disruption (~300 sites/sec)
python 5_designinnerprimers.py \
  --input 4.allguidesandgenomicseqs.csv \
  --output 5.allwormguidesinternalprimers.csv \
  --flush-every 100

# Step 6: Design outer homology arm primers (6 sites/min)
python 6_designouterprimers.py \
  --input 5.allwormguidesinternalprimers.csv \
  --output 6.allwormguidesbothprimers.csv \
  --genome wbcel235/ce11.fa \
  --flush-every 100

# Step 7: Design genotyping/initial amplification primers
python 7_designgenotypingprimers.py \
  --input 6.allwormguidesbothprimers.csv \
  --output 7.finalfile.csv \
  --genome wbcel235/ce11.fa \
  --flush-every 100
```

## Pipeline Details

### Script 1: Protein Feature Identification
For each gene in the input file, this script queries Ensembl to identify the canonical transcript and its associated UniProt ID. If the canonical transcript has no UniProt ID, the script will search for UniProt IDs associated with other, non-canonical transcripts. As a last resort, the script searches UniProt for any entries associated with the gene. In these cases, the UniProt entry and associated transcript with the longest amino acid sequence is chosen. The script then collects all the relevant protein feature annotations with their amino acid coordinates from the UniProt entry: lipidation sites, initiator methionines, signal peptides, propeptides, transit peptides, transmembrane regions, and chain regions. It also searches for CAAX, PTS1, PTS2, ER KDEL, and ER di-lysine motifs based on sequence.

### Script 2: Insertion Site Identification
This script maps the amino acid coordinates of the start and end of the chain region (the mature protein after e.g. signal peptides have been cleaved off) back to genomic coordinates (the N- and C-terminal insertion sites) via the Ensembl coding sequence information. In cases where an N-terminal or C-terminal residue is lipidated (e.g. N-myristoyl glycine or S-palmitoyl cysteine), the chain region is modified to exclude them, such that the insertion site is "inside" these lipidated residues. The insertion position is defined as the last base before the insertion site, i.e. if the knock-in is desired between bases 100 and 101, the position will be 100 (Note that before/after is from the gene's 5'>3' direction, not necessarily the genomic coordinates.)

### Script 3: Guide RNA Design and Selection
From this point onwards, the remaining scripts run entirely locally. Using a local genome FASTA file, this script extracts 50bp upstream and 50bp downstream of each insertion position, and scans them for suitable guide sequences. Using a defined PAM sequence (NGG by default), it searches for adjacent protospacer sequences (excluding those with extreme GC content, or homopolymer runs of more than 6 bases). Candidates are then screened for off-targets elsewhere in the genome using the PAM-aware FlashFry software (a FlashFry database must be created first). Where possible, only guides with no off-targets are chosen, and the closer to the insertion site, the better. Guides are classified as either suitable for in vivo transcription or in vitro transcription based on the absence/presence of 4+ thymines in a row (in vitro transcription by T7 RNA Polymerase can often read through this stop signal). If an in vitro transcription guide is located closer to the insertion site than the nearest in vivo guide, both will be returned. Otherwise, only the best in vivo guide will be.

### Script 4: Genomic Sequence Extraction
This script extracts 100bp upstream and 100bp downstream of each insertion position, and carefully formats them such that bases on the distal side of the insertion site to the mature protein coding sequence are lower case, and mature protein coding bases are upper case. This formatting is necessary for identifying codon reading frames in the next processing step.

### Script 5: Inner Primer Design with PAM Disruption
Using the local genomic region extracted in the previous step, this script designs Primer 2 and Primer 3 sequences. These are defined as the 30-35bp immediately upstream (Primer 2) and downstream (Primer 3) of the insertion site, using the gene's 5'>3' direction. Critically, silent mutations are introduced at this stage to prevent Cas9 from cutting the vector or the repaired allele. If the guide sequence spans the insertion site, such that fewer than 18 bases of the 23bp guide (protospacer plus PAM) are on one side of the insertion site, no mutations are needed, as these two sides of the guide will be split by the knock-in sequence. However, if the guide is only or primarily on the left or right side of the insertion site, then silent mutations need to be introduced into Primers 2 or 3, respectively. The script recognises the codon reading frame, and attempts to make a silent mutation that disrupts the PAM site. If successful, it introduces an additional 3 silent mutations to the 3' end of the guide sequence. If not, it introduces a total of 5 silent mutations to the 3' end of the guide sequence. All these numbers are adjustable. If the guide sequence is too far away from the insertion site for the 35bp length of the primer to reach it in order to introduce the silent mutations while maintaining 15 3' bases intact to ensure good binding, and additional primer, designated 2A or 3A is designed that overlaps with Primer 2 (2B) or 3 (3B) by 20bp.

### Script 6: Outer Primer Design
This script designs Primers 1 and 4, which pair with Primers 2 and 3, respectively, to amplify the homology arms on either side of the insertion site. It works by extracting genomic sequences and searching for primers that would amplify homology arms of the target size (500-1000bp by default), which have similar annealing temperatures to Primer 2/3. It prefers primers with a 3' GC-clamp, and also performs a BLAT search to ensure that the primers are specific. If no suitable primers can be found, it progressively searches larger and larger regions up to 10kb. In our experience, the length of the homology arms makes little difference to knock-in efficiency, provided they are at least 500bp. Primers with off-targets are only accepted as a last resort, and preferentially those with distant off-targets that are less likely to interfere.

### Script 7: Genotyping Primer Design
Finally, this script designs a set of primers to amplify the whole span of the genome of interest, which can be used for genotyping the successful insertion, and/or for amplifying the region in an initial PCR so the homology arms can be amplified from this fragment rather than the whole genome, where primer overhangs and long sequences could make the reaction harder to optimize. It finds the span from the start of Primer 1 to the end of Primer 4 (the target region) plus an additional flanking 1kb, and searches for suitable and specific (BLAT search) primers that amplify the target region. If none are found, it progressively widens the flanks up to a maximum (e.g. 5kb). If still none are found, the primer properties (length, Tm) are relaxed. Again, primers with off-targets are only accepted as a last resort, and preferentially those with distant off-targets that are less likely to interfere.


### Bonus Script: Visualizing the Inner Primer Design



## Performance Tips

### For Large Gene Lists
All processing scripts support `--flush-every` for incremental writing:
```bash
python 3_getguidesequences.py \
  --input-sites sites.csv \
  --output guides.csv \
  --genome genome.fa \
  --flush-every 100  # Write every 100 rows
```

### Resuming Interrupted Runs
By default, scripts resume from where they left off:
```bash
# Run interrupted at gene 500/1000
python 3_getguidesequences.py \
  --input-sites sites.csv \
  --output guides.csv \
  --genome genome.fa
# Automatically resumes from gene 501

# Force restart from the beginning:
python 3_getguidesequences.py \
  --input-sites sites.csv \
  --output guides.csv \
  --genome genome.fa \
  --no-resume
```

## Citation
If you use this pipeline in your research, please cite:
- xyz


**Note**: This pipeline was developed for *C. elegans* but the core logic should work for other organisms with appropriate genome files. The main species-specific component is the signal peptide prediction in script 1, which can be easily modified. Similarly, if you want to use your own list of insertion positions, you can skip to Step 3 and the pipeline should still work as intended, provided the insertion positions are all in-frame, otherwise scripts 4 and 5 will not correctly identify the reading frames for the silent mutations.
