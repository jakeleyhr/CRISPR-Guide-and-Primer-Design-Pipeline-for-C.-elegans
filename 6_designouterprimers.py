#!/usr/bin/env python3
"""
Homology Arm Primer Design Pipeline - Step 3: Outer Primer Design

Designs outer PCR primers (P1 and P4) that pair with inner primers to amplify
complete homology arms for CRISPR knock-in repair templates. Uses Primer3 for
optimal primer design with automatic off-target checking via BLAT.

Input:
  - CSV file from script 5 containing inner primers (P2 and P3) with genomic
    coordinates and sequences
    (requires: chrom, pos, strand, primer_2_seq_genomic, primer_2_locus,
     primer_3_seq_genomic, primer_3_locus)

Output:
  - CSV file with all input columns plus outer primer information:
      • primer_1_seq: Upstream outer primer (pairs with P2)
      • primer_4_seq: Downstream outer primer (pairs with P4)
      • primer_1_locus, primer_4_locus: Genomic coordinates
      • primer_1_tm, primer_4_tm: Melting temperatures
      • primer_1_gc, primer_4_gc: GC content percentages
      • upstream_arm_length, downstream_arm_length: Final amplicon sizes
      • Off-target descriptions

Primer Nomenclature (Gene-Centric):
  - P1 + P2 amplify gene UPSTREAM homology arm (500-1000bp target)
  - P3 + P4 amplify gene DOWNSTREAM homology arm (500-1000bp target)
  
  For minus strand genes: gene upstream = genomic RIGHT side
  (Primer names always follow gene orientation, not genomic coordinates)

Design Strategy:
  1. Extract flanking sequences (up to e.g. 1kb each side of insertion)
  2. Locate inner primers (P2, P3) within extracted sequences
  3. Design outer primers (P1, P4) using Primer3:
     - Target product size: 500-1000bp initially
     - Optimal Tm and length tries to match P1/P2
     - Prefentially chooses primer sequences with 3' GC-clamp
  4. Verify primer uniqueness with BLAT (optional but recommended)
  5. Annotate off-target hits and categorize by location
  6. Progressively falls back to target homology arms of up to 10kb if needed

Off-Target Detection (BLAT):
  - Checks for primer binding elsewhere in genome
  - Categories:
    * same_chrom_near: Within 10kb of expected location (bad)
    * same_chrom_far: Same chromosome, >10kb away (fine if an initial 
     amplifcation of region is performed prior to using these primers)
    * diff_chrom: Different chromosome
  - Falls back to best available primer if all have off-targets
  - Can be disabled with --skip-offtarget-check for speed

Required arguments:
  --input PATH       Input CSV from script 5 with inner primers
  --output PATH      Output CSV path with outer primers
  --genome PATH      Genome FASTA file

Optional arguments:
  --flush-every N              Write results every N rows for progress tracking (default: 10)
  --no-resume                  Overwrite output instead of resuming from checkpoint
  --skip-offtarget-check       Skip off-target search using BLAT (much faster but less safe)
  --min-arm N                  Minimum homology arm length in bp (default: 500)
  --max-arm N                  Maximum homology arm length in bp (default: 1000)
  --primer-min-size N          Minimum primer length in bp (default: 22)
  --primer-max-size N          Maximum primer length in bp (default: 36)
  --primer-min-tm N            Minimum primer Tm in °C (default: 50.0)
  --primer-max-tm N            Maximum primer Tm in °C (default: 72.0)
  --blat-min-identity N        Minimum % identity for BLAT hits (default: 90)
  --blat-min-score N           Minimum alignment score for BLAT (default: 20)

Example:
  python 6_designouterprimers.py \
    --input 5.allwormguidesinternalprimers.csv \
    --output 6.allwormguidesbothprimers.csv \
    --genome wbcel235/ce11.fa \
    --flush-every 100

Recommended: install BLAT (for off-target checking)
  - mkdir -p ~/bin # Create bin directory if it doesn't exist
  - cd ~/bin
  - curl -o blat http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/blat/blat # Download BLAT
  - chmod +x blat # Make it executable
  - echo 'export PATH="$HOME/bin:$PATH"' >> ~/.zshrc # Add to PATH permanently
  - blat # Test
"""


import argparse, sys, os, subprocess, tempfile, traceback
import pandas as pd
from collections import defaultdict

try:
    import primer3
except ImportError:
    print("ERROR: primer3-py not installed! Install with: pip install primer3-py")
    sys.exit(1)

try:
    from Bio import SeqIO
except ImportError:
    print("ERROR: biopython not installed! Install with: pip install biopython")
    sys.exit(1)

def load_genome(fasta_path):
    """Load genome sequences from FASTA file"""
    print(f"Loading genome from {fasta_path}...")
    genome = {}
    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            genome[record.id] = str(record.seq).upper()
        print(f"  Loaded {len(genome)} sequences")
        return genome
    except Exception as e:
        print(f"ERROR loading FASTA: {e}")
        sys.exit(1)

def extract_flanks(genome, chrom, position, strand, flank, max_arm):
    """
    Extract upstream and downstream flanking sequences based on gene strand
    
    For PLUS strand:
      - Upstream (P1+P2): pos - flank to pos (then RC so P2 is at start)
      - Downstream (P3+P4): pos to pos + flank (P3 at start)
    
    For MINUS strand:
      - Upstream (P1+P2): pos to pos + flank (P2 at start)
      - Downstream (P3+P4): pos - flank to pos (then RC so P3 is at start)
    
    Returns: (upstream_seq, downstream_seq, upstream_start, downstream_start)
    """
    if chrom not in genome:
        return None, None, None, None
    
    seq = genome[chrom]
    
    if strand == '+':
        # Upstream (gene left): extract left, then RC
        upstream_start = max(0, position - flank)
        upstream_seq_forward = seq[upstream_start:position]
        upstream_seq = rc(upstream_seq_forward)
        
        # Downstream (gene right): extract right
        downstream_start = position + 1
        downstream_end = min(len(seq), position + flank)
        downstream_seq = seq[downstream_start:downstream_end]
        
        print(f"  DEBUG extract_flanks (PLUS):")
        print(f"    Upstream: genomic {upstream_start+1}-{position} (RC'd), seq len={len(upstream_seq)}")
        print(f"    Downstream: genomic {downstream_start+1}-{downstream_end} (forward), seq len={len(downstream_seq)}")
        print(f"    Upstream first 50bp: {upstream_seq[:50]}")
        print(f"    Downstream first 50bp: {downstream_seq[:50]}")
        
    else:  # strand == '-'
        # Upstream (gene left, genomic right): extract right
        upstream_start = position
        upstream_end = min(len(seq), position + flank)
        upstream_seq = seq[upstream_start:upstream_end]
        
        # Downstream (gene right, genomic left): extract left, then RC
        downstream_start = max(0, position - flank)
        downstream_seq_forward = seq[downstream_start:position-1]
        downstream_seq = rc(downstream_seq_forward)
        
        print(f"  DEBUG extract_flanks (MINUS):")
        print(f"    Upstream: genomic {upstream_start+1}-{upstream_end} (forward), seq len={len(upstream_seq)}")
        print(f"    Downstream: genomic {downstream_start+1}-{position} (RC'd), seq len={len(downstream_seq)}")
        print(f"    Upstream first 50bp: {upstream_seq[:50]}")
        print(f"    Downstream first 50bp: {downstream_seq[:50]}")
    
    return upstream_seq, downstream_seq, upstream_start, downstream_start

def rc(seq):
    """Reverse complement"""
    comp = {'A':'T', 'T':'A', 'G':'C', 'C':'G',
            'a':'t', 't':'a', 'g':'c', 'c':'g',
            'N':'N', 'n':'n'}
    return ''.join([comp.get(b, b) for b in seq[::-1]])

def find_primer_flexible(primer, seq, name="primer"):
    """
    Find primer in sequence, trying multiple strategies:
    1. Exact match (forward)
    2. Exact match (RC)
    3. Fuzzy match (forward, up to 5 mismatches)
    4. Fuzzy match (RC, up to 5 mismatches)
    5. Truncated match at boundaries (handles primers extending beyond extraction)
    
    Returns: (position, length, orientation) or (None, None, None)
    """
    # Try exact forward
    pos = seq.find(primer)
    if pos != -1:
        return pos, len(primer), 'forward'
    
    # Try exact RC
    primer_rc = rc(primer)
    pos = seq.find(primer_rc)
    if pos != -1:
        return pos, len(primer_rc), 'RC'
    
    # Try fuzzy forward
    pos, mismatches = fuzzy_find(primer, seq, max_mismatches=5)
    if pos is not None:
        print(f"    Note: {name} found with {mismatches} mismatches (forward)")
        return pos, len(primer), 'forward_fuzzy'
    
    # Try fuzzy RC
    pos, mismatches = fuzzy_find(primer_rc, seq, max_mismatches=5)
    if pos is not None:
        print(f"    Note: {name} found with {mismatches} mismatches (RC)")
        return pos, len(primer_rc), 'RC_fuzzy'
    
    # Try truncated at end (forward) - handles primers extending beyond boundary
    for truncate in range(1, min(6, len(primer))):
        truncated = primer[:-truncate]
        if len(seq) >= len(truncated) and seq[-len(truncated):] == truncated:
            print(f"    Note: {name} found truncated by {truncate}bp at sequence end (forward)")
            return len(seq) - len(truncated), len(truncated), 'forward_truncated'
    
    # Try truncated at end (RC)
    for truncate in range(1, min(6, len(primer_rc))):
        truncated = primer_rc[:-truncate]
        if len(seq) >= len(truncated) and seq[-len(truncated):] == truncated:
            print(f"    Note: {name} found truncated by {truncate}bp at sequence end (RC)")
            return len(seq) - len(truncated), len(truncated), 'RC_truncated'
    
    # Try truncated at start (forward)
    for truncate in range(1, min(6, len(primer))):
        truncated = primer[truncate:]
        if seq[:len(truncated)] == truncated:
            print(f"    Note: {name} found truncated by {truncate}bp at sequence start (forward)")
            return 0, len(truncated), 'forward_truncated'
    
    # Try truncated at start (RC)
    for truncate in range(1, min(6, len(primer_rc))):
        truncated = primer_rc[truncate:]
        if seq[:len(truncated)] == truncated:
            print(f"    Note: {name} found truncated by {truncate}bp at sequence start (RC)")
            return 0, len(truncated), 'RC_truncated'
    
    return None, None, None

def fuzzy_find(primer, seq, max_mismatches=5):
    """Find primer allowing up to max_mismatches mismatches"""
    best_pos = None
    best_mismatches = max_mismatches + 1
    
    for i in range(len(seq) - len(primer) + 1):
        mismatches = sum(1 for j in range(len(primer)) if seq[i+j] != primer[j])
        if mismatches < best_mismatches:
            best_mismatches = mismatches
            best_pos = i
            if mismatches == 0:
                break
    
    if best_mismatches <= max_mismatches:
        return best_pos, best_mismatches
    return None, None

def parse_blat_psl(psl_path):
    """Parse BLAT PSL - simplified version"""
    hits = defaultdict(list)
    if not os.path.exists(psl_path):
        return hits
    
    with open(psl_path, 'r') as f:
        for line in f:
            if line.startswith('psLayout') or line.startswith('match') or line.startswith('---'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 17:
                continue
            try:
                matches = int(parts[0])
                mismatches = int(parts[1])
                qName = parts[9]
                tName = parts[13]
                tStart = int(parts[15])
                tEnd = int(parts[16])
                
                if matches + mismatches == 0:
                    continue
                identity = 100.0 * matches / (matches + mismatches)
                
                hits[qName].append({
                    'chrom': tName,
                    'start': tStart + 1,
                    'end': tEnd,
                    'identity': identity,
                    'score': matches
                })
            except:
                continue
    return hits

def check_offtargets_categorized(hits, primer_id, expected_chrom, expected_start, expected_end, cut_position, tolerance=5):
    """Categorize off-targets by location"""
    if primer_id not in hits:
        return 'clean', 0, "No hits"
    
    good_hits = [h for h in hits[primer_id] if h['identity'] >= 90 and h['score'] >= 18]
    if not good_hits:
        return 'clean', 0, "No quality hits"
    
    offtargets = []
    for hit in good_hits:
        is_expected = (hit['chrom'] == expected_chrom and 
                      abs(hit['start'] - expected_start) <= tolerance and 
                      abs(hit['end'] - expected_end) <= tolerance)
        if not is_expected:
            offtargets.append(hit)
    
    if not offtargets:
        return 'clean', 0, "Only on-target"
    
    # Categorize
    other_chrom = [h for h in offtargets if h['chrom'] != expected_chrom]
    if other_chrom:
        return 'other_chrom', len(other_chrom), f"{len(other_chrom)} hits on other chromosomes"
    
    same_chrom_far = [h for h in offtargets if abs(h['start'] - cut_position) > 10000]
    if same_chrom_far:
        return 'same_chrom', len(same_chrom_far), f"{len(same_chrom_far)} hits >10kb away"
    
    return 'near_position', len(offtargets), f"{len(offtargets)} hits within 10kb"

# prefer primers that end in G or C at the 3' end
def _gc_end_pref(seq):
    return 0 if seq[-1].upper() in ("G","C") else 1

def design_outer_primer_for_p1p2_arm(seq, p2_primer, seq_start, strand, cut_position, args, num_return=50):
    """
    Design P1 candidates to pair with P2 for gene upstream (P1+P2) homology arm
    Returns list of primer candidates, best first
    """
    print(f"  Searching for P2...")
    
    # Calculate P2 Tm and constraints
    p2_tm_result = primer3.calc_tm(p2_primer)
    p2_len = len(p2_primer)
    print(f"    P2 length: {p2_len}bp, Tm: {p2_tm_result:.1f}°C")
    
    # Constrain P1 size to match P2 within allowed range
    p1_min_size = args.primer_min_size
    p1_max_size = args.primer_max_size
    p1_opt_size = max(args.primer_min_size, min(args.primer_max_size, p2_len))
    
    # Constrain P1 Tm to match P2, with clamping to allowed range
    p1_opt_tm = p2_tm_result
    p1_min_tm = args.primer_min_tm
    p1_max_tm = args.primer_max_tm
    
    # Clamp optimal Tm to be within allowed range
    if p2_tm_result < args.primer_min_tm:
        p1_opt_tm = args.primer_min_tm
    elif p2_tm_result > args.primer_max_tm:
        p1_opt_tm = args.primer_max_tm
    
    # Ensure min < opt < max (primer3 requirement)
    if p1_min_tm >= p1_opt_tm:
        p1_min_tm = p1_opt_tm - 1.0
    if p1_max_tm <= p1_opt_tm:
        p1_max_tm = p1_opt_tm + 1.0
    
    print(f"    P1 constraints: {p1_min_size}-{p1_max_size}bp, Tm: {p1_min_tm:.1f}-{p1_max_tm:.1f}°C (opt: {p1_opt_tm:.1f})")
    
    # Find P2 in sequence
    p2_pos, p2_len_found, p2_orientation = find_primer_flexible(p2_primer, seq, "P2")
    
    if p2_pos is None:
        print(f"    ✗ P2 not found in sequence")
        return []
    
    p2_end = p2_pos + p2_len_found
    print(f"    ✓ P2 found at {p2_pos}-{p2_end} ({p2_orientation})")
    
    # Design P1 candidates
    try:
        search_end = min(len(seq), p2_pos + args.max_arm)
        search_start = search_end - 500
        
        p3_args_custom = {
            'PRIMER_MIN_SIZE': p1_min_size,
            'PRIMER_MAX_SIZE': p1_max_size,
            'PRIMER_OPT_SIZE': p1_opt_size,
            'PRIMER_MIN_TM': p1_min_tm,
            'PRIMER_MAX_TM': p1_max_tm,
            'PRIMER_OPT_TM': p1_opt_tm,
            'PRIMER_MIN_GC': 40.0,
            'PRIMER_MAX_GC': 60.0,
            'PRIMER_PRODUCT_SIZE_RANGE': [[args.min_arm, args.max_arm]],
            'PRIMER_NUM_RETURN': num_return,
        }
        
        res = primer3.bindings.design_primers({
            'SEQUENCE_TEMPLATE': seq,
            'SEQUENCE_INCLUDED_REGION': [search_start, search_end - search_start],
            'PRIMER_PICK_LEFT_PRIMER': 0,
            'PRIMER_PICK_RIGHT_PRIMER': 1,
            'PRIMER_PICK_INTERNAL_OLIGO': 0,
        }, p3_args_custom)
        
        num_returned = res.get('PRIMER_RIGHT_NUM_RETURNED', 0)
        if num_returned == 0:
            return []
        
        # Extract all candidates
        candidates = []
        for i in range(num_returned):
            try:
                seq_p1 = res[f'PRIMER_RIGHT_{i}_SEQUENCE']
                start = res[f'PRIMER_RIGHT_{i}'][0]
                length = res[f'PRIMER_RIGHT_{i}'][1]
                tm = res[f'PRIMER_RIGHT_{i}_TM']
                gc = res[f'PRIMER_RIGHT_{i}_GC_PERCENT']
                
                # Calculate genomic coordinates
                if strand == '+':
                    p1_start_genomic = cut_position - start
                    p1_end_genomic = cut_position - (start - length + 1)
                else:
                    p1_start_genomic = seq_start + (start - length + 1) + 1
                    p1_end_genomic = seq_start + start + 1
                
                candidates.append({
                    'primer_1_seq': seq_p1,
                    'primer_1_start': p1_start_genomic,
                    'primer_1_end': p1_end_genomic,
                    'primer_1_tm': round(tm, 1),
                    'primer_1_gc': round(gc, 1),
                    'primer_1_notes': ''
                })
            except KeyError:
                continue
        
        return candidates
        
    except Exception as e:
        return {
            'primer_1_seq': 'ERROR',
            'primer_1_start': None,
            'primer_1_end': None,
            'primer_1_tm': None,
            'primer_1_gc': None,
            'primer_1_notes': f'Exception: {str(e)}'
        }

def design_outer_primer_for_p3p4_arm(seq, p3_primer, seq_start, strand, cut_position, args, num_return=50):
    """
    Design P4 candidates to pair with P3 for gene downstream (P3+P4) homology arm
    Constrains P4 to match P3's length and Tm within allowed ranges
    Returns list of primer candidates, best first
    """
    print(f"  Searching for P3...")
    
    # Calculate P3's Tm using Primer3
    p3_tm_result = primer3.calc_tm(p3_primer)
    p3_len = len(p3_primer)
    print(f"    P3 length: {p3_len}bp, Tm: {p3_tm_result:.1f}°C")
    
    # Constrain P4 size to match P3 within allowed range
    p4_min_size = args.primer_min_size
    p4_max_size = args.primer_max_size
    p4_opt_size = max(args.primer_min_size, min(args.primer_max_size, p3_len))
    
    # Constrain P4 Tm to match P3, with clamping to allowed range
    p4_opt_tm = p3_tm_result
    p4_min_tm = args.primer_min_tm
    p4_max_tm = args.primer_max_tm
    
    # Clamp optimal Tm to be within allowed range
    if p3_tm_result < args.primer_min_tm:
        p4_opt_tm = args.primer_min_tm
    elif p3_tm_result > args.primer_max_tm:
        p4_opt_tm = args.primer_max_tm
    
    # Ensure min < opt < max (primer3 requirement)
    if p4_min_tm >= p4_opt_tm:
        p4_min_tm = p4_opt_tm - 1.0
    if p4_max_tm <= p4_opt_tm:
        p4_max_tm = p4_opt_tm + 1.0
    
    print(f"    P4 constraints: {p4_min_size}-{p4_max_size}bp, Tm: {p4_min_tm:.1f}-{p4_max_tm:.1f}°C (opt: {p4_opt_tm:.1f})")
    
    # Find P3 in sequence
    p3_pos, p3_len_found, p3_orientation = find_primer_flexible(p3_primer, seq, "P3")
    if p3_pos is None:
        print(f"    ✗ P3 not found in sequence (length: {len(seq)}bp)")
        print(f"    P3 sequence: {p3_primer}")
        print(f"    First 100bp of seq: {seq[:100]}")
        return []
    
    p3_end = p3_pos + p3_len_found
    print(f"    ✓ P3 found at {p3_pos}-{p3_end} ({p3_orientation})")
    
    # Design P4 to pair with P3
    try:
        search_end = min(len(seq), p3_pos + args.max_arm)
        search_start = search_end - 500
        
        if search_end <= search_start:
            return []
        
        print(f"    Searching for P4 in region {search_start}-{search_end} (product target: {args.min_arm}-{args.max_arm}bp)")
        
        p3_args_custom = {
            'PRIMER_MIN_SIZE': p4_min_size,
            'PRIMER_MAX_SIZE': p4_max_size,
            'PRIMER_OPT_SIZE': p4_opt_size,
            'PRIMER_MIN_TM': p4_min_tm,
            'PRIMER_MAX_TM': p4_max_tm,
            'PRIMER_OPT_TM': p4_opt_tm,
            'PRIMER_MIN_GC': 40.0,
            'PRIMER_MAX_GC': 60.0,
            'PRIMER_PRODUCT_SIZE_RANGE': [[args.min_arm, args.max_arm]],
            'PRIMER_NUM_RETURN': num_return,
        }
        
        res = primer3.bindings.design_primers({
            'SEQUENCE_TEMPLATE': seq,
            'SEQUENCE_INCLUDED_REGION': [search_start, search_end - search_start],
            'PRIMER_PICK_LEFT_PRIMER': 0,
            'PRIMER_PICK_RIGHT_PRIMER': 1,
            'PRIMER_PICK_INTERNAL_OLIGO': 0,
        }, p3_args_custom)

        num_returned = res.get('PRIMER_RIGHT_NUM_RETURNED', 0)
        if num_returned == 0:
            return []
        
        # Extract all candidates
        candidates = []
        for i in range(num_returned):
            try:
                seq_p4 = res[f'PRIMER_RIGHT_{i}_SEQUENCE']
                start = res[f'PRIMER_RIGHT_{i}'][0]
                length = res[f'PRIMER_RIGHT_{i}'][1]
                tm = res[f'PRIMER_RIGHT_{i}_TM']
                gc = res[f'PRIMER_RIGHT_{i}_GC_PERCENT']

                # Calculate genomic coordinates
                if strand == '+':
                    p4_start_genomic = seq_start + (start - length + 1) + 1
                    p4_end_genomic = seq_start + start + 1
                else:
                    p4_end_genomic = cut_position - (start - length + 1)
                    p4_start_genomic = cut_position - start
                
                candidates.append({
                    'primer_4_seq': seq_p4,
                    'primer_4_start': p4_start_genomic,
                    'primer_4_end': p4_end_genomic,
                    'primer_4_tm': round(tm, 1),
                    'primer_4_gc': round(gc, 1),
                    'primer_4_notes': ''
                })
            except KeyError:
                continue
        
        return candidates
        
    except Exception as e:
        return {
            'primer_4_seq': 'ERROR',
            'primer_4_start': None,
            'primer_4_end': None,
            'primer_4_tm': None,
            'primer_4_gc': None,
            'primer_4_notes': f'Exception: {str(e)}'
        }

def parse_locus_coords(locus_str):
    """
    Parse locus string like 'I:5107814-5107843 (-)' to get start and end coordinates.
    Returns (start, end) as integers, or (None, None) if parsing fails.
    """
    try:
        if not locus_str or pd.isna(locus_str):
            return None, None
        parts = locus_str.split(':')
        if len(parts) < 2:
            return None, None
        coords = parts[1].split(' ')[0]
        start_str, end_str = coords.split('-')
        return int(start_str), int(end_str)
    except:
        return None, None

def process(row, n, genome, args, genome_fasta_path=None):
    """
    Process one row - design P1 and P4
    
    Primer names follow GENE coordinates:
    - P1+P2 = gene upstream arm
    - P3+P4 = gene downstream arm
    
    For minus strand: gene upstream = genomic RIGHT
    """
    print(f"\n{'='*70}")
    print(f"Row {n}: {row.get('gene', '?')} ({row.get('site_type', '?')})")
    print(f"{'='*70}")
    
    # Get fields
    chrom = row.get('chrom', '')
    position = row.get('pos', None)
    p2 = row.get('primer_2_seq_genomic', '').upper()
    p3 = row.get('primer_3_seq_genomic', '').upper()
    strand = row.get('strand', '+')
    
    # Get P2 and P3 coordinates from locus strings
    p2_locus = row.get('primer_2_locus', '')
    p3_locus = row.get('primer_3_locus', '')
    p2_start, p2_end = parse_locus_coords(p2_locus)
    p3_start, p3_end = parse_locus_coords(p3_locus)

    print(f"  P2 from CSV: {p2}")
    print(f"  P2 locus: {p2_locus} → coords: {p2_start}-{p2_end}")
    print(f"  P3 from CSV: {p3}")
    print(f"  P3 locus: {p3_locus} → coords: {p3_start}-{p3_end}")
    
    # Validation
    if not chrom or position is None or not p2 or not p3:
        missing = []
        if not chrom: missing.append('chrom')
        if position is None: missing.append('pos')
        if not p2: missing.append('primer_2_seq_genomic')
        if not p3: missing.append('primer_3_seq_genomic')
        err = f"Missing: {', '.join(missing)}"
        return {
            'primer_1_seq': 'ERROR', 'primer_1_locus': None,
            'primer_1_tm': None, 'primer_1_gc': None, 'upstream_arm_length': None,
            'primer_1_notes': err,
            'primer_4_seq': 'ERROR', 'primer_4_locus': None,
            'primer_4_tm': None, 'primer_4_gc': None, 'downstream_arm_length': None,
            'primer_4_notes': err
        }
    
    try:
        position = int(position)
    except:
        err = f'Invalid position: {position}'
        return {
            'primer_1_seq': 'ERROR', 'primer_1_notes': err,
            'primer_4_seq': 'ERROR', 'primer_4_notes': err
        }
    
    # Extract genomic flanks
    upstream_seq, downstream_seq, upstream_start, downstream_start = extract_flanks(
        genome, chrom, position, strand, args.max_arm, args.max_arm
    )
    if upstream_seq is None:
        err = f'Chromosome {chrom} not found'
        return {
            'primer_1_seq': 'ERROR', 'primer_1_notes': err,
            'primer_4_seq': 'ERROR', 'primer_4_notes': err
        }
    
    print(f"Chr {chrom}:{position} (strand {strand})")
    print(f"Upstream arm (P1+P2): {len(upstream_seq)}bp")
    print(f"Downstream arm (P3+P4): {len(downstream_seq)}bp")
    
    # P1+P2 always use upstream, P3+P4 always use downstream
    p1p2_seq, p1p2_start = upstream_seq, upstream_start
    p3p4_seq, p3p4_start = downstream_seq, downstream_start
    
    # Design primers with retry logic and off-target checking
    print(f"\nDesigning P1 (pairs with P2):")
    r1 = None
    best_fallback_diff_chrom = None
    best_fallback_same_chrom = None
    best_fallback_near = None
    total_p1_tested = 0
    total_p1_clean = 0
    total_p1_offtargets = 0
    
    # Try all arm sizes from max_arm to max_arm_fallback
    for attempt_arm_size in range(args.max_arm, args.max_arm_fallback + 1, 500):
        if attempt_arm_size > args.max_arm:
            print(f"  Retry with MAX_ARM={attempt_arm_size}bp...")
            flanks = extract_flanks(genome, chrom, position, strand, attempt_arm_size + 100, attempt_arm_size)
            p1p2_seq, p1p2_start = flanks[0], flanks[2]  # upstream

        candidates = design_outer_primer_for_p1p2_arm(p1p2_seq, p2, p1p2_start, strand, position, args, num_return=50)
        
        # Sort by GC end preference, then score
        for c in candidates:
            c['_sort_score'] = c.get('score', 0)
        candidates.sort(key=lambda c: (_gc_end_pref(c['primer_1_seq']), -c['_sort_score']))
        for c in candidates:
            c.pop('_sort_score', None)
        
        if not candidates:
            continue
        
        print(f"    Testing {len(candidates)} P1 candidates...")
        
        if genome_fasta_path:
            # Batch BLAT search for all candidates
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
                batch_fasta = f.name
                for idx, candidate in enumerate(candidates):
                    f.write(f">P1_row{n}_cand{idx}\n{candidate['primer_1_seq']}\n")
            
            with tempfile.NamedTemporaryFile(mode='w', suffix='.psl', delete=False) as f:
                batch_psl = f.name
            
            # Run BLAT once for all candidates
            cmd = ['blat', genome_fasta_path, batch_fasta, batch_psl, '-stepSize=5', 
                   '-repMatch=1000000', f'-minScore={args.blat_min_score}', 
                   f'-minIdentity={args.blat_min_identity}', '-minMatch=1']
            result = subprocess.run(cmd, capture_output=True, timeout=300)
            
            if result.returncode == 0:
                all_hits = parse_blat_psl(batch_psl)
                
                for idx, candidate in enumerate(candidates):
                    total_p1_tested += 1
                    primer_id = f"P1_row{n}_cand{idx}"
                    
                    offtarget_hits = []
                    if primer_id in all_hits:
                        for hit in all_hits[primer_id][:5]:
                            if hit['identity'] >= args.blat_min_identity and hit['score'] >= args.blat_min_score:
                                is_expected = (hit['chrom'] == chrom and 
                                             abs(hit['start'] - candidate['primer_1_start']) <= 5 and 
                                             abs(hit['end'] - candidate['primer_1_end']) <= 5)
                                if not is_expected:
                                    offtarget_hits.append(hit)
                    
                    if not offtarget_hits:
                        total_p1_clean += 1
                        r1 = candidate
                        r1['primer_1_notes'] = 'No off-targets'
                        print(f"  ✓ Clean P1 found at MAX_ARM={attempt_arm_size}bp (candidate #{idx+1}, {total_p1_clean} clean / {total_p1_tested} tested)")
                        break
                    else:
                        total_p1_offtargets += 1
                        
                        # Categorize off-targets
                        other_chrom = [h for h in offtarget_hits if h['chrom'] != chrom]
                        same_chrom_far = [h for h in offtarget_hits if h['chrom'] == chrom and abs(h['start'] - position) > 10000]
                        same_chrom_near = [h for h in offtarget_hits if h['chrom'] == chrom and abs(h['start'] - position) <= 10000]
                        
                        # Categorize by WORST type present
                        if same_chrom_near:
                            category = 'same_chrom_near'
                        elif same_chrom_far:
                            category = 'same_chrom_far'
                        elif other_chrom:
                            category = 'diff_chrom'
                        else:
                            category = 'unknown'
                        
                        print(f"    Candidate #{idx+1} ({category}): {candidate['primer_1_seq']}")
                        
                        # Save fallbacks in priority order: diff_chrom > same_chrom_far > same_chrom_near
                        if category == 'diff_chrom' and not best_fallback_diff_chrom:
                            best_fallback_diff_chrom = candidate.copy()
                            best_fallback_diff_chrom['primer_1_notes'] = f'Off-targets on different chromosome(s): {len(other_chrom)} hit(s)'
                        elif category == 'same_chrom_far' and not best_fallback_same_chrom:
                            best_fallback_same_chrom = candidate.copy()
                            best_fallback_same_chrom['primer_1_notes'] = f'Off-targets on same chromosome >10kb away: {len(same_chrom_far)} hit(s)'
                        elif category == 'same_chrom_near' and not best_fallback_near:
                            best_fallback_near = candidate.copy()
                            best_fallback_near['primer_1_notes'] = f'Off-targets within 10kb: {len(same_chrom_near)} hit(s)'
            
            os.unlink(batch_fasta)
            os.unlink(batch_psl)
        else:
            r1 = candidates[0]
            r1['primer_1_notes'] = 'No off-target check performed'
            break
        
        # If we found a clean primer, stop searching
        if r1:
            break
    
    # If no clean primer found after trying all sizes, use best fallback
    if not r1:
        if best_fallback_diff_chrom:
            r1 = best_fallback_diff_chrom
            print(f"  ! Using fallback P1 (different chromosome off-targets) - tested {total_p1_tested} candidates, all had off-targets")
        elif best_fallback_same_chrom:
            r1 = best_fallback_same_chrom
            print(f"  ! Using fallback P1 (same chrom >10kb off-targets) - tested {total_p1_tested} candidates, all had off-targets")
        elif best_fallback_near:
            r1 = best_fallback_near
            print(f"  ! Using fallback P1 (same chrom <10kb off-targets - LAST RESORT) - tested {total_p1_tested} candidates, all had off-targets")
        else:
            r1 = {
                'primer_1_seq': 'ERROR',
                'primer_1_start': None,
                'primer_1_end': None,
                'primer_1_locus': None,
                'primer_1_tm': None,
                'primer_1_gc': None,
                'upstream_arm_length': None,
                'primer_1_notes': f'No suitable primer found (tested {total_p1_tested} candidates, all had off-targets)'
            }
    
    print(f"\nDesigning P4 (pairs with P3):")
    r4 = None
    best_fallback_diff_chrom = None
    best_fallback_same_chrom = None
    best_fallback_near = None
    total_p4_tested = 0
    total_p4_clean = 0
    total_p4_offtargets = 0
    
    # Try all arm sizes from max_arm to max_arm_fallback
    for attempt_arm_size in range(args.max_arm, args.max_arm_fallback + 1, 500):
        if attempt_arm_size > args.max_arm:
            print(f"  Retry with MAX_ARM={attempt_arm_size}bp...")
            flanks = extract_flanks(genome, chrom, position, strand, attempt_arm_size + 100, attempt_arm_size)
            p3p4_seq, p3p4_start = flanks[1], flanks[3]  # Downstream
        
        candidates = design_outer_primer_for_p3p4_arm(p3p4_seq, p3, p3p4_start, strand, position, args, num_return=50)
        
        for c in candidates:
            c['_sort_score'] = c.get('score', 0)
        candidates.sort(key=lambda c: (_gc_end_pref(c['primer_4_seq']), -c['_sort_score']))
        for c in candidates:
            c.pop('_sort_score', None)
        
        if not candidates:
            continue
        
        print(f"    Testing {len(candidates)} P4 candidates...")
        
        if genome_fasta_path:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
                batch_fasta = f.name
                for idx, candidate in enumerate(candidates):
                    f.write(f">P4_row{n}_cand{idx}\n{candidate['primer_4_seq']}\n")
            
            with tempfile.NamedTemporaryFile(mode='w', suffix='.psl', delete=False) as f:
                batch_psl = f.name
            
            cmd = ['blat', genome_fasta_path, batch_fasta, batch_psl, '-stepSize=5', 
                   '-repMatch=1000000', f'-minScore={args.blat_min_score}', 
                   f'-minIdentity={args.blat_min_identity}', '-minMatch=1']
            result = subprocess.run(cmd, capture_output=True, timeout=300)
            
            if result.returncode == 0:
                all_hits = parse_blat_psl(batch_psl)
                
                for idx, candidate in enumerate(candidates):
                    total_p4_tested += 1
                    primer_id = f"P4_row{n}_cand{idx}"
                    
                    offtarget_hits = []
                    if primer_id in all_hits:
                        for hit in all_hits[primer_id][:5]:
                            if hit['identity'] >= args.blat_min_identity and hit['score'] >= args.blat_min_score:
                                is_expected = (hit['chrom'] == chrom and 
                                            abs(hit['start'] - candidate['primer_4_start']) <= 5 and 
                                            abs(hit['end'] - candidate['primer_4_end']) <= 5)
                                if not is_expected:
                                    offtarget_hits.append(hit)
                    
                    if not offtarget_hits:
                        total_p4_clean += 1
                        r4 = candidate
                        r4['primer_4_notes'] = 'No off-targets'
                        print(f"  ✓ Clean P4 found at MAX_ARM={attempt_arm_size}bp (candidate #{idx+1}, {total_p4_clean} clean / {total_p4_tested} tested)")
                        break
                    else:
                        total_p4_offtargets += 1
                        
                        # Categorize off-targets
                        other_chrom = [h for h in offtarget_hits if h['chrom'] != chrom]
                        same_chrom_far = [h for h in offtarget_hits if h['chrom'] == chrom and abs(h['start'] - position) > 10000]
                        same_chrom_near = [h for h in offtarget_hits if h['chrom'] == chrom and abs(h['start'] - position) <= 10000]
                        
                        # Categorize by WORST type present
                        if same_chrom_near:
                            category = 'same_chrom_near'
                        elif same_chrom_far:
                            category = 'same_chrom_far'
                        elif other_chrom:
                            category = 'diff_chrom'
                        else:
                            category = 'unknown'
                        
                        # Save fallbacks in priority order: diff_chrom > same_chrom_far > same_chrom_near
                        if category == 'diff_chrom' and not best_fallback_diff_chrom:
                            best_fallback_diff_chrom = candidate.copy()
                            best_fallback_diff_chrom['primer_4_notes'] = f'Off-targets on different chromosome(s): {len(other_chrom)} hit(s)'
                        elif category == 'same_chrom_far' and not best_fallback_same_chrom:
                            best_fallback_same_chrom = candidate.copy()
                            best_fallback_same_chrom['primer_4_notes'] = f'Off-targets on same chromosome >10kb away: {len(same_chrom_far)} hit(s)'
                        elif category == 'same_chrom_near' and not best_fallback_near:
                            best_fallback_near = candidate.copy()
                            best_fallback_near['primer_4_notes'] = f'Off-targets within 10kb: {len(same_chrom_near)} hit(s)'
            
            os.unlink(batch_fasta)
            os.unlink(batch_psl)
        else:
            r4 = candidates[0]
            r4['primer_4_notes'] = 'No off-target check performed'
            break
        
        # If we found a clean primer, stop searching
        if r4:
            break
    
    # If no clean primer found after trying all sizes, use best fallback
    if not r4:
        if best_fallback_diff_chrom:
            r4 = best_fallback_diff_chrom
            print(f"  ! Using fallback P4 (different chromosome off-targets) - tested {total_p4_tested} candidates, all had off-targets")
        elif best_fallback_same_chrom:
            r4 = best_fallback_same_chrom
            print(f"  ! Using fallback P4 (same chrom >10kb off-targets) - tested {total_p4_tested} candidates, all had off-targets")
        elif best_fallback_near:
            r4 = best_fallback_near
            print(f"  ! Using fallback P4 (same chrom <10kb off-targets - LAST RESORT) - tested {total_p4_tested} candidates, all had off-targets")
        else:
            r4 = {
                'primer_4_seq': 'ERROR',
                'primer_4_start': None,
                'primer_4_end': None,
                'primer_4_locus': None,
                'primer_4_tm': None,
                'primer_4_gc': None,
                'downstream_arm_length': None,
                'primer_4_notes': f'No suitable primer found (tested {total_p4_tested} candidates, all had off-targets)'
            }
    
    # Summary
    if r1['primer_1_seq'] != 'ERROR':
        print(f"✓ P1: {r1['primer_1_seq']}")
        if r1.get('primer_1_notes'):
            print(f"  Note: {r1['primer_1_notes']}")
    else:
        print(f"✗ P1: {r1['primer_1_notes']}")
    
    if r4['primer_4_seq'] != 'ERROR':
        print(f"✓ P4: {r4['primer_4_seq']}")
        if r4.get('primer_4_notes'):
            print(f"  Note: {r4['primer_4_notes']}")
    else:
        print(f"✗ P4: {r4['primer_4_notes']}")
    
    # Calculate arm lengths from genomic coordinates
    upstream_arm_length = None
    downstream_arm_length = None
    
    if (r1.get('primer_1_start') is not None and r1.get('primer_1_end') is not None 
        and p2_start is not None and p2_end is not None):
        leftmost = min(r1['primer_1_start'], p2_start)
        rightmost = max(r1['primer_1_end'], p2_end)
        upstream_arm_length = rightmost - leftmost + 1
        print(f"  Upstream arm: {upstream_arm_length}bp")
    
    if (r4.get('primer_4_start') is not None and r4.get('primer_4_end') is not None 
        and p3_start is not None and p3_end is not None):
        leftmost = min(r4['primer_4_start'], p3_start)
        rightmost = max(r4['primer_4_end'], p3_end)
        downstream_arm_length = rightmost - leftmost + (2 if strand == '-' else 1)
        print(f"  Downstream arm: {downstream_arm_length}bp")
    
    r1['upstream_arm_length'] = upstream_arm_length
    r4['downstream_arm_length'] = downstream_arm_length            
    
    # Create locus strings
    if r1.get('primer_1_start') is not None and r1.get('primer_1_end') is not None:
        r1['primer_1_locus'] = f"{chrom}:{r1['primer_1_start']}-{r1['primer_1_end']} ({strand})"
    else:
        r1['primer_1_locus'] = None
    
    if r4.get('primer_4_start') is not None and r4.get('primer_4_end') is not None:
        opposite_strand = '-' if strand == '+' else '+'
        r4['primer_4_locus'] = f"{chrom}:{r4['primer_4_start']}-{r4['primer_4_end']} ({opposite_strand})"
    else:
        r4['primer_4_locus'] = None
    
    # Remove start/end columns
    r1_filtered = {k: v for k, v in r1.items() if k not in ['primer_1_start', 'primer_1_end']}
    r4_filtered = {k: v for k, v in r4.items() if k not in ['primer_4_start', 'primer_4_end']}
    
    return {**r1_filtered, **r4_filtered}

def main():
    parser = argparse.ArgumentParser(
        description='Design outer primers (P1, P4) with checkpoint/resume',
        epilog="Use --flush-every N to save progress every N rows"
    )
    parser.add_argument('--input', required=True, help='Input CSV')
    parser.add_argument('--output', required=True, help='Output CSV')
    parser.add_argument('--genome', help='Genome FASTA')
    parser.add_argument('--flush-every', type=int, default=10, 
                       help='Save progress every N rows (default: 10)')
    parser.add_argument('--no-resume', action='store_true',
                       help='Overwrite existing output instead of resuming from checkpoint')
    parser.add_argument('--skip-offtarget-check', action='store_true',
                       help='Skip BLAT off-target checking (faster but less safe)')
    
    # Arm length parameters
    parser.add_argument('--min-arm', type=int, default=500,
                       help='Minimum homology arm length in bp (default: 500)')
    parser.add_argument('--max-arm', type=int, default=1000,
                       help='Maximum homology arm length in bp (default: 1000)')
    parser.add_argument('--max-arm-fallback', type=int, default=10000,
                   help='Maximum homology arm length to try before using primers with off-targets (default: 10000)')

    
    # Primer design parameters
    parser.add_argument('--primer-min-size', type=int, default=22,
                       help='Minimum primer length in bp (default: 22)')
    parser.add_argument('--primer-max-size', type=int, default=36,
                       help='Maximum primer length in bp (default: 36)')
    parser.add_argument('--primer-min-tm', type=float, default=50.0,
                       help='Minimum primer Tm in °C (default: 50.0)')
    parser.add_argument('--primer-max-tm', type=float, default=72.0,
                       help='Maximum primer Tm in °C (default: 72.0)')
    
    # BLAT parameters
    parser.add_argument('--blat-min-identity', type=int, default=90,
                       help='Minimum %% identity for BLAT hits (default: 90)')
    parser.add_argument('--blat-min-score', type=int, default=20,
                       help='Minimum alignment score for BLAT (default: 20)')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.genome):
        print(f"ERROR: Genome FASTA not found: {args.genome}")
        sys.exit(1)
    
    genome = load_genome(args.genome)

    # Check BLAT availability
    genome_fasta_for_blat = None
    if not args.skip_offtarget_check:
        try:
            subprocess.run(['blat'], capture_output=True, text=True)
            genome_fasta_for_blat = args.genome
            print("✓ BLAT found - will check for off-targets")
        except FileNotFoundError:
            print("WARNING: BLAT not found - skipping off-target checks")
            print("  Install BLAT or use --skip-offtarget-check to suppress this warning")
    else:
        print("Skipping off-target checks (--skip-offtarget-check)")
    
    print(f"\nReading {args.input}...")
    df = pd.read_csv(args.input)
    print(f"Loaded {len(df)} rows")
    
    print(f"\nConfiguration:")
    print(f"  Homology arm length: {args.min_arm}-{args.max_arm} bp (will try up to {args.max_arm_fallback} bp before using primers with off-targets)")
    print(f"  Primer size range: {args.primer_min_size}-{args.primer_max_size} bp")
    print(f"  Primer Tm range: {args.primer_min_tm}-{args.primer_max_tm}°C")
    print(f"  BLAT parameters: min identity={args.blat_min_identity}%, min score={args.blat_min_score}")
    print(f"  Flush every: {args.flush_every} rows")
    
    # Check for existing output and resume
    start_idx = 0
    if args.no_resume or not os.path.exists(args.output):
        print(f"\n[start] Creating new output file: {args.output}")
    else:
        print(f"\n[resume] Found existing output: {args.output}")
        try:
            existing = pd.read_csv(args.output)
            if 'primer_1_seq' in existing.columns:
                completed = existing['primer_1_seq'].notna() | existing['primer_4_seq'].notna()
                if completed.any():
                    start_idx = completed[::-1].idxmax() + 1
                    print(f"[resume] Resuming from row {start_idx + 1}")
        except Exception as e:
            print(f"[resume] Could not read existing output: {e}")
    
    # Pre-allocate results
    results = [None] * len(df)
    
    # Load existing results if resuming
    if start_idx > 0:
        try:
            existing = pd.read_csv(args.output)
            for i in range(start_idx):
                results[i] = {
                    col: existing.loc[i, col] if col in existing.columns else None
                    for col in ['primer_1_seq', 'primer_1_locus',
                               'primer_1_tm', 'primer_1_gc', 'upstream_arm_length', 'primer_1_notes',
                               'primer_4_seq', 'primer_4_locus',
                               'primer_4_tm', 'primer_4_gc', 'downstream_arm_length', 'primer_4_notes']
                }
        except:
            pass
    
    print(f"\nProcessing rows {start_idx + 1} to {len(df)}...")
    print(f"Saving every {args.flush_every} rows\n")
    
    processed = 0
    for idx in range(start_idx, len(df)):
        try:
            results[idx] = process(df.iloc[idx], idx+1, genome, args, genome_fasta_for_blat)
            processed += 1
            
            # Periodic save
            if processed % args.flush_every == 0:
                print(f"\n{'='*70}")
                print(f"CHECKPOINT: {processed} rows ({idx+1}/{len(df)})")
                print(f"{'='*70}")
                
                # Fill missing results with empty dicts
                save_results = []
                for r in results:
                    if r is not None:
                        save_results.append(r)
                    else:
                        save_results.append({
                            'primer_1_seq': None, 'primer_1_locus': None,
                            'primer_1_tm': None, 'primer_1_gc': None, 'upstream_arm_length': None,
                            'primer_1_notes': None,
                            'primer_4_seq': None, 'primer_4_locus': None,
                            'primer_4_tm': None, 'primer_4_gc': None, 'downstream_arm_length': None,
                            'primer_4_notes': None
                        })
                
                out = pd.concat([df, pd.DataFrame(save_results)], axis=1)
                out.to_csv(args.output, index=False)
                print(f"✓ Saved to {args.output}")
                
                # Stats
                completed = [r for r in results if r is not None]
                success = sum(1 for r in completed if r['primer_1_seq'] != 'ERROR' and r['primer_4_seq'] != 'ERROR')
                print(f"Success: {success}/{len(completed)} ({100*success/len(completed):.1f}%)\n")
                
        except Exception as e:
            print(f"\nERROR row {idx+1}: {e}")
            traceback.print_exc()
            results[idx] = {
                'primer_1_seq': 'ERROR', 'primer_1_notes': str(e),
                'primer_4_seq': 'ERROR', 'primer_4_notes': str(e)
            }
    
    # Final save
    print(f"\n{'='*70}")
    print(f"FINAL SAVE")
    print(f"{'='*70}")
    
    save_results = []
    for r in results:
        if r is not None:
            save_results.append(r)
        else:
            save_results.append({
                'primer_1_seq': None, 'primer_1_locus': None,
                'primer_1_tm': None, 'primer_1_gc': None, 'upstream_arm_length': None,
                'primer_1_notes': None,
                'primer_4_seq': None, 'primer_4_locus': None,
                'primer_4_tm': None, 'primer_4_gc': None, 'downstream_arm_length': None,
                'primer_4_notes': None
            })
    
    out = pd.concat([df, pd.DataFrame(save_results)], axis=1)
    out.to_csv(args.output, index=False)
    print(f"✓ Saved to {args.output}")
    
    # Final stats
    completed = [r for r in results if r is not None]
    success = sum(1 for r in completed if r['primer_1_seq'] != 'ERROR' and r['primer_4_seq'] != 'ERROR')
    partial = sum(1 for r in completed if (r['primer_1_seq'] != 'ERROR') != (r['primer_4_seq'] != 'ERROR'))
    
    print(f"\n{'='*70}")
    print(f"FINAL SUMMARY")
    print(f"{'='*70}")
    print(f"Total: {len(completed)}")
    print(f"Both primers: {success} ({100*success/len(completed):.1f}%)")
    print(f"One primer: {partial}")
    print(f"Both failed: {len(completed) - success - partial}")
    print(f"{'='*70}\n")

if __name__ == '__main__':
    main()