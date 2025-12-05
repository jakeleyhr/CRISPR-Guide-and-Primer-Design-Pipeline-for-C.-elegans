#!/usr/bin/env python3
"""
Homology Arm Primer Design Pipeline - Step 3: Outer Primer Design
Designs outer PCR primers (P1 and P4) that pair with inner primers to amplify
complete homology arms for CRISPR knock-in repair templates. Uses batched 
processing with Primer3 for optimal primer design and single BLAT call for 
off-target checking across all candidates.

Input:
  - CSV file from script 5 containing inner primers (P2 and P3) with genomic
    coordinates and sequences
    (requires: chrom, pos, strand, primer_2_seq_genomic, primer_2_locus,
     primer_3_seq_genomic, primer_3_locus)

Output:
  - CSV file with all input columns plus outer primer information:
      • primer_1_seq: Upstream outer primer (pairs with P2)
      • primer_4_seq: Downstream outer primer (pairs with P3)
      • primer_1_locus, primer_4_locus: Genomic coordinates
      • primer_1_tm, primer_4_tm: Melting temperatures
      • primer_1_gc, primer_4_gc: GC content percentages
      • upstream_arm_length, downstream_arm_length: Final amplicon sizes
      • Off-target annotations with rescue status

Primer Nomenclature (Gene-centric):
  - P1 + P2 amplify gene UPSTREAM homology arm (500-1000bp target)
  - P3 + P4 amplify gene DOWNSTREAM homology arm (500-1000bp target)
  For minus strand genes: gene upstream = genomic RIGHT side
  (Primer names always follow gene orientation, not genomic coordinates)

Batched Design Strategy:
  1. Design ALL primer candidates for all rows at initial flank size
     - Extracts flanking sequences around insertion sites
     - Locates inner primers (P2, P3) within sequences
     - Generates 50 P1 candidates and 50 P4 candidates per row (configurable)
     - Target product size: 500-1000bp (configurable)
     - Primers designed with Primer3 to match P2/P3 Tm when possible
     - Preferentially selects primers with 3' GC-clamp
  
  2. Check ALL primers in one BLAT call
     - Creates single FASTA with all P1 and P4 candidates
     - Single BLAT run checks thousands of primers at once
     - Parses results to identify clean primers (no off-targets)
  
  3. Select best primer pair for each row
     - Prioritizes clean primers (no off-targets)
     - Falls back to primers with acceptable off-targets if needed
  
  4. Progressive flank widening for failed rows (500bp increments)
     - Rows with off-targets try wider flanks (up to 10kb)
     - Independent rescue: P1 and P4 can succeed at different flank widths
  
  5. Fallback for remaining failures
     - Uses best available primers even if off-targets remain
     - Annotates with detailed off-target information

Off-Target Detection (BLAT):
  - Single batched BLAT call checks all primers simultaneously
  - Categories:
    * same_chrom_near: Within 10kb of expected location (concerning)
    * same_chrom_far: Same chromosome, >10kb away (acceptable with pre-amplification)
    * diff_chrom: Different chromosome (acceptable)
  - Progressive rescue tries wider flanks before accepting off-targets
  - Can be disabled with --skip-offtarget-check (not recommended)

Required arguments:
  --input PATH       Input CSV from script 5 with inner primers (P2, P3)
  --output PATH      Output CSV path with outer primers (P1, P4)
  --genome PATH      Genome FASTA file (e.g., ce11.fa for C. elegans)

Optional arguments:
  --batch-size N               Process N rows at a time (default: 100, all at once)
  --flush-every N              Write results every N rows for progress tracking (default: 10)
  --no-resume                  Overwrite output instead of resuming from checkpoint
  --skip-offtarget-check       Skip BLAT off-target check (faster but not recommended)
  
  Homology arm parameters:
  --min-arm N                  Minimum homology arm length in bp (default: 500)
  --max-arm N                  Initial maximum arm length in bp (default: 1000)
  --max-arm-fallback N         Maximum arm length for rescue attempts (default: 10000)
  
  Primer design parameters:
  --primer-min-size N          Minimum primer length in bp (default: 22)
  --primer-max-size N          Maximum primer length in bp (default: 36)
  --primer-min-tm N            Minimum primer Tm in °C (default: 50.0)
  --primer-max-tm N            Maximum primer Tm in °C (default: 72.0)
  
  BLAT parameters:
  --blat-min-identity N        Minimum % identity for BLAT hits (default: 90)
  --blat-min-score N           Minimum alignment score for BLAT (default: 20)

Example:
  python 6_designouterprimersbatch.py \
    --input 5.allwormguidesinternalprimers.csv \
    --output 6.allwormguidesbothprimers.csv \
    --genome ~/Desktop/wbcel235/ce11.fa \
    --batch-size 100

Recommended: Install BLAT for off-target checking
  - mkdir -p ~/bin
  - cd ~/bin
  - curl -o blat http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/blat/blat
  - chmod +x blat
  - echo 'export PATH="$HOME/bin:$PATH"' >> ~/.zshrc
  - blat  # Test

"""

import argparse, sys, os, subprocess, tempfile, traceback, time
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

def extract_flanks(genome, chrom, position, strand, flank):
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
    
    # Safety margin needed because P2/P3 primer coordinates have edge case issues
    # This causes Primer3 to find slightly different primers than without margin,
    # but ensures all primers can be reliably found
    safety_margin = 3
    
    if strand == '+':
        # Upstream (gene left): extract left with margin, then RC
        upstream_start = max(0, position - flank - safety_margin)
        upstream_seq_forward = seq[upstream_start:position + safety_margin]
        upstream_seq = rc(upstream_seq_forward)
        
        # Downstream (gene right): extract right with margin
        downstream_start = position + 1 - safety_margin
        downstream_end = min(len(seq), position + flank + safety_margin)
        downstream_seq = seq[downstream_start:downstream_end]
        
    else:  # strand == '-'
        # Upstream (gene left, genomic right): extract right with margin
        upstream_start = position - safety_margin
        upstream_end = min(len(seq), position + flank + safety_margin)
        upstream_seq = seq[upstream_start:upstream_end]
        
        # Downstream (gene right, genomic left): extract left with margin, then RC
        downstream_start = max(0, position - flank - safety_margin)
        downstream_seq_forward = seq[downstream_start:position + safety_margin]
        downstream_seq = rc(downstream_seq_forward)
    
    return upstream_seq, downstream_seq, upstream_start, downstream_start

def rc(seq):
    """Reverse complement"""
    comp = {'A':'T', 'T':'A', 'G':'C', 'C':'G',
            'a':'t', 't':'a', 'g':'c', 'c':'g',
            'N':'N', 'n':'n'}
    return ''.join([comp.get(b, b) for b in seq[::-1]])

def find_primer_flexible(primer, seq, name="primer"):
    """Find primer in sequence, trying multiple strategies"""
    # Fast path: exact match forward
    pos = seq.find(primer)
    if pos != -1:
        return pos, len(primer), 'forward'
    
    # Fast path: exact match RC
    primer_rc = rc(primer)
    pos = seq.find(primer_rc)
    if pos != -1:
        return pos, len(primer_rc), 'RC'
    
    # SLOW: Only do fuzzy finding if exact matches failed
    # This is expensive - scanning entire sequence position by position
    pos = fuzzy_find(primer, seq, max_mismatches=5)
    if pos is not None:
        return pos, len(primer), 'forward_fuzzy'
    
    pos = fuzzy_find(primer_rc, seq, max_mismatches=5)
    if pos is not None:
        return pos, len(primer_rc), 'RC_fuzzy'
    
    # Try truncated matches at end
    for truncate in range(1, min(6, len(primer))):
        truncated = primer[:-truncate]
        if len(seq) >= len(truncated) and seq[-len(truncated):] == truncated:
            return len(seq) - len(truncated), len(truncated), 'forward_truncated'
    
    return None, None, None

def fuzzy_find(primer, seq, max_mismatches=5):
    """Find primer allowing up to max_mismatches mismatches"""
    best_pos = None
    best_mismatches = max_mismatches + 1
    primer_len = len(primer)
    seq_len = len(seq)
    
    # Early termination if we find a perfect or near-perfect match
    for i in range(seq_len - primer_len + 1):
        mismatches = 0
        # Count mismatches with early exit
        for j in range(primer_len):
            if seq[i+j] != primer[j]:
                mismatches += 1
                # Early exit if too many mismatches already
                if mismatches > max_mismatches:
                    break
        
        if mismatches < best_mismatches:
            best_mismatches = mismatches
            best_pos = i
            # Perfect match found - stop searching
            if mismatches == 0:
                break
            # Very good match - stop searching
            if mismatches <= 1:
                break
    
    if best_mismatches <= max_mismatches:
        return best_pos
    return None, None

def parse_blat_psl(psl_path):
    """Parse BLAT PSL output"""
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

def _gc_end_pref(seq):
    """Prefer primers that end in G or C at the 3' end"""
    return 0 if seq[-1].upper() in ("G","C") else 1

def design_outer_primer_for_p1p2_arm(seq, p2_primer, seq_start, strand, cut_position, args, target_arm_size, num_return=50, debug=False, safety_margin=3):
    """Design P1 candidates to pair with P2 for gene upstream (P1+P2) homology arm
    
    Args:
        safety_margin: Number of bases added to extraction for tolerance (default 3)
    
    Returns: (candidates_list, error_message)
        - candidates_list: list of primer dicts (empty if failed)
        - error_message: None if success, string if failed
    """
    p2_tm_result = primer3.calc_tm(p2_primer)
    p2_len = len(p2_primer)
    
    p1_min_size = args.primer_min_size
    p1_max_size = args.primer_max_size
    p1_opt_size = max(args.primer_min_size, min(args.primer_max_size, p2_len))
    
    p1_opt_tm = p2_tm_result
    p1_min_tm = args.primer_min_tm
    p1_max_tm = args.primer_max_tm
    
    if p2_tm_result < args.primer_min_tm:
        p1_opt_tm = args.primer_min_tm
    elif p2_tm_result > args.primer_max_tm:
        p1_opt_tm = args.primer_max_tm
    
    if p1_min_tm >= p1_opt_tm:
        p1_min_tm = p1_opt_tm - 1.0
    if p1_max_tm <= p1_opt_tm:
        p1_max_tm = p1_opt_tm + 1.0
    
    p2_pos, p2_len_found, p2_orientation = find_primer_flexible(p2_primer, seq, "P2")
    
    if p2_pos is None:
        if debug:
            print(f"      P2 NOT FOUND in upstream sequence")
            print(f"      P2 primer: {p2_primer}")
            print(f"      Tried forward, RC, fuzzy, and truncated matches")
        return [], "P2 not found in upstream sequence"
    
    if debug:
        print(f"      P2 found at position {p2_pos} ({p2_orientation})")
    
    try:
        search_end = len(seq)
        search_start = max(0, search_end - 500)
        
        p3_args_custom = {
            'PRIMER_MIN_SIZE': p1_min_size,
            'PRIMER_MAX_SIZE': p1_max_size,
            'PRIMER_OPT_SIZE': p1_opt_size,
            'PRIMER_MIN_TM': p1_min_tm,
            'PRIMER_MAX_TM': p1_max_tm,
            'PRIMER_OPT_TM': p1_opt_tm,
            'PRIMER_MIN_GC': 40.0,  
            'PRIMER_MAX_GC': 60.0,  
            'PRIMER_PRODUCT_SIZE_RANGE': [[args.min_arm, target_arm_size]],
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
            return [], "Primer3 found no valid P1 primers"
        
        candidates = []
        for i in range(num_returned):
            try:
                seq_p1 = res[f'PRIMER_RIGHT_{i}_SEQUENCE']
                start = res[f'PRIMER_RIGHT_{i}'][0]
                length = res[f'PRIMER_RIGHT_{i}'][1]
                tm = res[f'PRIMER_RIGHT_{i}_TM']
                gc = res[f'PRIMER_RIGHT_{i}_GC_PERCENT']
                
                if strand == '+':
                    # Plus strand uses cut_position directly
                    # Adjust start to account for safety margin in extraction
                    adjusted_start = start - safety_margin
                    p1_start_genomic = cut_position - adjusted_start
                    p1_end_genomic = cut_position - (adjusted_start - length + 1)
                else:
                    # Minus strand uses seq_start directly (already includes margin)
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
        
        return candidates, None
        
    except Exception as e:
        return [], f"Primer3 exception: {str(e)}"

def design_outer_primer_for_p3p4_arm(seq, p3_primer, seq_start, strand, cut_position, args, target_arm_size, num_return=50, debug=False, safety_margin=3):
    """Design P4 candidates to pair with P3 for gene downstream (P3+P4) homology arm
    
    Args:
        safety_margin: Number of bases added to extraction for tolerance (default 3)
    
    Returns: (candidates_list, error_message)
        - candidates_list: list of primer dicts (empty if failed)
        - error_message: None if success, string if failed
    """
    p3_tm_result = primer3.calc_tm(p3_primer)
    p3_len = len(p3_primer)
    
    p4_min_size = args.primer_min_size
    p4_max_size = args.primer_max_size
    p4_opt_size = max(args.primer_min_size, min(args.primer_max_size, p3_len))
    
    p4_opt_tm = p3_tm_result
    p4_min_tm = args.primer_min_tm
    p4_max_tm = args.primer_max_tm
    
    if p3_tm_result < args.primer_min_tm:
        p4_opt_tm = args.primer_min_tm
    elif p3_tm_result > args.primer_max_tm:
        p4_opt_tm = args.primer_max_tm
    
    if p4_min_tm >= p4_opt_tm:
        p4_min_tm = p4_opt_tm - 1.0
    if p4_max_tm <= p4_opt_tm:
        p4_max_tm = p4_opt_tm + 1.0
    
    p3_pos, p3_len_found, p3_orientation = find_primer_flexible(p3_primer, seq, "P3")
    if p3_pos is None:
        if debug:
            print(f"      P3 NOT FOUND in downstream sequence")
            print(f"      P3 primer: {p3_primer}")
            print(f"      Tried forward, RC, fuzzy, and truncated matches")
        return [], "P3 not found in downstream sequence"
    
    if debug:
        print(f"      P3 found at position {p3_pos} ({p3_orientation})")
    
    try:
        search_end = len(seq)
        search_start = max(0, search_end - 500)
        
        if search_end <= search_start:
            return [], "Downstream sequence too short"
        
        p3_args_custom = {
            'PRIMER_MIN_SIZE': p4_min_size,
            'PRIMER_MAX_SIZE': p4_max_size,
            'PRIMER_OPT_SIZE': p4_opt_size,
            'PRIMER_MIN_TM': p4_min_tm,
            'PRIMER_MAX_TM': p4_max_tm,
            'PRIMER_OPT_TM': p4_opt_tm,
            'PRIMER_MIN_GC': 40.0,
            'PRIMER_MAX_GC': 60.0,
            'PRIMER_PRODUCT_SIZE_RANGE': [[args.min_arm, target_arm_size]],
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
            return [], "Primer3 found no valid P4 primers"
        
        candidates = []
        for i in range(num_returned):
            try:
                seq_p4 = res[f'PRIMER_RIGHT_{i}_SEQUENCE']
                start = res[f'PRIMER_RIGHT_{i}'][0]
                length = res[f'PRIMER_RIGHT_{i}'][1]
                tm = res[f'PRIMER_RIGHT_{i}_TM']
                gc = res[f'PRIMER_RIGHT_{i}_GC_PERCENT']

                if strand == '+':
                    # Plus strand uses seq_start directly (already includes margin)
                    p4_start_genomic = seq_start + (start - length + 1) + 1
                    p4_end_genomic = seq_start + start + 1
                else:
                    # Minus strand uses cut_position directly
                    # Adjust start to account for safety margin in extraction
                    adjusted_start = start - safety_margin
                    p4_end_genomic = cut_position - (adjusted_start - length + 1)
                    p4_start_genomic = cut_position - adjusted_start
                
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
        
        return candidates, None
        
    except Exception as e:
        return [], f"Primer3 exception: {str(e)}"

def parse_locus_coords(locus_str):
    """Parse locus string like 'I:5107814-5107843 (-)' to get start and end coordinates"""
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

def get_gene_name(row, row_idx):
    """Safely extract gene name from row"""
    gene = row.get('gene')
    if gene is None or (isinstance(gene, float) and pd.isna(gene)):
        return f'row{row_idx+1}'
    return str(gene)

def run_blat_batch(primer_fasta_path, genome_fasta_path, output_psl_path, min_identity, min_score):
    """Run BLAT to search primers against genome"""
    cmd = [
        'blat',
        genome_fasta_path,
        primer_fasta_path,
        output_psl_path,
        '-stepSize=5',
        '-repMatch=1000000',
        f'-minScore={min_score}',
        f'-minIdentity={min_identity}',
        # Removed -minMatch=1 for better performance (uses BLAT default)
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        return result.returncode == 0
    except:
        return False

def categorize_offtargets(offtarget_hits, expected_chrom, cut_position):
    """Categorize off-targets by location"""
    if not offtarget_hits:
        return None, 0, None
    
    other_chrom = [h for h in offtarget_hits if h['chrom'] != expected_chrom]
    same_chrom_far = [h for h in offtarget_hits if h['chrom'] == expected_chrom and abs(h['start'] - cut_position) > 10000]
    same_chrom_near = [h for h in offtarget_hits if h['chrom'] == expected_chrom and abs(h['start'] - cut_position) <= 10000]
    
    if same_chrom_near:
        return 'same_chrom_near', len(same_chrom_near), f'{len(same_chrom_near)} hits within 10kb'
    elif same_chrom_far:
        return 'same_chrom_far', len(same_chrom_far), f'{len(same_chrom_far)} hits >10kb away'
    elif other_chrom:
        return 'diff_chrom', len(other_chrom), f'{len(other_chrom)} hits on other chromosomes'
    else:
        return 'unknown', len(offtarget_hits), f'{len(offtarget_hits)} unclassified hits'

def check_offtargets(hits, primer_id, expected_chrom, expected_start, expected_end, cut_position, args, tolerance=5):
    """Check if a primer has off-target hits"""
    if primer_id not in hits:
        return False, 0, "No hits", None
    
    good_hits = [h for h in hits[primer_id] 
                 if h['identity'] >= args.blat_min_identity and h['score'] >= args.blat_min_score]
    
    if not good_hits:
        return False, 0, "No quality hits", None
    
    offtargets = []
    for hit in good_hits:
        is_expected = (hit['chrom'] == expected_chrom and 
                      abs(hit['start'] - expected_start) <= tolerance and 
                      abs(hit['end'] - expected_end) <= tolerance)
        if not is_expected:
            offtargets.append(hit)
    
    if not offtargets:
        return False, 0, "Only on-target", None
    
    category, count, description = categorize_offtargets(offtargets, expected_chrom, cut_position)
    return True, count, description, category

def create_result_from_candidates(p1_cand, p4_cand, chrom, strand, p2_start, p2_end, p3_start, p3_end, p1_notes, p4_notes):
    """Create result dict from selected candidates"""
    upstream_arm_length = None
    downstream_arm_length = None
    
    if (p1_cand['primer_1_start'] is not None and p1_cand['primer_1_end'] is not None 
        and p2_start is not None and p2_end is not None):
        leftmost = min(p1_cand['primer_1_start'], p2_start)
        rightmost = max(p1_cand['primer_1_end'], p2_end)
        upstream_arm_length = rightmost - leftmost + 1
    
    if (p4_cand['primer_4_start'] is not None and p4_cand['primer_4_end'] is not None 
        and p3_start is not None and p3_end is not None):
        leftmost = min(p4_cand['primer_4_start'], p3_start)
        rightmost = max(p4_cand['primer_4_end'], p3_end)
        downstream_arm_length = rightmost - leftmost + (2 if strand == '-' else 1)
    
    p1_locus = None
    p4_locus = None
    
    if p1_cand['primer_1_start'] is not None and p1_cand['primer_1_end'] is not None:
        p1_locus = f"{chrom}:{p1_cand['primer_1_start']}-{p1_cand['primer_1_end']} ({strand})"
    
    if p4_cand['primer_4_start'] is not None and p4_cand['primer_4_end'] is not None:
        opposite_strand = '-' if strand == '+' else '+'
        p4_locus = f"{chrom}:{p4_cand['primer_4_start']}-{p4_cand['primer_4_end']} ({opposite_strand})"
    
    return {
        'primer_1_seq': p1_cand['primer_1_seq'],
        'primer_1_locus': p1_locus,
        'primer_1_tm': p1_cand['primer_1_tm'],
        'primer_1_gc': p1_cand['primer_1_gc'],
        'upstream_arm_length': upstream_arm_length,
        'primer_1_notes': p1_notes,
        'primer_4_seq': p4_cand['primer_4_seq'],
        'primer_4_locus': p4_locus,
        'primer_4_tm': p4_cand['primer_4_tm'],
        'primer_4_gc': p4_cand['primer_4_gc'],
        'downstream_arm_length': downstream_arm_length,
        'primer_4_notes': p4_notes
    }

def process_batch(batch_rows, genome, genome_fasta_path, batch_start_idx, args):
    """
    OPTIMIZED: Process a batch of rows by designing ALL primers first,
    then running ONE BLAT call for everything.
    
    Includes progressive flank widening: tries initial flank, then increases
    by 500bp increments up to max_arm_fallback if no clean primers found.
    """
    print(f"\n{'='*70}")
    print(f"BATCH: Processing rows {batch_start_idx+1} to {batch_start_idx+len(batch_rows)}")
    print(f"{'='*70}")
    
    batch_results = []
    
    # Step 1: Design ALL primer candidates for all rows at initial flank size
    import time
    step1_start = time.time()
    print(f"\nStep 1: Designing primers for {len(batch_rows)} rows at {args.max_arm}bp flank...")
    
    primer_designs = []
    errors = []
    
    primer3_time = 0
    extraction_time = 0
    
    for idx, (row_idx, row) in enumerate(batch_rows):
        if (idx + 1) % 10 == 0 or idx == 0:
            elapsed = time.time() - step1_start
            rate = (idx + 1) / elapsed if elapsed > 0 else 0
            print(f"  Progress: {idx + 1}/{len(batch_rows)} rows ({rate:.1f} rows/sec, {elapsed:.1f}s elapsed)...", end='\r')
        
        chrom = row.get('chrom', '')
        position = row.get('pos', None)
        
        # Handle NaN values in primer sequences
        p2_raw = row.get('primer_2_seq_genomic', '')
        p3_raw = row.get('primer_3_seq_genomic', '')
        p2 = str(p2_raw).upper() if pd.notna(p2_raw) else ''
        p3 = str(p3_raw).upper() if pd.notna(p3_raw) else ''
        
        strand = row.get('strand', '+')
        p2_locus = row.get('primer_2_locus', '')
        p3_locus = row.get('primer_3_locus', '')
        
        gene_name = get_gene_name(row, row_idx)
        
        # Debug: print first row details
        if idx == 0:
            print(f"\n  DEBUG first row (row {row_idx+1}, {gene_name}):")
            print(f"    chrom: '{chrom}' (empty={not chrom})")
            print(f"    position: {position} (None={position is None})")
            print(f"    strand: {strand}")
            print(f"    p2_raw: {p2_raw} (type: {type(p2_raw).__name__})")
            print(f"    p3_raw: {p3_raw} (type: {type(p3_raw).__name__})")
            print(f"    p2: '{p2[:20] if p2 else ''}...' (empty={not p2})")
            print(f"    p3: '{p3[:20] if p3 else ''}...' (empty={not p3})")
            
            # Check if there's a primer_2_seq column (template seq vs genomic seq)
            p2_template = row.get('primer_2_seq', '')
            p3_template = row.get('primer_3_seq', '')
            print(f"    primer_2_seq (template): {p2_template[:30] if p2_template and pd.notna(p2_template) else 'N/A'}...")
            print(f"    primer_3_seq (template): {p3_template[:30] if p3_template and pd.notna(p3_template) else 'N/A'}...")
            print(f"    primer_2_locus: {p2_locus}")
            print(f"    primer_3_locus: {p3_locus}")
        
        # Skip rows without inner primers
        if not p2 or not p3:
            # Don't add to errors or results - just skip silently
            continue
        
        if not chrom or position is None:
            missing = []
            if not chrom: missing.append('chrom')
            if position is None: missing.append('pos')
            errors.append(f"Row {row_idx+1} ({gene_name}): Missing required fields: {', '.join(missing)}")
            result = {
                'primer_1_seq': 'ERROR', 'primer_1_locus': None,
                'primer_1_tm': None, 'primer_1_gc': None, 'upstream_arm_length': None,
                'primer_1_notes': 'Missing required fields',
                'primer_4_seq': 'ERROR', 'primer_4_locus': None,
                'primer_4_tm': None, 'primer_4_gc': None, 'downstream_arm_length': None,
                'primer_4_notes': 'Missing required fields'
            }
            batch_results.append((row_idx, result))
            continue
        
        try:
            position = int(position)
        except:
            errors.append(f"Row {row_idx+1} ({gene_name}): Invalid position")
            result = {
                'primer_1_seq': 'ERROR', 'primer_1_locus': None,
                'primer_1_tm': None, 'primer_1_gc': None, 'upstream_arm_length': None,
                'primer_1_notes': 'Invalid position',
                'primer_4_seq': 'ERROR', 'primer_4_locus': None,
                'primer_4_tm': None, 'primer_4_gc': None, 'downstream_arm_length': None,
                'primer_4_notes': 'Invalid position'
            }
            batch_results.append((row_idx, result))
            continue
        
        # DEBUG: Print detailed info for first row
        if idx == 0:
            print(f"\n  === DEBUG ROW {row_idx+1} ({gene_name}) ===")
            print(f"    CSV position: {position} (type: {type(position)})")
            print(f"    Strand: {strand}")
            print(f"    P2 locus: {p2_locus}")
            print(f"    P3 locus: {p3_locus}")
            print(f"    P2 seq: {p2[:20]}...")
            print(f"    P3 seq: {p3[:20]}...")
        
        ext_start = time.time()
        upstream_seq, downstream_seq, upstream_start, downstream_start = extract_flanks(
            genome, chrom, position, strand, args.max_arm
        )
        extraction_time += time.time() - ext_start
        
        if idx == 0:
            print(f"    After extract_flanks:")
            print(f"      upstream_start={upstream_start}, len={len(upstream_seq) if upstream_seq else 0}")
            print(f"      downstream_start={downstream_start}, len={len(downstream_seq) if downstream_seq else 0}")
            if upstream_seq:
                print(f"      Upstream first 30bp: {upstream_seq[:30]}")
            if downstream_seq:
                print(f"      Downstream first 30bp: {downstream_seq[:30]}")
        
        if upstream_seq is None:
            errors.append(f"Row {row_idx+1} ({gene_name}): Chromosome not found")
            result = {
                'primer_1_seq': 'ERROR', 'primer_1_locus': None,
                'primer_1_tm': None, 'primer_1_gc': None, 'upstream_arm_length': None,
                'primer_1_notes': 'Chromosome not found',
                'primer_4_seq': 'ERROR', 'primer_4_locus': None,
                'primer_4_tm': None, 'primer_4_gc': None, 'downstream_arm_length': None,
                'primer_4_notes': 'Chromosome not found'
            }
            batch_results.append((row_idx, result))
            continue
        
        # Debug: print sequence lengths for first row
        if idx == 0:
            print(f"    upstream_seq length: {len(upstream_seq)}bp")
            print(f"    upstream_seq first 50bp: {upstream_seq[:50]}")
            print(f"    upstream_seq last 50bp: {upstream_seq[-50:]}")
            print(f"    downstream_seq length: {len(downstream_seq)}bp")
            print(f"    downstream_seq first 50bp: {downstream_seq[:50]}")
            print(f"    downstream_seq last 50bp: {downstream_seq[-50:]}")
            
            # Try to find P2 in the genomic region manually
            print(f"    Checking P2 manually:")
            if p2 in upstream_seq:
                pos = upstream_seq.find(p2)
                print(f"      P2 found FORWARD at position {pos} in upstream")
            elif rc(p2) in upstream_seq:
                pos = upstream_seq.find(rc(p2))
                print(f"      P2 found RC at position {pos} in upstream")
            else:
                # Check if it's just at the boundary
                if upstream_seq.startswith(p2[1:]):
                    print(f"      P2 missing first base - found from position 1: {upstream_seq[:len(p2)]}")
                else:
                    print(f"      P2 NOT found in upstream at all")
                
            # Try to find P3
            if p3 in upstream_seq:
                print(f"    P3 found FORWARD in upstream (wrong place!)")
            if rc(p3) in upstream_seq:
                print(f"    P3 found RC in upstream (wrong place!)")
            if p3 in downstream_seq:
                print(f"    P3 found FORWARD in downstream")
            if rc(p3) in downstream_seq:
                print(f"    P3 found RC in downstream")
        
        # Design P1 candidates
        p3_start = time.time()
        p1_candidates, p1_error = design_outer_primer_for_p1p2_arm(
            upstream_seq, p2, upstream_start, strand, position, args, 
            target_arm_size=args.max_arm, num_return=50, debug=(idx == 0)
        )
        primer3_time += time.time() - p3_start
        
        # Debug: print design results for first row
        if idx == 0 or not p1_candidates:
            print(f"    P1 candidates: {len(p1_candidates)}")
            if not p1_candidates:
                print(f"    P1 design failed: {p1_error}")
                print(f"    DEBUG: Row {row_idx+1} ({gene_name})")
                print(f"      Position: {position}, Strand: {strand}")
                print(f"      P2 locus: {p2_locus}")
                print(f"      P2 sequence: {p2[:30]}...")
                print(f"      Upstream start: {upstream_start}, length: {len(upstream_seq)}bp")
                print(f"      Upstream first 50bp: {upstream_seq[:50]}")
                print(f"      Looking for P2 in upstream:")
                if p2 in upstream_seq:
                    print(f"        Found at position {upstream_seq.find(p2)}")
                elif rc(p2) in upstream_seq:
                    print(f"        Found RC at position {upstream_seq.find(rc(p2))}")
                else:
                    print(f"        NOT FOUND (neither forward nor RC)")
        
        # Design P4 candidates
        p3_start = time.time()
        p4_candidates, p4_error = design_outer_primer_for_p3p4_arm(
            downstream_seq, p3, downstream_start, strand, position, args,
            target_arm_size=args.max_arm, num_return=50, debug=(idx == 0)
        )
        primer3_time += time.time() - p3_start
        
        # Debug: print design results for first row
        if idx == 0 or not p4_candidates:
            print(f"    P4 candidates: {len(p4_candidates)}")
            if not p4_candidates:
                print(f"    P4 design failed: {p4_error}")
                print(f"    DEBUG: Row {row_idx+1} ({gene_name})")
                print(f"      Position: {position}, Strand: {strand}")
                print(f"      P3 locus: {p3_locus}")
                print(f"      P3 sequence: {p3[:30]}...")
                print(f"      Downstream start: {downstream_start}, length: {len(downstream_seq)}bp")
                print(f"      Downstream first 50bp: {downstream_seq[:50]}")
                print(f"      Looking for P3 in downstream:")
                if p3 in downstream_seq:
                    print(f"        Found at position {downstream_seq.find(p3)}")
                elif rc(p3) in downstream_seq:
                    print(f"        Found RC at position {downstream_seq.find(rc(p3))}")
                else:
                    print(f"        NOT FOUND (neither forward nor RC)")
        
        if not p1_candidates or not p4_candidates:
            error_msg = []
            if not p1_candidates:
                error_msg.append(f"P1: {p1_error}")
            if not p4_candidates:
                error_msg.append(f"P4: {p4_error}")
            error_str = "; ".join(error_msg)
            errors.append(f"Row {row_idx+1} ({gene_name}): {error_str}")
            result = {
                'primer_1_seq': 'ERROR' if not p1_candidates else None,
                'primer_1_locus': None,
                'primer_1_tm': None,
                'primer_1_gc': None,
                'upstream_arm_length': None,
                'primer_1_notes': p1_error if not p1_candidates else None,
                'primer_4_seq': 'ERROR' if not p4_candidates else None,
                'primer_4_locus': None,
                'primer_4_tm': None,
                'primer_4_gc': None,
                'downstream_arm_length': None,
                'primer_4_notes': p4_error if not p4_candidates else None
            }
            batch_results.append((row_idx, result))
            continue
        
        # Sort by GC-clamp preference
        p1_candidates.sort(key=lambda c: _gc_end_pref(c['primer_1_seq']))
        p4_candidates.sort(key=lambda c: _gc_end_pref(c['primer_4_seq']))
        
        p2_start, p2_end = parse_locus_coords(p2_locus)
        p3_start, p3_end = parse_locus_coords(p3_locus)
        
        primer_designs.append({
            'row_idx': row_idx,
            'row': row,
            'p1_candidates': p1_candidates,
            'p4_candidates': p4_candidates,
            'chrom': chrom,
            'strand': strand,
            'position': position,
            'p2': p2,
            'p3': p3,
            'p2_start': p2_start,
            'p2_end': p2_end,
            'p3_start': p3_start,
            'p3_end': p3_end
        })
    
    step1_time = time.time() - step1_start
    print(f"\n  Step 1 completed in {step1_time:.1f}s")
    print(f"    - Average: {step1_time/len(batch_rows):.2f}s per row")
    print(f"    - Designed primers for {len(primer_designs)} rows")
    print(f"    - Timing breakdown:")
    print(f"      • Sequence extraction: {extraction_time:.2f}s ({extraction_time/step1_time*100:.1f}%)")
    print(f"      • Primer3 calls: {primer3_time:.2f}s ({primer3_time/step1_time*100:.1f}%)")
    other_time = step1_time - extraction_time - primer3_time
    print(f"      • Other (find_primer, etc): {other_time:.2f}s ({other_time/step1_time*100:.1f}%)")
    
    if errors:
        print(f"  Errors during design: {len(errors)}")
        # Print first few errors for debugging
        for error in errors[:5]:
            print(f"    - {error}")
        if len(errors) > 5:
            print(f"    ... and {len(errors)-5} more errors")
    
    if not primer_designs:
        return batch_results
    
    # Step 2: Create mega-FASTA with ALL primer candidates
    total_primers = sum(len(d['p1_candidates']) + len(d['p4_candidates']) 
                       for d in primer_designs)
    print(f"\nStep 2: Preparing to check {total_primers} primers in ONE BLAT call...")
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        primer_fasta = f.name
        for design in primer_designs:
            row_idx = design['row_idx']
            
            for cand_idx, cand in enumerate(design['p1_candidates']):
                primer_id = f"row{row_idx}_P1_cand{cand_idx}"
                f.write(f">{primer_id}\n{cand['primer_1_seq']}\n")
            
            for cand_idx, cand in enumerate(design['p4_candidates']):
                primer_id = f"row{row_idx}_P4_cand{cand_idx}"
                f.write(f">{primer_id}\n{cand['primer_4_seq']}\n")
    
    # Step 3: Run ONE BLAT call
    with tempfile.NamedTemporaryFile(mode='w', suffix='.psl', delete=False) as f:
        psl_output = f.name
    
    print(f"  Running BLAT on {total_primers} primers...")
    blat_start = time.time()
    blat_success = run_blat_batch(primer_fasta, genome_fasta_path, psl_output,
                                   args.blat_min_identity, args.blat_min_score)
    blat_time = time.time() - blat_start
    print(f"  BLAT completed in {blat_time:.1f}s ({total_primers/blat_time:.0f} primers/sec)")
    
    if not blat_success:
        print("  WARNING: BLAT failed, using first candidate without off-target check")
        for design in primer_designs:
            result = create_result_from_candidates(
                design['p1_candidates'][0], design['p4_candidates'][0],
                design['chrom'], design['strand'],
                design['p2_start'], design['p2_end'],
                design['p3_start'], design['p3_end'],
                'BLAT failed - no off-target check', 'BLAT failed - no off-target check'
            )
            batch_results.append((design['row_idx'], result))
        
        os.unlink(primer_fasta)
        os.unlink(psl_output)
        return batch_results
    
    # Step 4: Parse BLAT results
    blat_hits = parse_blat_psl(psl_output)
    print(f"  BLAT found hits for {len(blat_hits)} primers")
    
    # Step 5: Select best clean candidates for each row
    print(f"\nStep 3: Selecting best clean primer pair for each row...")
    
    clean_count = 0
    offtarget_count = 0
    
    for design in primer_designs:
        row_idx = design['row_idx']
        gene_name = get_gene_name(design['row'], row_idx)
        position = design['position']
        
        # Track best alternatives by category for fallback
        best_p1_clean = None
        best_p1_diff_chrom = None
        best_p1_same_chrom = None
        best_p1_near = None
        
        best_p4_clean = None
        best_p4_diff_chrom = None
        best_p4_same_chrom = None
        best_p4_near = None
        
        # Find best P1
        for cand_idx, cand in enumerate(design['p1_candidates']):
            primer_id = f"row{row_idx}_P1_cand{cand_idx}"
            has_ot, num_ot, details, category = check_offtargets(
                blat_hits, primer_id, design['chrom'],
                cand['primer_1_start'], cand['primer_1_end'], position, args
            )
            
            if not has_ot:
                best_p1_clean = (cand, 'No off-targets' if cand_idx == 0 else f'Clean (alternative #{cand_idx+1})')
                break
            else:
                # Save best fallback by category
                if category == 'diff_chrom' and not best_p1_diff_chrom:
                    best_p1_diff_chrom = (cand, f'Off-targets on different chromosome(s): {details}')
                elif category == 'same_chrom_far' and not best_p1_same_chrom:
                    best_p1_same_chrom = (cand, f'Off-targets on same chromosome >10kb away: {details}')
                elif category == 'same_chrom_near' and not best_p1_near:
                    best_p1_near = (cand, f'Off-targets within 10kb: {details}')
        
        # Select best P1 with priority: clean > diff_chrom > same_chrom > near
        if best_p1_clean:
            best_p1, p1_notes = best_p1_clean
        elif best_p1_diff_chrom:
            best_p1, p1_notes = best_p1_diff_chrom
        elif best_p1_same_chrom:
            best_p1, p1_notes = best_p1_same_chrom
        elif best_p1_near:
            best_p1, p1_notes = best_p1_near
        else:
            best_p1 = design['p1_candidates'][0]
            p1_notes = 'Off-targets present (best available)'
        
        # Find best P4
        for cand_idx, cand in enumerate(design['p4_candidates']):
            primer_id = f"row{row_idx}_P4_cand{cand_idx}"
            has_ot, num_ot, details, category = check_offtargets(
                blat_hits, primer_id, design['chrom'],
                cand['primer_4_start'], cand['primer_4_end'], position, args
            )
            
            if not has_ot:
                best_p4_clean = (cand, 'No off-targets' if cand_idx == 0 else f'Clean (alternative #{cand_idx+1})')
                break
            else:
                # Save best fallback by category
                if category == 'diff_chrom' and not best_p4_diff_chrom:
                    best_p4_diff_chrom = (cand, f'Off-targets on different chromosome(s): {details}')
                elif category == 'same_chrom_far' and not best_p4_same_chrom:
                    best_p4_same_chrom = (cand, f'Off-targets on same chromosome >10kb away: {details}')
                elif category == 'same_chrom_near' and not best_p4_near:
                    best_p4_near = (cand, f'Off-targets within 10kb: {details}')
        
        # Select best P4 with priority: clean > diff_chrom > same_chrom > near
        if best_p4_clean:
            best_p4, p4_notes = best_p4_clean
        elif best_p4_diff_chrom:
            best_p4, p4_notes = best_p4_diff_chrom
        elif best_p4_same_chrom:
            best_p4, p4_notes = best_p4_same_chrom
        elif best_p4_near:
            best_p4, p4_notes = best_p4_near
        else:
            best_p4 = design['p4_candidates'][0]
            p4_notes = 'Off-targets present (best available)'
        
        # Track stats
        p1_is_clean = 'No off-targets' in p1_notes or 'Clean' in p1_notes
        p4_is_clean = 'No off-targets' in p4_notes or 'Clean' in p4_notes
        
        if p1_is_clean and p4_is_clean:
            clean_count += 1
            print(f"  Row {row_idx+1:3d} ({gene_name:20s}): ✓ Clean")
        else:
            offtarget_count += 1
            print(f"  Row {row_idx+1:3d} ({gene_name:20s}): ⚠ Has off-targets")
        
        # Create result
        result = create_result_from_candidates(
            best_p1, best_p4,
            design['chrom'], design['strand'],
            design['p2_start'], design['p2_end'],
            design['p3_start'], design['p3_end'],
            p1_notes, p4_notes
        )
        batch_results.append((row_idx, result))
    
    print(f"\n  Summary: {clean_count} clean, {offtarget_count} with off-targets")
    
    # Step 6: Progressive flank widening for rows that couldn't find clean primers
    rows_needing_wider_flanks = []
    
    for design in primer_designs:
        row_idx = design['row_idx']
        # Check if this row was saved with off-target notes
        saved_result = None
        for saved_row_idx, saved_res in batch_results:
            if saved_row_idx == row_idx:
                saved_result = saved_res
                break
        
        if saved_result:
            p1_has_ot = saved_result.get('primer_1_notes', '') and 'Off-targets' in saved_result.get('primer_1_notes', '')
            p4_has_ot = saved_result.get('primer_4_notes', '') and 'Off-targets' in saved_result.get('primer_4_notes', '')
            if p1_has_ot or p4_has_ot:
                # Track which specific primers need rescue
                rows_needing_wider_flanks.append((design, saved_result, p1_has_ot, p4_has_ot))
    
    if rows_needing_wider_flanks and args.max_arm < args.max_arm_fallback:
        # Count how many need each primer
        total_need_p1 = sum(1 for _, _, p1_ot, _ in rows_needing_wider_flanks if p1_ot)
        total_need_p4 = sum(1 for _, _, _, p4_ot in rows_needing_wider_flanks if p4_ot)
        
        print(f"\n{'='*70}")
        print(f"STEP 4: WIDER FLANK RESCUE")
        print(f"{'='*70}")
        print(f"rows needing rescue: {len(rows_needing_wider_flanks)}")
        print(f"  • P1 needs rescue: {total_need_p1} rows")
        print(f"  • P4 needs rescue: {total_need_p4} rows")
        print(f"  • Both need rescue: {total_need_p1 + total_need_p4 - len(rows_needing_wider_flanks)} rows")
        print(f"\nTrying flanks: {args.max_arm + 500}bp to {args.max_arm_fallback}bp (500bp steps)")
        print(f"{'='*70}")
        
        widening_rescued = 0
        remaining_rows = rows_needing_wider_flanks
        
        for test_flank in range(args.max_arm + 500, args.max_arm_fallback + 1, 500):
            if not remaining_rows:
                break
            
            # Count needs for this iteration
            need_p1 = sum(1 for _, _, p1_ot, _ in remaining_rows if p1_ot)
            need_p4 = sum(1 for _, _, _, p4_ot in remaining_rows if p4_ot)
            
            print(f"\n  ┌─ Flank width: {test_flank}bp ─────────────────────")
            print(f"  │ Remaining: {len(remaining_rows)} rows (P1: {need_p1}, P4: {need_p4})")
            
            # Process all rows at this flank size
            rescued_this_flank = []
            flank_designs = []
            for design, old_result, needs_p1_rescue, needs_p4_rescue in remaining_rows:
                    row_idx = design['row_idx']
                    chrom = design['chrom']
                    position = design['position']
                    strand = design['strand']
                    p2 = design['p2']
                    p3 = design['p3']
                    
                    wider_upstream, wider_downstream, wider_up_start, wider_down_start = extract_flanks(
                        genome, chrom, position, strand, test_flank
                    )
                    
                    if wider_upstream is None:
                        continue
                    
                    # Only design primers that need rescue
                    p1_wider = None
                    p4_wider = None
                    
                    if needs_p1_rescue:
                        p1_wider, _ = design_outer_primer_for_p1p2_arm(
                            wider_upstream, p2, wider_up_start, strand, position, args,
                            target_arm_size=test_flank, num_return=50, debug=False
                        )
                    
                    if needs_p4_rescue:
                        p4_wider, _ = design_outer_primer_for_p3p4_arm(
                            wider_downstream, p3, wider_down_start, strand, position, args,
                            target_arm_size=test_flank, num_return=50, debug=False
                        )
                    
                    # Only add to flank_designs if at least one primer was designed
                    if p1_wider or p4_wider:
                        # Sort by GC-clamp preference
                        if p1_wider:
                            p1_wider.sort(key=lambda c: _gc_end_pref(c['primer_1_seq']))
                        if p4_wider:
                            p4_wider.sort(key=lambda c: _gc_end_pref(c['primer_4_seq']))
                        
                        flank_designs.append({
                            'design': design,
                            'old_result': old_result,
                            'row_idx': row_idx,
                            'p1_candidates': p1_wider if p1_wider else [],
                            'p4_candidates': p4_wider if p4_wider else [],
                            'chrom': chrom,
                            'position': position,
                            'needs_p1_rescue': needs_p1_rescue,
                            'needs_p4_rescue': needs_p4_rescue
                        })
                
            if not flank_designs:
                print(f"    No primers designed at {test_flank}bp")
                continue
            
            # ========================================
            # Process P1 primers separately
            # ========================================
            
            # Check if any rows need P1 rescue
            genes_needing_p1 = [fd for fd in flank_designs if fd.get('needs_p1_rescue', False)]
            
            best_p1_by_gene = {}  # row_idx -> (candidate, notes)
            
            if genes_needing_p1:
                # Create FASTA with only P1 candidates
                with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
                    p1_fasta = f.name
                    for flank_design in genes_needing_p1:
                        row_idx = flank_design['row_idx']
                        for cand_idx, cand in enumerate(flank_design['p1_candidates']):
                            f.write(f">row{row_idx}_flank{test_flank}_P1_{cand_idx}\n{cand['primer_1_seq']}\n")
            
                # Run BLAT for P1
                with tempfile.NamedTemporaryFile(mode='w', suffix='.psl', delete=False) as f:
                    p1_psl = f.name
                
                p1_count = sum(len(fd['p1_candidates']) for fd in genes_needing_p1)
                print(f"  │  Testing P1: {len(genes_needing_p1)} rows, {p1_count} primers")
                
                if not run_blat_batch(p1_fasta, genome_fasta_path, p1_psl,
                                     args.blat_min_identity, args.blat_min_score):
                    print(f"    BLAT failed for P1 at {test_flank}bp")
                    os.unlink(p1_fasta)
                    os.unlink(p1_psl)
                    continue
                
                p1_hits = parse_blat_psl(p1_psl)
                
                # Select best P1 for each gene that needs rescue
                for flank_design in genes_needing_p1:
                    row_idx = flank_design['row_idx']
                    position = flank_design['position']
                    
                    best_p1_clean = None
                    best_p1_fallback = None
                    
                    for cand_idx, cand in enumerate(flank_design['p1_candidates']):
                        primer_id = f"row{row_idx}_flank{test_flank}_P1_{cand_idx}"
                        has_ot, num_ot, details, category = check_offtargets(
                            p1_hits, primer_id, flank_design['chrom'],
                            cand['primer_1_start'], cand['primer_1_end'], position, args
                        )
                        
                        if not has_ot:
                            best_p1_clean = (cand, f'Wider flank rescue: {test_flank}bp')
                            break
                        elif category == 'diff_chrom' and not best_p1_fallback:
                            best_p1_fallback = (cand, f'Wider flank ({test_flank}bp) with off-targets on different chromosomes')
                    
                    if best_p1_clean:
                        best_p1_by_gene[row_idx] = best_p1_clean
                
                # Cleanup P1 files
                os.unlink(p1_fasta)
                os.unlink(p1_psl)
            else:
                print(f"  │ → P1: All clean (no rescue needed)")
            
            # ========================================
            # Process P4 primers separately
            # ========================================
            
            # Check if any rows need P4 rescue
            genes_needing_p4 = [fd for fd in flank_designs if fd.get('needs_p4_rescue', False)]
            
            best_p4_by_gene = {}  # row_idx -> (candidate, notes)
            
            if genes_needing_p4:
                # Create FASTA with only P4 candidates
                with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
                    p4_fasta = f.name
                    for flank_design in genes_needing_p4:
                        row_idx = flank_design['row_idx']
                        for cand_idx, cand in enumerate(flank_design['p4_candidates']):
                            f.write(f">row{row_idx}_flank{test_flank}_P4_{cand_idx}\n{cand['primer_4_seq']}\n")
                
                # Run BLAT for P4
                with tempfile.NamedTemporaryFile(mode='w', suffix='.psl', delete=False) as f:
                    p4_psl = f.name
                
                p4_count = sum(len(fd['p4_candidates']) for fd in genes_needing_p4)
                print(f"  │  Testing P4: {len(genes_needing_p4)} rows, {p4_count} primers")
                
                if not run_blat_batch(p4_fasta, genome_fasta_path, p4_psl,
                                     args.blat_min_identity, args.blat_min_score):
                    print(f"    BLAT failed for P4 at {test_flank}bp")
                    os.unlink(p4_fasta)
                    os.unlink(p4_psl)
                    continue
                
                p4_hits = parse_blat_psl(p4_psl)
                
                # Select best P4 for each gene that needs rescue
                for flank_design in genes_needing_p4:
                    row_idx = flank_design['row_idx']
                    position = flank_design['position']
                    
                    best_p4_clean = None
                    best_p4_fallback = None
                    
                    for cand_idx, cand in enumerate(flank_design['p4_candidates']):
                        primer_id = f"row{row_idx}_flank{test_flank}_P4_{cand_idx}"
                        has_ot, num_ot, details, category = check_offtargets(
                            p4_hits, primer_id, flank_design['chrom'],
                            cand['primer_4_start'], cand['primer_4_end'], position, args
                        )
                        
                        if not has_ot:
                            best_p4_clean = (cand, f'Wider flank rescue: {test_flank}bp')
                            break
                        elif category == 'diff_chrom' and not best_p4_fallback:
                            best_p4_fallback = (cand, f'Wider flank ({test_flank}bp) with off-targets on different chromosomes')
                    
                    if best_p4_clean:
                        best_p4_by_gene[row_idx] = best_p4_clean
            
                # Cleanup P4 files
                os.unlink(p4_fasta)
                os.unlink(p4_psl)
            else:
                print(f"  │ → P4: All clean (no rescue needed)")
            
            # ========================================
            # Combine results and update batch_results
            # ========================================
            
            for flank_design in flank_designs:
                row_idx = flank_design['row_idx']
                design = flank_design['design']
                gene_name = get_gene_name(design['row'], row_idx)
                
                # Get best P1 and P4 from the separate processing steps
                best_p1_clean = best_p1_by_gene.get(row_idx)
                best_p4_clean = best_p4_by_gene.get(row_idx)
            
                # Accept any clean primers found (independent rescue)
                # Get current result for this row to check which primers need rescue
                current_result = None
                for idx, res in batch_results:
                    if idx == row_idx:
                        current_result = res
                        break
            
                # Determine which primers to use (prioritize clean from wider flanks)
                use_p1 = None
                use_p1_notes = None
                use_p4 = None
                use_p4_notes = None
            
                # For P1: use clean from wider flank if available, else keep original
                if best_p1_clean:
                    p1_cand, p1_notes = best_p1_clean
                    use_p1 = p1_cand
                    use_p1_notes = p1_notes
                elif current_result:
                    # Keep original P1 (we'll extract from current_result below)
                    pass
            
                # For P4: use clean from wider flank if available, else keep original
                if best_p4_clean:
                    p4_cand, p4_notes = best_p4_clean
                    use_p4 = p4_cand
                    use_p4_notes = p4_notes
                elif current_result:
                    # Keep original P4 (we'll extract from current_result below)
                    pass
            
                # Only update if at least one primer improved to clean
                if best_p1_clean or best_p4_clean:
                    # If we need to keep an original primer, reconstruct candidate dict from current_result
                    if current_result and not use_p1:
                        # Reconstruct P1 candidate from current result
                        p1_start, p1_end = parse_locus_coords(current_result.get('primer_1_locus'))
                        use_p1 = {
                            'primer_1_seq': current_result['primer_1_seq'],
                            'primer_1_start': p1_start,
                            'primer_1_end': p1_end,
                            'primer_1_tm': current_result['primer_1_tm'],
                            'primer_1_gc': current_result['primer_1_gc']
                        }
                        use_p1_notes = current_result.get('primer_1_notes', '')
                
                    if current_result and not use_p4:
                        # Reconstruct P4 candidate from current result
                        p4_start, p4_end = parse_locus_coords(current_result.get('primer_4_locus'))
                        use_p4 = {
                            'primer_4_seq': current_result['primer_4_seq'],
                            'primer_4_start': p4_start,
                            'primer_4_end': p4_end,
                            'primer_4_tm': current_result['primer_4_tm'],
                            'primer_4_gc': current_result['primer_4_gc']
                        }
                        use_p4_notes = current_result.get('primer_4_notes', '')
                
                    # Create result with mix of new and/or original primers
                    if use_p1 and use_p4:
                        result = create_result_from_candidates(
                            use_p1, use_p4,
                            design['chrom'], design['strand'],
                            design['p2_start'], design['p2_end'],
                            design['p3_start'], design['p3_end'],
                            use_p1_notes, use_p4_notes
                        )
                    
                        # Update the batch_results (remove old result, add new one)
                        batch_results = [(idx, res) for idx, res in batch_results if idx != row_idx]
                        batch_results.append((row_idx, result))
                    
                        rescued_this_flank.append(design)
                        widening_rescued += 1
                    
                        # Show which primers were rescued
                        p1_symbol = "✓" if best_p1_clean else "·"
                        p4_symbol = "✓" if best_p4_clean else "·"
                        print(f"  │   [{p1_symbol}][{p4_symbol}] Row {row_idx+1:3d} {gene_name}")
            
            # Summary for this flank size
            if rescued_this_flank:
                # Count what was rescued
                p1_rescued = sum(1 for d in rescued_this_flank if d['row_idx'] in best_p1_by_gene)
                p4_rescued = sum(1 for d in rescued_this_flank if d['row_idx'] in best_p4_by_gene)
                print(f"  │")
                print(f"  └─ Rescued: {len(rescued_this_flank)} rows (P1: {p1_rescued}, P4: {p4_rescued})")
            else:
                print(f"  └─ No rows rescued at {test_flank}bp")
            
            # Update remaining_rows: remove fully rescued rows, update flags for partially rescued
            fully_rescued_designs = set()
            updated_flags = {}  # row_idx -> (new_p1_has_ot, new_p4_has_ot)
            
            for design in rescued_this_flank:
                row_idx = design['row_idx']
                # Check the updated result to see which primers are clean
                for idx, res in batch_results:
                    if idx == row_idx:
                        p1_clean = not (res.get('primer_1_notes', '') and 'Off-targets' in res.get('primer_1_notes', ''))
                        p4_clean = not (res.get('primer_4_notes', '') and 'Off-targets' in res.get('primer_4_notes', ''))
                        
                        if p1_clean and p4_clean:
                            # Both clean - remove from queue
                            fully_rescued_designs.add(row_idx)
                        else:
                            # Partially rescued - update flags
                            updated_flags[row_idx] = (not p1_clean, not p4_clean)
                        break
            
            # Rebuild remaining_rows with updated flags
            new_remaining_rows = []
            for design, old_result, p1_has_ot, p4_has_ot in remaining_rows:
                row_idx = design['row_idx']
                
                # Skip if fully rescued
                if row_idx in fully_rescued_designs:
                    continue
                
                # Update flags if partially rescued
                if row_idx in updated_flags:
                    p1_has_ot, p4_has_ot = updated_flags[row_idx]
                
                new_remaining_rows.append((design, old_result, p1_has_ot, p4_has_ot))
            
            remaining_rows = new_remaining_rows
        
        # Final summary
        print(f"\n{'='*70}")
        print(f"RESCUE COMPLETE")
        print(f"{'='*70}")
        print(f" Rescued: {widening_rescued} rows")
        
        if remaining_rows:
            still_need_p1 = sum(1 for _, _, p1_ot, _ in remaining_rows if p1_ot)
            still_need_p4 = sum(1 for _, _, _, p4_ot in remaining_rows if p4_ot)
            print(f"⚠ Remaining with off-targets: {len(remaining_rows)} rows")
            print(f"  • P1 still has off-targets: {still_need_p1}")
            print(f"  • P4 still has off-targets: {still_need_p4}")
            print(f"  Using best available primers (with off-target notes)")
        else:
            print(f" All rows rescued successfully!")
        print(f"{'='*70}")
    
    # Cleanup
    os.unlink(primer_fasta)
    os.unlink(psl_output)
    
    return batch_results


def main():
    parser = argparse.ArgumentParser(
        description='Design outer primers (P1, P4) with batch processing',
        epilog="Batch processing dramatically improves performance by running one BLAT call per batch"
    )
    parser.add_argument('--input', required=True, help='Input CSV')
    parser.add_argument('--output', required=True, help='Output CSV')
    parser.add_argument('--genome', required=True, help='Genome FASTA')
    parser.add_argument('--batch-size', type=int, default=100, 
                       help='Process N rows per batch (default: 100)')
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
                       help='Maximum homology arm length to try before giving up (default: 10000)')
    
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
    
    if not os.path.exists(args.input):
        print(f"ERROR: Input CSV not found: {args.input}")
        sys.exit(1)
    
    # Check BLAT availability
    genome_fasta_for_blat = None
    if not args.skip_offtarget_check:
        try:
            subprocess.run(['blat'], capture_output=True, text=True)
            genome_fasta_for_blat = args.genome
            print(" BLAT found - will check for off-targets")
        except FileNotFoundError:
            print("ERROR: BLAT not found in PATH")
            print("  Install BLAT or use --skip-offtarget-check to proceed without validation")
            sys.exit(1)
    else:
        print("WARNING: Skipping off-target checks (--skip-offtarget-check)")
    
    # Load genome
    genome = load_genome(args.genome)
    
    # Load input CSV
    print(f"\nReading {args.input}...")
    df = pd.read_csv(args.input)
    print(f"Loaded {len(df)} rows")
    
    # Validate that input has required primers from script 5
    required_cols = ['chrom', 'pos', 'strand', 'primer_2_seq_genomic', 'primer_3_seq_genomic']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"ERROR: Input CSV missing required columns: {missing_cols}")
        print(f"  This script requires output from script 5 (inner primers)")
        sys.exit(1)
    
    # Check how many rows have inner primers
    has_p2 = df['primer_2_seq_genomic'].notna() & (df['primer_2_seq_genomic'] != 'ERROR')
    has_p3 = df['primer_3_seq_genomic'].notna() & (df['primer_3_seq_genomic'] != 'ERROR')
    rows_with_primers = (has_p2 & has_p3).sum()
    
    print(f"  Rows with inner primers (P2 and P3): {rows_with_primers}/{len(df)}")
    if rows_with_primers == 0:
        print(f"ERROR: No rows have inner primers!")
        print(f"  This script requires primer_2_seq_genomic and primer_3_seq_genomic to be filled.")
        print(f"  Please run script 5 first to design inner primers.")
        sys.exit(1)
    elif rows_with_primers < len(df):
        print(f"  Will skip {len(df) - rows_with_primers} rows without inner primers.")
    
    print(f"\nConfiguration:")
    print(f"  Homology arm length: {args.min_arm}-{args.max_arm} bp")
    print(f"  Primer size range: {args.primer_min_size}-{args.primer_max_size} bp")
    print(f"  Primer Tm range: {args.primer_min_tm}-{args.primer_max_tm}°C")
    print(f"  BLAT parameters: min identity={args.blat_min_identity}%, min score={args.blat_min_score}")
    print(f"  Batch size: {args.batch_size} rows")
    
    # Check for existing output and resume
    start_idx = 0
    if args.no_resume or not os.path.exists(args.output):
        print(f"\n[start] Creating new output file: {args.output}")
    else:
        print(f"\n[resume] Found existing output: {args.output}")
        try:
            existing = pd.read_csv(args.output)
            if 'primer_1_seq' in existing.columns:
                # Find the last row that has BOTH primer_1 and primer_4 successfully designed
                completed = (existing['primer_1_seq'].notna() & 
                           (existing['primer_1_seq'] != 'ERROR') &
                           existing['primer_4_seq'].notna() & 
                           (existing['primer_4_seq'] != 'ERROR'))
                if completed.any():
                    # Find last completed row
                    last_completed = completed[::-1].idxmax()
                    start_idx = last_completed + 1
                    print(f"[resume] Last completed row: {last_completed + 1}")
                    print(f"[resume] Resuming from row {start_idx + 1}")
                else:
                    print(f"[resume] No completed rows found, starting from beginning")
            else:
                print(f"[resume] No primer_1_seq column found, starting from beginning")
        except Exception as e:
            print(f"[resume] Could not read existing output: {e}")
            print(f"[resume] Starting from beginning")
    
    # Pre-allocate results
    all_results = [None] * len(df)
    
    # Load existing results if resuming
    if start_idx > 0:
        try:
            existing = pd.read_csv(args.output)
            for i in range(start_idx):
                all_results[i] = {
                    col: existing.loc[i, col] if col in existing.columns else None
                    for col in ['primer_1_seq', 'primer_1_locus',
                               'primer_1_tm', 'primer_1_gc', 'upstream_arm_length', 'primer_1_notes',
                               'primer_4_seq', 'primer_4_locus',
                               'primer_4_tm', 'primer_4_gc', 'downstream_arm_length', 'primer_4_notes']
                }
        except:
            pass
    
    print(f"\nProcessing rows {start_idx + 1} to {len(df)}...")
    print(f"Saving every {args.batch_size} rows\n")
    
    # Process in batches
    idx = start_idx
    while idx < len(df):
        batch_end = min(idx + args.batch_size, len(df))
        batch_rows = [(i, df.iloc[i]) for i in range(idx, batch_end)]
        
        t0 = time.perf_counter() # start timer

        # Process batch
        if args.skip_offtarget_check:
            # Simple processing without off-target check
            batch_results = []
            for row_idx, row in batch_rows:
                chrom = row.get('chrom', '')
                position = row.get('pos', None)
                
                # Handle NaN values in primer sequences
                p2_raw = row.get('primer_2_seq_genomic', '')
                p3_raw = row.get('primer_3_seq_genomic', '')
                p2 = str(p2_raw).upper() if pd.notna(p2_raw) else ''
                p3 = str(p3_raw).upper() if pd.notna(p3_raw) else ''
                
                strand = row.get('strand', '+')
                p2_locus = row.get('primer_2_locus', '')
                p3_locus = row.get('primer_3_locus', '')
                
                # Skip rows without inner primers
                if not p2 or not p3:
                    continue
                
                if not chrom or position is None:
                    result = {
                        'primer_1_seq': 'ERROR', 'primer_1_locus': None,
                        'primer_1_tm': None, 'primer_1_gc': None, 'upstream_arm_length': None,
                        'primer_1_notes': 'Missing required fields',
                        'primer_4_seq': 'ERROR', 'primer_4_locus': None,
                        'primer_4_tm': None, 'primer_4_gc': None, 'downstream_arm_length': None,
                        'primer_4_notes': 'Missing required fields'
                    }
                else:
                    try:
                        position = int(position)
                        upstream_seq, downstream_seq, upstream_start, downstream_start = extract_flanks(
                            genome, chrom, position, strand, args.max_arm
                        )
                        
                        if upstream_seq is None:
                            result = {
                                'primer_1_seq': 'ERROR', 'primer_1_notes': 'Chromosome not found',
                                'primer_4_seq': 'ERROR', 'primer_4_notes': 'Chromosome not found'
                            }
                        else:
                            p1_cands, _ = design_outer_primer_for_p1p2_arm(
                                upstream_seq, p2, upstream_start, strand, position, args,
                                target_arm_size=args.max_arm, num_return=1
                            )
                            p4_cands, _ = design_outer_primer_for_p3p4_arm(
                                downstream_seq, p3, downstream_start, strand, position, args,
                                target_arm_size=args.max_arm, num_return=1
                            )
                            
                            if p1_cands and p4_cands:
                                p2_start, p2_end = parse_locus_coords(p2_locus)
                                p3_start, p3_end = parse_locus_coords(p3_locus)
                                result = create_result_from_candidates(
                                    p1_cands[0], p4_cands[0], chrom, strand,
                                    p2_start, p2_end, p3_start, p3_end,
                                    'No off-target check performed', 'No off-target check performed'
                                )
                            else:
                                result = {
                                    'primer_1_seq': 'ERROR', 'primer_1_notes': 'Primer3 failed',
                                    'primer_4_seq': 'ERROR', 'primer_4_notes': 'Primer3 failed'
                                }
                    except Exception as e:
                        result = {
                            'primer_1_seq': 'ERROR', 'primer_1_notes': str(e),
                            'primer_4_seq': 'ERROR', 'primer_4_notes': str(e)
                        }
                
                batch_results.append((row_idx, result))
        else:
            # Process with off-target checking
            batch_results = process_batch(
                batch_rows, genome, genome_fasta_for_blat, idx, args
            )
        
        # Store results
        for row_idx, result in batch_results:
            all_results[row_idx] = result
        
        # Save checkpoint
        print(f"\n{'='*70}")
        print(f"CHECKPOINT: Saving results for rows {idx+1}-{batch_end}")
        print(f"{'='*70}")

        dt = time.perf_counter() - t0
        rate = (args.batch_size / dt) if dt > 0 else float('inf')
        print(f"[timing] Elapsed {dt:.2f}s | {args.batch_size} rows | {rate:.2f} rows/s", file=sys.stderr)
        
        # Fill missing results with empty dicts
        save_results = []
        for r in all_results:
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
        print(f" Saved to {args.output}")
        
        # Stats
        completed = [r for r in all_results if r is not None]
        success = sum(1 for r in completed 
                    if r['primer_1_seq'] not in ['ERROR', None] 
                    and r['primer_4_seq'] not in ['ERROR', None])
        print(f"Success: {success}/{len(completed)} ({100*success/len(completed):.1f}%)\n")
        
        idx = batch_end
    
    # Final summary
    print(f"\n{'='*70}")
    print(f"FINAL SUMMARY")
    print(f"{'='*70}")
    completed = [r for r in all_results if r is not None]
    success = sum(1 for r in completed 
                if r['primer_1_seq'] not in ['ERROR', None] 
                and r['primer_4_seq'] not in ['ERROR', None])
    partial = sum(1 for r in completed 
                 if (r['primer_1_seq'] not in ['ERROR', None]) != (r['primer_4_seq'] not in ['ERROR', None]))
    
    print(f"Total: {len(completed)}")
    print(f"Both primers: {success} ({100*success/len(completed):.1f}%)")
    print(f"One primer: {partial}")
    print(f"Both failed: {len(completed) - success - partial}")
    print(f"{'='*70}\n")

if __name__ == '__main__':
    main()