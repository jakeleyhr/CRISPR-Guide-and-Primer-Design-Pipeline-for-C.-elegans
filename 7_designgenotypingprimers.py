#!/usr/bin/env python3
"""
Homology Arm Primer Design Pipeline - Step 4: Genotyping Primer Design

Designs genotyping primers to verify successful CRISPR knock-in by amplifying
the entire region of interest (ROI) spanning both homology arms. Uses Primer3
for optimal primer design with automatic off-target checking via BLAT.

Input:
  - CSV file from script 6 containing outer primers (P1 and P4) with genomic
    coordinates that define the boundaries of the homology arm region
    (requires: primer_1_locus, primer_4_locus)

Output:
  - CSV file with all input columns plus genotyping primer information:
      • genotyping_primer_left_seq: Forward primer upstream of ROI
      • genotyping_primer_right_seq: Reverse primer downstream of ROI
      • genotyping_primer_left_locus, genotyping_primer_right_locus: Genomic coordinates
      • genotyping_primer_left_tm, genotyping_primer_right_tm: Melting temperatures
      • genotyping_primer_left_gc, genotyping_primer_right_gc: GC content percentages
      • amplicon_size: Expected PCR product size
      • genotyping_primers_notes: Off-target descriptions and design notes

Primer Design Strategy:
  1. Extract flanking sequences (1000bp each side of ROI defined by P1-P4)
  2. Design standard primers using Primer3:
     - Target primer size: 20-30bp (optimal 25bp)
     - Optimal Tm: ~58°C (range 54-66°C)
     - Preferentially chooses primer sequences with 3' GC-clamp
     - Targets ROI region for amplification
  3. Batch off-target checking: ALL alternatives checked in ONE BLAT call
  4. Select best clean primer pair from pre-checked alternatives (prioritizes GC-clamp)

Required arguments:
  --input PATH       Input CSV from script 6 with outer primers (P1 and P4)
  --output PATH      Output CSV path with genotyping primers
  --genome PATH      Genome FASTA file

Optional arguments:
  General:
    --batch-size N               Process and write results every N rows for progress tracking (default: 100)
    --no-resume                  Overwrite output instead of resuming from checkpoint
    --skip-offtarget-check       Skip off-target search using BLAT (much faster but less safe)
  
  Flank and rescue strategy:
    --initial-flank N            Initial flanking size in bp (default: 1000)
    --max-flank N                Maximum flanking size in bp before desperation mode (default: 5000)
    --num-alternatives N         Number of primer alternatives to generate (default: 50)
    --max-offtarget-attempts N   Maximum alternatives to try if off-targets found (default: 50)
    --max-rescue-batches N       Number of batches for mix-and-match rescue (default: 1)
    --num-single-primers N       New single primers to generate for matching (default: 1000)
    --num-single-primers-search N  Search space for single primer generation (default: 10000)
  
  BLAT parameters:
    --blat-min-identity N        Minimum % identity for BLAT hits (default: 90)
    --blat-min-score N           Minimum alignment score for BLAT (default: 20)
  
  Primer3 parameters:
    --primer-opt-size N          Optimal primer length (default: 25)
    --primer-min-size N          Minimum primer length (default: 20)
    --primer-max-size N          Maximum primer length (default: 30)
    --primer-opt-tm N            Optimal primer Tm in °C (default: 58.0)
    --primer-min-tm N            Minimum primer Tm in °C (default: 54.0)
    --primer-max-tm N            Maximum primer Tm in °C (default: 66.0)
    --primer-min-gc N            Minimum primer GC% (default: 35.0)
    --primer-max-gc N            Maximum primer GC% (default: 65.0)

Example:
  python 7_designgenotypingprimers.py \
    --input 6_guides_with_all_primers.csv \
    --output 7_guides_with_genotyping_primers.csv \
    --genome wbcel235/ce11.fa \
    --batch-size 50

  With custom rescue parameters for difficult targets:
  python 7_designgenotypingprimers.py \
    --input input.csv \
    --output output.csv \
    --genome genome.fa \
    --num-alternatives 100 \
    --initial-flank 2000 \
    --max-flank 8000

  Recommended: install BLAT (for off-target checking)
    - mkdir -p ~/bin  # Create bin directory if it doesn't exist
    - cd ~/bin
    - curl -o blat http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/blat/blat
    - chmod +x blat  # Make it executable
    - echo 'export PATH="$HOME/bin:$PATH"' >> ~/.zshrc  # Add to PATH permanently
    - blat  # Test
"""


import argparse, sys, os, re, subprocess, tempfile
import pandas as pd
from collections import defaultdict

try:
    import primer3
except ImportError:
    print("ERROR: primer3-py not installed! Install with: pip install primer3-py")
    sys.exit(1)

try:
    from Bio import SeqIO
    #from Bio.Seq import Seq
except ImportError:
    print("ERROR: biopython not installed! Install with: biopython")
    sys.exit(1)
    

def build_p3_args(args, num_return=None):
    """Build Primer3 arguments dictionary from command-line args"""
    if num_return is None:
        num_return = args.num_alternatives
    
    return {
        'PRIMER_OPT_SIZE': args.primer_opt_size,
        'PRIMER_MIN_SIZE': args.primer_min_size,
        'PRIMER_MAX_SIZE': args.primer_max_size,
        'PRIMER_OPT_TM': args.primer_opt_tm,
        'PRIMER_MIN_TM': args.primer_min_tm,
        'PRIMER_MAX_TM': args.primer_max_tm,
        'PRIMER_MIN_GC': args.primer_min_gc,
        'PRIMER_MAX_GC': args.primer_max_gc,
        'PRIMER_NUM_RETURN': num_return,
    }

def check_blat_available():
    """Check if BLAT is available in PATH"""
    try:
        result = subprocess.run(['blat'], capture_output=True, text=True)
        return True
    except FileNotFoundError:
        return False

def _gc_end_pref(seq):
    """Prefer primers that end in G or C at the 3' end (returns 0 for G/C, 1 for A/T)"""
    return 0 if seq[-1].upper() in ("G", "C") else 1

def select_best_offtarget_alternative(all_alternatives_with_ot):
    """
    Select the alternative with the fewest off-targets.
    Returns (primer_pair, notes) tuple
    """
    if not all_alternatives_with_ot:
        return None, None
    
    # Sort by: total off-targets (ascending), then by alternative index (prefer earlier)
    sorted_alts = sorted(all_alternatives_with_ot, key=lambda x: (x['total_ot'], x['alt_idx']))
    best = sorted_alts[0]
    
    notes = f'Best alternative with off-targets (#{best["alt_idx"]+1}): left={best["left_details"]} ({best["left_num_ot"]} OT), right={best["right_details"]} ({best["right_num_ot"]} OT)'
    
    return best['primer_pair'], notes

def get_gene_name(row, row_idx):
    """
    Safely extract gene name from row, handling NaN and non-string values.
    Returns string gene name or default.
    """
    gene = row.get('gene')
    if gene is None or (isinstance(gene, float) and pd.isna(gene)):
        return f'row{row_idx+1}'
    return str(gene)

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

def parse_locus(locus_str):
    """
    Parse locus string like 'I:5107080-5107099 (+)' or 'IV:9246353-9246387 (+)'
    Returns: (chrom, start, end, strand) or (None, None, None, None)
    """
    if pd.isna(locus_str) or locus_str is None or locus_str == '':
        return None, None, None, None
    
    # Pattern: CHROM:START-END (STRAND)
    pattern = r'([IVX]+):(\d+)-(\d+)\s*\(([+-])\)'
    match = re.match(pattern, str(locus_str))
    
    if match:
        chrom, start, end, strand = match.groups()
        return chrom, int(start), int(end), strand
    
    return None, None, None, None

def get_roi_span(primer_1_locus, primer_4_locus):
    """
    Get the min and max coordinates from primer_1_locus and primer_4_locus
    Returns: (chrom, roi_start, roi_end, strand) or (None, None, None, None)
    """
    # Parse both loci
    chrom1, start1, end1, strand1 = parse_locus(primer_1_locus)
    chrom4, start4, end4, strand4 = parse_locus(primer_4_locus)
    
    # Validate
    if chrom1 is None or chrom4 is None:
        return None, None, None, None
    
    if chrom1 != chrom4:
        print(f"  WARNING: P1 and P4 on different chromosomes: {chrom1} vs {chrom4}")
        return None, None, None, None
    
    # Get min and max coordinates
    all_coords = [start1, end1, start4, end4]
    roi_start = min(all_coords)
    roi_end = max(all_coords)
    
    # Use gene strand from primer_1
    return chrom1, roi_start, roi_end, strand1

def extract_flanking_sequence(genome, chrom, roi_start, roi_end, flank):
    """
    Extract sequence with flanking regions around ROI
    Returns: (full_seq, upstream_start, downstream_end) or (None, None, None)
    """
    if chrom not in genome:
        print(f"  ERROR: Chromosome {chrom} not in genome")
        return None, None, None
    
    seq = genome[chrom]
    chrom_len = len(seq)
    
    # Calculate extraction coordinates (1-based, inclusive)
    extract_start = max(1, roi_start - flank)
    extract_end = min(chrom_len, roi_end + flank)
    
    # Convert to 0-based for Python slicing
    extract_start_0 = extract_start - 1
    extract_end_0 = extract_end
    
    full_seq = seq[extract_start_0:extract_end_0]
    
    return full_seq, extract_start, extract_end

def design_verification_primers_multi(genome, chrom, roi_start, roi_end, gene_strand, flank_size, num_return=50):
    """
    Design multiple pairs of standard primers in the flanking regions to amplify the ROI
    
    Args:
        num_return: Number of primer alternatives to return (default: 50, override with args.num_alternatives)
    
    Returns list of dicts, one per primer pair alternative
    """
    results = []
    
    # Extract sequence with flanks
    full_seq, extract_start, extract_end = extract_flanking_sequence(
        genome, chrom, roi_start, roi_end, flank_size
    )
    
    if full_seq is None:
        return []
    
    # Calculate where ROI is within our extracted sequence (0-based)
    roi_start_in_seq = roi_start - extract_start
    roi_end_in_seq = roi_end - extract_start
    roi_len = roi_end - roi_start + 1
    
    # Define target region for primer pair (the ROI itself)
    target = [[roi_start_in_seq, roi_len]]
    
    # Set product size range
    min_product = roi_len
    max_product = roi_len + 2 * flank_size
    
    # Update primer3 args with product size range
    p3_args = P3_ARGS.copy()
    p3_args['PRIMER_PRODUCT_SIZE_RANGE'] = [[min_product, max_product]]
    p3_args['PRIMER_NUM_RETURN'] = num_return
    
    # Run Primer3
    try:
        primer3_result = primer3.design_primers(
            {
                'SEQUENCE_ID': f'{chrom}:{extract_start}-{extract_end}',
                'SEQUENCE_TEMPLATE': full_seq,
                'SEQUENCE_TARGET': target,
            },
            p3_args
        )
    except Exception as e:
        return []
    
    # Check for primer pairs
    num_pairs = primer3_result.get('PRIMER_PAIR_NUM_RETURNED', 0)
    
    if num_pairs == 0:
        return []
    
    # Extract all primer pairs
    for pair_idx in range(num_pairs):
        try:
            # Get left primer info
            left_seq = primer3_result[f'PRIMER_LEFT_{pair_idx}_SEQUENCE']
            left_start = primer3_result[f'PRIMER_LEFT_{pair_idx}'][0]
            left_tm = primer3_result[f'PRIMER_LEFT_{pair_idx}_TM']
            left_gc = primer3_result[f'PRIMER_LEFT_{pair_idx}_GC_PERCENT']
            
            # Get right primer info
            right_seq = primer3_result[f'PRIMER_RIGHT_{pair_idx}_SEQUENCE']
            right_pos = primer3_result[f'PRIMER_RIGHT_{pair_idx}'][0]
            right_tm = primer3_result[f'PRIMER_RIGHT_{pair_idx}_TM']
            right_gc = primer3_result[f'PRIMER_RIGHT_{pair_idx}_GC_PERCENT']
            
            # Product size
            amplicon_size = primer3_result[f'PRIMER_PAIR_{pair_idx}_PRODUCT_SIZE']
            
            # Calculate genomic coordinates (1-based, inclusive)
            left_start_genomic = extract_start + left_start
            left_end_genomic = left_start_genomic + len(left_seq) - 1
            
            right_end_genomic = extract_start + right_pos
            right_start_genomic = right_end_genomic - len(right_seq) + 1
            
            # Format locus strings
            left_locus = f"{chrom}:{left_start_genomic}-{left_end_genomic} ({gene_strand})"
            right_locus = f"{chrom}:{right_start_genomic}-{right_end_genomic} ({gene_strand})"
            
            results.append({
                'left_seq': left_seq,
                'left_locus': left_locus,
                'left_tm': round(left_tm, 1),
                'left_gc': round(left_gc, 1),
                'left_start_in_extract': left_start,
                'right_seq': right_seq,
                'right_locus': right_locus,
                'right_tm': round(right_tm, 1),
                'right_gc': round(right_gc, 1),
                'right_start_in_extract': right_start_genomic - extract_start,
                'amplicon_size': amplicon_size,
                'extract_start': extract_start,
                'chrom': chrom
            })
        except KeyError:
            continue
    
    # Sort by GC-clamp preference: primers ending in G/C preferred, then by Primer3's order
    results.sort(key=lambda r: (_gc_end_pref(r['left_seq']), _gc_end_pref(r['right_seq'])))
    
    return results

def run_blat_batch(primer_fasta_path, genome_fasta_path, output_psl_path, min_identity=None, min_score=None):
    """
    Run BLAT to search primers against genome
    Returns True if successful, False otherwise
    """
    if min_identity is None:
        min_identity = BLAT_MIN_IDENTITY
    if min_score is None:
        min_score = BLAT_MIN_SCORE
        
    cmd = [
        'blat',
        genome_fasta_path,
        primer_fasta_path,
        output_psl_path,
        '-stepSize=5',
        '-repMatch=1000000',
        f'-minScore={min_score}',
        f'-minIdentity={min_identity}',
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            print(f"    BLAT error: {result.stderr}")
            return False
        return True
    except subprocess.TimeoutExpired:
        print(f"    BLAT timeout after 300s")
        return False
    except Exception as e:
        print(f"    BLAT exception: {e}")
        return False

def parse_blat_psl(psl_path):
    """
    Parse BLAT PSL output to find off-target hits
    
    Returns dict: {primer_id: [hit_dicts, ...]}
    """
    hits = defaultdict(list)
    
    if not os.path.exists(psl_path):
        return hits
    
    with open(psl_path, 'r') as f:
        for line in f:
            # Skip header lines
            if line.startswith('psLayout') or line.startswith('match') or line.startswith('---'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 17:
                continue
            
            try:
                matches = int(parts[0])
                mismatches = int(parts[1])
                qName = parts[9]
                qSize = int(parts[10])
                tName = parts[13]
                tStart = int(parts[15])
                tEnd = int(parts[16])
                strand = parts[8]
                
                total_aligned = matches + mismatches
                if total_aligned == 0:
                    continue
                
                identity = 100.0 * matches / total_aligned
                score = matches
                
                hits[qName].append({
                    'chrom': tName,
                    'start': tStart + 1,
                    'end': tEnd,
                    'strand': strand,
                    'identity': identity,
                    'score': score,
                    'matches': matches,
                    'query_size': qSize
                })
                
            except (ValueError, IndexError):
                continue
    
    return hits

def check_offtargets(hits, primer_id, expected_chrom, expected_start, expected_end, tolerance=5):
    """
    Check if a primer has off-target hits
    
    Returns: (has_offtarget, num_offtargets, details)
    """
    if primer_id not in hits:
        return False, 0, "No BLAT hits found"
    
    primer_hits = hits[primer_id]
    
    # Filter for high-quality hits
    good_hits = [h for h in primer_hits 
                 if h['identity'] >= BLAT_MIN_IDENTITY and h['score'] >= BLAT_MIN_SCORE]
    
    if not good_hits:
        return False, 0, "No high-quality hits"
    
    # Check if hits are only at the expected location
    offtargets = []
    for hit in good_hits:
        is_expected = (
            hit['chrom'] == expected_chrom and
            abs(hit['start'] - expected_start) <= tolerance and
            abs(hit['end'] - expected_end) <= tolerance
        )
        
        if not is_expected:
            offtargets.append(hit)
    
    if offtargets:
        details = []
        for ot in offtargets[:3]:
            details.append(f"{ot['chrom']}:{ot['start']}-{ot['end']} ({ot['identity']:.1f}%)")
        detail_str = "; ".join(details)
        if len(offtargets) > 3:
            detail_str += f" and {len(offtargets)-3} more"
        
        return True, len(offtargets), detail_str
    
    return False, 0, "Only on-target hit found"

def format_offtarget_details_for_log(hits, primer_id, primer_seq):
    """
    Format detailed off-target information for failure logging
    Returns dict with all off-target locations
    """
    if primer_id not in hits:
        return {
            'primer_seq': primer_seq,
            'num_offtargets': 0,
            'offtarget_locations': 'No hits found',
            'offtarget_identities': '',
            'offtarget_scores': ''
        }
    
    primer_hits = hits[primer_id]
    
    # Filter for high-quality hits
    good_hits = [h for h in primer_hits 
                 if h['identity'] >= BLAT_MIN_IDENTITY and h['score'] >= BLAT_MIN_SCORE]
    
    if not good_hits:
        return {
            'primer_seq': primer_seq,
            'num_offtargets': 0,
            'offtarget_locations': 'No high-quality hits',
            'offtarget_identities': '',
            'offtarget_scores': ''
        }
    
    # Format all hits
    locations = []
    identities = []
    scores = []
    
    for hit in good_hits:
        locations.append(f"{hit['chrom']}:{hit['start']}-{hit['end']}")
        identities.append(f"{hit['identity']:.1f}%")
        scores.append(str(hit['score']))
    
    return {
        'primer_seq': primer_seq,
        'num_offtargets': len(good_hits),
        'offtarget_locations': '; '.join(locations),
        'offtarget_identities': '; '.join(identities),
        'offtarget_scores': '; '.join(scores)
    }


def design_single_direction_primers(genome, chrom, roi_start, roi_end, gene_strand, 
                                     direction='left', target_region_start=None, 
                                     target_region_end=None, num_return=5):
    """
    Design primers in only one direction (left or right) within a specific region
    
    Returns: list of primer dicts with 'seq', 'locus', 'tm', 'gc', 'start', 'end', 'chrom'
    """
    results = []
    
    # If no target region specified, use flanking regions
    if target_region_start is None or target_region_end is None:
        if direction == 'left':
            target_region_start = roi_start - FLANK_SIZE
            target_region_end = roi_start - 1
        else:
            target_region_start = roi_end + 1
            target_region_end = roi_end + FLANK_SIZE
    
    if chrom not in genome:
        return []
    
    seq = genome[chrom]
    chrom_len = len(seq)
    
    extract_start = max(1, target_region_start)
    extract_end = min(chrom_len, target_region_end)
    
    if extract_start >= extract_end:
        return []
    
    extract_start_0 = extract_start - 1
    extract_end_0 = extract_end
    target_seq = seq[extract_start_0:extract_end_0]
    
    p3_args = P3_ARGS.copy()
    p3_args['PRIMER_NUM_RETURN'] = num_return
    
    included_region = [[0, len(target_seq)]]
    
    try:
        if direction == 'left':
            p3_args['PRIMER_PICK_LEFT_PRIMER'] = 1
            p3_args['PRIMER_PICK_RIGHT_PRIMER'] = 0
            p3_args['PRIMER_PICK_INTERNAL_OLIGO'] = 0
        else:
            p3_args['PRIMER_PICK_LEFT_PRIMER'] = 0
            p3_args['PRIMER_PICK_RIGHT_PRIMER'] = 1
            p3_args['PRIMER_PICK_INTERNAL_OLIGO'] = 0
        
        primer3_result = primer3.design_primers(
            {
                'SEQUENCE_ID': f'{chrom}:{extract_start}-{extract_end}',
                'SEQUENCE_TEMPLATE': target_seq,
                'SEQUENCE_INCLUDED_REGION': included_region,
            },
            p3_args
        )
    except Exception as e:
        return []
    
    # Extract primers
    if direction == 'left':
        num_primers = primer3_result.get('PRIMER_LEFT_NUM_RETURNED', 0)
        for i in range(num_primers):
            try:
                seq_primer = primer3_result[f'PRIMER_LEFT_{i}_SEQUENCE']
                pos = primer3_result[f'PRIMER_LEFT_{i}'][0]
                tm = primer3_result[f'PRIMER_LEFT_{i}_TM']
                gc = primer3_result[f'PRIMER_LEFT_{i}_GC_PERCENT']
                
                start_genomic = extract_start + pos
                end_genomic = start_genomic + len(seq_primer) - 1
                locus = f"{chrom}:{start_genomic}-{end_genomic} ({gene_strand})"
                
                results.append({
                    'seq': seq_primer,
                    'locus': locus,
                    'tm': round(tm, 1),
                    'gc': round(gc, 1),
                    'start': start_genomic,
                    'end': end_genomic,
                    'chrom': chrom
                })
            except KeyError:
                continue
    else:
        num_primers = primer3_result.get('PRIMER_RIGHT_NUM_RETURNED', 0)
        for i in range(num_primers):
            try:
                seq_primer = primer3_result[f'PRIMER_RIGHT_{i}_SEQUENCE']
                pos = primer3_result[f'PRIMER_RIGHT_{i}'][0]
                tm = primer3_result[f'PRIMER_RIGHT_{i}_TM']
                gc = primer3_result[f'PRIMER_RIGHT_{i}_GC_PERCENT']
                
                end_genomic = extract_start + pos
                start_genomic = end_genomic - len(seq_primer) + 1
                locus = f"{chrom}:{start_genomic}-{end_genomic} ({gene_strand})"
                
                results.append({
                    'seq': seq_primer,
                    'locus': locus,
                    'tm': round(tm, 1),
                    'gc': round(gc, 1),
                    'start': start_genomic,
                    'end': end_genomic,
                    'chrom': chrom
                })
            except KeyError:
                continue
    
    # Sort by GC-clamp preference: primers ending in G/C preferred
    results.sort(key=lambda r: _gc_end_pref(r['seq']))
    
    return results

def process_batch_with_offtarget_check(rows_data, genome, genome_fasta_path, batch_start_idx, 
                                        initial_flank, max_flank, args):
    """
    OPTIMIZED: Process a batch of rows by checking ALL primer alternatives upfront in ONE BLAT call
    
    Progressive flank widening: Tries initial_flank, then increases by 1kb increments up to max_flank
    before attempting desperation mode with relaxed parameters.
    
    Args:
        args: Parsed command-line arguments containing rescue parameters
    
    Returns: (batch_results, failure_details)
        batch_results: list of result dicts, one per row
        failure_details: list of dicts with detailed off-target info for failed primers
    """
    print(f"\n{'='*70}")
    print(f"BATCH: Processing rows {batch_start_idx+1} to {batch_start_idx+len(rows_data)}")
    print(f"{'='*70}")
    
    batch_results = []
    failure_details = []
    
    # Step 1: Design ALL primer alternatives for all rows
    print(f"\nStep 1: Designing primers for {len(rows_data)} rows...")
    
    primer_designs = []
    total_primer_count = 0
    errors = []
    
    for idx, (row_idx, row) in enumerate(rows_data):
        # Progress indicator every 10 rows
        if (idx + 1) % 10 == 0 or idx == 0:
            print(f"  Progress: {idx + 1}/{len(rows_data)} rows processed...", end='\r')
        
        chrom, roi_start, roi_end, gene_strand = get_roi_span(
            row.get('primer_1_locus'), 
            row.get('primer_4_locus')
        )
        
        gene_name = get_gene_name(row, row_idx)
        
        if chrom is None:
            errors.append(f"Row {row_idx+1} ({gene_name}): Could not parse primer loci")
            result = {
                'genotyping_primer_left_seq': 'ERROR',
                'genotyping_primer_left_locus': None,
                'genotyping_primer_left_tm': None,
                'genotyping_primer_left_gc': None,
                'genotyping_primer_right_seq': 'ERROR',
                'genotyping_primer_right_locus': None,
                'genotyping_primer_right_tm': None,
                'genotyping_primer_right_gc': None,
                'amplicon_size': None,
                'genotyping_primers_notes': 'Could not parse primer loci'
            }
            batch_results.append((row_idx, result))
            continue
        
        alternatives = design_verification_primers_multi(
            genome, chrom, roi_start, roi_end, gene_strand, initial_flank,
            num_return=args.num_alternatives
        )
        
        if not alternatives:
            errors.append(f"Row {row_idx+1} ({gene_name}): Primer3 failed to design primers")
            result = {
                'genotyping_primer_left_seq': 'ERROR',
                'genotyping_primer_left_locus': None,
                'genotyping_primer_left_tm': None,
                'genotyping_primer_left_gc': None,
                'genotyping_primer_right_seq': 'ERROR',
                'genotyping_primer_right_locus': None,
                'genotyping_primer_right_tm': None,
                'genotyping_primer_right_gc': None,
                'amplicon_size': None,
                'genotyping_primers_notes': 'Primer3 failed to design primers'
            }
            batch_results.append((row_idx, result))
            continue
        
        primer_designs.append({
            'row_idx': row_idx,
            'row': row,
            'alternatives': alternatives,
            'chrom': chrom,
            'roi_start': roi_start,
            'roi_end': roi_end,
            'gene_strand': gene_strand
        })
        total_primer_count += len(alternatives) * 2  # left + right for each alternative
    
    print(f"\r  Designed primers for {len(primer_designs)}/{len(rows_data)} rows                    ")
    print(f"  Total primer alternatives: {total_primer_count} primers")
    
    if errors:
        print(f"  Errors during design: {len(errors)}")
        for error in errors:
            print(f"    - {error}")
    
    if not primer_designs:
        return batch_results, failure_details
    
    # Step 2: Create mega-FASTA with ALL primer alternatives
    print(f"\nStep 2: Preparing to check {total_primer_count} primers in ONE BLAT call...")
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        primer_fasta = f.name
        for design in primer_designs:
            row_idx = design['row_idx']
            alternatives = design['alternatives']
            
            # Write all alternatives for this row
            for alt_idx, primer_pair in enumerate(alternatives):
                primer_id_left = f"row{row_idx}_left_alt{alt_idx}"
                f.write(f">{primer_id_left}\n{primer_pair['left_seq']}\n")
                
                primer_id_right = f"row{row_idx}_right_alt{alt_idx}"
                f.write(f">{primer_id_right}\n{primer_pair['right_seq']}\n")
    
    # Step 3: Run ONE BLAT call for everything
    with tempfile.NamedTemporaryFile(mode='w', suffix='.psl', delete=False) as f:
        psl_output = f.name
    
    print(f"  Running BLAT...")
    blat_success = run_blat_batch(primer_fasta, genome_fasta_path, psl_output)
    
    if not blat_success:
        print("  WARNING: BLAT failed, using first alternative without off-target check")
        for design in primer_designs:
            primer_pair = design['alternatives'][0]
            result = {
                'genotyping_primer_left_seq': primer_pair['left_seq'],
                'genotyping_primer_left_locus': primer_pair['left_locus'],
                'genotyping_primer_left_tm': primer_pair['left_tm'],
                'genotyping_primer_left_gc': primer_pair['left_gc'],
                'genotyping_primer_right_seq': primer_pair['right_seq'],
                'genotyping_primer_right_locus': primer_pair['right_locus'],
                'genotyping_primer_right_tm': primer_pair['right_tm'],
                'genotyping_primer_right_gc': primer_pair['right_gc'],
                'amplicon_size': primer_pair['amplicon_size'],
                'genotyping_primers_notes': 'BLAT failed - no off-target check performed'
            }
            batch_results.append((design['row_idx'], result))
        
        os.unlink(primer_fasta)
        os.unlink(psl_output)
        return batch_results, failure_details
    
    # Step 4: Parse BLAT results once
    blat_hits = parse_blat_psl(psl_output)
    print(f"  BLAT found hits for {len(blat_hits)} primers")
    
    # Step 5: For each row, find the best clean alternative
    print(f"\nStep 3: Selecting best clean primer pair for each row...")
    
    clean_count = 0
    offtarget_count = 0
    rows_needing_rescue = []  # Track rows that need rescue attempts
    
    for design in primer_designs:
        row_idx = design['row_idx']
        alternatives = design['alternatives']
        chrom = design['chrom']
        gene_name = get_gene_name(design['row'], row_idx)
        
        # Try alternatives in order (already sorted by GC-clamp preference)
        best_pair = None
        best_notes = None
        best_alt_idx = None
        
        # Also track which individual primers are clean for potential mix-and-match
        good_left_primers = []
        good_right_primers = []
        
        # Track ALL alternatives with their off-target counts for best fallback selection
        all_alternatives_with_ot = []
        
        for alt_idx, primer_pair in enumerate(alternatives):
            left_chrom, left_start, left_end, _ = parse_locus(primer_pair['left_locus'])
            right_chrom, right_start, right_end, _ = parse_locus(primer_pair['right_locus'])
            
            primer_id_left = f"row{row_idx}_left_alt{alt_idx}"
            primer_id_right = f"row{row_idx}_right_alt{alt_idx}"
            
            # Check off-targets for this alternative
            left_has_offtarget, left_num_ot, left_details = check_offtargets(
                blat_hits, primer_id_left, left_chrom, left_start, left_end
            )
            
            right_has_offtarget, right_num_ot, right_details = check_offtargets(
                blat_hits, primer_id_right, right_chrom, right_start, right_end
            )
            
            # Track this alternative with off-target info
            total_ot = left_num_ot + right_num_ot
            all_alternatives_with_ot.append({
                'primer_pair': primer_pair,
                'alt_idx': alt_idx,
                'left_num_ot': left_num_ot,
                'right_num_ot': right_num_ot,
                'total_ot': total_ot,
                'left_details': left_details,
                'right_details': right_details
            })
            
            # Track individual clean primers for mix-and-match
            if not left_has_offtarget:
                good_left_primers.append({
                    'seq': primer_pair['left_seq'],
                    'locus': primer_pair['left_locus'],
                    'tm': primer_pair['left_tm'],
                    'gc': primer_pair['left_gc'],
                    'start': left_start,
                    'end': left_end,
                    'chrom': chrom,
                    'alt_idx': alt_idx
                })
            
            if not right_has_offtarget:
                good_right_primers.append({
                    'seq': primer_pair['right_seq'],
                    'locus': primer_pair['right_locus'],
                    'tm': primer_pair['right_tm'],
                    'gc': primer_pair['right_gc'],
                    'start': right_start,
                    'end': right_end,
                    'chrom': chrom,
                    'alt_idx': alt_idx
                })
            
            if not left_has_offtarget and not right_has_offtarget:
                # Found clean pair!
                best_pair = primer_pair
                best_alt_idx = alt_idx
                if alt_idx == 0:
                    best_notes = 'Clean (first choice)'
                else:
                    best_notes = f'Clean (alternative #{alt_idx + 1})'
                clean_count += 1
                break
        
        # If no clean pair found, select the one with fewest off-targets
        if best_pair is None:
            # Sort by total off-targets (ascending), then by alt_idx to prefer earlier alternatives
            all_alternatives_with_ot.sort(key=lambda x: (x['total_ot'], x['alt_idx']))
            best_alt_info = all_alternatives_with_ot[0]
            best_pair = best_alt_info['primer_pair']
            best_alt_idx = best_alt_info['alt_idx']
            best_notes = f'Off-targets: left={best_alt_info["left_details"]}, right={best_alt_info["right_details"]}'
        
        # Print per-row status
        if 'Off-targets' in best_notes:
            offtarget_count += 1
            print(f"  Row {row_idx+1:3d} ({gene_name:20s}): ⚠ All {len(alternatives)} pairs have off-targets, will try rescue...")
            
            # Track for rescue attempts - include all alternatives with off-target info
            rows_needing_rescue.append({
                'design': design,
                'good_left_primers': good_left_primers,
                'good_right_primers': good_right_primers,
                'first_pair': best_pair,
                'first_notes': best_notes,
                'all_alternatives_with_ot': all_alternatives_with_ot
            })
        else:
            if best_alt_idx == 0:
                print(f"  Row {row_idx+1:3d} ({gene_name:20s}): ✓ Clean (first choice)")
            else:
                print(f"  Row {row_idx+1:3d} ({gene_name:20s}): ✓ Clean (alternative #{best_alt_idx + 1}/{len(alternatives)})")
            
            result = {
                'genotyping_primer_left_seq': best_pair['left_seq'],
                'genotyping_primer_left_locus': best_pair['left_locus'],
                'genotyping_primer_left_tm': best_pair['left_tm'],
                'genotyping_primer_left_gc': best_pair['left_gc'],
                'genotyping_primer_right_seq': best_pair['right_seq'],
                'genotyping_primer_right_locus': best_pair['right_locus'],
                'genotyping_primer_right_tm': best_pair['right_tm'],
                'genotyping_primer_right_gc': best_pair['right_gc'],
                'amplicon_size': best_pair['amplicon_size'],
                'genotyping_primers_notes': best_notes
            }
            batch_results.append((row_idx, result))
    
    print(f"\n  Initial Results: {clean_count} clean, {offtarget_count} need rescue")
    
    # Step 4: Mix-and-match rescue for rows with off-targets
    rescue_count = 0
    widening_rescued = 0
    desperation_count = 0

    if rows_needing_rescue:
        print(f"\nStep 4: Mix-and-match rescue for {len(rows_needing_rescue)} rows...")
        print(f"  Attempting to pair clean left primers with clean right primers")
        
        rows_still_failing = []
        
        for rescue_data in rows_needing_rescue:
            design = rescue_data['design']
            row_idx = design['row_idx']
            gene_name = get_gene_name(design['row'], row_idx)
            good_left = rescue_data['good_left_primers']
            good_right = rescue_data['good_right_primers']
            
            print(f"  Row {row_idx+1:3d} ({gene_name:20s}): Found {len(good_left)} clean left, {len(good_right)} clean right")
            
            if good_left and good_right:
                # Find best pairing by Tm match and reasonable amplicon size
                best_pair = None
                best_score = float('inf')
                
                for left in good_left:
                    for right in good_right:
                        tm_diff = abs(left['tm'] - right['tm'])
                        amplicon_size = right['end'] - left['start'] + 1
                        
                        # Prefer reasonable amplicon sizes
                        size_penalty = 0
                        if amplicon_size > 3000:
                            size_penalty = (amplicon_size - 3000) / 1000
                        
                        score = tm_diff + size_penalty
                        
                        if score < best_score:
                            best_score = score
                            best_pair = (left, right, amplicon_size, tm_diff)
                
                if best_pair:
                    left, right, amplicon_size, tm_diff = best_pair
                    result = {
                        'genotyping_primer_left_seq': left['seq'],
                        'genotyping_primer_left_locus': left['locus'],
                        'genotyping_primer_left_tm': left['tm'],
                        'genotyping_primer_left_gc': left['gc'],
                        'genotyping_primer_right_seq': right['seq'],
                        'genotyping_primer_right_locus': right['locus'],
                        'genotyping_primer_right_tm': right['tm'],
                        'genotyping_primer_right_gc': right['gc'],
                        'amplicon_size': amplicon_size,
                        'genotyping_primers_notes': f'Mix-and-match rescue: paired alt#{left["alt_idx"]+1} left with alt#{right["alt_idx"]+1} right (ΔTm={tm_diff:.1f}°C)'
                    }
                    batch_results.append((row_idx, result))
                    rescue_count += 1
                    print(f"       ✓ RESCUED: Paired alt#{left['alt_idx']+1} left with alt#{right['alt_idx']+1} right")
                else:
                    rows_still_failing.append(rescue_data)
            else:
                rows_still_failing.append(rescue_data)
        
        rows_needing_rescue = rows_still_failing
        print(f"  Mix-and-match rescued {rescue_count} rows, {len(rows_needing_rescue)} still need help")
    
    # Step 5: Progressive flank widening for remaining failures
    if rows_needing_rescue:
        print(f"\nStep 5: PROGRESSIVE FLANK WIDENING for {len(rows_needing_rescue)} remaining rows...")
        print(f"  Trying wider flanks from {initial_flank+1000}bp to {max_flank}bp in 1kb increments...")
        print(f"  Using BATCHED approach: all genes at each flank size")
        
        widening_rescued = 0
        remaining_rows = rows_needing_rescue.copy()
        
        # Try each flank size incrementally
        for test_flank in range(initial_flank + 1000, max_flank + 1, 1000):
            if not remaining_rows:
                break
                
            print(f"\n  Testing {test_flank}bp flanks for {len(remaining_rows)} genes...")
            
            # Design primers for all remaining rows at this flank size
            flank_designs = []
            for rescue_data in remaining_rows:
                design = rescue_data['design']
                row_idx = design['row_idx']
                chrom = design['chrom']
                roi_start = design['roi_start']
                roi_end = design['roi_end']
                gene_strand = design['gene_strand']
                
                wider_alternatives = design_verification_primers_multi(
                    genome, chrom, roi_start, roi_end, gene_strand, test_flank,
                    num_return=args.num_alternatives
                )
                
                if wider_alternatives:
                    flank_designs.append({
                        'rescue_data': rescue_data,
                        'row_idx': row_idx,
                        'alternatives': wider_alternatives,
                        'chrom': chrom
                    })
            
            if not flank_designs:
                print(f"    No primers designed at {test_flank}bp")
                continue
            
            # Create mega-FASTA with ALL primers from ALL rows at this flank size
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
                wider_fasta = f.name
                for flank_design in flank_designs:
                    row_idx = flank_design['row_idx']
                    for alt_idx, primer_pair in enumerate(flank_design['alternatives']):
                        f.write(f">row{row_idx}_flank{test_flank}_left{alt_idx}\n{primer_pair['left_seq']}\n")
                        f.write(f">row{row_idx}_flank{test_flank}_right{alt_idx}\n{primer_pair['right_seq']}\n")
            
            # Run ONE BLAT for all primers at this flank size
            with tempfile.NamedTemporaryFile(mode='w', suffix='.psl', delete=False) as f:
                wider_psl = f.name
            
            print(f"    Running BLAT on {len(flank_designs)} genes × {args.num_alternatives} alternatives...")
            if not run_blat_batch(wider_fasta, genome_fasta_path, wider_psl):
                print(f"    BLAT failed at {test_flank}bp")
                os.unlink(wider_fasta)
                os.unlink(wider_psl)
                continue
            
            wider_hits = parse_blat_psl(wider_psl)
            
            # Check each row's alternatives and collect clean individual primers
            rescued_this_flank = []
            rows_still_failing_this_flank = []
            
            for flank_design in flank_designs:
                row_idx = flank_design['row_idx']
                rescue_data = flank_design['rescue_data']
                design = rescue_data['design']
                gene_name = get_gene_name(design['row'], row_idx)
                chrom = flank_design['chrom']
                
                # Track clean individual primers for potential mix-and-match
                good_left_primers = []
                good_right_primers = []
                
                # Try alternatives for this row
                rescued = False
                for alt_idx, primer_pair in enumerate(flank_design['alternatives']):
                    left_chrom, left_start, left_end, _ = parse_locus(primer_pair['left_locus'])
                    right_chrom, right_start, right_end, _ = parse_locus(primer_pair['right_locus'])
                    
                    left_id = f"row{row_idx}_flank{test_flank}_left{alt_idx}"
                    right_id = f"row{row_idx}_flank{test_flank}_right{alt_idx}"
                    
                    left_has_ot, _, _ = check_offtargets(wider_hits, left_id, left_chrom, left_start, left_end)
                    right_has_ot, _, _ = check_offtargets(wider_hits, right_id, right_chrom, right_start, right_end)
                    
                    # Track individual clean primers
                    if not left_has_ot:
                        good_left_primers.append({
                            'seq': primer_pair['left_seq'],
                            'locus': primer_pair['left_locus'],
                            'tm': primer_pair['left_tm'],
                            'gc': primer_pair['left_gc'],
                            'start': left_start,
                            'end': left_end,
                            'chrom': chrom,
                            'alt_idx': alt_idx
                        })
                    
                    if not right_has_ot:
                        good_right_primers.append({
                            'seq': primer_pair['right_seq'],
                            'locus': primer_pair['right_locus'],
                            'tm': primer_pair['right_tm'],
                            'gc': primer_pair['right_gc'],
                            'start': right_start,
                            'end': right_end,
                            'chrom': chrom,
                            'alt_idx': alt_idx
                        })
                    
                    if not left_has_ot and not right_has_ot:
                        # Success!
                        result = {
                            'genotyping_primer_left_seq': primer_pair['left_seq'],
                            'genotyping_primer_left_locus': primer_pair['left_locus'],
                            'genotyping_primer_left_tm': primer_pair['left_tm'],
                            'genotyping_primer_left_gc': primer_pair['left_gc'],
                            'genotyping_primer_right_seq': primer_pair['right_seq'],
                            'genotyping_primer_right_locus': primer_pair['right_locus'],
                            'genotyping_primer_right_tm': primer_pair['right_tm'],
                            'genotyping_primer_right_gc': primer_pair['right_gc'],
                            'amplicon_size': primer_pair['amplicon_size'],
                            'genotyping_primers_notes': f'Wider flank rescue: {test_flank}bp flanks (alternative #{alt_idx+1})'
                        }
                        batch_results.append((row_idx, result))
                        rescued = True
                        widening_rescued += 1
                        rescued_this_flank.append(rescue_data)
                        print(f"    Row {row_idx+1:4d} ({gene_name:20s}): ✓ Rescued")
                        break
                
                if not rescued:
                    # Save for potential mix-and-match
                    rows_still_failing_this_flank.append({
                        'rescue_data': rescue_data,
                        'good_left_primers': good_left_primers,
                        'good_right_primers': good_right_primers,
                        'flank_size': test_flank
                    })
            
            os.unlink(wider_fasta)
            os.unlink(wider_psl)
            
            if rescued_this_flank:
                print(f"    Rescued {len(rescued_this_flank)} genes with matching pairs at {test_flank}bp")
            
            # Try mix-and-match rescue for rows that didn't find matching pairs
            if rows_still_failing_this_flank:
                print(f"    Attempting mix-and-match for {len(rows_still_failing_this_flank)} remaining genes at {test_flank}bp...")
                
                mixmatch_rescued = 0
                for row_data in rows_still_failing_this_flank:
                    rescue_data = row_data['rescue_data']
                    design = rescue_data['design']
                    row_idx = design['row_idx']
                    gene_name = get_gene_name(design['row'], row_idx)
                    chrom = design['chrom']
                    roi_start = design['roi_start']
                    roi_end = design['roi_end']
                    gene_strand = design['gene_strand']
                    
                    good_left = row_data['good_left_primers']
                    good_right = row_data['good_right_primers']
                    
                    if not good_left or not good_right:
                        continue
                    
                    # Try to find compatible pairs
                    found_pair = False
                    for left_primer in good_left:
                        for right_primer in good_right:
                            # Check if primers are from different alternatives (avoid re-testing same pairs)
                            if left_primer['alt_idx'] == right_primer['alt_idx']:
                                continue
                            
                            # Calculate amplicon size
                            amplicon_size = right_primer['end'] - left_primer['start'] + 1
                            
                            # Check if reasonable size (within expected range)
                            roi_len = roi_end - roi_start + 1
                            min_expected = roi_len
                            max_expected = roi_len + 2 * test_flank
                            
                            if min_expected <= amplicon_size <= max_expected:
                                result = {
                                    'genotyping_primer_left_seq': left_primer['seq'],
                                    'genotyping_primer_left_locus': left_primer['locus'],
                                    'genotyping_primer_left_tm': left_primer['tm'],
                                    'genotyping_primer_left_gc': left_primer['gc'],
                                    'genotyping_primer_right_seq': right_primer['seq'],
                                    'genotyping_primer_right_locus': right_primer['locus'],
                                    'genotyping_primer_right_tm': right_primer['tm'],
                                    'genotyping_primer_right_gc': right_primer['gc'],
                                    'amplicon_size': amplicon_size,
                                    'genotyping_primers_notes': f'Mix-match at {test_flank}bp: left alt#{left_primer["alt_idx"]+1}, right alt#{right_primer["alt_idx"]+1}'
                                }
                                batch_results.append((row_idx, result))
                                rescued_this_flank.append(rescue_data)
                                widening_rescued += 1
                                mixmatch_rescued += 1
                                found_pair = True
                                print(f"    Row {row_idx+1:4d} ({gene_name:20s}): ✓ Mix-matched")
                                break
                        
                        if found_pair:
                            break
                
                if mixmatch_rescued:
                    print(f"    Mix-and-match rescued {mixmatch_rescued} additional genes at {test_flank}bp")
            
            # Remove rescued rows from remaining
            remaining_rows = [r for r in remaining_rows if r not in rescued_this_flank]
        
        rows_needing_rescue = remaining_rows
        print(f"  Wider flanks rescued {widening_rescued} rows, {len(rows_needing_rescue)} still need desperation mode")
    
    # Step 6: Desperation mode for remaining failures
    desperation_count = 0
    if rows_needing_rescue:
        print(f"\nStep 6: DESPERATION MODE for {len(rows_needing_rescue)} remaining rows...")
        print(f"  Using relaxed parameters: longer primers (18-35bp), wider Tm (52-68°C), larger search region")
        print(f"  Using BATCHED approach: all genes together")
        
        # Build desperation parameters
        desperation_args = P3_ARGS.copy()
        desperation_args.update({
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 35,
            'PRIMER_OPT_SIZE': 28,
            'PRIMER_MIN_TM': 52.0,
            'PRIMER_MAX_TM': 68.0,
            'PRIMER_MIN_GC': 30.0,
            'PRIMER_MAX_GC': 70.0,
            'PRIMER_NUM_RETURN': 10
        })
        
        desperation_flank = max_flank + 1000
        
        # Design desperation primers for all rows
        desperation_designs = []
        for rescue_data in rows_needing_rescue:
            design = rescue_data['design']
            row_idx = design['row_idx']
            chrom = design['chrom']
            roi_start = design['roi_start']
            roi_end = design['roi_end']
            gene_strand = design['gene_strand']
            
            full_seq, extract_start, extract_end = extract_flanking_sequence(
                genome, chrom, roi_start, roi_end, desperation_flank
            )
            
            if full_seq is None:
                # Can't extract sequence - use fallback immediately
                gene_name = get_gene_name(design['row'], row_idx)
                best_pair, ot_notes = select_best_offtarget_alternative(rescue_data['all_alternatives_with_ot'])
                if best_pair is None:
                    best_pair = rescue_data['first_pair']
                    ot_notes = rescue_data['first_notes']
                
                result = {
                    'genotyping_primer_left_seq': best_pair['left_seq'],
                    'genotyping_primer_left_locus': best_pair['left_locus'],
                    'genotyping_primer_left_tm': best_pair['left_tm'],
                    'genotyping_primer_left_gc': best_pair['left_gc'],
                    'genotyping_primer_right_seq': best_pair['right_seq'],
                    'genotyping_primer_right_locus': best_pair['right_locus'],
                    'genotyping_primer_right_tm': best_pair['right_tm'],
                    'genotyping_primer_right_gc': best_pair['right_gc'],
                    'amplicon_size': best_pair['amplicon_size'],
                    'genotyping_primers_notes': f'FAILED: Could not extract sequence. {ot_notes}'
                }
                batch_results.append((row_idx, result))
                failure_details.append({
                    'row_idx': row_idx,
                    'gene': gene_name,
                    'left_seq': best_pair['left_seq'],
                    'right_seq': best_pair['right_seq'],
                    'notes': ot_notes
                })
                print(f"  Row {row_idx+1:4d} ({gene_name:20s}): ✗ Could not extract sequence")
                continue
            
            # Calculate ROI position
            roi_start_in_seq = roi_start - extract_start
            roi_len = roi_end - roi_start + 1
            target = [[roi_start_in_seq, roi_len]]
            
            min_product = roi_len
            max_product = roi_len + 2 * desperation_flank
            desp_args = desperation_args.copy()
            desp_args['PRIMER_PRODUCT_SIZE_RANGE'] = [[min_product, max_product]]
            
            try:
                primer3_result = primer3.design_primers(
                    {
                        'SEQUENCE_ID': f'{chrom}:{extract_start}-{extract_end}',
                        'SEQUENCE_TEMPLATE': full_seq,
                        'SEQUENCE_TARGET': target,
                    },
                    desp_args
                )
                
                num_pairs = primer3_result.get('PRIMER_PAIR_NUM_RETURNED', 0)
                
                if num_pairs == 0:
                    # Primer3 failed - use fallback
                    continue
                
                # Extract primer pairs
                alternatives = []
                for pair_idx in range(num_pairs):
                    try:
                        left_seq = primer3_result[f'PRIMER_LEFT_{pair_idx}_SEQUENCE']
                        left_start = primer3_result[f'PRIMER_LEFT_{pair_idx}'][0]
                        left_tm = primer3_result[f'PRIMER_LEFT_{pair_idx}_TM']
                        left_gc = primer3_result[f'PRIMER_LEFT_{pair_idx}_GC_PERCENT']
                        
                        right_seq = primer3_result[f'PRIMER_RIGHT_{pair_idx}_SEQUENCE']
                        right_pos = primer3_result[f'PRIMER_RIGHT_{pair_idx}'][0]
                        right_tm = primer3_result[f'PRIMER_RIGHT_{pair_idx}_TM']
                        right_gc = primer3_result[f'PRIMER_RIGHT_{pair_idx}_GC_PERCENT']
                        
                        amplicon_size = primer3_result[f'PRIMER_PAIR_{pair_idx}_PRODUCT_SIZE']
                        
                        # Calculate coordinates
                        left_start_genomic = extract_start + left_start
                        left_end_genomic = left_start_genomic + len(left_seq) - 1
                        
                        right_end_genomic = extract_start + right_pos
                        right_start_genomic = right_end_genomic - len(right_seq) + 1
                        
                        left_locus = f"{chrom}:{left_start_genomic}-{left_end_genomic} ({gene_strand})"
                        right_locus = f"{chrom}:{right_start_genomic}-{right_end_genomic} ({gene_strand})"
                        
                        alternatives.append({
                            'left_seq': left_seq,
                            'left_locus': left_locus,
                            'left_tm': round(left_tm, 1),
                            'left_gc': round(left_gc, 1),
                            'right_seq': right_seq,
                            'right_locus': right_locus,
                            'right_tm': round(right_tm, 1),
                            'right_gc': round(right_gc, 1),
                            'amplicon_size': amplicon_size
                        })
                    except KeyError:
                        continue
                
                if alternatives:
                    desperation_designs.append({
                        'rescue_data': rescue_data,
                        'row_idx': row_idx,
                        'alternatives': alternatives,
                        'chrom': chrom
                    })
                    
            except Exception as e:
                # Primer3 exception - use fallback
                continue
        
        if not desperation_designs:
            print(f"  No desperation primers designed")
        else:
            # Create mega-FASTA with ALL desperation primers
            print(f"  Designed desperation primers for {len(desperation_designs)} genes")
            print(f"  Running BLAT on all desperation primers...")
            
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
                desp_fasta = f.name
                for desp_design in desperation_designs:
                    row_idx = desp_design['row_idx']
                    for alt_idx, primer_pair in enumerate(desp_design['alternatives']):
                        f.write(f">row{row_idx}_desp_left{alt_idx}\n{primer_pair['left_seq']}\n")
                        f.write(f">row{row_idx}_desp_right{alt_idx}\n{primer_pair['right_seq']}\n")
            
            with tempfile.NamedTemporaryFile(mode='w', suffix='.psl', delete=False) as f:
                desp_psl = f.name
            
            if run_blat_batch(desp_fasta, genome_fasta_path, desp_psl):
                desp_hits = parse_blat_psl(desp_psl)
                
                # Check each row's alternatives
                for desp_design in desperation_designs:
                    row_idx = desp_design['row_idx']
                    rescue_data = desp_design['rescue_data']
                    design = rescue_data['design']
                    gene_name = get_gene_name(design['row'], row_idx)
                    
                    rescued = False
                    for alt_idx, primer_pair in enumerate(desp_design['alternatives']):
                        left_chrom, left_start, left_end, _ = parse_locus(primer_pair['left_locus'])
                        right_chrom, right_start, right_end, _ = parse_locus(primer_pair['right_locus'])
                        
                        left_id = f"row{row_idx}_desp_left{alt_idx}"
                        right_id = f"row{row_idx}_desp_right{alt_idx}"
                        
                        left_has_ot, _, _ = check_offtargets(desp_hits, left_id, left_chrom, left_start, left_end)
                        right_has_ot, _, _ = check_offtargets(desp_hits, right_id, right_chrom, right_start, right_end)
                        
                        if not left_has_ot and not right_has_ot:
                            result = {
                                'genotyping_primer_left_seq': primer_pair['left_seq'],
                                'genotyping_primer_left_locus': primer_pair['left_locus'],
                                'genotyping_primer_left_tm': primer_pair['left_tm'],
                                'genotyping_primer_left_gc': primer_pair['left_gc'],
                                'genotyping_primer_right_seq': primer_pair['right_seq'],
                                'genotyping_primer_right_locus': primer_pair['right_locus'],
                                'genotyping_primer_right_tm': primer_pair['right_tm'],
                                'genotyping_primer_right_gc': primer_pair['right_gc'],
                                'amplicon_size': primer_pair['amplicon_size'],
                                'genotyping_primers_notes': f'Desperation rescue: longer primers ({len(primer_pair["left_seq"])}bp L, {len(primer_pair["right_seq"])}bp R), wider search'
                            }
                            batch_results.append((row_idx, result))
                            rescued = True
                            desperation_count += 1
                            print(f"  Row {row_idx+1:4d} ({gene_name:20s}): ✓ Desperation rescued")
                            break
                    
                    if not rescued:
                        # Use fallback
                        best_pair, ot_notes = select_best_offtarget_alternative(rescue_data['all_alternatives_with_ot'])
                        if best_pair is None:
                            best_pair = rescue_data['first_pair']
                            ot_notes = rescue_data['first_notes']
                        
                        result = {
                            'genotyping_primer_left_seq': best_pair['left_seq'],
                            'genotyping_primer_left_locus': best_pair['left_locus'],
                            'genotyping_primer_left_tm': best_pair['left_tm'],
                            'genotyping_primer_left_gc': best_pair['left_gc'],
                            'genotyping_primer_right_seq': best_pair['right_seq'],
                            'genotyping_primer_right_locus': best_pair['right_locus'],
                            'genotyping_primer_right_tm': best_pair['right_tm'],
                            'genotyping_primer_right_gc': best_pair['right_gc'],
                            'amplicon_size': best_pair['amplicon_size'],
                            'genotyping_primers_notes': f'FAILED: Desperation primers had off-targets. {ot_notes}'
                        }
                        batch_results.append((row_idx, result))
                        failure_details.append({
                            'row_idx': row_idx,
                            'gene': gene_name,
                            'left_seq': best_pair['left_seq'],
                            'right_seq': best_pair['right_seq'],
                            'notes': ot_notes
                        })
                        print(f"  Row {row_idx+1:4d} ({gene_name:20s}): ✗ Desperation primers had off-targets")
                
                os.unlink(desp_fasta)
                os.unlink(desp_psl)
            else:
                print(f"  BLAT failed for desperation primers")
                os.unlink(desp_fasta)
                os.unlink(desp_psl)
        
        # Handle any rows that didn't get desperation primers designed
        processed_row_ids = {d['row_idx'] for d in desperation_designs}
        for rescue_data in rows_needing_rescue:
            design = rescue_data['design']
            row_idx = design['row_idx']
            
            if row_idx not in processed_row_ids:
                gene_name = get_gene_name(design['row'], row_idx)
                best_pair, ot_notes = select_best_offtarget_alternative(rescue_data['all_alternatives_with_ot'])
                if best_pair is None:
                    best_pair = rescue_data['first_pair']
                    ot_notes = rescue_data['first_notes']
                
                result = {
                    'genotyping_primer_left_seq': best_pair['left_seq'],
                    'genotyping_primer_left_locus': best_pair['left_locus'],
                    'genotyping_primer_left_tm': best_pair['left_tm'],
                    'genotyping_primer_left_gc': best_pair['left_gc'],
                    'genotyping_primer_right_seq': best_pair['right_seq'],
                    'genotyping_primer_right_locus': best_pair['right_locus'],
                    'genotyping_primer_right_tm': best_pair['right_tm'],
                    'genotyping_primer_right_gc': best_pair['right_gc'],
                    'amplicon_size': best_pair['amplicon_size'],
                    'genotyping_primers_notes': f'FAILED: Desperation mode failed. {ot_notes}'
                }
                batch_results.append((row_idx, result))
                failure_details.append({
                    'row_idx': row_idx,
                    'gene': gene_name,
                    'left_seq': best_pair['left_seq'],
                    'right_seq': best_pair['right_seq'],
                    'notes': ot_notes
                })
        
        print(f"  Desperation rescued {desperation_count} rows, {len(rows_needing_rescue) - desperation_count} still failed")
    
    print(f"\n  Final Summary: {clean_count + rescue_count + widening_rescued + desperation_count} clean total")
    print(f"    - {clean_count} from initial alternatives")
    print(f"    - {rescue_count} from mix-and-match rescue")
    print(f"    - {widening_rescued} from progressive widening")
    print(f"    - {desperation_count} from desperation mode")
    print(f"    - {len(failure_details)} with unavoidable off-targets")
    
    # Cleanup
    os.unlink(primer_fasta)
    os.unlink(psl_output)
    
    return batch_results, failure_details


def main():
    parser = argparse.ArgumentParser(
        description='Design genotyping primers with off-target checking and progressive flank widening',
        epilog='Designs primers starting with specified flank size, progressively widening by 1kb increments before desperation mode'
    )
    parser.add_argument('--input', required=True, help='Input CSV (must have primer_1_locus and primer_4_locus columns)')
    parser.add_argument('--output', required=True, help='Output CSV')
    parser.add_argument('--genome', help='Genome FASTA')
    parser.add_argument('--batch-size', type=int, default=100, 
                       help='Process and save every N rows (default: 100)')
    parser.add_argument('--skip-offtarget-check', action='store_true',
                       help='Skip BLAT off-target checking (faster but less safe)')
    parser.add_argument('--failure-log', 
                       help='Output CSV for detailed off-target failures (default: <output>_failures.csv)')
    parser.add_argument('--no-resume', action='store_true',
                       help='Overwrite existing output instead of resuming from checkpoint')
    parser.add_argument('--initial-flank', type=int, default=1000,
                       help='Initial flanking size in bp (default: 1000)')
    parser.add_argument('--max-flank', type=int, default=5000,
                       help='Maximum flanking size in bp before desperation mode (default: 5000)')
    
    # Primer design parameters
    parser.add_argument('--primer-opt-size', type=int, default=25,
                       help='Optimal primer length in bp (default: 25)')
    parser.add_argument('--primer-min-size', type=int, default=20,
                       help='Minimum primer length in bp (default: 20)')
    parser.add_argument('--primer-max-size', type=int, default=30,
                       help='Maximum primer length in bp (default: 30)')
    parser.add_argument('--primer-opt-tm', type=float, default=58.0,
                       help='Optimal primer Tm in °C (default: 58.0)')
    parser.add_argument('--primer-min-tm', type=float, default=54.0,
                       help='Minimum primer Tm in °C (default: 54.0)')
    parser.add_argument('--primer-max-tm', type=float, default=66.0,
                       help='Maximum primer Tm in °C (default: 66.0)')
    parser.add_argument('--primer-min-gc', type=float, default=35.0,
                       help='Minimum primer GC%% (default: 35.0)')
    parser.add_argument('--primer-max-gc', type=float, default=65.0,
                       help='Maximum primer GC%% (default: 65.0)')
    
    # BLAT parameters
    parser.add_argument('--blat-min-identity', type=int, default=90,
                       help='Minimum %% identity for BLAT hits (default: 90)')
    parser.add_argument('--blat-min-score', type=int, default=20,
                       help='Minimum alignment score for BLAT (default: 20)')
    
    # Primer alternative and rescue parameters
    parser.add_argument('--num-alternatives', type=int, default=50,
                       help='Number of primer pair alternatives to generate from Primer3 (default: 50)')
    parser.add_argument('--max-offtarget-attempts', type=int, default=50,
                       help='Maximum alternative primer pairs to try if off-targets found (default: 50)')
    parser.add_argument('--max-rescue-batches', type=int, default=1,
                       help='Number of batches for mix-and-match rescue (default: 1)')
    parser.add_argument('--num-single-primers', type=int, default=1000,
                       help='Number of new single primers to generate for matching (default: 1000)')
    parser.add_argument('--num-single-primers-search', type=int, default=10000,
                       help='Search space for single primer generation (default: 10000)')
    
    args = parser.parse_args()
    
    # Validate flank parameters
    if args.initial_flank < 100:
        print(f"ERROR: initial-flank must be at least 100bp")
        sys.exit(1)
    if args.max_flank < args.initial_flank:
        print(f"ERROR: max-flank ({args.max_flank}) must be >= initial-flank ({args.initial_flank})")
        sys.exit(1)
    if args.max_flank > 20000:
        print(f"WARNING: max-flank > 20kb may result in very long PCR products")
    
    # Build global P3_ARGS from command-line arguments
    global P3_ARGS, BLAT_MIN_IDENTITY, BLAT_MIN_SCORE
    P3_ARGS = build_p3_args(args)
    BLAT_MIN_IDENTITY = args.blat_min_identity
    BLAT_MIN_SCORE = args.blat_min_score
    
    # Set default failure log name
    if args.failure_log is None:
        base_name = os.path.splitext(args.output)[0]
        args.failure_log = f"{base_name}_failures.csv"
    
    # Check inputs
    if not os.path.exists(args.genome):
        print(f"ERROR: FASTA not found: {args.genome}")
        sys.exit(1)
    
    if not os.path.exists(args.input):
        print(f"ERROR: Input CSV not found: {args.input}")
        sys.exit(1)
    
    # Check BLAT availability
    if not args.skip_offtarget_check:
        if not check_blat_available():
            print("ERROR: BLAT not found in PATH")
            print("  Install BLAT or use --skip-offtarget-check to proceed without validation")
            print("  BLAT download: http://hgdownload.soe.ucsc.edu/admin/exe/")
            sys.exit(1)
        print("✓ BLAT found")
    else:
        print("WARNING: Skipping off-target check (--skip-offtarget-check)")
    
    # Load genome
    genome = load_genome(args.genome)
    
    # Load input CSV
    print(f"\nReading {args.input}...")
    df = pd.read_csv(args.input)
    print(f"Loaded {len(df)} rows")
    
    # Validate required columns
    required_cols = ['primer_1_locus', 'primer_4_locus']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"ERROR: Missing required columns: {missing_cols}")
        sys.exit(1)
    
    # Check for existing output and resume
    start_idx = 0
    if os.path.exists(args.output):
        print(f"Found existing output: {args.output}")
        try:
            existing = pd.read_csv(args.output)
            if 'genotyping_primer_left_seq' in existing.columns:
                completed = existing['genotyping_primer_left_seq'].notna()
                if completed.any():
                    start_idx = completed[::-1].idxmax() + 1
                    print(f"Resuming from row {start_idx + 1}")
        except Exception as e:
            print(f"Could not read existing output: {e}")
    
    # Pre-allocate results
    all_results = [None] * len(df)
    
    # Load existing results if resuming
    if start_idx > 0:
        try:
            existing = pd.read_csv(args.output)
            for i in range(start_idx):
                all_results[i] = {
                    col: existing.loc[i, col] if col in existing.columns else None
                    for col in ['genotyping_primer_left_seq', 'genotyping_primer_left_locus',
                               'genotyping_primer_left_tm', 'genotyping_primer_left_gc',
                               'genotyping_primer_right_seq', 'genotyping_primer_right_locus',
                               'genotyping_primer_right_tm', 'genotyping_primer_right_gc',
                               'amplicon_size', 'genotyping_primers_notes']
                }
        except:
            pass
    
    print(f"\nProcessing rows {start_idx + 1} to {len(df)}...")
    print(f"Batch size: {args.batch_size} rows\n")
    
    # Track all failures for logging
    all_failures = []
    
    # Process in batches
    idx = start_idx
    while idx < len(df):
        batch_end = min(idx + args.batch_size, len(df))
        batch_rows = [(i, df.iloc[i]) for i in range(idx, batch_end)]
        
        # Process batch
        if args.skip_offtarget_check:
            # Simple processing without off-target check
            batch_results = []
            for row_idx, row in batch_rows:
                chrom, roi_start, roi_end, gene_strand = get_roi_span(
                    row.get('primer_1_locus'), 
                    row.get('primer_4_locus')
                )
                
                if chrom is None:
                    result = {
                        'genotyping_primer_left_seq': 'ERROR',
                        'genotyping_primer_left_locus': None,
                        'genotyping_primer_left_tm': None,
                        'genotyping_primer_left_gc': None,
                        'genotyping_primer_right_seq': 'ERROR',
                        'genotyping_primer_right_locus': None,
                        'genotyping_primer_right_tm': None,
                        'genotyping_primer_right_gc': None,
                        'amplicon_size': None,
                        'genotyping_primers_notes': 'Could not parse primer loci'
                    }
                else:
                    alternatives = design_verification_primers_multi(
                        genome, chrom, roi_start, roi_end, gene_strand, args.initial_flank, num_return=1
                    )
                    
                    if alternatives:
                        primer_pair = alternatives[0]
                        result = {
                            'genotyping_primer_left_seq': primer_pair['left_seq'],
                            'genotyping_primer_left_locus': primer_pair['left_locus'],
                            'genotyping_primer_left_tm': primer_pair['left_tm'],
                            'genotyping_primer_left_gc': primer_pair['left_gc'],
                            'genotyping_primer_right_seq': primer_pair['right_seq'],
                            'genotyping_primer_right_locus': primer_pair['right_locus'],
                            'genotyping_primer_right_tm': primer_pair['right_tm'],
                            'genotyping_primer_right_gc': primer_pair['right_gc'],
                            'amplicon_size': primer_pair['amplicon_size'],
                            'genotyping_primers_notes': 'No off-target check performed'
                        }
                    else:
                        result = {
                            'genotyping_primer_left_seq': 'ERROR',
                            'genotyping_primer_left_locus': None,
                            'genotyping_primer_left_tm': None,
                            'genotyping_primer_left_gc': None,
                            'genotyping_primer_right_seq': 'ERROR',
                            'genotyping_primer_right_locus': None,
                            'genotyping_primer_right_tm': None,
                            'genotyping_primer_right_gc': None,
                            'amplicon_size': None,
                            'genotyping_primers_notes': 'Primer3 failed'
                        }
                
                batch_results.append((row_idx, result))
            batch_failures = []  # No failures tracked in skip mode
        else:
            # Process with off-target checking
            batch_results, batch_failures = process_batch_with_offtarget_check(
                batch_rows, genome, args.genome, idx, args.initial_flank, args.max_flank, args
            )
        
        # Store results
        for row_idx, result in batch_results:
            all_results[row_idx] = result
        
        # Accumulate failures
        all_failures.extend(batch_failures)
        
        # Save checkpoint
        print(f"\n{'='*70}")
        print(f"CHECKPOINT: Saving results for rows {idx+1}-{batch_end}")
        print(f"{'='*70}")
        
        # Fill missing results with empty dicts
        save_results = []
        for r in all_results:
            if r is not None:
                save_results.append(r)
            else:
                save_results.append({
                    'genotyping_primer_left_seq': None,
                    'genotyping_primer_left_locus': None,
                    'genotyping_primer_left_tm': None,
                    'genotyping_primer_left_gc': None,
                    'genotyping_primer_right_seq': None,
                    'genotyping_primer_right_locus': None,
                    'genotyping_primer_right_tm': None,
                    'genotyping_primer_right_gc': None,
                    'amplicon_size': None,
                    'genotyping_primers_notes': None
                })
        
        out = pd.concat([df, pd.DataFrame(save_results)], axis=1)
        out.to_csv(args.output, index=False)
        print(f"✓ Saved to {args.output}")
        
        # Save failure log incrementally if there are failures
        if all_failures and not args.skip_offtarget_check:
            failures_df = pd.DataFrame(all_failures)
            # Append to existing file, write header only if file doesn't exist
            write_header = not os.path.exists(args.failure_log)
            failures_df.to_csv(args.failure_log, mode='a', header=write_header, index=False)
            print(f"✓ Updated failure log: {args.failure_log} ({len(all_failures)} failures)")
        
        # Stats
        completed = [r for r in all_results if r is not None]
        success = sum(1 for r in completed 
                    if r['genotyping_primer_left_seq'] not in ['ERROR', None] 
                    and r['genotyping_primer_right_seq'] not in ['ERROR', None])
        with_warnings = sum(1 for r in completed
                          if r.get('genotyping_primers_notes', '').startswith('WARNING:'))
        print(f"Success: {success}/{len(completed)} ({100*success/len(completed):.1f}%)")
        if with_warnings > 0:
            print(f"With off-target warnings: {with_warnings}")
        print()
        
        idx = batch_end
    
    # Final summary
    print(f"\n{'='*70}")
    print(f"FINAL SUMMARY")
    print(f"{'='*70}")
    completed = [r for r in all_results if r is not None]
    success = sum(1 for r in completed 
                if r['genotyping_primer_left_seq'] not in ['ERROR', None] 
                and r['genotyping_primer_right_seq'] not in ['ERROR', None])
    with_warnings = sum(1 for r in completed
                      if r.get('genotyping_primers_notes', '').startswith('WARNING:'))
    clean_success = success - with_warnings
    
    print(f"Total processed: {len(completed)}")
    print(f"Clean success (no off-targets): {clean_success} ({100*clean_success/len(completed):.1f}%)")
    print(f"Success with off-target warnings: {with_warnings}")
    print(f"Failed: {len(completed) - success}")
    if all_failures:
        print(f"\nDetailed failure log: {args.failure_log}")
    print(f"{'='*70}\n")

if __name__ == '__main__':
    main()