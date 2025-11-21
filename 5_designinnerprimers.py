#!/usr/bin/env python3
"""
Homology Arm Primer Design Pipeline - Step 2: Inner Primer Design with PAM/Guide Mutations

Designs internal PCR primers for CRISPR knock-in homology arms with intelligent
silent mutation placement to prevent guide RNA re-cutting after successful integration.

Input:
  - CSV file from script 4 containing genomic sequences and guide coordinates
    (requires: genomic_context, genomic_context_coords, guide_seq, guide_strand,
     site_type, strand, chrom, pos)

Output:
  - CSV file with all input columns plus primer sequences and mutation details:
      • primer_2_seq: Upstream internal primer (reverse complement)
      • primer_3_seq: Downstream internal primer (forward)
      • primer_2a_seq/primer_3a_seq: Bridge primers if mutations span primer regions
      • mutations_detail: Summary of silent mutations introduced
      • Genomic coordinates (locus) for all primers and guides

Mutation Strategy:
  - PAM disruption: If possible, mutate the NGG PAM motif (prevents Cas9 recognition)
  - Seed region mutations: If PAM is coding/immutable, introduce 3-5 silent mutations
    in the guide seed region (positions 8-12 from PAM, critical for targeting)
  - Silent mutations only: All mutations preserve amino acid sequence using
    synonymous codons
  - Adaptive primer design: If mutations fall within primer region, designs
    bridging primers (2a/2b or 3a/3b) to handle mutated templates

Decision Logic:
  1. Analyze guide binding: Determine if guide spans insertion site
  2. Check intact guide length: If ≥18bp of guide remains on either side of
     insertion, mutations are required to prevent re-cutting
  3. Determine mutation location:
     - If guide spans insertion: mutate the side with more intact guide sequence
     - If guide is one-sided: mutate that entire side
  4. Design primers:
     - No mutations needed: simple 30-35bp primers flanking insertion
     - Mutations in non-primer region: standard primers with mutated template
     - Mutations in primer region: bridge primer strategy (primer a + b pairs)

Primer Specifications:
  - Length: 30-35bp (optimized for Tm ~60°C)
  - primer_2: Upstream, binds to coding strand (reverse complement output)
  - primer_3: Downstream, binds to non-coding strand (forward output)
  - Bridge primers: When mutations overlap primer binding sites, creates two
    overlapping primers that together span the mutated region

Processing Steps:
  1. Parse genomic context and locate guide binding site
  2. Identify insertion point (case transition in genomic_context)
  3. Determine if mutations are needed based on intact guide length
  4. If mutations needed:
     a. Attempt PAM site disruption (if in non-coding region)
     b. If PAM immutable, introduce seed region mutations
     c. Verify all mutations are silent (preserve amino acids)
  5. Design primers accounting for mutation positions
  6. Calculate genomic coordinates for all primers

Optional arguments:
  --flush-every N                 Write results every N rows for progress tracking (default: 1000)
  --no-resume                     Overwrite output instead of resuming from checkpoint
  --min-intact-length N           Minimum guide length (bp) on left or right side to require mutations (default: 18)
  --num-mutations-without-pam N   Seed mutations when PAM cannot be mutated (default: 5)
  --num-mutations-with-pam N      Additional seed mutations when PAM is successfully mutated (default: 3)
  --min-primer-len N              Minimum primer length in bp (default: 30)
  --max-primer-len N              Maximum primer length in bp (default: 35)

Example (basic):
  python 5_designinnerprimers.py \
    --input 4.allguidesandgenomicseqs.csv \
    --output 5.allwormguidesinternalprimers.csv

Example (custom parameters):
  python 5_designinnerprimers.py \
    --input 4.allguidesandgenomicseqs.csv \
    --output 5.allwormguidesinternalprimers.csv \
    --min-intact-length 20 \
    --num-mutations-with-pam 4 \
    --min-primer-len 28 \
    --flush-every 5000

Required arguments:
  --input PATH      Input CSV from script 4 with genomic sequences
  --output PATH     Output CSV path with inner primer sequences

Example:
  python 5_designinnerprimers.py \
    -input 4_guides_with_sequences.csv \
    -output 5_guides_with_inner_primers.csv
"""

import argparse, re, os, traceback, time
import pandas as pd
from typing import Dict, List, Tuple, Optional

# Genetic code table
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Amino acid single-letter to three-letter and full name mapping
AA_NAMES = {
    'A': ('Ala', 'Alanine'),
    'C': ('Cys', 'Cysteine'),
    'D': ('Asp', 'Aspartate'),
    'E': ('Glu', 'Glutamate'),
    'F': ('Phe', 'Phenylalanine'),
    'G': ('Gly', 'Glycine'),
    'H': ('His', 'Histidine'),
    'I': ('Ile', 'Isoleucine'),
    'K': ('Lys', 'Lysine'),
    'L': ('Leu', 'Leucine'),
    'M': ('Met', 'Methionine'),
    'N': ('Asn', 'Asparagine'),
    'P': ('Pro', 'Proline'),
    'Q': ('Gln', 'Glutamine'),
    'R': ('Arg', 'Arginine'),
    'S': ('Ser', 'Serine'),
    'T': ('Thr', 'Threonine'),
    'V': ('Val', 'Valine'),
    'W': ('Trp', 'Tryptophan'),
    'Y': ('Tyr', 'Tyrosine'),
    '*': ('Ter', 'Stop'),
    'X': ('Xaa', 'Unknown')
}

# Reverse lookup: amino acid -> list of codons
AA_TO_CODONS = {}
for codon, aa in CODON_TABLE.items():
    if aa not in AA_TO_CODONS:
        AA_TO_CODONS[aa] = []
    AA_TO_CODONS[aa].append(codon)


def select_best_synonymous_codon(original_codon: str, synonymous_codons: List[str], position_to_change: Optional[int] = None) -> Optional[str]:
    """
    Select the best synonymous codon, preferring minimal nucleotide changes.
    If position_to_change is specified, prioritize changing that position.
    """
    if not synonymous_codons:
        return None
    
    best_codon = None
    min_changes = float('inf')
    
    for syn_codon in synonymous_codons:
        if syn_codon == original_codon:
            continue
        
        # Count number of nucleotide changes
        changes = sum(1 for i in range(3) if syn_codon[i] != original_codon[i])
        
        # If we want to change a specific position, prioritize codons that change it
        if position_to_change is not None:
            changes_target_position = syn_codon[position_to_change] != original_codon[position_to_change]
            if not changes_target_position:
                continue  # Skip if doesn't change the target position
            # Prefer codons with fewer total changes
            if changes < min_changes:
                min_changes = changes
                best_codon = syn_codon
        else:
            # Just prefer codons with fewer changes
            if changes < min_changes:
                min_changes = changes
                best_codon = syn_codon
    
    return best_codon


def translate_codon(codon: str) -> str:
    """Translate a codon to amino acid"""
    return CODON_TABLE.get(codon.upper(), 'X')


def get_aa_name(aa_letter: str, format='short') -> str:
    """
    Get amino acid name from single letter code.
    
    Args:
        aa_letter: Single letter amino acid code (e.g., 'F')
        format: 'short' for 3-letter (e.g., 'Phe'), 'long' for full name
    
    Returns:
        Amino acid name
    """
    if aa_letter not in AA_NAMES:
        return 'Xaa' if format == 'short' else 'Unknown'
    return AA_NAMES[aa_letter][0] if format == 'short' else AA_NAMES[aa_letter][1]


def reverse_complement(seq: str) -> str:
    """Return reverse complement, preserving case"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    return ''.join([complement.get(b, b) for b in seq[::-1]])


def calculate_tm(seq: str) -> float:
    """Calculate melting temperature using Wallace rule (simple approximation)"""
    seq_upper = seq.upper()
    gc_count = seq_upper.count('G') + seq_upper.count('C')
    at_count = seq_upper.count('A') + seq_upper.count('T')
    return 4 * gc_count + 2 * at_count


def find_case_transition(seq: str) -> Optional[int]:
    """Find the index where uppercase transitions to lowercase or vice versa"""
    for i in range(len(seq) - 1):
        if seq[i].isupper() != seq[i+1].isupper():
            return i + 1
    return None


def get_synonymous_codons(aa: str) -> List[str]:
    """Get list of synonymous codons for an amino acid"""
    return AA_TO_CODONS.get(aa, [])


def get_codon_overlap_with_guide(codon_idx: int, site_type: str, reading_frame_start: int,
                                   guide_start_in_context: int, guide_end_in_context: int,
                                   target_start_in_context: int) -> Tuple[bool, List[int]]:
    """
    Check if a codon overlaps with the guide binding region and which positions overlap.
    
    Args:
        codon_idx: Index of the codon in the codons list
        site_type: 'Nterm' or 'Cterm'
        reading_frame_start: Starting position of reading frame
        guide_start_in_context: Start position of guide binding in genomic_context
        guide_end_in_context: End position of guide binding in genomic_context
        target_start_in_context: Start position of target sequence in genomic_context
    
    Returns:
        Tuple of (has_overlap, list_of_positions_in_codon_that_overlap)
        - has_overlap: True if any base of the codon overlaps with guide
        - list_of_positions_in_codon: [0,1,2] indicating which positions overlap (e.g., [0,1] means first two bases)
    """
    # Calculate codon positions in genomic_context
    # Since we're working with full genomic context, codon position is simply:
    codon_start_in_context = reading_frame_start + codon_idx * 3
    
    # Check each position of the codon (0, 1, 2)
    positions_in_guide = []
    for pos_in_codon in range(3):
        pos_in_context = codon_start_in_context + pos_in_codon
        if guide_start_in_context <= pos_in_context < guide_end_in_context:
            positions_in_guide.append(pos_in_codon)
    
    has_overlap = len(positions_in_guide) > 0
    return has_overlap, positions_in_guide


def apply_partial_codon_mutation(original_codon: str, mutated_codon: str, 
                                   positions_to_change: List[int]) -> str:
    """
    Apply a codon mutation but only change bases at specified positions.
    Other positions keep their original bases.
    
    Args:
        original_codon: Original codon sequence (e.g., "ATG")
        mutated_codon: Desired mutated codon (e.g., "ATC")
        positions_to_change: List of positions to actually change (e.g., [2] means only change 3rd base)
    
    Returns:
        Partially mutated codon (e.g., "ATG" -> "ATC" with positions=[2] -> "ATC")
    """
    result = list(original_codon)
    for pos in positions_to_change:
        if pos < len(mutated_codon):
            result[pos] = mutated_codon[pos]
    return ''.join(result)


def extract_coding_sequence(seq: str) -> Tuple[str, List[int]]:
    """
    Extract all bases and their positions (treating all as coding since insertion is in-frame)
    Returns: (coding_sequence, list_of_positions_in_original)
    """
    coding_bases = []
    positions = []
    for i, base in enumerate(seq):
        coding_bases.append(base.upper())  # Convert to uppercase for consistency
        positions.append(i)
    return ''.join(coding_bases), positions


def find_subsequence_with_case(haystack: str, needle: str) -> Optional[int]:
    """Find subsequence ignoring case but preserving it"""
    haystack_upper = haystack.upper()
    needle_upper = needle.upper()
    pos = haystack_upper.find(needle_upper)
    return pos if pos != -1 else None


def process_guide(row: pd.Series, row_num: int, args) -> Dict:
    """Process a single guide row and design primers"""
    
    print(f"\n{'='*80}")
    print(f"Processing Row {row_num}: {row['gene']} - {row['site_type']}")
    print(f"{'='*80}")
    
    result = {
        'mutations_needed': False,
        'mutation_reason': '',
        'mutation_side': '',
        'num_mutations': 0,
        'mutations_detail': '',
        'original_sequence': '',
        'mutated_sequence': '',
        'aa_check': '',
        'primer_2_seq': '',
        'primer_2_seq_genomic': '',
        'primer_2_locus': '',
        'primer_2_has_mutations': False,
        'primer_2_unmutated_3p_bp': 0,
        'primer_2_category': '',
        'primer_2a_seq': '',
        'primer_2a_seq_genomic': '',
        'primer_2a_locus': '',
        'primer_2a_overlap_with_2b': '',
        'primer_3_seq': '',
        'primer_3_seq_genomic': '',
        'primer_3_locus': '',
        'primer_3_has_mutations': False,
        'primer_3_unmutated_3p_bp': 0,
        'primer_3_category': '',
        'primer_3a_seq': '',
        'primer_3a_seq_genomic': '',
        'primer_3a_locus': '',
        'primer_3a_overlap_with_3b': '',
        'full_guide_with_pam': '',
        'guide_split_info': '',
        'guide_binding_sequence': '',
        'guide_binding_start': '',
        'guide_binding_end': '',
        'guide_orientation': '',
        'guide_locus': '',
        'error': ''
    }
    
    try:
        # ============================================================
        # STEP 0: PARSE GENOMIC COORDINATES
        # ============================================================
        
        genomic_context_coords = row['genomic_context_coords']  # Format: "I:5107744-5107943"
        
        # Parse the coordinates
        match = re.match(r'(.+):(\d+)-(\d+)', genomic_context_coords)
        if match:
            context_chrom = match.group(1)
            context_start = int(match.group(2))
            context_end = int(match.group(3))
            print(f"Genomic context: {context_chrom}:{context_start}-{context_end}")
        else:
            result['error'] = "Could not parse genomic_context_coords"
            return result
        
        # ============================================================
        # STEP 1: FIND GUIDE SEQUENCE IN GENOMIC CONTEXT
        # ============================================================
        
        # The guide_seq column contains the full guide+PAM sequence in 5'-3' orientation
        # from the guide's perspective. We need to find where this binds in genomic_context.
        guide_seq = row['guide_seq']  # Full guide+PAM, 5'-3' from guide perspective
        site_type = row['site_type']
        genomic_context = row['genomic_context']
        
        # Get coordinate columns for verification
        guide_5p_coord = row.get('guide_5p_coord', 'N/A')
        guide_3p_coord = row.get('guide_3p_coord', 'N/A')
        pam_coords = row.get('pam_coords', 'N/A')
        cut_site_coord = row.get('cut_site_coord', 'N/A')
        
        print(f"\n--- STEP 1: Find Guide Binding Site ---")
        print(f"Guide sequence (5'-3' from guide perspective): {guide_seq}")
        print(f"Genomic context length: {len(genomic_context)}")
        print(f"Expected coordinates:")
        print(f"  guide_5p_coord: {guide_5p_coord}")
        print(f"  guide_3p_coord: {guide_3p_coord}")
        print(f"  pam_coords: {pam_coords}")
        print(f"  cut_site_coord: {cut_site_coord}")
        
        # Try to find guide_seq in genomic_context (case-insensitive)
        genomic_upper = genomic_context.upper()
        guide_upper = guide_seq.upper()
        
        # Try forward orientation first
        forward_pos = genomic_upper.find(guide_upper)
        
        if forward_pos != -1:
            # Found in forward orientation
            print(f"\n✓ Guide found in FORWARD orientation at position {forward_pos}")
            guide_binding_start = forward_pos
            guide_binding_end = forward_pos + len(guide_seq)
            guide_binding_sequence = genomic_context[guide_binding_start:guide_binding_end]
            guide_orientation = 'forward'
            full_guide = guide_seq
            
        else:
            # Try reverse complement
            guide_rc = reverse_complement(guide_seq)
            guide_rc_upper = guide_rc.upper()
            rc_pos = genomic_upper.find(guide_rc_upper)
            
            if rc_pos != -1:
                print(f"\n✓ Guide found in REVERSE COMPLEMENT orientation at position {rc_pos}")
                guide_binding_start = rc_pos
                guide_binding_end = rc_pos + len(guide_rc)
                guide_binding_sequence = genomic_context[guide_binding_start:guide_binding_end]
                guide_orientation = 'reverse_complement'
                full_guide = guide_rc
                
            else:
                result['error'] = "Guide sequence not found in genomic context (tried both orientations)"
                print(f"\n✗ ERROR: {result['error']}")
                print(f"Guide seq: {guide_seq}")
                print(f"Guide RC:  {guide_rc}")
                return result
        
        # Store guide binding information
        result['guide_binding_sequence'] = guide_binding_sequence
        result['guide_binding_start'] = guide_binding_start
        result['guide_binding_end'] = guide_binding_end
        result['guide_orientation'] = guide_orientation
        result['full_guide_with_pam'] = full_guide
        
        # Add genomic coordinates for guide
        guide_chrom = context_chrom
        gene_strand = row['strand']  # '+' or '-'
        
        if gene_strand == '+':
            # Plus strand: direct mapping
            guide_start = context_start + guide_binding_start
            guide_end = context_start + guide_binding_end - 1
            guide_strand = '+' if guide_orientation == 'forward' else '-'
        else:
            # Minus strand: coordinates are reversed
            # The leftmost position in context corresponds to the highest genomic coordinate
            guide_start = context_end - guide_binding_end + 1
            guide_end = context_end - guide_binding_start
            # Strand is flipped relative to orientation
            guide_strand = '-' if guide_orientation == 'forward' else '+'

        result['guide_locus'] = f"{guide_chrom}:{guide_start}-{guide_end} ({guide_strand})"
        
        print(f"Guide orientation: {guide_orientation}")
        print(f"Guide binding region: positions {guide_binding_start} to {guide_binding_end}")
        print(f"Guide binding sequence: {guide_binding_sequence}")
        print(f"  (showing case as in genomic_context)")
        print(f"Guide genomic coordinates: {result['guide_locus']}")
        
        # Extract just the guide portion and PAM portion for downstream analysis
        # The guide_seq should already include PAM, but let's identify which part is which
        # Typically: 20bp protospacer + 3bp PAM = 23bp total
        # We'll infer PAM location based on the guide_strand if available, otherwise use guide_orientation
        
        guide_strand = row.get('guide_strand', None)
        
        # If we have strand information, use it to determine PAM location
        if guide_strand == '+':
            # For + strand, PAM is at 3' end (right side) of guide_seq
            pam_side = 'right'
            print(f"Guide strand: + (PAM on right/3' end)")
        elif guide_strand == '-':
            # For - strand, PAM is at 3' end, but when written 5'-3' from guide perspective,
            # this appears on the left side
            pam_side = 'left'
            print(f"Guide strand: - (PAM on left when written 5'-3')")

        print(f"PAM side in full_guide: {pam_side}")
        
        # ============================================================
        # STEP 2: FIND INSERTION SITE
        # ============================================================
        
        print(f"\n--- STEP 2: Find Insertion Site ---")
        
        insertion_index = find_case_transition(genomic_context)
        if insertion_index is None:
            result['error'] = "No case transition found in genomic_context"
            print(f"ERROR: {result['error']}")
            return result
        
        print(f"Insertion index in genomic_context: {insertion_index}")
        
        left_seq = genomic_context[:insertion_index]
        right_seq = genomic_context[insertion_index:]
        
        print(f"Left sequence: {left_seq[-50:]}... (showing last 50bp)")
        print(f"Right sequence: ...{right_seq[:50]} (showing first 50bp)")
        
        # Check if guide spans insertion site using the binding coordinates we already found
        guide_spans_insertion = (guide_binding_start < insertion_index < guide_binding_end)
        
        if guide_spans_insertion:
            # Calculate where within the guide the insertion site falls
            guide_insertion_index = insertion_index - guide_binding_start
            print(f"✓ Guide SPANS insertion site at position {guide_insertion_index} within the guide")
            result['guide_split_info'] = f"Guide splits at position {guide_insertion_index}"
        else:
            guide_insertion_index = None
            if guide_binding_end <= insertion_index:
                print(f"Guide is entirely LEFT of insertion site")
                result['guide_split_info'] = "Guide entirely on left side"
            else:
                print(f"Guide is entirely RIGHT of insertion site")
                result['guide_split_info'] = "Guide entirely on right side"
        
        # ============================================================
        # STEP 3: DETERMINE IF MUTATIONS ARE NEEDED
        # ============================================================
        
        print(f"\n--- STEP 3: Check if Mutations Needed ---")
        
        if guide_strand == '+':
            # For + strand guide: 5' end is at guide_binding_start, 3' end (with PAM) is at guide_binding_end
            # Check how much of the guide is intact on each side of insertion
            
            guide_length = guide_binding_end - guide_binding_start
            left_intact_length = max(0, min(guide_binding_end, insertion_index) - guide_binding_start) # How much guide is left of insertion
            right_intact_length = max(0, guide_binding_end - max(guide_binding_start, insertion_index)) # How much guide is right of insertion

            print(f"Guide strand +: Full guide length = {guide_length}bp")
            print(f"  Left of insertion:  {left_intact_length}bp intact")
            print(f"  Right of insertion: {right_intact_length}bp intact")
            print(f"  Checking if either side has >= {args.min_intact_length}bp intact...")
            
            # If either side has >= args.min_intact_length bp of intact guide, mutations are needed
            if left_intact_length >= args.min_intact_length:
                result['mutations_needed'] = True
                result['mutation_side'] = 'left'
                result['mutation_reason'] = f"Left side has {left_intact_length}bp intact guide (>= {args.min_intact_length}bp threshold)"
                print(f"✗ Mutations ARE needed on LEFT side: {result['mutation_reason']}")
            elif right_intact_length >= args.min_intact_length:
                result['mutations_needed'] = True
                result['mutation_side'] = 'right'
                result['mutation_reason'] = f"Right side has {right_intact_length}bp intact guide (>= {args.min_intact_length}bp threshold)"
                print(f"✗ Mutations ARE needed on RIGHT side: {result['mutation_reason']}")
            else:
                result['mutations_needed'] = False
                result['mutation_reason'] = f"Neither side has >= {args.min_intact_length}bp intact (left: {left_intact_length}bp, right: {right_intact_length}bp)"
                print(f"✓ Mutations NOT needed: {result['mutation_reason']}")
        
        else:  # guide_strand == '-'
            # For - strand guide: 5' end is at guide_binding_end, 3' end (with PAM) is at guide_binding_start
            # Check how much of the guide is intact on each side of insertion
            
            guide_length = guide_binding_end - guide_binding_start
            left_intact_length = max(0, min(guide_binding_end, insertion_index) - guide_binding_start) # How much guide is left of insertion
            right_intact_length = max(0, guide_binding_end - max(guide_binding_start, insertion_index)) # How much guide is right of insertion

            print(f"Guide strand -: Full guide length = {guide_length}bp")
            print(f"  Left of insertion:  {left_intact_length}bp intact")
            print(f"  Right of insertion: {right_intact_length}bp intact")
            print(f"  Checking if either side has >= {args.min_intact_length}bp intact...")
            
            # If either side has >= args.min_intact_length bp of intact guide, mutations are needed
            if left_intact_length >= args.min_intact_length:
                result['mutations_needed'] = True
                result['mutation_side'] = 'left'
                result['mutation_reason'] = f"Left side has {left_intact_length}bp intact guide (>= {args.min_intact_length}bp threshold)"
                print(f"✗ Mutations ARE needed on LEFT side: {result['mutation_reason']}")
            elif right_intact_length >= args.min_intact_length:
                result['mutations_needed'] = True
                result['mutation_side'] = 'right'
                result['mutation_reason'] = f"Right side has {right_intact_length}bp intact guide (>= {args.min_intact_length}bp threshold)"
                print(f"✗ Mutations ARE needed on RIGHT side: {result['mutation_reason']}")
            else:
                result['mutations_needed'] = False
                result['mutation_reason'] = f"Neither side has >= {args.min_intact_length}bp intact (left: {left_intact_length}bp, right: {right_intact_length}bp)"
                print(f"✓ Mutations NOT needed: {result['mutation_reason']}")
        
        # ============================================================
        # STEP 4: SIMPLE CASE - NO MUTATIONS NEEDED
        # ============================================================
        
        if not result['mutations_needed']:
            print(f"\n--- STEP 4: Design Simple Primers (No Mutations) ---")
            
            # Primer 2: Left internal, reverse complement
            primer_2_template = left_seq[-args.min_primer_len:]
            result['primer_2_seq'] = reverse_complement(primer_2_template).upper()
            result['primer_2_seq_genomic'] = result['primer_2_seq']
            result['primer_2_has_mutations'] = False
            
            # Primer 3: Right internal, forward
            result['primer_3_seq'] = right_seq[:args.min_primer_len].upper()
            result['primer_3_seq_genomic'] = result['primer_3_seq']
            result['primer_3_has_mutations'] = False
            
            primer_2_chrom = context_chrom
            primer_3_chrom = context_chrom

            # Get the gene strand from input
            gene_strand = row['strand']  # '+' or '-'
            print(f"***Gene Strand: {gene_strand}***")

            # Calculate insertion point in genomic coordinates
            insertion_genomic = context_start + insertion_index

            if gene_strand == '+':
                # Plus strand: Primer 2 is left (lower coords), Primer 3 is right (higher coords)
                # Primer 2 (reverse primer, binds to minus strand)
                primer_2_start = insertion_genomic - args.min_primer_len
                primer_2_end = insertion_genomic - 1
                primer_2_strand = '-'
                
                # Primer 3 (forward primer, binds to plus strand)
                primer_3_start = insertion_genomic
                primer_3_end = insertion_genomic + args.min_primer_len - 1
                primer_3_strand = '+'

            else:  # gene_strand == '-'
                # Minus strand: Primer 2 is right (higher coords), Primer 3 is left (lower coords)
                # Primer 2 (forward primer relative to gene, but gene is on minus strand)
                primer_2_start = insertion_genomic
                primer_2_end = insertion_genomic + args.min_primer_len - 1
                primer_2_strand = '+'
                
                # Primer 3 (reverse primer relative to gene, but gene is on minus strand)
                primer_3_start = insertion_genomic - args.min_primer_len
                primer_3_end = insertion_genomic - 1
                primer_3_strand = '-'

            result['primer_2_locus'] = (f"{primer_2_chrom}:{primer_2_start}-{primer_2_end} ({primer_2_strand})")
            result['primer_3_locus'] = (f"{primer_3_chrom}:{primer_3_start}-{primer_3_end} ({primer_3_strand})")     
                                       
            print(f"Primer 2 (left, reverse):")
            print(f"  Template: {primer_2_template}")
            print(f"  Primer:   {result['primer_2_seq']}")
            print(f"  Coordinates: {result['primer_2_locus']}")
            
            print(f"Primer 3 (right, forward):")
            print(f"  Primer:   {result['primer_3_seq']}")
            print(f"  Coordinates: {result['primer_3_locus']}")     

            # Categorize primers - both have no mutations, so all bp are unmutated
            result['primer_2_unmutated_3p_bp'] = args.min_primer_len
            result['primer_2_category'] = 'sufficient_unmutated_3p'
            result['primer_3_unmutated_3p_bp'] = args.min_primer_len
            result['primer_3_category'] = 'sufficient_unmutated_3p'
            
            # Nested primers not needed when no mutations
            result['primer_2a_seq'] = 'N/A'
            result['primer_2a_seq_genomic'] = 'N/A'
            result['primer_2a_overlap_with_2b'] = 'N/A'
            result['primer_3a_seq'] = 'N/A'
            result['primer_3a_seq_genomic'] = 'N/A'
            result['primer_3a_overlap_with_3b'] = 'N/A'
            
            print(f"\nPrimer categorization (no mutations needed):")
            print(f"  Primer 2: {args.min_primer_len}bp unmutated at 3' end -> sufficient_unmutated_3p")
            print(f"  Primer 3: {args.min_primer_len}bp unmutated at 3' end -> sufficient_unmutated_3p")
            
            return result
        
        # ============================================================
        # STEP 5: COMPLEX CASE - MUTATIONS NEEDED
        # ============================================================
        
        print(f"\n--- STEP 5: Design Primers with Mutations ---")
        
        # We already know exactly which bases the guide binds to
        # Only these bases need to be considered for mutations
        print(f"Guide binding region in genomic_context: {guide_binding_start} to {guide_binding_end}")
        print(f"Guide binding sequence: {guide_binding_sequence}")
        
        # Define the mutation target region:
        # From insertion point to end of guide binding + 2bp (to ensure complete codons)
        # This ensures mutations ONLY happen within or immediately adjacent to guide binding
        
        if result['mutation_side'] == 'left':
            print(f"Target sequence is LEFT of insertion")
            # Mutation region: from some point to the left up to insertion
            # Ensure we include the guide binding region
            mutation_region_start_in_context = min(guide_binding_start, insertion_index - 50)
            mutation_region_end_in_context = insertion_index
            
            # Extract this region from genomic_context
            target_seq = genomic_context[mutation_region_start_in_context:mutation_region_end_in_context]
            target_start_in_context = mutation_region_start_in_context
            target_end_in_context = mutation_region_end_in_context
            
            print(f"Mutation region in context: {mutation_region_start_in_context} to {mutation_region_end_in_context}")
            
            # For original_sequence, use the FULL left sequence (0 to insertion_index)
            result['original_sequence'] = left_seq.upper()
            
        else:  # right side
            print(f"Target sequence is RIGHT of insertion")
            # Mutation region: from insertion to end of guide binding + 2bp
            # BUT if guide SPANS insertion, need to start from guide_binding_start
            mutation_region_start_in_context = min(guide_binding_start, insertion_index)
            mutation_region_end_in_context = min(guide_binding_end + 2, len(genomic_context))
            
            # Extract this region from genomic_context
            target_seq = genomic_context[mutation_region_start_in_context:mutation_region_end_in_context]
            target_start_in_context = mutation_region_start_in_context
            target_end_in_context = mutation_region_end_in_context
            
            print(f"Mutation region in context: {mutation_region_start_in_context} to {mutation_region_end_in_context}")
            
            # For original_sequence, use the FULL right sequence (insertion_index to end)
            result['original_sequence'] = right_seq.upper()
        
        # Calculate the overlap between guide binding and mutation target region
        overlap_start = max(guide_binding_start, target_start_in_context)
        overlap_end = min(guide_binding_end, target_end_in_context)
        
        if overlap_start >= overlap_end:
            result['error'] = "Guide binding does not overlap with mutation target side"
            print(f"ERROR: {result['error']}")
            return result
        
        print(f"Overlap region (where mutations allowed): {overlap_start} to {overlap_end}")
        overlap_in_guide = genomic_context[overlap_start:overlap_end]
        print(f"Overlap sequence: {overlap_in_guide}")
        print(f"*** MUTATIONS WILL ONLY BE MADE WITHIN THIS OVERLAP REGION ***")
        
        # Extract codons from ENTIRE genomic_context for proper reading frame
        # Then only mutate codons that fall within the mutation region
        print(f"\n--- Extract Coding Sequence for Mutations ---")
        
        # Get all bases from the guide binding sequence (treat all as coding)
        guide_bases_for_mutation = list(guide_binding_sequence)
        print(f"Guide bases (all treated as coding): {''.join(guide_bases_for_mutation)}")
        
        # Extract coding sequence from ENTIRE genomic_context
        coding_seq_full, coding_positions_full = extract_coding_sequence(genomic_context)
        print(f"Full genomic coding sequence length: {len(coding_seq_full)}")
        
        # Get codons from entire genomic_context with proper reading frame
        # For genomic_context with midpoint insertion, reading frame = (len // 2) % 3
        insertion_index = len(genomic_context) // 2
        reading_frame_start = insertion_index % 3
        
        codons = [coding_seq_full[i:i+3] for i in range(reading_frame_start, len(coding_seq_full), 3)]
        
        print(f"Site type: {site_type}")
        print(f"Insertion index in genomic_context: {insertion_index}")
        print(f"Reading frame start: {reading_frame_start}")
        print(f"Number of codons: {len(codons)}")
        print(f"Codons: {' '.join(codons[:10])}..." if len(codons) > 10 else f"Codons: {' '.join(codons)}")
        
        # Translate
        aa_original = ''.join([translate_codon(c) for c in codons if len(c) == 3])
        print(f"Original AA sequence: {aa_original}")
        
        # ============================================================
        # LOCATE PAM IN CODING SEQUENCE
        # ============================================================
        
        print(f"\n--- Locate PAM for Mutation ---")
        
        if guide_strand == '+':
            pam_seq = full_guide[-3:]
        else:
            pam_seq = full_guide[:3]
        
        print(f"PAM sequence: {pam_seq}")
        print(f"PAM side: {pam_side}")
        
        # ============================================================
        # DESIGN MUTATIONS
        # ============================================================
        
        print(f"\n--- Design Silent Mutations ---")
        
        mutations_list = []
        pam_mutated = False
        
        # ============================================================
        # CALCULATE PAM POSITION IN CODING SEQUENCE
        # ============================================================
        
        print(f"\n--- Calculate PAM Position in Full Genomic Context ---")
        
        # Since we're now working with the full genomic_context,
        # PAM position is simply based on guide_binding_start
        if guide_strand == '+':
            # PAM is last 3bp of guide
            pam_pos_in_coding = guide_binding_start + len(full_guide) - 3
        else:
            # PAM is first 3bp of guide  
            pam_pos_in_coding = guide_binding_start
        
        print(f"PAM position in full genomic_context: {pam_pos_in_coding}")
        
        # Verify this is correct
        if pam_pos_in_coding >= 0 and pam_pos_in_coding + 3 <= len(coding_seq_full):
            found_pam = coding_seq_full[pam_pos_in_coding:pam_pos_in_coding+3]
            print(f"Verification: PAM at calculated position = {found_pam}")
            if found_pam.upper() != pam_seq.upper():
                print(f"WARNING: PAM mismatch! Expected {pam_seq.upper()}, found {found_pam.upper()}")
                pam_pos_in_coding = None
        else:
            print(f"WARNING: PAM position out of bounds")
            pam_pos_in_coding = None
        
        # ============================================================
        # Now continue with mapping guide to codons
        # ============================================================
        
        # Now we can map guide positions to codon indices in the full context
        guide_positions = []
        for i in range(guide_binding_start, guide_binding_end):
            guide_positions.append(i)
        
        print(f"Guide overlaps with positions in genomic_context: {guide_positions[:10]}..." if len(guide_positions) > 10 else f"Guide overlaps with positions: {guide_positions}")
        
        # Convert genomic positions to codon indices
        # Codon index = (position - reading_frame_start) // 3
        guide_codon_indices = set()
        for pos in guide_positions:
            adjusted_pos = pos - reading_frame_start
            if adjusted_pos >= 0:
                codon_idx = adjusted_pos // 3
                if codon_idx < len(codons):
                    guide_codon_indices.add(codon_idx)
        
        print(f"Guide overlaps with codon indices: {sorted(guide_codon_indices)}")
        
        # Try to mutate PAM (the GG in NGG, not the N)
        # For minus strand, this means the CC in CCN
        print(f"\nAttempting PAM mutation:")
        print(f"Note: Only targeting the GG in NGG (or CC in CCN for minus strand)")
        
        pam_mutation_options = []
        
        if pam_pos_in_coding is not None:
            # Determine which PAM positions to target based on strand
            if guide_strand == '+':
                # PAM is NGG - target positions 1 and 2 (the two Gs)
                pam_positions_to_target = [1, 2]
                print(f"Plus strand guide: PAM = {pam_seq}, targeting positions 1 and 2 (GG)")
            else:
                # PAM is CCN (reverse complement of NGG) - target positions 0 and 1 (the two Cs)
                pam_positions_to_target = [0, 1]
                print(f"Minus strand guide: PAM = {pam_seq}, targeting positions 0 and 1 (CC)")
            
            # Try mutating each targeted PAM base and collect all options
            for i in pam_positions_to_target:
                if i >= len(pam_seq):
                    continue
                    
                pam_base = pam_seq[i]
                pam_base_pos = pam_pos_in_coding + i
                
                # Which codon? (using full genomic context coordinates)
                adjusted_pos = pam_base_pos - reading_frame_start
                if adjusted_pos < 0:
                    continue
                codon_index = adjusted_pos // 3
                position_in_codon = adjusted_pos % 3
                
                if codon_index >= len(codons):
                    continue
                
                # CHECK: Does this codon overlap with guide, and which positions?
                has_overlap, positions_in_guide = get_codon_overlap_with_guide(
                    codon_index, site_type, reading_frame_start,
                    guide_binding_start, guide_binding_end, 0  # target_start is 0 since we use full context
                )
                
                if not has_overlap:
                    print(f"  PAM base {i}: codon {codon_index} does NOT overlap with guide - skipping")
                    continue
                
                print(f"  PAM base {i}: codon {codon_index} overlaps guide at positions {positions_in_guide}")
                
                original_codon = codons[codon_index]
                if len(original_codon) < 3:
                    continue
                    
                original_aa = translate_codon(original_codon)
                
                print(f"  PAM base {i} ('{pam_base}'): codon {codon_index} = {original_codon} ({original_aa}), position {position_in_codon}")
                
                # Try to find synonymous codon with different base at this position
                synonymous = get_synonymous_codons(original_aa)
                print(f"    Synonymous codons for {original_aa}: {synonymous}")
                
                # Select best codon that changes this position with minimal total changes
                best_codon = select_best_synonymous_codon(original_codon, synonymous, position_to_change=position_in_codon)
                
                if best_codon:
                    # Check if the PAM position we want to change is actually in the guide
                    if position_in_codon not in positions_in_guide:
                        print(f"    → Skipping: target PAM position {position_in_codon} not in guide region {positions_in_guide}")
                        continue
                    
                    # Apply mutation only at positions within guide
                    final_codon = apply_partial_codon_mutation(original_codon, best_codon, positions_in_guide)
                    
                    # CRITICAL: Verify the partial mutation is still synonymous!
                    final_aa = translate_codon(final_codon)
                    if final_aa != original_aa:
                        print(f"    → Skipping: partial mutation creates non-synonymous codon!")
                        print(f"       {original_codon} ({original_aa}) -> {final_codon} ({final_aa})")
                        continue
                    
                    # Count how many nucleotide changes this requires (only in guide region)
                    num_changes = sum(1 for pos in positions_in_guide if best_codon[pos] != original_codon[pos])
                    
                    print(f"    → Option: {original_codon} -> {final_codon} ({num_changes} nt changes in guide region, positions {positions_in_guide})")
                    pam_mutation_options.append({
                        'codon_index': codon_index,
                        'original_codon': original_codon,
                        'mutated_codon': final_codon,
                        'amino_acid': original_aa,
                        'reason': 'PAM disruption',
                        'position_in_codon': position_in_codon,
                        'num_changes': num_changes,
                        'pam_base_index': i,
                        'positions_in_guide': positions_in_guide
                    })
            
            # Select the PAM mutation with minimal nucleotide changes
            if pam_mutation_options:
                best_pam_mutation = min(pam_mutation_options, key=lambda x: x['num_changes'])
                print(f"    ✓ Selected best PAM mutation: {best_pam_mutation['original_codon']} -> {best_pam_mutation['mutated_codon']} ({best_pam_mutation['num_changes']} changes)")
                mutations_list.append(best_pam_mutation)
                pam_mutated = True
        
        if pam_mutated:
            print(f"\n✓ PAM successfully mutated")
            print(f"Now searching entire guide for {args.num_mutations_with_pam} additional mutation(s), prioritizing positions closest to PAM...")
            
            # Identify all codons that overlap with the guide (excluding PAM codon if already mutated)
            candidate_codon_indices = set()
            for i in range(len(full_guide) - 3):  # Exclude PAM (last 3bp for +, first 3bp for -)
                if guide_strand == '+':
                    guide_pos = i  # 0 to len-3
                else:
                    guide_pos = i + 3  # 3 to len
                
                genomic_pos = guide_binding_start + guide_pos
                adjusted_pos = genomic_pos - reading_frame_start
                if adjusted_pos >= 0:
                    codon_idx = adjusted_pos // 3
                    if codon_idx < len(codons):
                        candidate_codon_indices.add(codon_idx)
            
            print(f"Candidate codon indices across full guide: {sorted(candidate_codon_indices)}")
            
            # Sort codons by distance from PAM (closer to PAM = higher priority)
            if guide_strand == '+':
                # PAM is at end, so higher indices are closer
                sorted_seed_codons = sorted(candidate_codon_indices, reverse=True)
            else:
                # PAM is at start, so lower indices are closer
                sorted_seed_codons = sorted(candidate_codon_indices)
            
            print(f"Codons sorted by distance from PAM (closest first): {sorted_seed_codons[:10]}..." if len(sorted_seed_codons) > 10 else f"Codons sorted: {sorted_seed_codons}")
            
            # Try to mutate codons in seed region that haven't been mutated yet
            # Avoid mutations where the changed nucleotides are adjacent
            mutated_codon_info = [(m['codon_index'], m['original_codon'], m['mutated_codon']) 
                                   for m in mutations_list]
            
            # Get positions of all mutated nucleotides
            mutated_nt_positions = set()
            for codon_idx, orig, mut in mutated_codon_info:
                for pos_in_codon in range(3):
                    if orig[pos_in_codon] != mut[pos_in_codon]:
                        # This nucleotide position is mutated (in genomic_context coordinates)
                        nt_pos = reading_frame_start + codon_idx * 3 + pos_in_codon
                        mutated_nt_positions.add(nt_pos)
            
            print(f"Already mutated nucleotide positions: {sorted(mutated_nt_positions)}")
            
            seed_mutations_added = 0
            for codon_idx in sorted_seed_codons:
                if seed_mutations_added >= args.num_mutations_with_pam:
                    break
                    
                if codon_idx in [m['codon_index'] for m in mutations_list]:
                    continue
                
                # CHECK: Does this codon overlap with guide?
                has_overlap, positions_in_guide = get_codon_overlap_with_guide(
                    codon_idx, site_type, reading_frame_start,
                    guide_binding_start, guide_binding_end, 0  # target_start is 0 for full context
                )
                
                if not has_overlap:
                    print(f"  Seed codon {codon_idx} does NOT overlap with guide - skipping")
                    continue
                
                original_codon = codons[codon_idx]
                if len(original_codon) < 3:
                    continue
                
                original_aa = translate_codon(original_codon)
                synonymous = get_synonymous_codons(original_aa)
                
                # Try all synonymous codons, sorted by number of changes (fewest first)
                synonymous_sorted = sorted(
                    [s for s in synonymous if s != original_codon],
                    key=lambda s: sum(1 for i in range(3) if s[i] != original_codon[i])
                )
                
                mutation_added = False
                for candidate_codon in synonymous_sorted:
                    # Apply mutation only at positions within guide
                    final_codon = apply_partial_codon_mutation(original_codon, candidate_codon, positions_in_guide)
                    
                    # CRITICAL: Verify the partial mutation is still synonymous!
                    final_aa = translate_codon(final_codon)
                    if final_aa != original_aa:
                        print(f"  Skipping option {original_codon} -> {final_codon} (partial mutation creates non-synonymous codon: {original_aa} -> {final_aa})")
                        continue
                    
                    # Check if this mutation would create adjacent nucleotide changes
                    # Only count changes that are actually in the guide region
                    new_mutation_positions = []
                    for pos_in_codon in positions_in_guide:
                        if candidate_codon[pos_in_codon] != original_codon[pos_in_codon]:
                            if site_type == 'Nterm':
                                nt_pos = codon_idx * 3 + pos_in_codon
                            else:
                                nt_pos = reading_frame_start + codon_idx * 3 + pos_in_codon
                            new_mutation_positions.append(nt_pos)
                    
                    # Check if any new mutation position is adjacent to existing ones
                    has_adjacent = False
                    for new_pos in new_mutation_positions:
                        if (new_pos - 1 in mutated_nt_positions) or (new_pos + 1 in mutated_nt_positions):
                            has_adjacent = True
                            break
                    
                    if has_adjacent:
                        print(f"  Skipping option {original_codon} -> {final_codon} (would create adjacent nucleotide mutations at positions {new_mutation_positions})")
                        continue
                    
                    # Found a valid mutation!
                    print(f"  Adding seed mutation: codon {codon_idx} {original_codon} -> {final_codon} ({original_aa}), positions {positions_in_guide}")
                    mutations_list.append({
                        'codon_index': codon_idx,
                        'original_codon': original_codon,
                        'mutated_codon': final_codon,
                        'amino_acid': original_aa,
                        'reason': f'Seed mutation ({seed_mutations_added+1} of {args.num_mutations_with_pam})',
                        'positions_in_guide': positions_in_guide
                    })
                    
                    # Update mutated_nt_positions for next iteration
                    for new_pos in new_mutation_positions:
                        mutated_nt_positions.add(new_pos)
                    
                    seed_mutations_added += 1
                    mutation_added = True
                    break  # Found a valid mutation, move to next codon
                
                if not mutation_added:
                    print(f"  No valid synonymous mutations found for codon {codon_idx} {original_codon} (all options create adjacency issues)")
        
        else:
            print(f"\n✗ PAM could not be mutated (e.g., in-frame Proline)")
            print(f"Attempting to add {args.num_mutations_without_pam} mutations across entire guide, prioritizing positions closest to PAM...")
            
            # Identify all codons that overlap with the guide (excluding PAM region)
            candidate_codon_indices = set()
            for i in range(len(full_guide) - 3):  # Exclude PAM
                if guide_strand == '+':
                    guide_pos = i
                else:
                    guide_pos = i + 3
                
                genomic_pos = guide_binding_start + guide_pos
                adjusted_pos = genomic_pos - reading_frame_start
                if adjusted_pos >= 0:
                    codon_idx = adjusted_pos // 3
                    if codon_idx < len(codons):
                        candidate_codon_indices.add(codon_idx)
            
            print(f"Candidate codon indices: {sorted(candidate_codon_indices)}")
            
            # Sort codons by distance from PAM (closer to PAM = higher priority)
            if guide_strand == '+':
                sorted_seed_codons_no_pam = sorted(candidate_codon_indices, reverse=True)
            else:
                sorted_seed_codons_no_pam = sorted(candidate_codon_indices)
            
            print(f"Codons sorted by distance from PAM (closest first): {sorted_seed_codons_no_pam[:10]}..." if len(sorted_seed_codons_no_pam) > 10 else f"Codons sorted: {sorted_seed_codons_no_pam}")
            
            # Get positions of all mutated nucleotides (from PAM mutation attempts, if any)
            mutated_nt_positions = set()
            for mut in mutations_list:
                codon_idx = mut['codon_index']
                orig = mut['original_codon']
                mutated = mut['mutated_codon']
                for pos_in_codon in range(3):
                    if orig[pos_in_codon] != mutated[pos_in_codon]:
                        if site_type == 'Nterm':
                            nt_pos = codon_idx * 3 + pos_in_codon
                        else:
                            nt_pos = reading_frame_start + codon_idx * 3 + pos_in_codon
                        mutated_nt_positions.add(nt_pos)
            
            print(f"Already mutated nucleotide positions: {sorted(mutated_nt_positions)}")
            
            for codon_idx in sorted_seed_codons_no_pam:
                if len(mutations_list) >= args.num_mutations_without_pam:
                    break
                
                # CHECK: Does this codon overlap with guide?
                has_overlap, positions_in_guide = get_codon_overlap_with_guide(
                    codon_idx, site_type, reading_frame_start,
                    guide_binding_start, guide_binding_end, 0  # target_start is 0 for full context
                )
                
                if not has_overlap:
                    print(f"  Seed codon {codon_idx} does NOT overlap with guide - skipping")
                    continue
                
                original_codon = codons[codon_idx]
                if len(original_codon) < 3:
                    continue
                
                original_aa = translate_codon(original_codon)
                synonymous = get_synonymous_codons(original_aa)
                
                # Try all synonymous codons, sorted by number of changes (fewest first)
                synonymous_sorted = sorted(
                    [s for s in synonymous if s != original_codon],
                    key=lambda s: sum(1 for i in range(3) if s[i] != original_codon[i])
                )
                
                mutation_added = False
                for candidate_codon in synonymous_sorted:
                    # Apply mutation only at positions within guide
                    final_codon = apply_partial_codon_mutation(original_codon, candidate_codon, positions_in_guide)
                    
                    # CRITICAL: Verify the partial mutation is still synonymous!
                    final_aa = translate_codon(final_codon)
                    if final_aa != original_aa:
                        print(f"  Skipping option {original_codon} -> {final_codon} (partial mutation creates non-synonymous codon: {original_aa} -> {final_aa})")
                        continue
                    
                    # Check if this mutation would create adjacent nucleotide changes
                    # Only count changes that are actually in the guide region
                    new_mutation_positions = []
                    for pos_in_codon in positions_in_guide:
                        if candidate_codon[pos_in_codon] != original_codon[pos_in_codon]:
                            if site_type == 'Nterm':
                                nt_pos = codon_idx * 3 + pos_in_codon
                            else:
                                nt_pos = reading_frame_start + codon_idx * 3 + pos_in_codon
                            new_mutation_positions.append(nt_pos)
                    
                    # Check if any new mutation position is adjacent to existing ones
                    has_adjacent = False
                    for new_pos in new_mutation_positions:
                        if (new_pos - 1 in mutated_nt_positions) or (new_pos + 1 in mutated_nt_positions):
                            has_adjacent = True
                            break
                    
                    if has_adjacent:
                        print(f"  Skipping option {original_codon} -> {final_codon} (would create adjacent nucleotide mutations at positions {new_mutation_positions})")
                        continue
                    
                    # Found a valid mutation!
                    print(f"  Adding seed mutation: codon {codon_idx} {original_codon} -> {final_codon} ({original_aa}), positions {positions_in_guide}")
                    mutations_list.append({
                        'codon_index': codon_idx,
                        'original_codon': original_codon,
                        'mutated_codon': final_codon,
                        'amino_acid': original_aa,
                        'reason': f'Seed mutation ({len(mutations_list)+1} of {args.num_mutations_without_pam})',
                        'positions_in_guide': positions_in_guide
                    })
                    
                    # Update mutated_nt_positions for next iteration
                    for new_pos in new_mutation_positions:
                        mutated_nt_positions.add(new_pos)
                    
                    mutation_added = True
                    break  # Found a valid mutation, move to next codon
                
                if not mutation_added:
                    print(f"  No valid synonymous mutations found for codon {codon_idx} {original_codon} (all options create adjacency issues)")
        
        result['num_mutations'] = len(mutations_list)
        print(f"\nTotal mutations designed: {result['num_mutations']}")
        
        # ============================================================
        # APPLY MUTATIONS
        # ============================================================
        
        print(f"\n--- Apply Mutations ---")
        
        mutated_coding_seq = coding_seq_full
        
        for mut in mutations_list:
            codon_idx = mut['codon_index']
            # Codon position in full genomic_context
            codon_start = reading_frame_start + codon_idx * 3
            
            if codon_start + 3 <= len(mutated_coding_seq):
                mutated_coding_seq = (
                    mutated_coding_seq[:codon_start] +
                    mut['mutated_codon'] +
                    mutated_coding_seq[codon_start+3:]
                )
                print(f"  Codon {codon_idx}: {mut['original_codon']} -> {mut['mutated_codon']} ({mut['amino_acid']})")
        
        # Verify AA sequence unchanged
        mutated_codons = [mutated_coding_seq[i:i+3] for i in range(reading_frame_start, len(mutated_coding_seq), 3)]
        aa_mutated = ''.join([translate_codon(c) for c in mutated_codons if len(c) == 3])
        
        print(f"\nOriginal AA: {aa_original}")
        print(f"Mutated AA:  {aa_mutated}")
        
        if aa_original == aa_mutated:
            print(f"✓ AA sequences match!")
            result['aa_check'] = "PASS"
        else:
            print(f"✗ WARNING: AA sequences differ!")
            result['aa_check'] = "FAIL"
        
        # Reconstruct FULL sequence with mutations applied
        # We need to map mutations back to the full left_seq or right_seq
        print(f"\n--- Reconstruct Full Sequence with Mutations ---")
        
        if result['mutation_side'] == 'left':
            # Work with full left_seq (positions 0 to insertion_index in genomic_context)
            full_seq = left_seq
            full_seq_start_in_context = 0
            full_seq_end_in_context = insertion_index
            print(f"Working with LEFT sequence (full_seq length: {len(full_seq)})")
        else:
            # Work with full right_seq (positions insertion_index to end in genomic_context)
            full_seq = right_seq
            full_seq_start_in_context = insertion_index
            full_seq_end_in_context = len(genomic_context)
            print(f"Working with RIGHT sequence (full_seq length: {len(full_seq)})")
        
        print(f"full_seq spans genomic_context positions: {full_seq_start_in_context} to {full_seq_end_in_context}")
        print(f"reading_frame_start: {reading_frame_start}")
        
        # Convert full_seq to uppercase for mutation application
        full_seq_upper = full_seq.upper()
        
        # Apply each mutation to the full sequence
        # The mutated_coding_seq already has all mutations in genomic_context coordinates
        # We just need to extract the relevant portion for this side
        mutated_full_seq = mutated_coding_seq[full_seq_start_in_context:full_seq_end_in_context]
        
        print(f"\nApplied mutations from mutated genomic context")
        print(f"Original full_seq:  {full_seq_upper[:50]}...")
        print(f"Mutated full_seq:   {mutated_full_seq[:50]}...")
        
        result['mutated_sequence'] = mutated_full_seq
        
        print(f"\nFull sequence with mutations (length {len(mutated_full_seq)}):")
        print(f"  Original: {full_seq_upper[:50]}...")
        print(f"  Mutated:  {mutated_full_seq[:50]}...")
        
        # Check if any mutations were actually applied
        if full_seq_upper == mutated_full_seq:
            print(f"   WARNING: No mutations were applied to full sequence!")
        else:
            print(f"  ✓ Mutations successfully applied to full sequence")
        
        # Format mutations detail with independent AA verification for both codons
        mutations_detail = []
        for mut in mutations_list:
            # Independently verify each codon against codon table
            original_aa = translate_codon(mut['original_codon'])
            mutated_aa = translate_codon(mut['mutated_codon'])
            
            # Get three-letter codes
            original_aa_short = get_aa_name(original_aa, format='short')
            mutated_aa_short = get_aa_name(mutated_aa, format='short')
            
            # Format: TTC(F/Phe)->TTT(F/Phe)
            mutations_detail.append(
                f"{mut['original_codon']}({original_aa}/{original_aa_short})->"
                f"{mut['mutated_codon']}({mutated_aa}/{mutated_aa_short})"
            )
        result['mutations_detail'] = '; '.join(mutations_detail)
        
        # ============================================================
        # DESIGN PRIMERS
        # ============================================================
        
        print(f"\n--- Design Primers ---")

        primer_2_chrom = context_chrom
        primer_3_chrom = context_chrom
        
        if result['mutation_side'] == 'left':
            # Primer 2 has mutations (left, reverse)
            # Take last args.max_primer_len bp from the FULL mutated sequence
            primer_2_template = mutated_full_seq[-args.max_primer_len:]
            primer_2_template_genomic = left_seq[-args.max_primer_len:].upper()
            # Mark mutations with lowercase
            primer_2_template_marked = ''
            for i in range(len(primer_2_template)):
                if primer_2_template[i].upper() != primer_2_template_genomic[i]:
                    primer_2_template_marked += primer_2_template[i].lower()
                else:
                    primer_2_template_marked += primer_2_template[i].upper()

            result['primer_2_seq'] = reverse_complement(primer_2_template_marked)
            result['primer_2_seq_genomic'] = reverse_complement(primer_2_template_genomic)
            result['primer_2_has_mutations'] = True

            # Primer 3 no mutations (right, forward)
            result['primer_3_seq']         = right_seq[:args.min_primer_len].upper()
            result['primer_3_seq_genomic'] = right_seq[:args.min_primer_len].upper()
            result['primer_3_has_mutations'] = False
            
            # Get the gene strand from input
            gene_strand = row['strand']  # '+' or '-'
            print(f"***Gene Strand: {gene_strand}***")

            # Calculate insertion point in genomic coordinates
            insertion_genomic = context_start + insertion_index

            if gene_strand == '+':
                # Plus strand: Primer 2 is left (lower coords), Primer 3 is right (higher coords)
                # Primer 2 (reverse primer, binds to minus strand)
                primer_2_start = insertion_genomic - args.max_primer_len
                primer_2_end = insertion_genomic - 1
                primer_2_strand = '-'
                
                # Primer 3 (forward primer, binds to plus strand)
                primer_3_start = insertion_genomic
                primer_3_end = insertion_genomic + args.min_primer_len - 1
                primer_3_strand = '+'

            else:  # gene_strand == '-'
                # Minus strand: Primer 2 is right (higher coords), Primer 3 is left (lower coords)
                # Primer 2 (forward primer relative to gene, but gene is on minus strand)
                primer_2_start = insertion_genomic
                primer_2_end = insertion_genomic + args.max_primer_len - 1
                primer_2_strand = '+'
                
                # Primer 3 (reverse primer relative to gene, but gene is on minus strand)
                primer_3_start = insertion_genomic - args.min_primer_len
                primer_3_end = insertion_genomic - 1
                primer_3_strand = '-'

            result['primer_2_locus'] = (f"{primer_2_chrom}:{primer_2_start}-{primer_2_end} ({primer_2_strand})")
            result['primer_3_locus'] = (f"{primer_3_chrom}:{primer_3_start}-{primer_3_end} ({primer_3_strand})")     
                 
            print(f"Primer 2 (left, reverse, WITH mutations):")
            print(f"  Template (last {args.max_primer_len}bp of mutated left_seq): {primer_2_template}")
            print(f"  Primer (mutated):   {result['primer_2_seq']}")
            print(f"  Primer (genomic):   {result['primer_2_seq_genomic']}")
            print(f"  Coordinates: {result['primer_3_locus']}")     

            print(f"Primer 3 (right, forward, NO mutations):")
            print(f"  Primer:   {result['primer_3_seq']}")
            print(f"  Coordinates: {result['primer_3_locus']}") 

        else:  # mutation_side == 'right'
            # Primer 2 no mutations (left, reverse)
            primer_2_template = left_seq[-args.min_primer_len:]
            result['primer_2_seq']         = reverse_complement(primer_2_template).upper()
            result['primer_2_seq_genomic'] = reverse_complement(primer_2_template).upper()
            result['primer_2_has_mutations'] = False

            # Primer 3 has mutations (right, forward)
            # Take first args.max_primer_len bp from the FULL mutated sequence
            primer_3_template = mutated_full_seq[:args.max_primer_len]
            primer_3_template_genomic = right_seq[:args.max_primer_len].upper()

            # Mark mutations with lowercase
            primer_3_marked = ''
            for i in range(len(primer_3_template)):
                if primer_3_template[i].upper() != primer_3_template_genomic[i]:
                    primer_3_marked += primer_3_template[i].lower()
                else:
                    primer_3_marked += primer_3_template[i].upper()

            result['primer_3_seq'] = primer_3_marked
            result['primer_3_seq_genomic'] = primer_3_template_genomic
            result['primer_3_has_mutations'] = True
            
            primer_2_chrom = context_chrom
            primer_3_chrom = context_chrom

            # Get the gene strand from input
            gene_strand = row['strand']  # '+' or '-'
            print(f"***Gene Strand: {gene_strand}***")

            # Calculate insertion point in genomic coordinates
            insertion_genomic = context_start + insertion_index

            if gene_strand == '+':
                # Plus strand: Primer 2 is left (lower coords), Primer 3 is right (higher coords)
                # Primer 2 (reverse primer, binds to minus strand)
                primer_2_start = insertion_genomic - args.min_primer_len
                primer_2_end = insertion_genomic - 1
                primer_2_strand = '-'
                
                # Primer 3 (forward primer, binds to plus strand)
                primer_3_start = insertion_genomic
                primer_3_end = insertion_genomic + args.min_primer_len - 1
                primer_3_strand = '+'

            else:  # gene_strand == '-'
                # Minus strand: Primer 2 is right (higher coords), Primer 3 is left (lower coords)
                # Primer 2 (forward primer relative to gene, but gene is on minus strand)
                primer_2_start = insertion_genomic
                primer_2_end = insertion_genomic + args.min_primer_len - 1
                primer_2_strand = '+'
                
                # Primer 3 (reverse primer relative to gene, but gene is on minus strand)
                primer_3_start = insertion_genomic - args.min_primer_len
                primer_3_end = insertion_genomic - 1
                primer_3_strand = '-'

            result['primer_2_locus'] = (f"{primer_2_chrom}:{primer_2_start}-{primer_2_end} ({primer_2_strand})")
            result['primer_3_locus'] = (f"{primer_3_chrom}:{primer_3_start}-{primer_3_end} ({primer_3_strand})")     
                 
            print(f"Primer 2 (left, reverse, NO mutations):")
            print(f"  Template (last {args.min_primer_len}bp of left_seq): {primer_2_template}")
            print(f"  Primer:   {result['primer_2_seq']}")
            print(f"  Coordinates: {result['primer_2_locus']}")
            
            print(f"Primer 3 (right, forward, WITH mutations):")
            print(f"  Primer (mutated, first {args.max_primer_len}bp of mutated right_seq): {result['primer_3_seq']}")
            print(f"  Primer (genomic, first {args.max_primer_len}bp of right_seq): {result['primer_3_seq_genomic']}")
            print(f"  Coordinates: {result['primer_3_locus']}")     

        # ============================================================
        # CATEGORIZE PRIMERS BY UNMUTATED 3' END
        # ============================================================
        
        print(f"\n--- Categorize Primers by 3' End Mutation Status ---")
        
        # Determine unmutated length at 3' end for each primer
        # For primers, the 3' end is the end that extends during PCR
        
        if result['mutation_side'] == 'left':
            # Primer 2 has mutations (left, reverse)
            # For reverse primer, 3' end corresponds to the START of the template (beginning of primer_2_template)
            # Check how many bp at the start of primer_2_template are unmutated
            primer_2_template_original = left_seq[-args.max_primer_len:].upper()
            primer_2_unmutated_3p_length = 0
            for i in range(min(len(primer_2_template), len(primer_2_template_original))):
                if primer_2_template[i].upper() == primer_2_template_original[i]:
                    primer_2_unmutated_3p_length += 1
                else:
                    break  # Stop at first mutation
            
            result['primer_2_unmutated_3p_bp'] = primer_2_unmutated_3p_length
            result['primer_2_category'] = 'sufficient_unmutated_3p' if primer_2_unmutated_3p_length >= 15 else 'insufficient_unmutated_3p'
            
            print(f"Primer 2 (reverse):")
            print(f"  Unmutated bp at 3' end: {primer_2_unmutated_3p_length}")
            print(f"  Category: {result['primer_2_category']}")
            
            # Primer 3 has no mutations (right, forward)
            # For forward primer, 3' end corresponds to the END of the sequence
            # All args.min_primer_len bp are unmutated
            result['primer_3_unmutated_3p_bp'] = args.min_primer_len
            result['primer_3_category'] = 'sufficient_unmutated_3p'
            
            print(f"Primer 3 (forward):")
            print(f"  Unmutated bp at 3' end: {args.min_primer_len} (all bp unmutated)")
            print(f"  Category: {result['primer_3_category']}")
        
        else:  # mutation_side == 'right'
            # Primer 2 has no mutations (left, reverse)
            # All args.min_primer_len bp are unmutated
            result['primer_2_unmutated_3p_bp'] = args.min_primer_len
            result['primer_2_category'] = 'sufficient_unmutated_3p'
            
            print(f"Primer 2 (reverse):")
            print(f"  Unmutated bp at 3' end: {args.min_primer_len} (all bp unmutated)")
            print(f"  Category: {result['primer_2_category']}")
            
            # Primer 3 has mutations (right, forward)
            # For forward primer, 3' end corresponds to the END of the sequence (last bp)
            # Check how many bp at the end are unmutated
            primer_3_original = right_seq[:args.max_primer_len].upper()
            primer_3_unmutated_3p_length = 0
            for i in range(1, min(len(result['primer_3_seq']), len(primer_3_original)) + 1):
                if result['primer_3_seq'][-i].upper() == primer_3_original[-i]:
                    primer_3_unmutated_3p_length += 1
                else:
                    break  # Stop at first mutation from the end
            
            result['primer_3_unmutated_3p_bp'] = primer_3_unmutated_3p_length
            result['primer_3_category'] = 'sufficient_unmutated_3p' if primer_3_unmutated_3p_length >= 15 else 'insufficient_unmutated_3p'
            
            print(f"Primer 3 (forward):")
            print(f"  Unmutated bp at 3' end: {primer_3_unmutated_3p_length}")
            print(f"  Category: {result['primer_3_category']}")
        
        # ============================================================
        # DESIGN NESTED PRIMER 2A (if needed)
        # ============================================================
        
        print(f"\n--- Design Nested Primer 2A (if insufficient 3' end) ---")

        # Design nested primer 2A only if primer 2 has insufficient unmutated 3' end OR is supposed to have mutations but doesn't
        if result['primer_2_category'] == 'insufficient_unmutated_3p' or (result['primer_2_has_mutations'] and result['primer_2_seq'] == result['primer_2_seq_genomic']):
            print(f"Primer 2 has insufficient unmutated 3' end ({result['primer_2_unmutated_3p_bp']}bp < 15bp)")
            print(f"Designing nested Primer 2A for first-round PCR...")
            
            if result['mutation_side'] == 'left':
                # Primer 2 is reverse primer on left side with mutations
                # Primer 2A design:
                # - Must have 20bp overlap with 3' end of Primer 2B (which is 5' end of template)
                # - Must have 15bp unmutated sequence at 3' end of the PRIMER (which is 5' end of template)
                
                # For reverse primer, 3' end of primer = 5' end of template (beginning)
                # Find the position of the first (leftmost) mutation in the mutated_full_seq
                first_mutation_pos = -1
                for i in range(len(mutated_full_seq)):
                    if i < len(left_seq) and mutated_full_seq[i].upper() != left_seq[i].upper():
                        first_mutation_pos = i
                        break  # We want the FIRST mutation from the left
                
                if first_mutation_pos == -1:
                    # No mutations found (shouldn't happen, but handle it)
                    print(f"  WARNING: No mutations detected in sequence")
                    extension = 15
                    total_length = args.max_primer_len
                else:
                    print(f"  First mutation position in mutated_full_seq: {first_mutation_pos}")
                    # For reverse primer:
                    # - Primer 2B uses last args.max_primer_len bp of mutated_full_seq (positions -args.max_primer_len to end)
                    # - Primer 2B's 3' end (primer) = 5' end of template = start of those last args.max_primer_len bp
                    # - That's at position len(mutated_full_seq) - args.max_primer_len
                    # - Primer 2A should start earlier to have 20bp overlap and 15bp unmutated at its 3' end
                    # - Primer 2A's 3' end (primer) = 5' end of template should be >= first_mutation_pos + 15
                    
                    primer_2b_template_start = len(mutated_full_seq) - args.max_primer_len
                    
                    # We want the START of Primer 2A template to be at least 15bp BEFORE first mutation
                    # So: primer_2a_start <= first_mutation_pos - 15
                    required_start = max(0, first_mutation_pos - 15)
                    
                    # Primer 2A should overlap 20bp with Primer 2B at the end
                    # So Primer 2A ends at len(mutated_full_seq) - 15 (15bp before the insertion)
                    primer_2a_end = len(mutated_full_seq) - 15
                    
                    # Calculate length
                    total_length = primer_2a_end - required_start
                    extension = len(mutated_full_seq) - primer_2a_end
                    
                    print(f"  Primer 2A template: positions {required_start} to {primer_2a_end} (length {total_length}bp)")
                    print(f"  This ensures {15}bp extension beyond Primer 2B and {primer_2a_end - first_mutation_pos}bp after first mutation")
                
                # Get sequence from mutated sequence for Primer 2A (extending further left)
                primer_2a_template = mutated_full_seq[-(total_length + extension):-extension]
                primer_2a_template_genomic = left_seq[-(total_length + extension):-extension].upper()
                
                # Mark mutations with lowercase
                primer_2a_template_marked = ''
                for i in range(len(primer_2a_template)):
                    if primer_2a_template[i].upper() != primer_2a_template_genomic[i]:
                        primer_2a_template_marked += primer_2a_template[i].lower()
                    else:
                        primer_2a_template_marked += primer_2a_template[i].upper()

                result['primer_2a_seq'] = reverse_complement(primer_2a_template_marked)
                result['primer_2a_seq_genomic'] = reverse_complement(primer_2a_template_genomic)
                result['primer_2a_overlap_with_2b'] = 20  # 20bp overlap with Primer 2B
                
                # Primer 2A genomic coordinates
                primer_2a_chrom = context_chrom

                if gene_strand == '+':
                    # Plus strand: Primer 2 is left (lower coords), Primer 3 is right (higher coords)
                    # Primer 2 (reverse primer, binds to minus strand)
                    primer_2a_start = insertion_genomic - total_length - extension
                    primer_2a_end = insertion_genomic - extension - 1
                    primer_2a_strand = '-'

                else:  # gene_strand == '-'
                    # Minus strand: Primer 2 is right (higher coords), Primer 3 is left (lower coords)
                    # Primer 2 (forward primer relative to gene, but gene is on minus strand)
                    primer_2a_start = insertion_genomic + extension
                    primer_2a_end = insertion_genomic + total_length + extension - 1
                    primer_2a_strand = '+'
            
                result['primer_2a_locus'] = (f"{primer_2a_chrom}:{primer_2a_start}-{primer_2a_end} ({primer_2a_strand})")

                print(f"Primer 2A (left, reverse, nested):")
                print(f"  Template (length {len(primer_2a_template)}bp, ending {extension}bp before insertion): {primer_2a_template}")
                print(f"  Primer (mutated):   {result['primer_2a_seq']}")
                print(f"  Primer (genomic):   {result['primer_2a_seq_genomic']}")
                print(f"  Length:   {len(result['primer_2a_seq'])}bp")
                print(f"  Overlap with Primer 2B: {result['primer_2a_overlap_with_2b']}bp")
                print(f"  Coordinates: {result['primer_2a_locus']}") 

                # Verify 15bp unmutated at 3' end of PRIMER (= 5' end of template)
                template_start_pos = len(mutated_full_seq) - len(primer_2a_template) - extension
                primer_2a_template_original = left_seq[template_start_pos:template_start_pos + len(primer_2a_template)].upper()
                primer_2a_unmutated_3p = 0
                for i in range(len(primer_2a_template)):
                    if primer_2a_template[i].upper() == primer_2a_template_original[i]:
                        primer_2a_unmutated_3p += 1
                    else:
                        break
                print(f"  Verified: {primer_2a_unmutated_3p}bp unmutated at 3' end (5' of template)")
                
                # Check overlap
                primer_2b_3p = result['primer_2_seq'][-20:]  # 3' end of primer 2B (last 20bp)
                primer_2a_5p = result['primer_2a_seq'][:20]  # 5' end of primer 2A (first 20bp)
                
                if primer_2b_3p == primer_2a_5p:
                    print(f"  ✓ Verified: 20bp overlap matches between Primer 2A and 2B")
                else:
                    print(f"    WARNING: Overlap mismatch!")
                    print(f"    Primer 2B 3' end: {primer_2b_3p}")
                    print(f"    Primer 2A 5' end: {primer_2a_5p}")
            
            else:  # mutation_side == 'right', but Primer 2 has no mutations in this case
                # This shouldn't happen since Primer 2 on right side has no mutations
                print(f"  Note: Primer 2 is on left side without mutations - nested primer not needed")
                result['primer_2a_seq'] = 'N/A'
                result['primer_2a_seq_genomic'] = 'N/A'
                result['primer_2a_overlap_with_2b'] = 'N/A'
        
        else:
            print(f"Primer 2 has sufficient unmutated 3' end - nested Primer 2A not needed")
            result['primer_2a_seq'] = 'N/A'
            result['primer_2a_seq_genomic'] = 'N/A'
            result['primer_2a_overlap_with_2b'] = 'N/A'
        
        # ============================================================
        # DESIGN NESTED PRIMER 3A (if needed)
        # ============================================================
        
        print(f"\n--- Design Nested Primer 3A (if insufficient 3' end) ---")
        
        # Design nested primer 3A only if primer 3 has insufficient unmutated 3' end OR is supposed to have mutations but doesn't
        if result['primer_3_category'] == 'insufficient_unmutated_3p' or (result['primer_3_has_mutations'] and result['primer_3_seq'] == result['primer_3_seq_genomic']):
            print(f"Primer 3 has insufficient unmutated 3' end ({result['primer_3_unmutated_3p_bp']}bp < 15bp)")
            print(f"Designing nested Primer 3A for first-round PCR...")
            
            if result['mutation_side'] == 'right':
                # Primer 3 is forward primer on right side with mutations
                # Primer 3A design:
                # - Must have 20bp overlap with 3' end of Primer 3B
                # - Must have 15bp unmutated sequence at 3' end (after last mutation)
                
                # Find the position of the last mutation in the mutated_full_seq
                # Compare original vs mutated to find rightmost (last) mutation position
                last_mutation_pos = -1
                for i in range(len(mutated_full_seq)):
                    if i < len(right_seq) and mutated_full_seq[i].upper() != right_seq[i].upper():
                        last_mutation_pos = i
                
                if last_mutation_pos == -1:
                    # No mutations found (shouldn't happen, but handle it)
                    print(f"  WARNING: No mutations detected in sequence")
                    offset = 15
                    new_length = args.max_primer_len
                else:
                    print(f"  Last mutation position in mutated_full_seq: {last_mutation_pos}")
                    # Calculate offset to ensure 15bp after last mutation
                    # Primer 3A needs to start early enough that it has:
                    # - 20bp overlap with Primer 3B's 3' end
                    # - last mutation is at least 20bp from Primer 3A's 3' end (15bp + 5bp buffer)
                    # 
                    # Primer 3B is first args.max_primer_len bp, so its 3' end is at position args.max_primer_len-1
                    # Primer 3A should have 20bp overlap, so Primer 3A's 5' end should be at position 15
                    # Primer 3A's 3' end will be at position 15 + args.max_primer_len - 1 = args.max_primer_len+14
                    # 
                    # We want: last_mutation_pos + 15 <= (Primer 3A's 3' end)
                    # So we need to extend Primer 3A such that its 3' end is at least last_mutation_pos + 15
                    
                    # Calculate required 3' end position
                    required_3p_end = last_mutation_pos + 15
                    
                    # Default design: start at position 15, length args.max_primer_len, so 3' end at args.max_primer_len+14
                    default_offset = 15
                    default_length = args.max_primer_len
                    default_3p_end = default_offset + default_length - 1
                    
                    if required_3p_end > default_3p_end:
                        # Need to extend Primer 3A
                        new_length = required_3p_end - default_offset + 1
                        print(f"  Extending Primer 3A from {default_length}bp to {new_length}bp to ensure 15bp after last mutation")
                    else:
                        new_length = default_length
                    
                    offset = default_offset
                
                # Get sequence from mutated_full_seq for Primer 3A
                primer_3a_template = mutated_full_seq[offset:offset+new_length]
                primer_3a_template_genomic = right_seq[offset:offset+new_length].upper()
                
                # Mark mutations with lowercase
                primer_3a_marked = ''
                for i in range(len(primer_3a_template)):
                    if primer_3a_template[i].upper() != primer_3a_template_genomic[i]:
                        primer_3a_marked += primer_3a_template[i].lower()
                    else:
                        primer_3a_marked += primer_3a_template[i].upper()

                result['primer_3a_seq'] = primer_3a_marked
                result['primer_3a_seq_genomic'] = primer_3a_template_genomic
                result['primer_3a_overlap_with_3b'] = 20  # 20bp overlap with Primer 3B

                 # Primer 3A genomic coordinates
                primer_3a_chrom = context_chrom

                if gene_strand == '+':
                    # Plus strand: Primer 3 is right (higher coords), Primer 2 is left (lower coords)
                    # Primer 3 (forward primer relative to gene)
                    primer_3a_start = insertion_genomic + offset
                    primer_3a_end = insertion_genomic + offset + new_length - 1
                    primer_3a_strand = '+'

                else:  # gene_strand == '-'
                    # Minus strand: Primer 3 is left (lower coords), Primer 2 is right (higher coords)
                    # Primer 3 (forward primer relative to gene, but gene is on minus strand)
                    primer_3a_start = insertion_genomic - offset - new_length
                    primer_3a_end = insertion_genomic - offset - 1
                    primer_3a_strand = '-'
                
                result['primer_3a_locus'] = (f"{primer_3a_chrom}:{primer_3a_start}-{primer_3a_end} ({primer_3a_strand})")
                
                print(f"Primer 3A (right, forward, nested):")
                print(f"  Template (starting at position {offset}, length {new_length}bp): {primer_3a_template}")
                print(f"  Primer (mutated):   {result['primer_3a_seq']}")
                print(f"  Primer (genomic):   {result['primer_3a_seq_genomic']}")
                print(f"  Length:   {len(result['primer_3a_seq'])}bp")
                print(f"  Overlap with Primer 3B: {result['primer_3a_overlap_with_3b']}bp")
                print(f"  Coordinates: {result['primer_3a_locus']}") 
                
                # Verify 15bp unmutated at 3' end
                primer_3a_original = right_seq[offset:offset+new_length].upper()
                primer_3a_unmutated_3p = 0
                for i in range(1, min(len(primer_3a_template), len(primer_3a_original)) + 1):
                    if primer_3a_template[-i].upper() == primer_3a_original[-i]:
                        primer_3a_unmutated_3p += 1
                    else:
                        break
                print(f"  Verified: {primer_3a_unmutated_3p}bp unmutated at 3' end")
                
                # Check overlap
                primer_3b_3p = result['primer_3_seq'][-20:]  # 3' end of primer 3B (last 20bp)
                primer_3a_5p = result['primer_3a_seq'][:20]  # 5' end of primer 3A (first 20bp)
                
                if primer_3b_3p == primer_3a_5p:
                    print(f"  ✓ Verified: 20bp overlap matches between Primer 3A and 3B")
                else:
                    print(f"    WARNING: Overlap mismatch!")
                    print(f"    Primer 3B 3' end: {primer_3b_3p}")
                    print(f"    Primer 3A 5' end: {primer_3a_5p}")
            
            else:  # mutation_side == 'left', but Primer 3 has no mutations in this case
                # This shouldn't happen since Primer 3 on left side has no mutations
                print(f"  Note: Primer 3 is on right side without mutations - nested primer not needed")
                result['primer_3a_seq'] = 'N/A'
                result['primer_3a_seq_genomic'] = 'N/A'
                result['primer_3a_overlap_with_3b'] = 'N/A'
        
        else:
            print(f"Primer 3 has sufficient unmutated 3' end - nested Primer 3A not needed")
            result['primer_3a_seq'] = 'N/A'
            result['primer_3a_seq_genomic'] = 'N/A'
            result['primer_3a_overlap_with_3b'] = 'N/A'
        
        # ============================================================
        # SUMMARY: GUIDE BINDING VERIFICATION
        # ============================================================
        
        print(f"\n{'='*80}")
        print(f"GUIDE BINDING SUMMARY")
        print(f"{'='*80}")
        print(f"Guide orientation: {result['guide_orientation']}")
        print(f"Guide binding positions in genomic_context: {result['guide_binding_start']} to {result['guide_binding_end']}")
        print(f"Guide binding sequence: {result['guide_binding_sequence']}")
        print(f"Guide genomic coordinates: {result['guide_locus']}")
        print(f"\nCoordinate Verification:")
        print(f"  Expected guide_5p_coord:  {guide_5p_coord}")
        print(f"  Expected guide_3p_coord:  {guide_3p_coord}")
        print(f"  Expected pam_coords:      {pam_coords}")
        print(f"  Expected cut_site_coord:  {cut_site_coord}")
        print(f"\n  Detected guide_start:     {guide_binding_start}")
        print(f"  Detected guide_end:       {guide_binding_end}")
        print(f"{'='*80}\n")
        
        return result
    
    except Exception as e:
        print(f"\nERROR processing row: {str(e)}")
        traceback.print_exc()
        result['error'] = str(e)
        return result

# ===============================================================================
# Main
# ===============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Primer design script with genomic coordinates and configurable mutation strategy.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
                Mutation Strategy Parameters:
                --min-intact-length: Threshold for requiring mutations (default: 18bp)
                                    If ≥18bp of guide remains intact on either side of insertion,
                                    mutations are needed to prevent guide re-cutting
                
                --num-mutations-with-pam: Number of seed mutations when PAM is mutated (default: 3)
                --num-mutations-without-pam: Number of seed mutations when PAM is immutable (default: 5)
                
                --min-primer-len: Minimum primer length in bp (default: 30)
                --max-primer-len: Maximum primer length in bp (default: 35)
                        """
    )
    
    # Required arguments
    parser.add_argument('--input', help='Input CSV file', required=True)
    parser.add_argument('--output', help='Output CSV file', required=True)
    
    # Progress tracking
    parser.add_argument('--flush-every', type=int, default=1000,
                        help='Write results to disk every N rows for progress tracking (default: 1000)')
    parser.add_argument('--no-resume', action='store_true',
                        help='Overwrite existing output instead of resuming from checkpoint')
    
    # Optional mutation strategy arguments
    parser.add_argument('--min-intact-length', type=int, default=18,
                        help='Minimum intact guide length (bp) that requires mutations (default: 18)')
    parser.add_argument('--num-mutations-with-pam', type=int, default=3,
                        help='Number of additional mutations when PAM is successfully mutated (default: 3)')
    parser.add_argument('--num-mutations-without-pam', type=int, default=5,
                        help='Number of mutations when PAM cannot be mutated (default: 5)')
    parser.add_argument('--min-primer-len', type=int, default=30,
                        help='Minimum primer length in bp (default: 30)')
    parser.add_argument('--max-primer-len', type=int, default=35,
                        help='Maximum primer length in bp (default: 35)')
    
    args = parser.parse_args()

    t0 = time.time()

    input_file = args.input
    output_file = args.output
    
    print(f"Reading input file: {input_file}")
    df = pd.read_csv(input_file)
    
    print(f"Found {len(df)} rows")
    print(f"\nConfiguration:")
    print(f"  Minimum intact guide length for mutations: {args.min_intact_length} bp")
    print(f"  Seed mutations (with PAM disruption): {args.num_mutations_with_pam}")
    print(f"  Seed mutations (without PAM disruption): {args.num_mutations_without_pam}")
    print(f"  Primer length range: {args.min_primer_len}-{args.max_primer_len} bp")
    print(f"  Flush every: {args.flush_every} rows")
    
    # Define result columns (must match the keys in process_guide result dict)
    result_columns = [
        'mutations_needed', 'mutation_reason', 'mutation_side', 'num_mutations',
        'mutations_detail', 'original_sequence', 'mutated_sequence', 'aa_check',
        'primer_2_seq', 'primer_2_seq_genomic', 'primer_2_locus', 
        'primer_2_has_mutations', 'primer_2_unmutated_3p_bp', 'primer_2_category',
        'primer_2a_seq', 'primer_2a_seq_genomic', 'primer_2a_locus', 'primer_2a_overlap_with_2b',
        'primer_3_seq', 'primer_3_seq_genomic', 'primer_3_locus',
        'primer_3_has_mutations', 'primer_3_unmutated_3p_bp', 'primer_3_category',
        'primer_3a_seq', 'primer_3a_seq_genomic', 'primer_3a_locus', 'primer_3a_overlap_with_3b',
        'full_guide_with_pam', 'guide_split_info', 'guide_binding_sequence',
        'guide_binding_start', 'guide_binding_end', 'guide_orientation', 'guide_locus', 'error'
    ]
    
    if args.no_resume or not os.path.exists(output_file):
        # Create new output file with header including both input and result columns
        print(f"[start] Creating new output file: {output_file}")
        # Create empty dataframe with all columns
        all_columns = list(df.columns) + result_columns
        empty_df = pd.DataFrame(columns=all_columns)
        empty_df.to_csv(output_file, index=False)
        processed_rows = set()
    else:
        # Resume mode: read existing output to determine what's been processed
        print(f"[resume] Found existing output file: {output_file}")
        try:
            existing_df = pd.read_csv(output_file)
            # Use row index from input to track what's been processed
            # We'll check if 'primer_2_seq' column exists and is not empty
            if 'primer_2_seq' in existing_df.columns:
                processed_rows = set(existing_df[existing_df['primer_2_seq'].notna()].index)
                print(f"[resume] Found {len(processed_rows)} already-processed rows")
            else:
                processed_rows = set()
                print(f"[resume] No processed rows found (missing primer_2_seq column)")
        except Exception as e:
            print(f"[resume] Error reading existing file: {e}")
            print(f"[resume] Starting fresh")
            processed_rows = set()
    
    # Process rows incrementally
    results_buffer = []
    rows_buffer = []
    flush_every = args.flush_every
    processed_count = 0
    skipped_count = 0
    
    for idx, row in df.iterrows():
        # Skip if already processed (resume mode)
        if not args.no_resume and idx in processed_rows:
            skipped_count += 1
            continue
        
        # Print skip summary before processing first new row
        if skipped_count > 0 and processed_count == 0:
            print(f"[resume] Skipped {skipped_count} already-processed row(s), continuing with new rows...\n")
        
        # Process this row
        result = process_guide(row, idx + 1, args)
        results_buffer.append(result)
        rows_buffer.append(row)
        processed_count += 1
        
        # Flush to disk every N rows
        if len(results_buffer) >= flush_every:
            print(f"\n[flush] Writing {len(results_buffer)} rows to disk (total processed: {processed_count + skipped_count}/{len(df)})...")
            
            # Combine input rows with results
            input_chunk = pd.DataFrame(rows_buffer)
            results_chunk = pd.DataFrame(results_buffer)
            output_chunk = pd.concat([input_chunk.reset_index(drop=True), results_chunk], axis=1)
            
            # Append to output file
            output_chunk.to_csv(output_file, mode='a', header=False, index=False)
            
            # Clear buffers
            results_buffer = []
            rows_buffer = []
    
    # Final flush for remaining rows
    if results_buffer:
        print(f"\n[flush] Writing final {len(results_buffer)} rows to disk...")
        input_chunk = pd.DataFrame(rows_buffer)
        results_chunk = pd.DataFrame(results_buffer)
        output_chunk = pd.concat([input_chunk.reset_index(drop=True), results_chunk], axis=1)
        output_chunk.to_csv(output_file, mode='a', header=False, index=False)
    
    print(f"\n{'='*80}")
    print(f"Processing complete!")
    print(f"  Total rows in input: {len(df)}")
    print(f"  Rows processed: {processed_count}")
    if skipped_count > 0:
        print(f"  Rows skipped (resume): {skipped_count}")
    print(f"Output saved to: {output_file}")
    # Timing
    elapsed = time.time() - t0
    rate = processed_count / elapsed if elapsed > 0 else 0
    print(f"[timing] Processed {processed_count} new rows in {elapsed:.2f}s ({rate:.2f} rows/s)")
    print(f"{'='*80}")

if __name__ == '__main__':
    main()