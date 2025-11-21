#!/usr/bin/env python3
"""
Homology Arm Primer Design Pipeline - Step 1: Genomic Sequence Extraction

Extracts genomic sequences flanking insertion sites for PCR primer design in
CRISPR knock-in experiments.

Input:
  - CSV file from script 3 containing guide sequences and insertion coordinates
    (requires: chrom, pos, strand, site_type)

Output:
  - CSV file with all input columns plus:
      • upstream_arm: Upstream genomic sequence (default 100bp)
      • downstream_arm: Downstream genomic sequence (default 100bp)
      • Case formatting indicates insertion boundary:
        - N-terminal: lowercase UPSTREAM (5' UTR), uppercase DOWNSTREAM (CDS with start codon)
        - C-terminal: uppercase UPSTREAM (CDS), lowercase DOWNSTREAM (stop codon + 3' UTR)

Coordinate System:
  - Genome positions are 1-based (first nucleotide = position 1)
  - C-terminal insertions (plus strand):
    * pos marks last coding base before stop codon
    * Upstream arm:  bases [pos-window, pos-1]   → UPPERCASE (coding sequence)
    * Downstream arm: bases [pos, pos+window-1]   → lowercase (stop codon + 3' UTR)
  - N-terminal insertions (plus strand):
    * pos marks last non-coding base before start codon
    * Upstream arm:  bases [pos-window+1, pos]   → lowercase (5' UTR)
    * Downstream arm: bases [pos+1, pos+window]   → UPPERCASE (start codon + CDS)
  - Minus strand genes:
    * Coordinates adjusted for reverse complement
    * After reverse complement, case boundaries align with start/stop codons

Example (C-terminal, plus strand, pos=999, window=100):
  Genomic context:
    Position: ...997 998 999 | 1000 1001 1002...
    Bases:    ... G   C   A  |  T    A    A  ...
                          ^
                  Last coding base (pos=999)
  
  Output sequences:
    upstream_arm:  GCAGCT...GCA (100bp, bases 900-999, UPPERCASE = coding)
    downstream_arm: taa...tgcatc (100bp, bases 1000-1099, lowercase = stop + 3'UTR)

Processing Steps:
  1. Read input CSV with insertion site coordinates
  2. For each site, extract window-sized sequences upstream and downstream
  3. Apply strand-aware reverse complement if gene is on minus strand
  4. Format with case to indicate coding boundaries
  5. Write output CSV with original data plus genomic sequences

Required arguments:
  --input PATH      Input CSV with genomic coordinates from script 3
  --output PATH     Output CSV path with added sequence columns
  --genome PATH     Path to local genome FASTA file (auto-creates .fai index)

Optional arguments:
  --window N        Bases to extract on each side of insertion (default: 100)
  --flush-every N   Write to disk every N rows for progress tracking (default: 1000)
  --no-resume       Overwrite existing output instead of resuming

Example:
  python 4_getgenomicsequencearoundinsertion.py \
    --input 3.allguideseqs.csv \
    --output 4.allguidesandgenomicseqs.csv \
    --genome wbcel235/ce11.fa \
    --window 100 \
    --flush-every 10000
"""

import argparse, sys, os, time
import sys
import pandas as pd

try:
    from pyfaidx import Fasta
except ImportError:
    print("ERROR: pyfaidx not installed. Install with: pip install pyfaidx", file=sys.stderr)
    sys.exit(1)

# ===============================================================================
# Local Genome Access
# ===============================================================================

class GenomeReader:
    """Handle local genome FASTA file access."""
    
    def __init__(self, fasta_path: str):
        """Initialize genome reader."""
        print(f"[genome] Loading genome from {fasta_path}...", file=sys.stderr)
        self.genome = Fasta(fasta_path, as_raw=True, sequence_always_upper=True)
        chroms = list(self.genome.keys())
        print(f"[genome] Loaded {len(chroms)} sequences", file=sys.stderr)
    
    def fetch_sequence(self, chrom: str, start: int, end: int, strand: str = '+') -> str:
        """
        Fetch genomic sequence.
        
        Args:
            chrom: Chromosome name
            start: Start position (1-based, inclusive)
            end: End position (1-based, inclusive)
            strand: '+' for forward, '-' for reverse
        
        Returns:
            DNA sequence (uppercase)
        """
        if start < 1:
            start = 1
        if end < start:
            end = start
        
        # Try to find chromosome
        chrom_key = self._find_chromosome(chrom)
        if not chrom_key:
            print(f"  [genome] WARNING: Chromosome not found: {chrom}", file=sys.stderr)
            return ""
        
        try:
            # pyfaidx uses 0-based half-open intervals [start, end)
            seq_obj = self.genome[chrom_key][start-1:end]
            
            # Get sequence string
            if isinstance(seq_obj, str):
                seq = seq_obj
            elif hasattr(seq_obj, 'seq'):
                seq = seq_obj.seq
            else:
                seq = str(seq_obj)
            
            # Reverse complement if on minus strand
            if strand == '-':
                seq = self._reverse_complement(seq)
            
            return seq.upper()
            
        except Exception as e:
            print(f"  [genome] ERROR fetching {chrom}:{start}-{end}: {e}", file=sys.stderr)
            return ""
    
    def _find_chromosome(self, chrom: str):
        """Find chromosome key in genome, handling different naming conventions."""
        # Direct match
        if chrom in self.genome:
            return chrom
        
        # Try variations
        variations = [
            chrom,
            chrom.upper(),
            chrom.lower(),
            f"chr{chrom}",
            f"Chr{chrom}",
            f"chromosome.{chrom}",
        ]
        
        for var in variations:
            if var in self.genome:
                return var
        
        # Try partial match
        chrom_lower = chrom.lower()
        for key in self.genome.keys():
            if chrom_lower in key.lower():
                return key
        
        return None
    
    @staticmethod
    def _reverse_complement(seq: str) -> str:
        """Return reverse complement of DNA sequence."""
        complement = str.maketrans("ACGTNacgtn", "TGCANtgcan")
        return seq.translate(complement)[::-1]
    
    def close(self):
        """Close the genome file."""
        if hasattr(self, 'genome'):
            self.genome.close()

# ===============================================================================
# Guide Coordinate Calculation
# ===============================================================================

def find_guide_in_context(genomic_context: str, guide_seq: str) -> dict:
    """
    Find guide sequence within genomic context and calculate coordinates.
    
    The guide_seq includes both the 20bp guide and 3bp PAM (23bp total).
    We need to locate this sequence in the genomic context and determine:
    - Guide coordinates (first 20bp)
    - PAM coordinates (last 3bp)
    - Cut site (3bp upstream of PAM)
    
    Args:
        genomic_context: 200bp sequence with mixed case
        guide_seq: Guide sequence including PAM (e.g., "TTCAGAATGAGCACAACACCTGG")
        site_type: 'Nterm' or 'Cterm'
        strand: '+' or '-' (gene strand)
    
    Returns:
        Dictionary with coordinates (0-based):
        - guide_5p: 5' end of guide (20bp targeting sequence only)
        - guide_3p: 3' end of guide
        - pam_start: Start of PAM
        - pam_end: End of PAM
        - cut_site: Position of cut site
        - guide_only: 20bp guide sequence without PAM
        - guide_strand: '+' or '-' indicating guide orientation
    """
    if not genomic_context or not guide_seq:
        return {
            'guide_5p': None,
            'guide_3p': None,
            'pam_start': None,
            'pam_end': None,
            'cut_site': None,
            'guide_only': '',
            'guide_strand': ''
        }
    
    # Convert context to uppercase for searching
    context_upper = genomic_context.upper()
    guide_upper = guide_seq.upper()
    
    # Try to find guide+PAM in context (forward orientation)
    pos = context_upper.find(guide_upper)
    found_as_rc = False
    
    # If not found, try reverse complement
    if pos == -1:
        def reverse_complement(seq):
            complement = str.maketrans("ACGTN", "TGCAN")
            return seq.translate(complement)[::-1]
        
        guide_rc = reverse_complement(guide_upper)
        pos = context_upper.find(guide_rc)
        
        if pos != -1:
            found_as_rc = True
        else:
            # Guide not found in either orientation
            return {
                'guide_5p': None,
                'guide_3p': None,
                'pam_start': None,
                'pam_end': None,
                'cut_site': None,
                'guide_only': '',
                'guide_strand': ''
            }
    
    # Calculate coordinates
    # Guide+PAM is 23bp: 20bp guide + 3bp PAM
    if found_as_rc:
        # Guide is on reverse strand relative to context
        # PAM is at the 5' end (upstream in guide orientation, which is UPSTREAM in context)
        pam_start = pos + 2
        pam_end = pos
        guide_start = pos + 3
        guide_end = pos + 22
        
        guide_strand = '-'
        guide_5p = guide_end  # 5' end is at higher coordinate
        guide_3p = guide_start  # 3' end is at lower coordinate
        
        # Cut site is 3bp upstream of PAM in guide strand orientation
        # which is 3bp to the UPSTREAM (lower coordinates) in context
        cut_site = pam_start + 3
        
    else:
        # Guide is on forward strand relative to context
        pam_start = pos + 20
        pam_end = pos + 22
        guide_start = pos
        guide_end = pos + 19
        
        guide_strand = '+'
        guide_5p = guide_start
        guide_3p = guide_end
        # Cut site is 3bp upstream of PAM
        cut_site = pam_start - 3

    # Extract guide sequence (20bp, no PAM)
    guide_only = context_upper[guide_start:guide_end + 1]
    
    return {
        'guide_5p': guide_5p,
        'guide_3p': guide_3p,
        'pam_start': pam_start,
        'pam_end': pam_end,
        'cut_site': cut_site,
        'guide_only': guide_only,
        'guide_strand': guide_strand
    }

# ===============================================================================
# Sequence Extraction
# ===============================================================================

def extract_insertion_context(genome: GenomeReader, chrom: str, pos: int, 
                              strand: str, site_type: str, 
                              window: int = 100) -> str:
    """
    Extract genomic sequence around insertion site with case formatting.
    
    For plus strand genes:
    - C-terminal: pos marks last coding base; stop codon starts at pos (in lowercase)
    - N-terminal: pos marks position before start; start codon at pos+1 (in uppercase)
    
    For minus strand genes:
    - The insertion position 'pos' represents a genomic coordinate
    - After reverse complementation, the case boundary falls at the correct position
      relative to the gene's 5'->3' orientation
    - For N-term: start codon at the beginning of uppercase
    - For C-term: stop codon at the end of uppercase (beginning of lowercase after RC)
    
    Args:
        genome: GenomeReader instance
        chrom: Chromosome name
        pos: Insertion position (1-based)
        strand: Gene strand ('+' or '-')
        site_type: 'N' for N-terminal, 'C' for C-terminal
        window: Bases to extract on each side (default 100)
    
    Returns:
        Formatted sequence with case indicating insertion position
        - N-terminal: lowercase upstream, UPPERCASE downstream
        - C-terminal: UPPERCASE upstream, lowercase downstream
    """
    
    print(f"  [extract] Insertion pos: {pos}, strand: {strand}, site_type: {site_type}", file=sys.stderr)
    
    if strand == '+':
        # Plus strand coordinate calculation
        
        site_type_upper = site_type.upper() if isinstance(site_type, str) else 'C'
        
        if site_type_upper == 'CTERM':
            # C-terminal: stop codon should be entirely in lowercase (downstream arm)
            # pos marks the last coding base, so upstream arm ends at pos-1, downstream starts at pos
            # This ensures the stop codon (starting at genomic pos after last coding base)
            # falls entirely in the lowercase region
            genomic_upstream_start = pos - window
            genomic_upstream_end = pos - 1
            genomic_downstream_start = pos
            genomic_downstream_end = pos + window - 1
        else:
            # N-terminal: start codon should be at start of uppercase (downstream arm)
            # pos marks position before start codon
            # Upstream arm (lowercase) ends at pos, downstream arm (uppercase) starts at pos+1
            genomic_upstream_start = pos - window + 1
            genomic_upstream_end = pos
            genomic_downstream_start = pos + 1
            genomic_downstream_end = pos + window
        
        print(f"  [extract] Plus strand ({site_type}) - genomic upstream: {genomic_upstream_start}-{genomic_upstream_end}, downstream: {genomic_downstream_start}-{genomic_downstream_end}", file=sys.stderr)
        
        upstream = genome.fetch_sequence(chrom, genomic_upstream_start, genomic_upstream_end, strand='+')
        downstream = genome.fetch_sequence(chrom, genomic_downstream_start, genomic_downstream_end, strand='+')
        
    else:
        # Minus strand:
        # 
        # For minus strand genes, the insertion position 'pos' from script 1
        # is calculated as codon_end_gen + 1 (for N-term) or codon_start_gen - 1 (for C-term)
        # 'pos' marks the BOUNDARY, but it's INCLUSIVE of the coding sequence
        #
        # After reverse complement, we need:
        #   upstream (5' in gene) comes from higher genomic coords
        #   downstream (3' in gene) comes from lower genomic coords
        #
        # CORRECT coordinate calculation for minus strand:
        # - The genomic region that becomes the upstream homology arm (gene's 5') should START at pos
        # - The genomic region that becomes the downstream homology arm (gene's 3') should END at pos-1
        
        # Fetch genomic regions (will be reverse complemented)
        genomic_upstream_start = pos  # This will become gene's 5' end (upstream) - INCLUDES pos
        genomic_upstream_end = pos + window - 1
        genomic_downstream_start = pos - window  # This will become gene's 3' end (downstream)
        genomic_downstream_end = pos - 1  # EXCLUDES pos
        
        print(f"  [extract] Minus strand - genomic coords for upstream arm (5'): {genomic_upstream_start}-{genomic_upstream_end}", file=sys.stderr)
        print(f"  [extract] Minus strand - genomic coords for downstream arm (3'): {genomic_downstream_start}-{genomic_downstream_end}", file=sys.stderr)
        
        upstream = genome.fetch_sequence(chrom, genomic_upstream_start, genomic_upstream_end, strand='-')
        downstream = genome.fetch_sequence(chrom, genomic_downstream_start, genomic_downstream_end, strand='-')
    
    if not upstream or not downstream:
        print(f"  [warn] Failed to fetch sequences", file=sys.stderr)
        return ""
    
    # Verify we got expected lengths
    if len(upstream) != window:
        print(f"  [warn] Upstream arm: {len(upstream)}bp vs expected {window}bp", 
              file=sys.stderr)
    if len(downstream) != window:
        print(f"  [warn] Downstream arm: {len(downstream)}bp vs expected {window}bp", 
              file=sys.stderr)
    
    # Format based on site type
    site_type_upper = site_type.upper() if isinstance(site_type, str) else 'C'
    
    print(f"  [extract] Site type: '{site_type}' -> '{site_type_upper}'", file=sys.stderr)

    if site_type_upper == 'NTERM':
        # N-terminal: lowercase UPSTREAM, uppercase DOWNSTREAM
        formatted = upstream.lower() + downstream.upper()
        print(f"  [extract] N-terminal formatting: lowercase(upstream) + UPPERCASE(downstream)", file=sys.stderr)
        print(f"  [extract]   Boundary: ...{formatted[97:100]}|{formatted[100:103]}... (start codon in uppercase)", file=sys.stderr)
    else:
        # C-terminal: uppercase UPSTREAM, lowercase DOWNSTREAM
        formatted = upstream.upper() + downstream.lower()
        print(f"  [extract] C-terminal formatting: UPPERCASE(upstream) + lowercase(downstream)", file=sys.stderr)
        print(f"  [extract]   Boundary: ...{formatted[97:100]}|{formatted[100:103]}... (stop codon in lowercase)", file=sys.stderr)
    
    print(f"  [extract] Final sequence: {len(formatted)}bp", file=sys.stderr)
    print(f"  [extract]   First 30bp: {formatted[:30]}", file=sys.stderr)
    print(f"  [extract]   Last 30bp: {formatted[-30:]}", file=sys.stderr)
    
    return formatted

# ===============================================================================
# Row Processing
# ===============================================================================

def process_row(genome: GenomeReader, row: pd.Series, window: int = 100) -> dict:
    """Process a single row to extract genomic context and find guide."""
    result = {
        'genomic_context': '',
        'genomic_context_coords': '',
        'guide_5p_coord': None,
        'guide_3p_coord': None,
        'pam_coords': '',
        'cut_site_coord': None,
        'guide_sequence': '',
        'guide_strand': '',
        'verify_guide_seq': '',
        'verify_pam_seq': '',
        'guide_position': '',
        'pam_distance_category': ''
    }
    
    try:
        # Extract required fields
        chrom = str(row.get('chrom', '')).strip()
        pos = row.get('pos')
        strand = str(row.get('strand', '+')).strip()
        site_type = str(row.get('site_type', 'C')).strip()
        guide_seq = str(row.get('guide_seq', '')).strip()
        
        # Validate
        if not chrom or pd.isna(pos):
            print(f"  [skip] Missing chrom or pos", file=sys.stderr)
            return result
        
        pos = int(pos)
        
        # Extract genomic context
        context = extract_insertion_context(
            genome, chrom, pos, strand, site_type, window=window
        )
        
        result["genomic_context"] = context
        
        # Calculate absolute genomic coordinates of the 200bp extracted region
        if context:
            if strand == '+':
                # Plus strand coordinates
                site_type_upper = site_type.upper() if isinstance(site_type, str) else 'C'
                if site_type_upper == 'CTERM':
                    context_start = pos - window
                    context_end = pos + window - 1
                else:
                    context_start = pos - window + 1
                    context_end = pos + window
            else:
                # Extraction uses: upstream from pos to pos+window-1, downstream from pos-window to pos-1
                # So the full region is pos-window to pos+window-1
                context_start = pos - window
                context_end = pos + window - 1
            result["genomic_context_coords"] = f"{chrom}:{context_start}-{context_end}"
        
        if context and guide_seq:
            # Find guide in context
            coords = find_guide_in_context(context, guide_seq)
            
            if coords['guide_5p'] is not None:
                result["guide_5p_coord"] = coords['guide_5p']
                result["guide_3p_coord"] = coords['guide_3p']
                result["pam_coords"] = f"{coords['pam_start']}-{coords['pam_end']}"
                result["cut_site_coord"] = coords['cut_site']
                result["guide_sequence"] = coords['guide_only']
                result["guide_strand"] = coords['guide_strand']
                
                # Extract guide and PAM for verification
                if coords['guide_strand'] == '+':
                    # Forward: extract upstream to downstream
                    verify_guide = context[coords['guide_5p']:coords['guide_3p'] + 1]
                    verify_pam = context[coords['pam_start']:coords['pam_end'] + 1]
                else:
                    # Reverse: extract downstream to upstream
                    verify_guide = context[coords['guide_3p']:coords['guide_5p'] + 1]
                    verify_pam = context[coords['pam_end']:coords['pam_start'] + 1]
                
                result["verify_guide_seq"] = verify_guide
                result["verify_pam_seq"] = verify_pam
                
                # Categorize guide position relative to insertion site
                # Insertion point is at position 'window' in the (2*window)bp frame
                insertion_pos = window
                
                # Get the span of the guide
                guide_min = min(coords['guide_5p'], coords['guide_3p'])
                guide_max = max(coords['guide_5p'], coords['guide_3p'])
                
                if guide_max < insertion_pos:
                    # Guide entirely on upstream side
                    guide_position = "left_only"
                elif guide_min >= insertion_pos:
                    # Guide entirely on downstream side
                    guide_position = "right_only"
                else:
                    # Guide spans the insertion point
                    guide_position = "overlap"
                
                result["guide_position"] = guide_position
                
                # Calculate PAM distance from insertion site
                pam_min = min(coords['pam_start'], coords['pam_end'])
                pam_max = max(coords['pam_start'], coords['pam_end'])
                
                # Distance is the closest point of PAM to insertion
                if pam_max < insertion_pos:
                    # PAM is on upstream side
                    pam_distance = insertion_pos - pam_max
                elif pam_min >= insertion_pos:
                    # PAM is on downstream side
                    pam_distance = pam_min - insertion_pos
                else:
                    # PAM spans insertion (shouldn't happen normally)
                    pam_distance = 0
                
                # Categorize based on 35bp threshold
                if pam_distance <= 35:
                    pam_distance_category = "<=35bp"
                else:
                    pam_distance_category = ">35bp"
                
                result["pam_distance_category"] = pam_distance_category
        
        # Log extraction
        gene = row.get('gene', 'unknown')
        print(f"  [{gene}] {chrom}:{pos} ({strand} strand, {site_type}) -> {len(context)}bp", 
              file=sys.stderr)
        if result["genomic_context_coords"]:
            print(f"    Context region: {result['genomic_context_coords']}", 
                  file=sys.stderr)
        if result["guide_sequence"]:
            print(f"    Guide: {result['guide_sequence']} (pos {result['guide_5p_coord']}-{result['guide_3p_coord']}, {result['guide_strand']} strand)", 
                  file=sys.stderr)
            print(f"    PAM: {result['pam_coords']}, Cut: {result['cut_site_coord']}", 
                  file=sys.stderr)
            print(f"    Position: {result['guide_position']}, PAM distance: {result['pam_distance_category']}", 
                  file=sys.stderr)
            print(f"    Verify - Guide: {result['verify_guide_seq']}, PAM: {result['verify_pam_seq']}", 
                  file=sys.stderr)
        
    except Exception as e:
        print(f"  [error] {type(e).__name__}: {e}", file=sys.stderr)
    
    return result

def _write_results_chunk(output_path: str, df_chunk: pd.DataFrame, new_columns: list):
    """Write a chunk of results to the output file (append mode)."""
    if df_chunk.empty:
        return
    
    # Check if file exists and has content (to determine if we need header)
    write_header = not (os.path.exists(output_path) and os.path.getsize(output_path) > 0)
    
    # Write to CSV in append mode
    df_chunk.to_csv(output_path, mode='a', header=write_header, index=False)


# ===============================================================================
# Main
# ===============================================================================

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Extract genomic context around insertion sites"
    )
    parser.add_argument("--input", required=True, help="Input CSV file")
    parser.add_argument("--output", required=True, help="Output CSV file")
    parser.add_argument("--genome", required=True, help="Genome FASTA file")
    parser.add_argument("--window", type=int, default=100, 
                       help="Window size (bp) on each side (default: 100)")
    parser.add_argument("--flush-every", type=int, default=1000,
                       help="Write results to disk every N rows (default: 1000)")
    parser.add_argument("--no-resume", action="store_true",
                       help="Do not resume; overwrite output file")
    
    args = parser.parse_args()
    
    t0 = time.time()
    
    print(f"[start] Extracting genomic contexts")
    print(f"[input] {args.input}")
    print(f"[output] {args.output}")
    print(f"[genome] {args.genome}")
    print(f"[flush-every] {args.flush_every} rows")
    
    # Load input
    try:
        df = pd.read_csv(args.input)
        print(f"[input] Loaded {len(df)} rows")
    except Exception as e:
        print(f"[error] Failed to read input: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Check required columns
    required = ['chrom', 'pos', 'strand', 'site_type']
    missing = [col for col in required if col not in df.columns]
    if missing:
        print(f"[error] Missing required columns: {missing}", file=sys.stderr)
        print(f"[error] Available columns: {list(df.columns)}", file=sys.stderr)
        sys.exit(1)
    
    # Initialize or clear output file
    out_path = args.output
    if args.no_resume and os.path.exists(out_path):
        print(f"[no-resume] Removing existing output file", file=sys.stderr)
        os.remove(out_path)
    
    # Determine which rows have already been processed (for resume functionality)
    processed_indices = set()
    if not args.no_resume and os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        try:
            existing_df = pd.read_csv(out_path)
            # Assuming the output preserves the order of the input, count how many rows exist
            processed_indices = set(range(len(existing_df)))
            print(f"[resume] Found {len(existing_df)} existing rows; will skip those", file=sys.stderr)
        except Exception as e:
            print(f"[resume] Could not read existing output: {e}", file=sys.stderr)
            print(f"[resume] Starting from scratch", file=sys.stderr)
    
    # Define the column order for output
    new_columns_ordered = [
        'genomic_context_coords',
        'genomic_context',
        'guide_5p_coord',
        'guide_3p_coord',
        'pam_coords',
        'cut_site_coord',
        'guide_sequence',
        'guide_strand',
        'verify_guide_seq',
        'verify_pam_seq',
        'guide_position',
        'pam_distance_category'
    ]
    
    # Load genome
    try:
        genome = GenomeReader(args.genome)
    except Exception as e:
        print(f"[error] Failed to load genome: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Process each row with incremental flushing
    print(f"\n[process] Extracting sequences...")
    buffer = []  # Buffer to accumulate rows before flushing
    total = len(df)
    processed_count = 0
    
    for idx, row in df.iterrows():
        # Skip if already processed (resume functionality)
        if idx in processed_indices:
            continue
        
        print(f"\n[{idx+1}/{total}]", end=" ", file=sys.stderr)
        
        # Process the row
        result = process_row(genome, row, window=args.window)
        
        # Create a copy of the original row and add new columns
        row_dict = row.to_dict()
        for col in new_columns_ordered:
            row_dict[col] = result.get(col, '')
        
        buffer.append(row_dict)
        processed_count += 1
        
        # Flush buffer to disk periodically
        if len(buffer) >= args.flush_every:
            print(f"\n[flush] Processed {processed_count} rows, writing {len(buffer)} to disk...", file=sys.stderr)
            buffer_df = pd.DataFrame(buffer)
            _write_results_chunk(out_path, buffer_df, new_columns_ordered)
            buffer.clear()
    
    # Final flush for remaining rows
    if buffer:
        print(f"\n[flush] Final flush: writing {len(buffer)} rows to disk...", file=sys.stderr)
        buffer_df = pd.DataFrame(buffer)
        _write_results_chunk(out_path, buffer_df, new_columns_ordered)
        buffer.clear()
    
    # Close genome
    genome.close()
    
    # Read back the complete output for summary statistics
    try:
        final_df = pd.read_csv(out_path)
        print(f"\n\n[done] Wrote {len(final_df)} rows to {args.output}")
        
        # Summary
        successful = (final_df['genomic_context'] != "").sum()
        guides_found = (final_df['guide_sequence'] != "").sum() if 'guide_seq' in final_df.columns else 0
        print(f"[summary] Successfully extracted: {successful}/{len(final_df)} sequences")
        if 'guide_seq' in final_df.columns:
            print(f"[summary] Guide sequences found: {guides_found}/{len(final_df)}")
        
        # Timing
        elapsed = time.time() - t0
        rate = processed_count / elapsed if elapsed > 0 else 0
        print(f"[timing] Processed {processed_count} new rows in {elapsed:.2f}s ({rate:.2f} rows/s)")
        
    except Exception as e:
        print(f"\n[warning] Could not read output file for summary: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()