#!/usr/bin/env python3
"""
Bonus Script - Inner Primer Alignment Visualization

Generates visual alignments showing genomic context with guide RNA and inner primers
(P2, P3) for CRISPR knock-in experiments. Creates one PNG image per site, displaying 
nucleotide-level alignments with highlighted silent mutations to primers and the PAM.

Input:
  - CSV file containing:
      • genomic_context: Nucleotide sequence with CDS (proximal to insertion site) in uppercase 
      • guide_seq: Guide RNA sequence (20bp + 3bp PAM)
      • primer_2_seq, primer_3_seq: Inner primer sequences
      • primer_2a_seq, primer_3a_seq: Optional nested primers
      • gene: Gene name for labeling
      • site_type: 'Nterm' or 'Cterm' for terminal targeting
      • mutations_detail: Codon changes (e.g., "TTC->TTT(F); AAT->AAC(N)")
      • primer_2_has_mutations, primer_3_has_mutations: Boolean flags

Output:
  - PNG images (one per site) showing:
      • Genomic context sequence (black text)
      • Guide RNA (20bp) with PAM site highlighted in red in genomic sequence
      • Inner primers P2 and P3 (blue and green respectively)
      • Nested primers P2A and P3A if present
      • Pink highlighting for mutations in primers
      • Codon change annotations above mutated primers

Required arguments:
  --input PATH   Input CSV file with guide and primer data
  --outputdir PATH   Directory where PNG images will be saved

Optional arguments:
 --overwrite    Overwrite any already existing images in the output directory

Example:
  python 8.innerprimeralignments.py \
    --input 5.allwormguidesinternalprimers.csv \
    --outputdir alignments/

Output files named: {gene}_{site_type}_{guide_num}.png
  e.g., "unc-54_Nterm_1.png", "dpy-10_Cterm_2.png"
"""

import os, sys, argparse, csv
import matplotlib.pyplot as plt

# ============================================================
# Color customization
# ============================================================
PAM_COLOR = '#FF6B6B'          # Red background for PAM site
MISMATCH_COLOR = '#FFB6C1'     # Pink background for mismatches

# Text colors for sequences
GUIDE_TEXT_COLOR = '#000000'     # Guide RNA letter color
PRIMER2_TEXT_COLOR = '#0000FF'   # Primer 2 letter color
PRIMER3_TEXT_COLOR = '#006000'   # Primer 3 letter color
GENOMIC_TEXT_COLOR = '#000000'   # Genomic context letter color
# ============================================================

def find_sequence_in_context(guide, context, allow_mismatches=False, min_match_ratio=0.7):
    """Find where a sequence appears in the genomic context."""
    guide_upper = guide.upper()
    context_upper = context.upper()
    
    if not allow_mismatches:
        pos = context_upper.find(guide_upper)
        if pos != -1:
            return pos, pos + len(guide), 'forward'
        
        rev_comp = reverse_complement(guide_upper)
        pos = context_upper.find(rev_comp)
        if pos != -1:
            return pos, pos + len(rev_comp), 'reverse'
    else:
        best_match_pos = -1
        best_match_count = 0
        best_orientation = None
        min_matches = int(len(guide) * min_match_ratio)
        
        for i in range(len(context) - len(guide) + 1):
            matches = sum(1 for j in range(len(guide)) 
                         if context_upper[i+j] == guide_upper[j])
            if matches > best_match_count and matches >= min_matches:
                best_match_count = matches
                best_match_pos = i
                best_orientation = 'forward'
        
        rev_comp = reverse_complement(guide_upper)
        for i in range(len(context) - len(guide) + 1):
            matches = sum(1 for j in range(len(guide)) 
                         if context_upper[i+j] == rev_comp[j])
            if matches > best_match_count and matches >= min_matches:
                best_match_count = matches
                best_match_pos = i
                best_orientation = 'reverse'
        
        if best_match_pos != -1:
            return best_match_pos, best_match_pos + len(guide), best_orientation
    
    return None, None, None

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement.get(base.upper(), base) for base in seq[::-1])

def draw_sequence(ax, sequence, y_pos, start_pos, genomic_context, label, 
                  check_mismatches=False, label_position='above', text_color='black'):
    """Draw a sequence with appropriate coloring and return mismatch positions.
    
    label_position options: 'above', 'below', 'left', 'right'
    """
    mismatch_positions = []
    
    for j, base in enumerate(sequence):
        i = start_pos + j
        bg_color = 'none'
        
        if check_mismatches:
            if i < len(genomic_context) and base.upper() != genomic_context[i].upper():
                bg_color = MISMATCH_COLOR
                mismatch_positions.append(i)
        
        if bg_color != 'none':
            ax.text(i + 0.5, y_pos, base.upper(), ha='center', va='center', 
                   fontsize=9, fontweight='bold', fontfamily='monospace', color=text_color,
                   bbox=dict(boxstyle='square,pad=0.25', facecolor=bg_color, edgecolor='none'))
        else:
            ax.text(i + 0.5, y_pos, base.upper(), ha='center', va='center', 
                   fontsize=9, fontweight='bold', fontfamily='monospace', color=text_color)
    
    if label_position == 'above':
        label_y = y_pos + 0.25
        label_x = start_pos + len(sequence)/2
        ha = 'center'
    elif label_position == 'below':
        label_y = y_pos - 0.25
        label_x = start_pos + len(sequence)/2
        ha = 'center'
    elif label_position == 'left':
        label_y = y_pos
        label_x = start_pos - 1.0
        ha = 'right'
    elif label_position == 'right':
        label_y = y_pos
        label_x = start_pos + len(sequence) + 1.0
        ha = 'left'
    else:
        label_y = y_pos + 0.25
        label_x = start_pos + len(sequence)/2
        ha = 'center'
    
    ax.text(label_x, label_y, label, ha=ha, va='center', 
            fontsize=9, fontweight='bold')
    
    return mismatch_positions

def highlight_pam_in_genomic(ax, pam_start, pam_end, genomic_context, y_pos):
    """Highlight PAM site in the genomic context with red background."""
    for i in range(pam_start, pam_end):
        if i < len(genomic_context):
            base = genomic_context[i]
            ax.text(i + 0.5, y_pos, base, ha='center', va='center', 
                   fontsize=9, fontweight='bold', fontfamily='monospace', 
                   color=GENOMIC_TEXT_COLOR,
                   bbox=dict(boxstyle='square,pad=0.25', facecolor=PAM_COLOR, edgecolor='none'))

def annotate_codon_changes_from_csv(ax, primer_start, primer_end, annotation_y, mutations_detail_str):
    """Annotate codon changes using the mutations_detail field from the CSV."""
    if not mutations_detail_str or mutations_detail_str.strip() == '':
        return
    
    mutations = [m.strip() for m in mutations_detail_str.split(';')]
    primer_length = primer_end - primer_start
    primer_center = primer_start + primer_length / 2
    num_mutations = len(mutations)
    
    if num_mutations == 0:
        return
    
    row_spacing = 0.15
    
    for idx, mutation in enumerate(mutations):
        if '->' not in mutation:
            continue
            
        if num_mutations == 1:
            x_pos = primer_center
        else:
            spacing = primer_length / (num_mutations + 1)
            x_pos = primer_start + spacing * (idx + 1)
        
        y_pos = annotation_y + (idx * row_spacing)
        
        ax.text(x_pos, y_pos, mutation, ha='center', va='center',
               fontsize=7, color='darkred', fontweight='bold')

def create_alignment_image(row, output_path, title=None):
    """Create an alignment visualization for a single row."""
    genomic_context = row['genomic_context']
    guide_seq_with_pam = row['guide_seq']
    primer_2 = row['primer_2_seq']
    primer_3 = row['primer_3_seq']
    gene = row['gene']
    
    # Remove last 3 bases (PAM) from guide sequence for display
    guide_seq = guide_seq_with_pam[:-3] if len(guide_seq_with_pam) > 3 else guide_seq_with_pam
    
    primer_2_has_mutations = row.get('primer_2_has_mutations', 'False').lower() == 'true'
    primer_3_has_mutations = row.get('primer_3_has_mutations', 'False').lower() == 'true'
    
    primer_2a = row.get('primer_2a_seq', 'N/A')
    primer_3a = row.get('primer_3a_seq', 'N/A')
    has_primer_2a = primer_2a != 'N/A' and primer_2a != ''
    has_primer_3a = primer_3a != 'N/A' and primer_3a != ''
    
    # Find guide sequence (without PAM) in genomic context
    guide_start, guide_end, guide_orient = find_sequence_in_context(guide_seq, genomic_context, 
                                                                     allow_mismatches=False)
    
    # Determine PAM location based on orientation
    if guide_start is not None:
        if guide_orient == 'forward':
            # PAM is after the guide (3' end)
            pam_start = guide_end
            pam_end = guide_end + 3
        else:
            # PAM is before the guide (which is the 3' end in reverse)
            pam_start = guide_start - 3
            pam_end = guide_start
    else:
        pam_start = None
        pam_end = None
    
    primer2_start, primer2_end, primer2_orient = find_sequence_in_context(primer_2, genomic_context, 
                                                                           allow_mismatches=primer_2_has_mutations)
    primer3_start, primer3_end, primer3_orient = find_sequence_in_context(primer_3, genomic_context, 
                                                                           allow_mismatches=primer_3_has_mutations)
    
    primer2a_start, _, primer2a_orient = None, None, None
    primer3a_start, _, primer3a_orient = None, None, None
    
    if has_primer_2a:
        primer2a_start, _, primer2a_orient = find_sequence_in_context(primer_2a, genomic_context, 
                                                                                  allow_mismatches=True)
    if has_primer_3a:
        primer3a_start, _, primer3a_orient = find_sequence_in_context(primer_3a, genomic_context, 
                                                                                  allow_mismatches=True)
    
    num_primer_rows = 2
    if has_primer_2a:
        num_primer_rows += 1
    if has_primer_3a:
        num_primer_rows += 1
    
    fig_height = 2.5 + (num_primer_rows * 0.3)
    
    _, ax = plt.subplots(figsize=(20, fig_height))
    ax.set_xlim(-1, len(genomic_context) + 1)
    ax.set_ylim(-0.5, 1.2)
    ax.axis('off')
    
    if has_primer_2a:
        primer2a_y = 0.2
        primer2_y = 0.4
    else:
        primer2a_y = None
        primer2_y = 0.2
    
    if has_primer_3a:
        primer3a_y = 0.2
        primer3_y = 0.4
    else:
        primer3a_y = None
        primer3_y = 0.2
    
    genomic_y = 0.0
    guide_y = -0.2
    
    if has_primer_3a and primer3a_start is not None:
        display_seq = reverse_complement(primer_3a) if primer3a_orient == 'reverse' else primer_3a
        draw_sequence(ax, display_seq, primer3a_y, primer3a_start, genomic_context, 
                     'Primer 3A', check_mismatches=True, label_position='right',
                     text_color=PRIMER3_TEXT_COLOR)
    
    if primer3_start is not None:
        label = 'Primer 3B' if has_primer_3a else 'Primer 3'
        display_seq = reverse_complement(primer_3) if primer3_orient == 'reverse' else primer_3
        primer3_mismatches = draw_sequence(ax, display_seq, primer3_y, primer3_start, genomic_context, 
                     label, check_mismatches=True, label_position='right',
                     text_color=PRIMER3_TEXT_COLOR)
        if primer3_mismatches and row.get('mutations_detail'):
            annotate_codon_changes_from_csv(ax, primer3_start, primer3_end, 
                                          primer3_y + 0.15, row.get('mutations_detail', ''))
    
    if has_primer_2a and primer2a_start is not None:
        display_seq = reverse_complement(primer_2a) if primer2a_orient == 'reverse' else primer_2a
        draw_sequence(ax, display_seq, primer2a_y, primer2a_start, genomic_context, 
                     'Primer 2A', check_mismatches=True, label_position='left',
                     text_color=PRIMER2_TEXT_COLOR)
    
    if primer2_start is not None:
        label = 'Primer 2B' if has_primer_2a else 'Primer 2'
        display_seq = reverse_complement(primer_2) if primer2_orient == 'reverse' else primer_2
        primer2_mismatches = draw_sequence(ax, display_seq, primer2_y, primer2_start, genomic_context, 
                     label, check_mismatches=True, label_position='left',
                     text_color=PRIMER2_TEXT_COLOR)
        if primer2_mismatches and row.get('mutations_detail'):
            annotate_codon_changes_from_csv(ax, primer2_start, primer2_end, 
                                          primer2_y + 0.15, row.get('mutations_detail', ''))
    
    # Draw genomic context (will be overwritten at PAM positions)
    for i, base in enumerate(genomic_context):
        ax.text(i + 0.5, genomic_y, base, ha='center', va='center', 
                fontsize=9, fontweight='bold', fontfamily='monospace', color=GENOMIC_TEXT_COLOR)
    
    # Highlight PAM site in genomic context
    if pam_start is not None and pam_end is not None:
        highlight_pam_in_genomic(ax, pam_start, pam_end, genomic_context, genomic_y)
    
    # Draw Guide RNA (without PAM)
    if guide_start is not None:
        display_seq = reverse_complement(guide_seq) if guide_orient == 'reverse' else guide_seq
        
        # Label position depends on strand
        if guide_orient == 'forward':
            label_pos = 'left'
            pam_label_x = guide_end + 1.5  # PAM label to the right
        else:
            label_pos = 'right'
            pam_label_x = guide_start - 1.5  # PAM label to the left
        
        draw_sequence(ax, display_seq, guide_y, guide_start, genomic_context, 
                     'Guide', label_position=label_pos, text_color=GUIDE_TEXT_COLOR)
        
        # Add PAM label (positioned slightly higher than guide)
        ax.text(pam_label_x, guide_y + 0.07, 'PAM', ha='center', va='center',
                fontsize=9, fontweight='bold', color='#CC0000') 
    
    if title is None:
        title = f'{gene} - Nucleotide Alignment'
    
    plt.title(title, fontsize=14, fontweight='bold', pad=20)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(
        description='Generate nucleotide alignment visualizations from CSV data.'
    )
    parser.add_argument('--input', required=True, help='input CSV file')
    parser.add_argument('--outputdir', required=True, help='Directory where PNG images will be saved')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite existing images instead of skipping them')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' not found!")
        sys.exit(1)
    
    os.makedirs(args.outputdir, exist_ok=True)
    
    guide_counter = {}
    skipped_count = 0
    created_count = 0
    
    with open(args.input, 'r') as f:
        reader = csv.DictReader(f)
        for idx, row in enumerate(reader):
            if not all([row.get('genomic_context'), row.get('guide_seq'), 
                       row.get('primer_2_seq'), row.get('primer_3_seq')]):
                print(f"Skipping row {idx + 1}: missing required fields")
                continue
            
            gene = row.get('gene', f'row_{idx}')
            site_type = row.get('site_type', 'unknown')
            
            if site_type.lower() == 'nterm':
                site_label = 'N-terminal'
            elif site_type.lower() == 'cterm':
                site_label = 'C-terminal'
            else:
                site_label = site_type
            
            key = f"{gene}_{site_type}"
            if key not in guide_counter:
                guide_counter[key] = 0
            guide_counter[key] += 1
            guide_num = guide_counter[key]
            
            filename = f"{gene}_{site_type}_{guide_num}.png"
            title = f"{gene} {site_label} {guide_num}"
            output_path = os.path.join(args.outputdir, filename)
            
            # Check if file already exists
            if os.path.exists(output_path) and not args.overwrite:
                print(f"Skipping {filename} (already exists)")
                skipped_count += 1
                continue
            
            try:
                create_alignment_image(row, output_path, title)
                print(f"Created: {output_path}")
                created_count += 1
            except Exception as e:
                print(f"Error processing row {idx + 1} ({gene}): {e}")
    
    print(f"\n{'='*60}")
    print(f"Summary:")
    print(f"  Created: {created_count} images")
    print(f"  Skipped: {skipped_count} images (already existed)")
    print(f"  All alignments saved to: {args.outputdir}")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()