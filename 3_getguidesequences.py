#!/usr/bin/env python3
"""
sgRNA Design and Selection Pipeline for CRISPR-Mediated Protein Tagging

Identifies and ranks optimal sgRNA sequences for N- and C-terminal knock-in
insertion sites, with off-target filtering and deduplication.

Input:
  - CSV file from script 2 containing insertion site coordinates (requires:
    WBGene, Site, Chromosome, Position, Strand)

Output:
  - CSV file with up to --max-sgrnas sgRNAs per insertion site, containing:
      • Gene identifiers and site metadata
      • sgRNA sequence and PAM
      • Distance from insertion site (bp)
      • Off-target metrics (FlashFry hits and sequences)
      • Ranking information and selection rationale
      • FlashFry scores and genomic coordinates

Processing Strategy:
  1. Sequence extraction: Retrieves genomic context around each insertion site
     (default ±50bp window) from local reference genome using strand-aware coordinates.

  2. sgRNA scanning: Uses FlashFry to identify all possible sgRNA sequences
     matching the specified PAM pattern (default: SpCas9 NGG PAM).

  3. Filtering: Removes sgRNAs with:
     - On-target scores below threshold (--min-score)
     - Homopolymer runs (≥4 identical bases)
     - Extreme GC content (<25% or >80%)

  4. Ranking: Distance-first ranking prioritizes sgRNAs closest to insertion site

  5. Off-target analysis: FlashFry identifies off-target hits with PAM-aware mismatch detection (more accurate than BLAST).

Off-Target Detection:
  - Uses FlashFry for PAM-aware off-target detection
  - Configurable max mismatches in protospacer (--flashfry-max-mismatches, default: 3)
  - Reports hit count and off-target sequences for each sgRNA

Required arguments:
  --input PATH      Input CSV with genomic coordinates from script 3
  --output PATH     Output CSV path with added sequence columns
  --genome PATH     Path to local genome FASTA file (auto-creates .fai index)

Optional arguments:
  --max-sgrnas N                Maximum sgRNAs per site to output (default: 1)
  --pam SEQUENCE                PAM sequence in IUPAC notation (default: NGG for SpCas9). If you use FlashFry to 
                                 identify off-targets, ensure that the FlashFry db is made using the same PAM sequence
  --guide-lens N [N ...]        Guide protospacer lengths to search (default: 20). Can specify multiple: 19 20 21
  --window-bp N                 Half-window around insertion site for sequence extraction (default: 50bp = ±50bp window)
  --max-cut-distance N          Maximum distance between Cas9 cut site and insertion site (default: 50bp)
  --min-gc FRACTION             Minimum GC content fraction (default: 0.25 = 25%)
  --max-gc FRACTION             Maximum GC content fraction (default: 0.80 = 80%)
  --max-candidates N            Internal cap on candidates before off-target testing (default: 50)
  --offtarget-mode MODE         Off-target detection: none|flashfry (default: none)
  --flashfry-db PATH            Path to FlashFry database (required if offtarget-mode=flashfry)
  --flashfry-jar PATH           Path to FlashFry JAR file (default: FlashFry-assembly-1.15.jar)
  --flashfry-max-mismatches N   Max mismatches in protospacer for off-target search (default: 3, stringent)
  --allow-offtargets            Collect all guides up to max-sgrnas, even those with off-targets
  --log-output PATH             Optional detailed log CSV tracking all candidates and rejection reasons
  --flush-every N               Write results to disk every N input sites (default: 10)
  --no-resume                   Overwrite existing output instead of resuming from checkpoint

Example:
  python 3_getguidesequences.py \
    --input 2.insertion_sites.csv \
    --output 3.allguideseqs.csv \
    --genome wbcel235/ce11.fa \
    --max-sgrnas 1 \
    --offtarget-mode flashfry \
    --flashfry-db flashfry_cel_db/ce11_spcas9ngg_db \
    --flashfry-max-mismatches 3 \
    --flush-every 100
"""

from __future__ import annotations
import argparse, os, re, subprocess, sys, time, tempfile, traceback
import pandas as pd
from pyfaidx import Fasta
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple, Union

# ===============================================================================
# Constants
# ===============================================================================

HOMOPOLYMER7_RE = re.compile(r"(A{7,}|C{7,}|G{7,}|T{7,})", re.I)  # Reject ≥7 homopolymers
T4PLUS_RE       = re.compile(r"(T{4,})", re.I)                    # Identify ≥4 T-runs

IUPAC = {
    "A": "[A]","C": "[C]","G": "[G]","T": "[T]",
    "R": "[AG]","Y": "[CT]","S": "[GC]","W": "[AT]","K": "[GT]","M": "[AC]",
    "B": "[CGT]","D": "[AGT]","H": "[ACT]","V": "[ACG]","N": "[ACGT]",
}
IUPAC_COMP = {"A":"T","C":"G","G":"C","T":"A","R":"Y","Y":"R","S":"S","W":"W","K":"M","M":"K","B":"V","D":"H","H":"D","V":"B","N":"N"}

# ===============================================================================
# Local Genome Access
# ===============================================================================

def fetch_region_seq(fasta: Fasta, chrom: str, start: int, end: int) -> str:
    """
    Fetch sequence from local genome FASTA file using pyfaidx.

    Args:
        fasta: pyfaidx Fasta object
        chrom: chromosome name (e.g., 'I', 'II', 'X' or 'chrI', 'chrII')
        start: 1-based start position
        end: 1-based end position (inclusive)

    Returns Uppercase DNA sequence string
    """
    if start < 1:
        start = 1
    if end < start:
        end = start

    def _resolve_chrom(fasta_obj, qchrom):
        # Try exact matches first
        for cand in (qchrom, f"chr{qchrom}" if not qchrom.lower().startswith("chr") else qchrom[3:]):
            try:
                # quick membership check (pyfaidx supports dict-like lookup)
                if cand in fasta_obj:
                    return cand
            except Exception:
                pass
        # Fallback: case-insensitive or suffix match against available contigs
        q_low = qchrom.lower().lstrip("chr")
        for key in fasta_obj.keys():
            if key.lower() == q_low or key.lower().endswith(q_low):
                return key
        return None

    resolved = _resolve_chrom(fasta, chrom)
    if resolved is None:
        raise RuntimeError(f"Chromosome '{chrom}' not found in genome file (tried variants like 'chr{chrom}').")

    try:
        # pyfaidx uses 1-based coordinates; slice via [start-1:end]
        seq = fasta[resolved][start - 1 : end].seq
        return seq.upper()
    except Exception as e:
        raise RuntimeError(f"Failed to fetch {resolved}:{start}-{end}: {e}")



# ===============================================================================
# Sequence Helpers
# ===============================================================================-

def revcomp(seq: str) -> str:
    tb = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(tb)[::-1]

def iupac_to_regex(pam: str) -> str:
    parts = []
    for ch in pam.upper():
        parts.append(IUPAC.get(ch, "[ACGT]"))
    return "".join(parts)

def iupac_revcomp(pam: str) -> str:
    return "".join(IUPAC_COMP.get(ch.upper(), "N") for ch in pam)[::-1]

# ===============================================================================
# Data
# ===============================================================================

@dataclass
class Guide:
    site_name: str
    chrom: str
    strand: str  # '+' or '-'
    protospacer: str
    pam_seq: str
    guide_len: int
    cut_genomic: int
    distance_to_insertion: int
    gc: float
    on_target_score: Optional[float] = None
    off_target_hits: Optional[Union[int, str]] = None  # int for count, str for status messages
    off_target_details: Optional[List[Dict[str, str]]] = None
    transcription_type: Optional[str] = None  # 'in_vivo', 'in_vitro_only', or None

@dataclass
class LogEntry:
    """Tracks each candidate guide and why it was accepted/rejected"""
    wbgene: str
    gene_name: str
    transcript: str
    site_type: str
    chrom: str
    pos: int
    strand: str
    guide_seq: str
    guide_strand: str
    cut_genomic: int
    distance_to_insertion: int
    gc: float
    on_target_score: Optional[float]
    off_target_hits: Optional[Union[int, str]]  # int for count, str for status messages
    status: str  # 'accepted', 'rejected', 'filtered_out'
    rejection_reason: Optional[str]
    final_rank: Optional[int]  # 1, 2, 3, etc. for accepted guides

# ===============================================================================
# Scoring and Off-Targets
# ===============================================================================

def count_mismatches(seq1: str, seq2: str) -> int:
    """Count the number of mismatches between two sequences of equal length."""
    if len(seq1) != len(seq2):
        return -1
    return sum(1 for a, b in zip(seq1, seq2) if a != b)

def heuristic_on_target_score(protospacer: str) -> float:
    seq = protospacer.upper()
    L = len(seq)
    if L == 0 or any(ch not in "ACGT" for ch in seq):
        return 0.0
    gc = (seq.count("G") + seq.count("C")) / L
    score = 60.0 + 40.0 * (1.0 - abs(gc - 0.5) / 0.5)
    if HOMOPOLYMER7_RE.search(seq):
        score -= 20.0
    if "TTTT" in seq or "AAAA" in seq:
        score -= 5.0
    return max(0.0, min(100.0, score))

def guide_query_with_pam(g: Guide) -> str:
    return (g.protospacer + g.pam_seq).upper()

def flashfry_offtargets(query_seq: str, db_path: str, *, max_mismatches: int = 3,
                        flashfry_jar: str = "FlashFry-assembly-1.15.jar",
                        ) -> Tuple[Optional[int], List[Dict[str, str]]]:
    """
    Find off-targets using FlashFry (PAM-aware, more accurate than BLAST)
    
    Args:
        query_seq: Guide sequence (protospacer + PAM, e.g., 23bp for 20bp + NGG)
        db_path: Path to FlashFry database
        max_mismatches: Max mismatches in protospacer only (default 3 for stringency)
        flashfry_jar: Path to FlashFry JAR file
    
    Returns:
        Tuple of (number of TRUE off-target sites (excludes on-target), list of off-target details)
        Each off-target detail dict contains: 'sequence', 'mismatches'
        NOTE: 'mismatches' counts mismatches in PROTOSPACER ONLY (excludes PAM region)
        Returns (None, []) on error
    """

    with tempfile.NamedTemporaryFile('w', delete=False, suffix='.fasta') as fh:
        fasta_path = fh.name
        fh.write(">guide\n")
        fh.write(query_seq + "\n")
    
    output_path = fasta_path + ".output"
    
    try:
        discover_cmd = [
            'java', '-Xmx4g', '-jar', flashfry_jar,
            'discover',
            '--database', db_path,
            '--fasta', fasta_path,
            '--output', output_path,
            '--maxMismatch', str(max_mismatches)
        ]
        
        result = subprocess.run(discover_cmd, capture_output=True, text=True, check=False, timeout=60)
        
        if result.returncode != 0:
            stderr = result.stderr.strip()
            last = stderr.splitlines()[-1] if stderr else "unknown error"
            print(f"    [flashfry] error: {last}")
            return None, []
        
        if not os.path.exists(output_path):
            return None, []
        
        # Parse FlashFry output using otCount and offTargets columns
        off_target_details = []
        total_hits = 0
        
        with open(output_path) as f:
            lines = f.readlines()
            if len(lines) <= 1:  # Only header or empty
                return 0, []
            
            # Parse header to find column indices
            header = lines[0].strip().split('\t')
            try:
                ot_count_idx = header.index('otCount')  # FlashFry's total hit count
                off_targets_idx = header.index('offTargets')  # Off-target sequences
            except ValueError as e:
                print(f"    [flashfry] warning: could not parse header: {e}")
                return None, []
            
            # Parse the first data line (should be our query guide)
            if len(lines) > 1:
                line = lines[1].strip()
                if line:
                    parts = line.split('\t')
                    
                    # Get total hit count from FlashFry
                    try:
                        total_hits = int(parts[ot_count_idx])
                    except (IndexError, ValueError):
                        total_hits = 0
                    
                    # Parse off-target sequences from offTargets column
                    # Format: SEQUENCE_X_Y,SEQUENCE_X_Y where X and Y are FlashFry metrics
                    if off_targets_idx < len(parts):
                        off_targets_str = parts[off_targets_idx]
                        
                        if off_targets_str and off_targets_str not in ['NONE', 'NA', '']:
                            for ot_entry in off_targets_str.split(','):
                                if not ot_entry.strip():
                                    continue
                                
                                # Parse format: SEQUENCE_X_Y
                                ot_parts = ot_entry.strip().rsplit('_', 2)
                                if len(ot_parts) >= 3:
                                    sequence = ot_parts[0]
                                    
                                    # Count mismatches in PROTOSPACER ONLY (not PAM)
                                    # The PAM (last 3bp) is needed for Cas9 recognition, but mismatches
                                    # in the first position (N in NGG) are biologically less significant.
                                    # FlashFry searches based on protospacer mismatches, so we should
                                    # report the same metric for consistency.
                                    protospacer_len = len(query_seq) - 3  # Assume 3bp PAM (NGG)
                                    query_protospacer = query_seq[:protospacer_len]
                                    ot_protospacer = sequence[:protospacer_len]
                                    
                                    actual_mismatches = count_mismatches(query_protospacer, ot_protospacer)
                                    
                                    # Only include true off-targets (not the on-target with 0 mismatches)
                                    if actual_mismatches > 0:
                                        off_target_details.append({
                                            'sequence': sequence,
                                            'mismatches': actual_mismatches
                                        })
        
        # Calculate true off-target count (excluding on-target)
        true_off_targets = len(off_target_details)
        
        print(f"    [flashfry] total genome hits={total_hits}, true off-targets={true_off_targets}")
        if off_target_details:
            print(f"    [flashfry] off-target sequences:")
            for ot in off_target_details[:5]:  # Show first 5
                print(f"      {ot['sequence']} ({ot['mismatches']} mismatches)")
            if len(off_target_details) > 5:
                print(f"      ... and {len(off_target_details) - 5} more")
    
        return true_off_targets, off_target_details
        
    except subprocess.TimeoutExpired:
        print(f"    [flashfry] timeout after 60 seconds")
        return None, []
    except Exception as e:
        print(f"    [flashfry] exception: {e}")
        return None, []
        
    finally:
        for path in [fasta_path, output_path]:
            try:
                if os.path.exists(path):
                    os.unlink(path)
            except:
                pass

def flashfry_offtargets_batch(guides: List[Guide], db_path: str, *, max_mismatches: int = 3,
                              flashfry_jar: str = "FlashFry-assembly-1.15.jar",
                              ) -> Dict[str, Tuple[Optional[int], List[Dict[str, str]]]]:
    """
    Find off-targets for multiple guides in a single FlashFry run (much faster!)
    
    Args:
        guides: List of Guide objects to check
        db_path: Path to FlashFry database
        max_mismatches: Max mismatches in protospacer only (default 3)
        flashfry_jar: Path to FlashFry JAR file
    
    Returns:
        Dictionary mapping guide_sequence -> (off_target_count, off_target_details)
        Returns empty dict on error
    """
    if not guides:
        return {}
    
    # Create a single FASTA with all guides
    with tempfile.NamedTemporaryFile('w', delete=False, suffix='.fasta') as fh:
        fasta_path = fh.name
        for idx, g in enumerate(guides):
            guide_seq = guide_query_with_pam(g)
            fh.write(f">guide_{idx}\n")
            fh.write(guide_seq + "\n")
    
    output_path = fasta_path + ".output"
    
    try:
        discover_cmd = [
            'java', '-Xmx4g', '-jar', flashfry_jar,
            'discover',
            '--database', db_path,
            '--fasta', fasta_path,
            '--output', output_path,
            '--maxMismatch', str(max_mismatches)
        ]
        
        result = subprocess.run(discover_cmd, capture_output=True, text=True, check=False, timeout=300)
        
        if result.returncode != 0:
            stderr = result.stderr.strip()
            last = stderr.splitlines()[-1] if stderr else "unknown error"
            print(f"    [flashfry] batch error: {last}")
            return {}
        
        if not os.path.exists(output_path):
            return {}
        
        # Parse FlashFry output for all guides
        results = {}
        
        with open(output_path) as f:
            lines = f.readlines()
            if len(lines) <= 1:  # Only header or empty
                return {guide_query_with_pam(g): (0, []) for g in guides}
            
            # Parse header
            header = lines[0].strip().split('\t')
            try:
                contig_idx = header.index('contig')  # Guide identifier
                ot_count_idx = header.index('otCount')
                off_targets_idx = header.index('offTargets')
            except ValueError as e:
                print(f"    [flashfry] batch warning: could not parse header: {e}")
                return {}
            
            # Parse each guide's results
            for line in lines[1:]:
                if not line.strip():
                    continue
                    
                parts = line.strip().split('\t')
                
                # Get guide index and sequence
                try:
                    guide_id = parts[contig_idx]
                    guide_idx = int(guide_id.split('_')[1])
                    g = guides[guide_idx]
                    query_seq = guide_query_with_pam(g)
                except (IndexError, ValueError):
                    continue
                
                # Get total hit count
                try:
                    total_hits = int(parts[ot_count_idx])
                except (IndexError, ValueError):
                    total_hits = 0
                
                # Parse off-target sequences
                off_target_details = []
                if off_targets_idx < len(parts):
                    off_targets_str = parts[off_targets_idx]
                    
                    if off_targets_str and off_targets_str not in ['NONE', 'NA', '']:
                        for ot_entry in off_targets_str.split(','):
                            if not ot_entry.strip():
                                continue
                            
                            ot_parts = ot_entry.strip().rsplit('_', 2)
                            if len(ot_parts) >= 3:
                                sequence = ot_parts[0]
                                
                                # Count mismatches in protospacer only
                                protospacer_len = len(query_seq) - 3
                                query_protospacer = query_seq[:protospacer_len]
                                ot_protospacer = sequence[:protospacer_len]
                                
                                actual_mismatches = count_mismatches(query_protospacer, ot_protospacer)
                                
                                # Only include true off-targets
                                if actual_mismatches > 0:
                                    off_target_details.append({
                                        'sequence': sequence,
                                        'mismatches': actual_mismatches
                                    })
                
                # Store results for this guide
                true_off_targets = len(off_target_details)
                results[query_seq] = (true_off_targets, off_target_details)
        
        print(f"    [flashfry] batch processed {len(results)}/{len(guides)} guides")
        
        return results
        
    except subprocess.TimeoutExpired:
        print(f"    [flashfry] batch timeout after 300 seconds")
        return {}
    except Exception as e:
        print(f"    [flashfry] batch exception: {e}")
        return {}
        
    finally:
        for path in [fasta_path, output_path]:
            try:
                if os.path.exists(path):
                    os.unlink(path)
            except:
                pass

# ===============================================================================
# Scan for Guides
# ===============================================================================

def scan_guides_in_window(seq_plus: str, chrom: str, window_start: int, insertion_pos: int,
                          pam: str, guide_lens: Iterable[int], max_cut_distance: int) -> List[Guide]:
    seq = seq_plus.upper()
    pam_regex = iupac_to_regex(pam)
    pam_rev_regex = iupac_to_regex(iupac_revcomp(pam))

    guides: List[Guide] = []
    for L in guide_lens:
        # Fwd: protospacer(L) + PAM
        fwd = re.compile(rf"(?=(?P<prot>[ACGT]{{{L}}})(?P<pam>{pam_regex}))")
        for m in fwd.finditer(seq):
            idx0 = m.start()
            prot = m.group("prot")
            pam_seq = m.group("pam")
            # Exclude GGG PAM and NGGG (extra G immediately after PAM)
            if len(pam_seq) == 3:
                # GGG PAM (N=G)
                if pam_seq == "GGG":
                    continue
                # NGG with a G right after the PAM (i.e., prot + NGGG)
                next_i = idx0 + L + 3
                if pam_seq.endswith("GG") and next_i < len(seq) and seq[next_i] == "G":
                    continue
            cut0 = idx0 + L - 3  # SpCas9
            cut_genomic = window_start + cut0
            dist = abs(cut_genomic - insertion_pos)
            if dist <= max_cut_distance and "N" not in prot:
                guides.append(Guide(
                    site_name="", chrom=chrom, strand="+",
                    protospacer=prot, pam_seq=pam_seq, guide_len=L,
                    cut_genomic=int(cut_genomic), distance_to_insertion=int(dist),
                    gc=(prot.count("G")+prot.count("C"))/L,
                ))
        # Rev: revPAM + prot(L) on + strand; report protospacer in - orientation
        rev = re.compile(rf"(?=(?P<pam>{pam_rev_regex})(?P<prot>[ACGT]{{{L}}}))")
        for m in rev.finditer(seq):
            idx0 = m.start()
            prot_plus = m.group("prot")
            pam_plus = m.group("pam")
            prot = revcomp(prot_plus)
            pam_seq = revcomp(pam_plus)
            # Exclude GGG PAM and NGGG in canonical orientation on the minus strand.
            # Here idx0 points to the start of pam_plus on the + strand.
            # The base immediately AFTER the canonical PAM corresponds to the base
            # immediately BEFORE pam_plus on the + strand (i.e., seq[idx0 - 1]).
            if len(pam_seq) == 3:
                # GGG PAM (N=G) in canonical orientation
                if pam_seq == "GGG":
                    continue
                # NGG with a G right after the PAM (canonical). That base on the + strand is 'C'
                # just before pam_plus (because rc(G) == C).
                prev_i = idx0 - 1
                if pam_seq.endswith("GG") and prev_i >= 0 and seq[prev_i] == "C":
                    continue
            cut0 = idx0 + len(pam_plus) + 3
            cut_genomic = window_start + cut0
            dist = abs(cut_genomic - insertion_pos)
            if dist <= max_cut_distance and "N" not in prot:
                guides.append(Guide(
                    site_name="", chrom=chrom, strand="-",
                    protospacer=prot, pam_seq=pam_seq, guide_len=L,
                    cut_genomic=int(cut_genomic), distance_to_insertion=int(dist),
                    gc=(prot.count("G")+prot.count("C"))/L,
                ))
    return guides

# ===============================================================================
# Write/Resume Helpers
# ===============================================================================

def _write_rows(out_path: str, rows: list, fieldnames: list):
    if not rows:
        return
    df = pd.DataFrame(rows)
    # enforce columns and order
    for c in fieldnames:
        if c not in df.columns:
            df[c] = None
    df = df[fieldnames]
    header = not (os.path.exists(out_path) and os.path.getsize(out_path) > 0)
    df.to_csv(out_path, mode="a", header=header, index=False)

def _write_log_rows(log_path: str, rows: list):
    """Write log rows incrementally with consistent column order"""
    if not rows:
        return
    df = pd.DataFrame(rows)
    log_cols = ['wbgene', 'gene_name', 'transcript', 'site_type', 'chrom', 'pos', 'strand',
                'guide_seq', 'guide_strand', 'cut_genomic', 'distance_to_insertion', 
                'gc', 'on_target_score', 'off_target_hits',
                'status', 'rejection_reason', 'final_rank']
    # enforce columns and order
    for c in log_cols:
        if c not in df.columns:
            df[c] = None
    df = df[log_cols]
    header = not (os.path.exists(log_path) and os.path.getsize(log_path) > 0)
    df.to_csv(log_path, mode="a", header=header, index=False)

def _existing_done_site_keys(out_path: str) -> set[tuple]:
    """
    Return set of (WBGene, input_transcript, site_type, chrom, pos, strand) keys
    already present in the output (i.e., this site has been processed).
    """
    done = set()
    try:
        if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
            df = pd.read_csv(
                out_path,
                usecols=["WBGene","input_transcript","site_type","chrom","pos","strand"],
                dtype=str,
                engine="python"
            )
            df = df.dropna(subset=["WBGene","input_transcript","site_type","chrom","pos","strand"])
            # normalize whitespace / BOMs
            for c in ["WBGene","input_transcript","site_type","chrom","strand"]:
                df[c] = df[c].astype(str).str.strip().str.lstrip("\ufeff")
            # pos to int-safe string
            def _safe_pos(x):
                try:
                    return str(int(float(x)))
                except Exception:
                    return None
            df["pos"] = df["pos"].map(_safe_pos)
            df = df.dropna(subset=["pos"])
            for _, r in df.iterrows():
                done.add((r["WBGene"], r["input_transcript"], r["site_type"], r["chrom"], r["pos"], r["strand"]))
    except Exception:
        pass
    return done

def _dedupe_file(out_path: str):
    """
    Remove duplicate rows from output file, keeping the last occurrence of each unique site.
    Deduplicates based on (WBGene, input_transcript, site_type, chrom, pos, strand) combination.
    """
    try:
        if not os.path.exists(out_path) or os.path.getsize(out_path) == 0:
            return
        
        df = pd.read_csv(out_path, dtype=str, engine="python", on_bad_lines="skip")
        
        # Check for required columns
        required_cols = ["WBGene", "input_transcript", "site_type", "chrom", "pos", "strand"]
        if not all(col in df.columns for col in required_cols):
            print(f"[dedupe][warn] missing required columns for deduplication", file=sys.stderr)
            return
        
        # Normalize whitespace and BOMs
        for c in ["WBGene", "input_transcript", "site_type", "chrom", "strand"]:
            df[c] = df[c].astype(str).str.strip().str.lstrip("\ufeff")
        
        # Normalize pos to int-safe string
        def _safe_pos(x):
            try:
                return str(int(float(x)))
            except Exception:
                return None
        df["pos"] = df["pos"].map(_safe_pos)
        
        # Remove rows with missing key values
        df = df.dropna(subset=required_cols)
        df = df[df["WBGene"] != ""]
        
        original_count = len(df)
        
        # Remove duplicates, keeping the last occurrence (most recent)
        df = df.drop_duplicates(subset=required_cols + ["guide_seq"], keep="last")
        
        deduped_count = len(df)
        
        # Write back to file
        df.to_csv(out_path, index=False)
        
        if original_count != deduped_count:
            print(f"[dedupe] removed {original_count - deduped_count} duplicate row(s), kept {deduped_count} unique site(s)")
        else:
            print(f"[dedupe] no duplicates found, {deduped_count} unique site(s)")
            
    except Exception as e:
        print(f"[dedupe][warn] deduplication failed: {e}", file=sys.stderr)

def _dedupe_log_file(log_path: str):
    """
    Remove duplicate log entries, keeping the last occurrence of each unique guide candidate.
    Deduplicates based on (wbgene, site_type, guide_seq) combination.
    """
    try:
        if not os.path.exists(log_path) or os.path.getsize(log_path) == 0:
            return
        
        df = pd.read_csv(log_path, dtype=str, engine="python", on_bad_lines="skip")
        
        # Check for required columns
        required_cols = ["wbgene", "site_type", "guide_seq"]
        if not all(col in df.columns for col in required_cols):
            print(f"[dedupe_log][warn] missing required columns for log deduplication", file=sys.stderr)
            return
        
        # Normalize whitespace
        for c in required_cols:
            if c in df.columns:
                df[c] = df[c].astype(str).str.strip().str.lstrip("\ufeff")
        
        # Remove rows with missing key values
        df = df.dropna(subset=required_cols)
        df = df[df["wbgene"] != ""]
        
        original_count = len(df)
        
        # Remove duplicates, keeping the last occurrence (most recent)
        df = df.drop_duplicates(subset=required_cols, keep="last")
        
        deduped_count = len(df)
        
        # Write back to file
        df.to_csv(log_path, index=False)
        
        if original_count != deduped_count:
            print(f"[dedupe_log] removed {original_count - deduped_count} duplicate log entry(ies), kept {deduped_count} unique guide(s)")
            
    except Exception as e:
        print(f"[dedupe_log][warn] log deduplication failed: {e}", file=sys.stderr)

# ===============================================================================
# Main
# ===============================================================================
def main():
    t0 = time.perf_counter()
    ap = argparse.ArgumentParser(description="One-stop CSV of tagging sites + sgRNAs (standalone)")
    ap.add_argument('--input', required=False, help="CSV from make_insertion_sites_chainaware.py")
    ap.add_argument('--output', required=True, help="Output CSV path")
    ap.add_argument('--log-output', default=None, help="Optional detailed log CSV tracking all candidates and rejection reasons")
    ap.add_argument('--flush-every', type=int, default=10, help='Append+flush every N input rows (default: 10)')
    ap.add_argument('--no-resume', action='store_true', help='Do not resume; overwrite output')
    ap.add_argument('--max-sgrnas', type=int, default=1, help="Max sgRNAs per site to include")

    # Primary scanning
    ap.add_argument('--genome', required=True, help='Path to local genome FASTA file (e.g., genome.fa)')
    ap.add_argument('--pam', default='NGG', help='PAM sequence (IUPAC)')
    ap.add_argument('--guide-lens', nargs='+', type=int, default=[20], help='Guide lengths (e.g., 19 20 21 22)')
    ap.add_argument('--window-bp', type=int, default=50, help='Half-window around insertion for sequence fetch')
    ap.add_argument('--max-cut-distance', type=int, default=50, help='Max |cut-insertion| allowed (bp)')
    ap.add_argument('--min-gc', type=float, default=0.25, help='Minimum GC fraction')
    ap.add_argument('--max-gc', type=float, default=0.80, help='Maximum GC fraction')
    ap.add_argument('--max-candidates', type=int, default=50, help='Internal cap before selection')

    # Off-targets with FlashFry (PAM-aware, more accurate than BLAST)
    ap.add_argument('--offtarget-mode', choices=['none', 'flashfry'], default='none',
                    help='Off-target detection (flashfry=PAM-aware and stringent)')
    ap.add_argument('--flashfry-db', default=None, help='Path to FlashFry database')
    ap.add_argument('--flashfry-jar', default='FlashFry-assembly-1.15.jar',
                    help='Path to FlashFry JAR file (default: FlashFry-assembly-1.15.jar in current dir)')
    ap.add_argument('--flashfry-max-mismatches', type=int, default=3,
                    help='Max mismatches in protospacer only (3=stringent, 4=moderate, default: 3)')
    ap.add_argument('--allow-offtargets', action='store_true',
                    help='Allow guides with off-targets; return top N ranked by fewest off-targets (default: only return guides with 0 off-targets)')

    args = ap.parse_args()
    if not args.input:
        ap.error('You must provide --input-sites (or --input)')

    # Read input
    df = pd.read_csv(args.input)
    print("[start] Opening genome FASTA file and reading input …")
    
    # Open genome FASTA file
    if not os.path.exists(args.genome):
        print(f"[error] Genome file not found: {args.genome}")
        sys.exit(1)
    
    # pyfaidx will automatically create the .fai index if needed
    genome_fasta = Fasta(args.genome)
    print(f"[start] Genome loaded from {args.genome}")
    print(f"[start] Available chromosomes: {', '.join(list(genome_fasta.keys())[:10])}" + 
          (f" ... and {len(genome_fasta.keys()) - 10} more" if len(genome_fasta.keys()) > 10 else ""))
    print(f"[start] Loaded {len(df)} sites from {args.input}")

    # Prepare output for incremental write/resume
    cols_for_write = [
        'gene','WBGene','input_transcript','uniprot','protein_len','mode','chain_reason',
        'chain_aa_start','chain_aa_end','chain_adjusted_start','chain_adjusted_end','kept_aa_len',
        'warning', 'site_type','chrom','pos','strand',
        'guide_seq','guide_len','cut_genomic','distance_to_insertion','gc','on_target_score','off_target_hits','off_target_locations','transcription_type'
    ]
    out_path = args.output
    # initialize/overwrite per --no-resume
    if args.no_resume and os.path.exists(out_path):
        open(out_path, "w").close()
    # create header if file is empty/missing
    if not os.path.exists(out_path) or os.path.getsize(out_path) == 0:
        pd.DataFrame(columns=cols_for_write).to_csv(out_path, index=False)
    
    # Initialize log file if specified
    if args.log_output:
        if args.no_resume or not os.path.exists(args.log_output) or os.path.getsize(args.log_output) == 0:
            log_cols = ['wbgene', 'gene_name', 'transcript', 'site_type', 'chrom', 'pos', 'strand',
                        'guide_seq', 'guide_strand', 'cut_genomic', 'distance_to_insertion', 
                        'gc', 'on_target_score', 'off_target_hits',
                        'status', 'rejection_reason', 'final_rank']
            pd.DataFrame(columns=log_cols).to_csv(args.log_output, index=False)
        elif not args.no_resume:
            # Deduplicate log file at startup if resuming
            _dedupe_log_file(args.log_output)
    
    # determine already-completed site keys
    done_keys = set() if args.no_resume else _existing_done_site_keys(out_path)
    
    # Deduplicate output file at startup to ensure clean state
    if not args.no_resume:
        _dedupe_file(out_path)
    
    buffer: List[Dict[str, object]] = []
    log_buffer: List[Dict[str, object]] = []  # Detailed log of all candidates
    flush_every = max(1, int(args.flush_every))

    print(f"[start] scanning {len(df)} sites (pam={args.pam}, lenses={args.guide_lens}, window=±{args.window_bp}bp, max_cut_distance={args.max_cut_distance})")

    # Config prints
    print(f"[cfg] GC={args.min_gc:.2f}-{args.max_gc:.2f}, max_candidates={args.max_candidates}")
    print(f"[cfg] offtarget_mode={args.offtarget_mode}, flashfry_db={args.flashfry_db or '(none)'}")

    # FlashFry setup
    flashfry_disabled_reason = None
    if args.offtarget_mode == 'flashfry':
        if not args.flashfry_db:
            flashfry_disabled_reason = 'no --flashfry-db provided'
        elif not os.path.exists(args.flashfry_db):
            flashfry_disabled_reason = f'database not found: {args.flashfry_db}'
        else:
            print(f"[cfg] FlashFry enabled: db={args.flashfry_db}, max_mismatches={args.flashfry_max_mismatches}")
            print(f"[cfg] FlashFry JAR: {args.flashfry_jar}")
    
    if args.offtarget_mode == 'flashfry' and flashfry_disabled_reason:
        print(f"[warn] FlashFry disabled: {flashfry_disabled_reason}", file=sys.stderr)

    # Validate columns
    req_cols = ['gene','WBGene','input_transcript','uniprot','protein_len','site_type','mode','chain_reason','chain_aa_start','chain_aa_end','Chain_adjusted_start','Chain_adjusted_end','chrom','pos','strand']
    missing = [c for c in req_cols if c not in df.columns]
    if missing:
        print(f"[error] Input CSV missing required columns: {missing}")
        sys.exit(1)

    rows_out_total = 0  # incremental writes; count rows written
    skipped_count = 0  # track skipped rows for summary
    skip_summary_printed = False  # flag to print summary once

    for i, row in df.iterrows():
        try:
            # Resume-skip logic
            key = (
                str(row.get('WBGene','') or '').strip(),
                str(row.get('input_transcript','') or '').strip(),
                str(row.get('site_type','') or '').strip(),
                str(row.get('chrom','') or '').strip(),
                str(int(row.get('pos'))) if pd.notna(row.get('pos')) else '',
                str(row.get('strand','') or '').strip()
            )
            if (not args.no_resume) and all(key) and key in done_keys:
                skipped_count += 1
                continue

            # Print skip summary before processing first new row
            if not skip_summary_printed and skipped_count > 0:
                print(f"\n[resume] Skipped {skipped_count} already-processed site(s), continuing with new sites...\n")
                skip_summary_printed = True

            gene = row.get('gene') or row.get('WBGene')
            wbgene = row.get('WBGene')
            transcript = row.get('input_transcript')
            site_type = str(row['site_type'])
            chrom = str(row['chrom'])
            pos = int(row['pos'])
            strand = str(row['strand']).strip()
            name = f"{gene}-{site_type}"
            mode = row.get('mode', '')
            chain_start = row.get('chain_aa_start')
            chain_end = row.get('chain_aa_end')

            print(f"\n=== [{i+1}/{len(df)}] {wbgene} ===")
            print(f"  [gene] name={gene} transcript={transcript}")
            print(f"  [{site_type}] mode={mode} chrom={chrom}:{pos} strand={strand} chain_aa={chain_start}-{chain_end}")
        
        except KeyError as e:
            # Handle missing required columns
            print(f"\n{'='*80}")
            print(f"ERROR: Row {i+1} is missing required column: {str(e)}")
            print(f"{'='*80}")
            print("Available row data:")
            for col in row.index:
                print(f"  {col}: {row[col]}")
            print(f"\nWriting error placeholder and skipping to next row...")
            
            # Create placeholder row with available data and error marker
            error_row = {
                'gene': row.get('gene'),
                'WBGene': row.get('WBGene'),
                'input_transcript': row.get('input_transcript'),
                'uniprot': row.get('uniprot'),
                'protein_len': row.get('protein_len'),
                'mode': row.get('mode'),
                'chain_reason': row.get('chain_reason'),
                'chain_aa_start': row.get('chain_aa_start'),
                'chain_aa_end': row.get('chain_aa_end'),
                'kept_aa_len': (int(row.get('chain_aa_end')) - int(row.get('chain_aa_start')) + 1
                                if pd.notna(row.get('chain_aa_end')) and pd.notna(row.get('chain_aa_start')) else None),
                'warning': f"ERROR_MISSING_COLUMN: {str(e)}",
                'site_type': row.get('site_type'),
                'chrom': row.get('chrom'),
                'pos': row.get('pos'),
                'strand': row.get('strand'),
                'guide_seq': None,
                'guide_len': None,
                'cut_genomic': None,
                'distance_to_insertion': None,
                'gc': None,
                'on_target_score': None,
                'off_target_hits': None,
            }
            buffer.append(error_row)
            
            # Force flush to write error row immediately
            if buffer:
                tmp = pd.DataFrame(buffer)
                tmp.to_csv(out_path, mode='a', header=False, index=False)
                rows_out_total += len(tmp)
                buffer.clear()
            continue
        
        except Exception as e:
            # Handle any other unexpected errors
            print(f"\n{'='*80}")
            print(f"UNEXPECTED ERROR processing Row {i+1}: {str(e)}")
            print(f"{'='*80}")
            traceback.print_exc()
            print("Available row data:")
            for col in row.index:
                print(f"  {col}: {row[col]}")
            print(f"\nWriting error placeholder and skipping to next row...")
            
            # Create placeholder row with available data and error marker
            error_row = {
                'gene': row.get('gene'),
                'WBGene': row.get('WBGene'),
                'input_transcript': row.get('input_transcript'),
                'uniprot': row.get('uniprot'),
                'protein_len': row.get('protein_len'),
                'mode': row.get('mode'),
                'chain_reason': row.get('chain_reason'),
                'chain_aa_start': row.get('chain_aa_start'),
                'chain_aa_end': row.get('chain_aa_end'),
                'kept_aa_len': (int(row.get('chain_aa_end')) - int(row.get('chain_aa_start')) + 1
                                if pd.notna(row.get('chain_aa_end')) and pd.notna(row.get('chain_aa_start')) else None),
                'warning': f"ERROR_PROCESSING: {str(e)}",
                'site_type': row.get('site_type'),
                'chrom': row.get('chrom'),
                'pos': row.get('pos'),
                'strand': row.get('strand'),
                'guide_seq': None,
                'guide_len': None,
                'cut_genomic': None,
                'distance_to_insertion': None,
                'gc': None,
                'on_target_score': None,
                'off_target_hits': None,
            }
            buffer.append(error_row)
            
            # Force flush to write error row immediately
            if buffer:
                tmp = pd.DataFrame(buffer)
                tmp.to_csv(out_path, mode='a', header=False, index=False)
                rows_out_total += len(tmp)
                buffer.clear()
            continue

        # Initialize variables that will be used later (prevents UnboundLocalError)
        tested: List[Guide] = []
        final_guides: List[Guide] = []

        # Sequence window
        window_start = max(1, pos - args.window_bp)
        window_end   = pos + args.window_bp
        try:
            seq_plus_strand = fetch_region_seq(genome_fasta, chrom, window_start, window_end)
            print(f"  [seq] fetched {len(seq_plus_strand)} bp for window {chrom}:{window_start}-{window_end}")
        except Exception as e:
            print(f"  [seq] fetch failed: {e}")
            # placeholder row
            buffer.append({
                'gene': row.get('gene'), 'WBGene': row.get('WBGene'),
                'input_transcript': row.get('input_transcript'), 'uniprot': row.get('uniprot'),
                'protein_len': row.get('protein_len'), 'mode': row.get('mode'), 'chain_reason': row.get('chain_reason'),
                'chain_aa_start': row.get('chain_aa_start'), 'chain_aa_end': row.get('chain_aa_end'),
                'kept_aa_len': (int(row.get('chain_aa_end')) - int(row.get('chain_aa_start')) + 1
                                if pd.notna(row.get('chain_aa_end')) and pd.notna(row.get('chain_aa_start')) else None), 'warning': row.get('warning'), 'site_type': site_type, 'chrom': chrom, 'pos': pos, 'strand': strand,
                'guide_seq': None, 'guide_len': None, 'cut_genomic': None, 'distance_to_insertion': None,
                'gc': None, 'on_target_score': None, 'off_target_hits': None,
            })
            continue

        # Scan
        print(f"  [scan] scanning for PAMs/guides in {args.window_bp}bp window...")
        raw_guides = scan_guides_in_window(seq_plus_strand, chrom, window_start, pos,
                                           args.pam, args.guide_lens, args.max_cut_distance)
        print(f"  [scan] found {len(raw_guides)} raw candidates")

        # Filter + score - separate in vivo and in vitro guides
        kept_in_vivo: List[Guide] = []  # No 4+ homopolymers (good for in vivo)
        kept_in_vitro: List[Guide] = []  # Allows 4-6bp homopolymers (in vitro only)
        rejected_count = {'gc_filter': 0, 'homopolymer_7plus': 0, 'distance': 0}
        
        for g in raw_guides:
            g.site_name = name
            rejection_reason = None
            
            # Check GC content
            if not (args.min_gc <= g.gc <= args.max_gc):
                rejection_reason = f"GC out of range ({g.gc:.2f}, need {args.min_gc:.2f}-{args.max_gc:.2f})"
                rejected_count['gc_filter'] += 1
            # Reject ≥7
            elif HOMOPOLYMER7_RE.search(g.protospacer):
                rejection_reason = "Contains homopolymer run (7+ identical bases)"
                rejected_count['homopolymer_7plus'] += 1

            # T-run ≥4 → in vitro only
            elif T4PLUS_RE.search(g.protospacer):
                g.on_target_score = heuristic_on_target_score(g.protospacer)
                g.transcription_type = 'in_vitro_only'
                kept_in_vitro.append(g)
                rejection_reason = None

            # Everything else → in vivo
            else:
                g.on_target_score = heuristic_on_target_score(g.protospacer)
                g.transcription_type = 'in_vivo'
                kept_in_vivo.append(g)
                rejection_reason = None
            
            # Log every candidate if log output is enabled
            if args.log_output:
                log_buffer.append({
                    'wbgene': wbgene,
                    'gene_name': gene,
                    'transcript': transcript,
                    'site_type': site_type,
                    'chrom': chrom,
                    'pos': pos,
                    'strand': strand,
                    'guide_seq': g.protospacer + g.pam_seq,
                    'guide_strand': g.strand,
                    'cut_genomic': g.cut_genomic,
                    'distance_to_insertion': g.distance_to_insertion,
                    'gc': round(g.gc, 3),
                    'on_target_score': round(g.on_target_score, 3) if g.on_target_score else None,
                    'off_target_hits': None,
                    'status': 'rejected' if rejection_reason else 'passed_initial_filter',
                    'rejection_reason': rejection_reason,
                    'final_rank': None
                })
        
        print(f"  [filter] in_vivo: {len(kept_in_vivo)}, in_vitro: {len(kept_in_vitro)}, rejected: {len(raw_guides) - len(kept_in_vivo) - len(kept_in_vitro)} (GC={rejected_count['gc_filter']}, homopolymer_7+={rejected_count['homopolymer_7plus']})")

        # Rank in vivo guides (allow up to 5bp homopolymers; disallow T-run ≥4)
        kept_in_vivo = [g for g in kept_in_vivo if g.on_target_score is not None]
        kept_in_vivo = sorted(kept_in_vivo, key=lambda x: (x.distance_to_insertion, -(x.on_target_score or -1.0)))[:args.max_candidates]
        
        # Rank in vitro guides (4-6bp homopolymers allowed)
        kept_in_vitro = [g for g in kept_in_vitro if g.on_target_score is not None]
        kept_in_vitro = sorted(kept_in_vitro, key=lambda x: (x.distance_to_insertion, -(x.on_target_score or -1.0)))[:args.max_candidates]
        
        # Determine which guides to keep for selection
        # 1) If NO in vivo guide is found, test ALL in vitro guides (up to args.max_candidates cap applied earlier).
        # 2) If an in vivo guide IS found, test ONLY in vitro guides that are CLOSER than the best in vivo distance.
        kept = []
        if kept_in_vivo:
            kept.extend(kept_in_vivo)  # keep all eligible in vivo (already sorted/truncated)
            best_in_vivo_dist = kept_in_vivo[0].distance_to_insertion
            print(f"  [rank] in_vivo guides: {len(kept_in_vivo)} (closest: {best_in_vivo_dist}bp, score: {kept_in_vivo[0].on_target_score:.1f})")
        else:
            best_in_vivo_dist = float('inf')
            print(f"  [rank] in_vivo guides: 0")
        
        if kept_in_vitro:
            closest_in_vitro_dist = kept_in_vitro[0].distance_to_insertion
            print(f"  [rank] in_vitro guides: {len(kept_in_vitro)} (closest: {closest_in_vitro_dist}bp, score: {kept_in_vitro[0].on_target_score:.1f})")
            if best_in_vivo_dist == float('inf'):
                # No in vivo: include ALL in vitro
                kept.extend(kept_in_vitro)
                print("  [rank] including all in_vitro guides since no in_vivo guide was found")
            else:
                # In vivo present: include ALL in vitro that are closer than best in vivo
                closer_in_vitro = [g for g in kept_in_vitro if g.distance_to_insertion < best_in_vivo_dist]
                kept.extend(closer_in_vitro)
                if closer_in_vitro:
                    print(f"  [rank] including {len(closer_in_vitro)} in_vitro guide(s) closer than best in_vivo ({best_in_vivo_dist}bp)")
                else:
                    print(f"  [rank] no in_vitro guides closer than best in_vivo ({best_in_vivo_dist}bp)")
        else:
            print(f"  [rank] in_vitro guides: 0")

        # Off-targets (optional) - test guides until finding ones with 0 off-targets
        
        # Separate in_vivo and in_vitro guides for proper selection
        in_vivo_from_kept = [g for g in kept if g.transcription_type == 'in_vivo']
        in_vitro_from_kept = [g for g in kept if g.transcription_type == 'in_vitro_only']
        
        if kept and args.offtarget_mode == 'flashfry' and not flashfry_disabled_reason:
            # If NO in_vivo guide, test ALL in_vitro guides (override max_candidates cap)
            # If there IS an in_vivo guide, test kept (in_vivo + closer in_vitro), capped by max_candidates
            if len(in_vivo_from_kept) == 0 and len(in_vitro_from_kept) > 0:
                guides_to_test = in_vitro_from_kept  # all in vitro
                max_to_test = len(guides_to_test)
            else:
                guides_to_test = kept
                max_to_test = min(len(guides_to_test), args.max_candidates)
            
            print(f"  [flashfry] testing all {max_to_test} candidate guides for off-targets (max_mismatches={args.flashfry_max_mismatches})")
            print(f"  [flashfry] using BATCH MODE for faster processing...")
            
            zero_hits: List[Guide] = []
            first_error = None
            
            # Use batch processing - much faster than sequential!
            try:
                batch_results = flashfry_offtargets_batch(
                    guides_to_test[:max_to_test],
                    args.flashfry_db,
                    max_mismatches=args.flashfry_max_mismatches,
                    flashfry_jar=args.flashfry_jar
                )
                
                # Process results
                for idx_b, g in enumerate(guides_to_test[:max_to_test], 1):
                    q = guide_query_with_pam(g)
                    
                    if q in batch_results:
                        hits, off_target_details = batch_results[q]
                        g.off_target_hits = hits
                        g.off_target_details = off_target_details if off_target_details else None
                        tested.append(g)
                        
                        if hits == 0:
                            zero_hits.append(g)
                            print(f"  [flashfry] guide {idx_b}/{max_to_test}: {g.protospacer}{g.pam_seq} has 0 off-targets")
                        elif hits is not None:
                            print(f"  [flashfry] guide {idx_b}/{max_to_test}: {g.protospacer}{g.pam_seq} has {hits} off-target(s)")
                    else:
                        # Failed to get results for this guide
                        g.off_target_hits = None
                        g.off_target_details = None
                        tested.append(g)
                        print(f"  [flashfry] guide {idx_b}/{max_to_test}: {g.protospacer}{g.pam_seq} failed to check")
                        
            except Exception as e:
                first_error = str(e)
                print(f"  [flashfry] batch error: {first_error}")
                print(f"  [flashfry] falling back to sequential processing...")
                
                # Fallback to sequential processing if batch fails
                for idx_b, g in enumerate(guides_to_test[:max_to_test], 1):
                    q = guide_query_with_pam(g)
                    hits = None
                    off_target_details = []
                    try:
                        hits, off_target_details = flashfry_offtargets(
                            q, args.flashfry_db,
                            max_mismatches=args.flashfry_max_mismatches,
                            flashfry_jar=args.flashfry_jar
                        )
                    except Exception as e:
                        if not first_error:
                            first_error = str(e)
                            print(f"  [flashfry] error: {first_error}; continuing to test remaining guides")
                        hits = None
                        
                    g.off_target_hits = hits
                    g.off_target_details = off_target_details if off_target_details else None
                    tested.append(g)
                    
                    if hits == 0:
                        zero_hits.append(g)
                        print(f"  [flashfry] guide {idx_b}/{max_to_test}: {g.protospacer}{g.pam_seq} has 0 off-targets")
                    elif hits is not None:
                        print(f"  [flashfry] guide {idx_b}/{max_to_test}: {g.protospacer}{g.pam_seq} has {hits} off-target(s)")
                    else:
                        print(f"  [flashfry] guide {idx_b}/{max_to_test}: {g.protospacer}{g.pam_seq} failed to check")
                
                # Break only on fatal error
                if first_error and hits is None:
                    print(f"  [flashfry] fatal error, stopping FlashFry checks")
                    break
            
            # After FlashFry, select guides prioritizing zero off-targets
            tested_in_vivo = [g for g in tested if g.transcription_type == 'in_vivo']
            tested_in_vitro = [g for g in tested if g.transcription_type == 'in_vitro_only']
            
            # Rank all tested guides by off-targets
            def _rk(g: Guide):
                off = g.off_target_hits
                # Handle string values (status messages) by treating as very high number
                if isinstance(off, str):
                    off = 10**9
                elif off is None:
                    off = 10**9
                else:
                    off = int(off)
                return (off, g.distance_to_insertion, -(g.on_target_score or -1.0))
            
            tested_in_vivo_sorted = sorted(tested_in_vivo, key=_rk)
            tested_in_vitro_sorted = sorted(tested_in_vitro, key=_rk)
            
            # Check if we have guides with 0 off-targets (numeric zero, not string)
            zero_in_vivo = [g for g in tested_in_vivo_sorted if isinstance(g.off_target_hits, int) and g.off_target_hits == 0]
            zero_in_vitro = [g for g in tested_in_vitro_sorted if isinstance(g.off_target_hits, int) and g.off_target_hits == 0]
            has_zero_hits = bool(zero_in_vivo or zero_in_vitro)
            
            # Select guides based on --allow-offtargets flag
            if not args.allow_offtargets and has_zero_hits:
                # Default behavior: only return guides with 0 off-targets
                if zero_in_vivo:
                    # Add up to max_sgrnas in_vivo guides with 0 off-targets
                    final_guides.extend(zero_in_vivo[:args.max_sgrnas])
                    # ALWAYS add closer in_vitro guides (even if exceeds max_sgrnas)
                    if zero_in_vitro:
                        best_in_vivo_dist = zero_in_vivo[0].distance_to_insertion
                        closer_in_vitro = [g for g in zero_in_vitro if g.distance_to_insertion < best_in_vivo_dist]
                        if closer_in_vitro:
                            final_guides.extend(closer_in_vitro[:args.max_sgrnas])
                elif zero_in_vitro:
                    # No in_vivo found, use in_vitro
                    final_guides.extend(zero_in_vitro[:args.max_sgrnas])
                
                print(f"  [flashfry] selected {len(final_guides)} guide(s) with 0 off-targets (tested {len(tested)}/{max_to_test} candidates)")
            else:
                # --allow-offtargets: return top N ranked by fewest off-targets
                if not args.allow_offtargets:
                    print(f"  [flashfry] warning: no guides with 0 off-targets found (tested {len(tested)}/{max_to_test}), selecting guide(s) with fewest off-targets")
                else:
                    print(f"  [flashfry] selecting top {args.max_sgrnas} guide(s) ranked by fewest off-targets (tested {len(tested)}/{max_to_test} candidates)")
                
                if tested_in_vivo_sorted:
                    # Add up to max_sgrnas in_vivo guides
                    final_guides.extend(tested_in_vivo_sorted[:args.max_sgrnas])
                    # ALWAYS add closer in_vitro guides (even if exceeds max_sgrnas)
                    if tested_in_vitro_sorted:
                        best_in_vivo_dist = tested_in_vivo_sorted[0].distance_to_insertion
                        closer_in_vitro = [g for g in tested_in_vitro_sorted if g.distance_to_insertion < best_in_vivo_dist]
                        if closer_in_vitro:
                            final_guides.extend(closer_in_vitro[:args.max_sgrnas])
                elif tested_in_vitro_sorted:
                    # No in_vivo found, use in_vitro
                    final_guides.extend(tested_in_vitro_sorted[:args.max_sgrnas])
        
        # Set appropriate off_target_hits value for guides not tested with FlashFry
        elif kept:
            # Determine status message based on configuration
            if args.offtarget_mode == 'none':
                off_target_status = "not checked"
            elif args.offtarget_mode == 'flashfry' and flashfry_disabled_reason:
                off_target_status = "not checked - no database found"
            else:
                off_target_status = "not checked"
            
            # Set status for all guides
            for g in kept:
                if g.off_target_hits is None:
                    g.off_target_hits = off_target_status
            
            # No FlashFry: select best in_vivo, plus closer in_vitro if applicable
            if in_vivo_from_kept:
                # Add up to max_sgrnas in_vivo guides
                final_guides.extend(in_vivo_from_kept[:args.max_sgrnas])
                # ALWAYS add closer in_vitro guides (even if exceeds max_sgrnas)
                if in_vitro_from_kept:
                    best_in_vivo_dist = in_vivo_from_kept[0].distance_to_insertion
                    closer_in_vitro = [g for g in in_vitro_from_kept if g.distance_to_insertion < best_in_vivo_dist]
                    if closer_in_vitro:
                        final_guides.extend(closer_in_vitro[:args.max_sgrnas])
            elif in_vitro_from_kept:
                # No in_vivo guides found, use in_vitro
                final_guides.extend(in_vitro_from_kept[:args.max_sgrnas])

        # Emit rows (or placeholder if none)
        site_rows: List[Dict[str, object]] = []
        if not final_guides:
            site_rows.append({
                'gene': row.get('gene'), 'WBGene': row.get('WBGene'),
                'input_transcript': row.get('input_transcript'),
                'uniprot': row.get('uniprot'),
                'protein_len': row.get('protein_len'),
                'mode': row.get('mode'),
                'chain_reason': row.get('chain_reason'),
                'chain_aa_start': row.get('chain_aa_start'),
                'chain_aa_end': row.get('chain_aa_end'),
                'kept_aa_len': (int(row.get('chain_aa_end')) - int(row.get('chain_aa_start')) + 1
                                if pd.notna(row.get('chain_aa_end')) and pd.notna(row.get('chain_aa_start')) else None), 'warning': row.get('warning'), 'site_type': site_type, 'chrom': chrom, 'pos': pos, 'strand': strand,
                'guide_seq': None, 'guide_len': None, 'cut_genomic': None, 'distance_to_insertion': None,
                'gc': None, 'on_target_score': None, 'off_target_hits': None,
                'off_target_locations': None,
                'transcription_type': None,
            })
        else:
            for g in final_guides:
                # Format off-target details as a semicolon-separated string
                off_target_locations = None
                if g.off_target_details:
                    # Format: sequence(Xmm) where X is number of mismatches
                    locations = []
                    for ot in g.off_target_details:
                        mm = ot.get('mismatches', '?')
                        locations.append(f"{ot['sequence']}({mm}mm)")
                    off_target_locations = "; ".join(locations)
                
                site_rows.append({
                    'gene': row.get('gene'), 
                    'WBGene': row.get('WBGene'),
                    'input_transcript': row.get('input_transcript'), 
                    'uniprot': row.get('uniprot'),
                    'protein_len': row.get('protein_len'), 
                    'mode': row.get('mode'),
                    'chain_reason': row.get('chain_reason'),
                    'chain_aa_start': row.get('chain_aa_start'), 
                    'chain_aa_end': row.get('chain_aa_end'),
                    'chain_adjusted_start': row.get('Chain_adjusted_start'), 
                    'chain_adjusted_end': row.get('Chain_adjusted_end'),
                    'kept_aa_len': (int(row.get('chain_aa_end')) - int(row.get('chain_aa_start')) + 1
                                    if pd.notna(row.get('chain_aa_end')) and pd.notna(row.get('chain_aa_start')) else None), 'warning': row.get('warning'), 'site_type': site_type, 'chrom': chrom, 'pos': pos, 'strand': strand,
                    'guide_seq': g.protospacer + g.pam_seq, 'guide_len': g.guide_len,
                    'cut_genomic': g.cut_genomic, 'distance_to_insertion': g.distance_to_insertion,
                    'gc': round(g.gc, 3),
                    'on_target_score': round(g.on_target_score or 0.0, 3),
                    'off_target_hits': g.off_target_hits,
                    'off_target_locations': off_target_locations,
                    'transcription_type': g.transcription_type,
                })
        
        # Update log with off-target counts for all tested guides
        if args.log_output and tested:
            for g in tested:
                guide_seq = g.protospacer + g.pam_seq
                # Find this guide in log and update its off-target count
                for log_entry in reversed(log_buffer):
                    if (log_entry['guide_seq'] == guide_seq and 
                        log_entry['wbgene'] == wbgene and 
                        log_entry['site_type'] == site_type):
                        log_entry['off_target_hits'] = g.off_target_hits
                        break
        
        # Update log with final ranks for accepted guides
        if args.log_output and final_guides:
            for rank, g in enumerate(final_guides, 1):
                guide_seq = g.protospacer + g.pam_seq
                # Find this guide in log and update its rank
                for log_entry in reversed(log_buffer):
                    if (log_entry['guide_seq'] == guide_seq and 
                        log_entry['wbgene'] == wbgene and 
                        log_entry['site_type'] == site_type):
                        log_entry['status'] = 'accepted'
                        log_entry['final_rank'] = rank
                        # off_target_hits already set above
                        break

        # Print summary for this site
        if final_guides:
            print(f"  [done] selected {len(final_guides)} guide(s):")
            for idx, g in enumerate(final_guides, 1):
                off_str = f"off={g.off_target_hits}" if g.off_target_hits is not None else "off=N/A"
                tx_type = g.transcription_type or "unknown"
                print(f"    #{idx}: {g.protospacer}{g.pam_seq} (dist={g.distance_to_insertion}bp, score={g.on_target_score:.1f}, {off_str}, {tx_type})")
                # Print off-target details if any exist
                if g.off_target_details and len(g.off_target_details) > 0:
                    print(f"         Off-target sequences found:")
                    for ot_idx, ot in enumerate(g.off_target_details[:10], 1):  # Show up to 10
                        mm = ot.get('mismatches', '?')
                        print(f"           {ot_idx}. {ot['sequence']} ({mm} mismatches)")
                    if len(g.off_target_details) > 10:
                        print(f"           ... and {len(g.off_target_details) - 10} more off-targets")
        else:
            print(f"  [done] NO guides found - wrote placeholder row")
        
        buffer.extend(site_rows)
        rows_out_total += len(site_rows)
        # flush checkpoint
        if ((i+1) % flush_every) == 0:
            print(f"[flush] {i+1}/{len(df)} input rows; flushing {len(buffer)} output rows...")
            _write_rows(out_path, buffer, cols_for_write)
            buffer.clear()
            _dedupe_file(out_path)  # Deduplicate after each flush
            if args.log_output and log_buffer:
                print(f"[flush] flushing {len(log_buffer)} log entries...")
                _write_log_rows(args.log_output, log_buffer)
                log_buffer.clear()
                _dedupe_log_file(args.log_output)  # Deduplicate log after each flush

    # Final flush
    if buffer:
        print(f"[flush] final: flushing {len(buffer)} output rows...")
        _write_rows(out_path, buffer, cols_for_write)
        buffer.clear()
        _dedupe_file(out_path)  # Deduplicate after final flush
    if args.log_output and log_buffer:
        print(f"[flush] final: flushing {len(log_buffer)} log entries...")
        _write_log_rows(args.log_output, log_buffer)
        log_buffer.clear()
        _dedupe_log_file(args.log_output)  # Deduplicate log after final flush
    
    print(f"[done] Incremental write complete. Wrote ~{rows_out_total} rows (see output file)")
    
    # If all rows were skipped and summary wasn't printed yet
    if skipped_count > 0 and not skip_summary_printed:
        print(f"[resume] All {skipped_count} site(s) were already processed - nothing new to do")
    
    # Print log summary if log was written
    if args.log_output:
        print(f"[log] Detailed log saved to: {args.log_output}")
        # Read back to get summary stats
        if os.path.exists(args.log_output) and os.path.getsize(args.log_output) > 0:
            try:
                log_df = pd.read_csv(args.log_output)
                total_candidates = len(log_df)
                accepted = sum(log_df['status'] == 'accepted')
                rejected = sum(log_df['status'] == 'rejected')
                print(f"[log] Summary: {total_candidates} total candidates, {accepted} accepted, {rejected} rejected")
            except Exception:
                pass  # If we can't read it back, ignore

    # timing
    dt = time.perf_counter() - t0
    sites_processed = len(df) - skipped_count
    print(f"\n[timing] Elapsed {dt:.2f}s | {sites_processed}/{len(df)} sites processed")
    
    # Close genome FASTA file
    genome_fasta.close()
    print("[cleanup] Closed genome FASTA file")

if __name__ == "__main__":
    main()