#!/usr/bin/env python3
"""
Genomic Coordinate Calculator for N/C-Terminal CRISPR Knock-In Insertion Sites

Computes precise genomic coordinates for N- and C-terminal tag insertion sites,
accounting for post-translational processing and mature protein boundaries.

Input:
  - CSV file from script 1 containing WBGene IDs, canonical transcripts, and protein 
    features (minimally requires: WBGene, canonical_transcript)

Output:
  - CSV file with two rows per gene (Nterm and Cterm), containing:
      • Gene identifiers and protein metadata
      • Site type (Nterm/Cterm) and processing mode (chain/full)
      • Chain boundaries (amino acid coordinates of mature protein)
      • Genomic coordinates (chromosome, position, strand)
      • Processing rationale (why this insertion site was chosen)
      • Features summary and validation warnings

Processing Strategy:
  1. Chain mode (when PTM detected): Uses Chain_regions from UniProt to define
     the mature protein after signal peptide cleavage, propeptide removal, or methionine
     excision. Insertion sites are computed relative to the mature chain boundaries.
     If a terminal lipidation site(s) is present, the insertion is made 5aa proximal to
     the terminus to avoid disrupting the lipidation site.
     
  2. Full mode (no PTMs): Uses complete translation boundaries from start to
     stop codon when no post-translational modifications are detected.
     
  3. C-terminal override: If UniProt chain is shorter than Ensembl translation
     but no C-terminal processing features exist, extends C-terminus to Ensembl length
     to capture full coding sequence.

Coordinate Mapping:
  - Uses Ensembl /map/cds endpoint for strand-aware codon-to-genome mapping
  - N-terminal: Position immediately before first codon of mature protein
  - C-terminal: Position immediately before stop codon of mature protein
  - Handles splice junctions, introns, and strand orientation automatically
  - Provides fallback to cDNA mapping if CDS mapping fails

Features:
  - If interrupted, automatically resumes from where it left off by skipping gene-transcript
     pairs already present in the output file
  - Adjusts chain boundaries to exclude C-terminal lipidation sites
  - Flags uncertain coordinates from UniProt annotations
  - Additionally predicts CAAX motifs, PTS signals, ER retention motifs based on sequence
  - Detailed progress logging with per-gene summaries

Required arguments:
  --input PATH      Input CSV with protein features from script 1
  --output PATH     Output CSV path with added genomic coordinates
  
Optional arguments:
  --flush-every N     Write results every N genes (default: 10)
  --no-resume         Overwrite existing output instead of resuming
  --max-rows N        Process only first N rows (for testing)

Example:
  python 2_getinsertionsites.py \
    --input 1.protein_features.csv \
    --output 2.insertion_sites.csv \
    --flush-every 100
"""

from __future__ import annotations
import argparse, re, os, sys, time, math, requests
import pandas as pd
from typing import Dict, List, Optional, Tuple
from requests.adapters import HTTPAdapter, Retry


ENSEMBL_REST = "https://rest.ensembl.org"
SPECIES = "caenorhabditis_elegans"  # WBcel235

# --- helpers for parsing AA coordinate ranges ---
_feature_range_regexes = [
    re.compile(r"\b(\d+)\s*[_\-—]\s*(\d+)\b"),        # e.g. 1_20 or 1-20
    re.compile(r"\[\s*(\d+)\s*[-—]\s*(\d+)\s*\]"),   # e.g. [1-20]
]

def build_http_session() -> requests.Session:
    session = requests.Session()
    retries = Retry(
        total=6,
        backoff_factor=0.5,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET", "POST"),
    )
    session.mount("https://", HTTPAdapter(max_retries=retries))
    session.headers.update({
        "User-Agent": "celegans-chain-aware-sites/2.1 (contact: you@example.com)",
        "Accept": "application/json",
    })
    return session

def get_json(session: requests.Session, path: str, params: Dict | None = None) -> Optional[dict]:
    url = f"{ENSEMBL_REST}{path}"
    resp = session.get(url, params=params or {}, timeout=45)
    if resp.status_code == 200:
        return resp.json()
    try:
        msg = resp.json()
        print(f"[http] {path} -> {resp.status_code}: {msg}", file=sys.stderr)
    except Exception:
        print(f"[http] {path} -> {resp.status_code}", file=sys.stderr)
    return None

def lookup_id(session: requests.Session, stable_id: str, expand: int = 1) -> Optional[dict]:
    return get_json(session, f"/lookup/id/{stable_id}", {"expand": expand})

def map_cdna_to_genomic(session: requests.Session, transcript_id: str, cdna_pos: int) -> Optional[Tuple[str,int,int]]:
    if cdna_pos is None or cdna_pos < 1:
        return None
    j = get_json(session, f"/map/cdna/{transcript_id}/{cdna_pos}..{cdna_pos}")
    if not j or not isinstance(j.get("mappings"), list) or not j["mappings"]:
        return None
    # Prefer chromosome mappings; accept nested or flat
    best = None
    for m in j.get("mappings", []):
        mm = m.get("mapped") if isinstance(m.get("mapped"), dict) else (m if isinstance(m, dict) else None)
        if not isinstance(mm, dict):
            continue
        if mm.get("coord_system") == "chromosome":
            return (mm.get("seq_region_name"), mm.get("start"), mm.get("strand"))
        if best is None:
            best = (mm.get("seq_region_name"), mm.get("start"), mm.get("strand"))
    if best:
        return best
    m0 = j.get("mappings", [{}])[0]
    mm = m0.get("mapped") if isinstance(m0.get("mapped"), dict) else (m0 if isinstance(m0, dict) else None)
    if isinstance(mm, dict) and all(k in mm for k in ("seq_region_name","start","strand")):
        return (mm.get("seq_region_name"), mm.get("start"), mm.get("strand"))
    return None

# ===============================================================================
# CDS Mapping Helpers
# ===============================================================================

def map_cds_to_genome(session: requests.Session, transcript_id: str,
                      cds_start_nt: int, cds_end_nt: int) -> Optional[List[dict]]:
    j = get_json(session, f"/map/cds/{transcript_id}/{cds_start_nt}..{cds_end_nt}")
    if not j or "mappings" not in j or not j["mappings"]:
        return None
    return j["mappings"]

def insertion_from_codon_blocks(mapping_blocks: List[dict], which: str) -> Tuple[Optional[str], Optional[int], Optional[str]]:
    """
    Compute insertion coordinate from mapped codon blocks.
      N: just before start codon; C: just before stop (immediately after last coding codon).
    """
    if not mapping_blocks:
        return None, None, None
    mb0 = mapping_blocks[0]
    chrom = mb0.get("seq_region_name")
    strand = int(mb0.get("strand", 0) or 0)
    codon_start_gen = min(int(m["start"]) for m in mapping_blocks if "start" in m)
    codon_end_gen   = max(int(m["end"])   for m in mapping_blocks if "end" in m)
    if which == 'N':
        pos = (codon_start_gen - 1) if strand == 1 else (codon_end_gen + 1)
    else:
        pos = (codon_end_gen + 1) if strand == 1 else (codon_start_gen)
    if pos is not None and pos < 1:
        pos = 1
    strand_symbol = "+" if strand == 1 else "-" if strand == -1 else None
    return chrom, pos, strand_symbol

# ===============================================================================
# Parsing helpers
# ===============================================================================

_rng_re = re.compile(r"^\s*(\d+)\s*[_-]\s*(\d+)\b")

def get_gene_name(row: 'pd.Series') -> str | None:
    # Attempt to read a symbol from common columns if present
    for key in ["gene", "Gene", "locus", "public_name", "Public_name", "Gene name", "gene_name", "GeneName", "WB_symbol", "Symbol"]:
        try:
            val = row.get(key)
        except Exception:
            val = None
        if val is not None:
            s = str(val).strip()
            if s and s.lower() != "nan":
                return s
    return None

def fetch_gene_name(session: requests.Session, transcript_obj: dict) -> str | None:
    """
    Resolve a short gene symbol (e.g., 'fbl-1') from a transcript lookup object.
      1) Use transcript_obj['Parent'] to lookup the gene; prefer 'display_name' or 'external_name'.
      2) Fallback to transcript_obj.get('display_name')/'external_name' if they look like a gene symbol.
    """
    try:
        gene_id = transcript_obj.get("Parent")
        if gene_id:
            gj = lookup_id(session, gene_id, expand=0)
            if gj:
                for key in ("display_name", "external_name", "external_id"):
                    name = gj.get(key)
                    if name and isinstance(name, str) and name.strip():
                        return name.strip()
        for key in ("display_name", "external_name"):
            name = transcript_obj.get(key)
            if name and isinstance(name, str) and name.strip():
                return name.strip()
    except Exception:
        pass
    return None

def first_range(s: str | float | None) -> Optional[Tuple[int,int]]:
    if s is None or (isinstance(s, float) and math.isnan(s)):
        return None
    s = str(s)
    s = re.sub(r'(?<=\d)\s*\(\s*UNSURE\s*\)', '', s, flags=re.I)
    m = _rng_re.match(s)
    if not m:
        return None
    a, b = int(m.group(1)), int(m.group(2))
    if a < 1 or b < a:
        return None
    return (a, b)

def any_nonempty(*vals) -> bool:
    for v in vals:
        if v is None:
            continue
        if isinstance(v, float) and math.isnan(v):
            continue
        if str(v).strip() != "":
            return True
    return False

def build_chain_reason_pair(row: 'pd.Series', use_chain: bool, aa_start: int | None, aa_end: int | None, aa_len: int | None) -> tuple[str, str]:
    """
    Return (n_reason, c_reason) *based on where the chain actually trims*:
      - If aa_start > 1 -> N-terminus trimmed -> put reasons on N row
      - If aa_end < aa_len -> C-terminus trimmed -> put reasons on C row
    Feature text comes from columns when present; otherwise a generic message.
    """
    if not use_chain or not aa_len or not aa_start or not aa_end:
        return "", ""

    n_reason_parts = []
    c_reason_parts = []

    # N-side trimmed?
    if aa_start > 1:
        if any_nonempty(row.get("SP_regions")):
            n_reason_parts.append("N-terminal signal peptide cleaved off")
        if any_nonempty(row.get("Transit_peptides")):
            n_reason_parts.append("N-terminal transit peptide cleaved off")
        if any_nonempty(row.get("Propeptides")):
            n_reason_parts.append("N-terminal propeptide removed")
        if any_nonempty(row.get("iMET_events")):
            n_reason_parts.append("Initiator methionine removed")
        if not n_reason_parts:
            n_reason_parts.append("N-terminus cleaved (chain start > 1)")

    # C-side trimmed?
    if aa_end < aa_len:
        if any_nonempty(row.get("Propeptides")):
            c_reason_parts.append("C-terminal propeptide removed")
        if any_nonempty(row.get("Transit_peptides")):
            c_reason_parts.append("C-terminal transit peptide cleaved off")
        # We do NOT use CAAX here (not cleaved; prenylation motif)
        if not c_reason_parts:
            c_reason_parts.append("C-terminus cleaved (chain end < protein length)")

    return "; ".join(n_reason_parts), "; ".join(c_reason_parts)

def _detect_unknown_terminals(s: object) -> list[str]:
    """Detect strings like '1_?(UNKNOWN)' or '?(UNKNOWN)_360' and return detail messages."""
    if s is None:
        return []
    try:
        s_str = str(s).strip()
    except Exception:
        return []
    if not s_str or s_str.lower() == "nan":
        return []

    details = []
    # unknown start like '?(UNKNOWN)_360'
    if re.search(r'^\s*\?\s*(?:\(\s*UNKNOWN\s*\))?\s*[_\-—]\s*\d+', s_str, flags=re.I):
        details.append(f"unknown N-terminal boundary: '{s_str}'")
    # unknown end like '1_?(UNKNOWN)'
    if re.search(r'\d+\s*[_\-—]\s*\?\s*(?:\(\s*UNKNOWN\s*\))?\b', s_str, flags=re.I):
        details.append(f"unknown C-terminal boundary: '{s_str}'")
    # generic unknown
    if not details and re.search(r'UNKNOWN|\?\(', s_str, flags=re.I):
        details.append(f"unknown coordinate in: '{s_str}'")

    return details


def _extract_ranges_from_text(s) -> List[Tuple[int, int]]:
    """Return all (start,end) AA ranges found in a free-text field."""
    out: List[Tuple[int, int]] = []
    if s is None:
        return out
    
    try:
        s = str(s)
    except Exception:
        return out
    if not s or s.strip().lower() == "nan":
        return out
    for rx in _feature_range_regexes:
        for a, b in rx.findall(s):
            try:
                a_i, b_i = int(a), int(b)
            except Exception:
                continue
            if a_i >= 1 and b_i >= a_i:
                out.append((a_i, b_i))
    return out

def _dedupe_warning_field(rec: dict) -> None:
    """Remove duplicate tokens from rec['warning'] (pipe-separated), preserving order."""
    w = rec.get("warning")
    if not w:
        return
    toks = [t.strip() for t in str(w).split("|") if t and str(t).strip()]
    seen = set()
    out = []
    for t in toks:
        if t not in seen:
            seen.add(t)
            out.append(t)
    rec["warning"] = " | ".join(out)

def _merge_warning(prev, new):
    """
    Merge warning strings like "A | B" with deduplication (order-preserving).
    Accepts None/empty; returns a single string (possibly empty).
    """
    tokens = []
    for blob in (prev, new):
        if not blob:
            continue
        tokens.extend(t.strip() for t in str(blob).split("|") if t and str(t).strip())
    seen = set()
    out = []
    for t in tokens:
        if t not in seen:
            seen.add(t)
            out.append(t)
    return " | ".join(out)


# ===============================================================================
# Core Logic
# ===============================================================================

def compute_sites_for_row(session: requests.Session, row: pd.Series) -> List[Dict]:
    """Return two dicts (Nterm, Cterm) with genomic coordinates; strand-aware via /map/cds.
    Chain mode uses Chain_regions; reasons printed on the trimmed side only.
    Includes C-terminus override to Ensembl aa_len when no C-side feature is present but UniProt is shorter.
    """
    gene_name = get_gene_name(row)
    wbgene = row.get("WBGene")
    t_in = row.get("canonical_transcript") or row.get("transcript") or row.get("Transcript") or None
    uniprot = row.get("uniprot")

    use_chain = any_nonempty(row.get("iMET_events"), row.get("SP_regions"), row.get("Propeptides"), row.get("Transit_peptides"))
    chain_rng = first_range(row.get("Chain_regions")) if use_chain else None

    tj = lookup_id(session, str(t_in).strip(), expand=1) if t_in else None
    if not tj or tj.get("object_type") != "Transcript":
        rows = []
        for st in ("Nterm","Cterm"):
            n_reason, c_reason = build_chain_reason_pair(row, True if chain_rng else False, chain_rng[0] if chain_rng else None, chain_rng[1] if chain_rng else None, int(row.get("protein_len")) if row.get("protein_len")==row.get("protein_len") else None)
            rows.append({
                "gene": gene_name,
                "WBGene": wbgene,
                "input_transcript": t_in,
                "uniprot": uniprot,
                "protein_len": row.get("protein_len"),
                "site_type": st,
                "mode": "chain" if chain_rng else "full",
                "chain_reason": n_reason if st=="Nterm" else (c_reason if st=="Cterm" else ""),
                "chain_aa_start": chain_rng[0] if chain_rng else None,
                "chain_aa_end": chain_rng[1] if chain_rng else None,
                "chrom": None, "pos": None, "strand": None,
                "boundary_cdna": None, "before_cdna": None,
                "boundary_chrom": None, "boundary_pos": None, "boundary_strand": None,
                "source": "lookup_failed",
                "warning": "Could not resolve transcript"
            })
        return rows

    if not gene_name:
        gene_name = fetch_gene_name(session, tj) or gene_name

    tx_id = tj["id"]
    tr = tj.get("Translation")
    if not tr:
        rows = []
        for st in ("Nterm","Cterm"):
            n_reason, c_reason = build_chain_reason_pair(row, True if chain_rng else False, chain_rng[0] if chain_rng else None, chain_rng[1] if chain_rng else None, int(row.get("protein_len")) if row.get("protein_len")==row.get("protein_len") else None)
            rows.append({
                "gene": gene_name,
                "WBGene": wbgene,
                "input_transcript": t_in,
                "uniprot": uniprot,
                "protein_len": row.get("protein_len"),
                "site_type": st,
                "mode": "chain" if chain_rng else "full",
                "chain_reason": n_reason if st=="Nterm" else (c_reason if st=="Cterm" else ""),
                "chain_aa_start": chain_rng[0] if chain_rng else None,
                "chain_aa_end": chain_rng[1] if chain_rng else None,
                "chrom": None, "pos": None, "strand": None,
                "boundary_cdna": None, "before_cdna": None,
                "boundary_chrom": None, "boundary_pos": None, "boundary_strand": None,
                "source": "no_translation",
                "warning": "Transcript lacks a translation"
            })
        return rows

    # AA / CDS lengths
    try:
        aa_len = tr.get("length")
        if aa_len is None:
            aa_len = (int(tr["end"]) - int(tr["start"]) + 1) // 3
        aa_len = int(aa_len)
        cds_len_nt = aa_len * 3
    except Exception:
        aa_len = None
        cds_len_nt = None

    # Choose chain vs full AA span - ALWAYS use original chain first
    if chain_rng and aa_len:
        aa_start, aa_end = chain_rng
        aa_start = max(1, min(int(aa_start), aa_len))
        aa_end   = max(1, min(int(aa_end),   aa_len))
        if aa_end < aa_start:
            aa_start, aa_end = aa_end, aa_start
        mode = "chain"
    else:
        aa_start, aa_end = 1, aa_len if aa_len else 1
        mode = "full"

    # Store the original chain values (these will go in chain_aa_start/end columns)
    original_chain_start = aa_start
    original_chain_end = aa_end

    # --- Calculate adjusted coordinates based on lipidation ---
    try:
        _chain_start, _chain_end = (aa_start, aa_end)
    except Exception:
        _chain_start, _chain_end = (None, None)

    adjusted_start = None
    adjusted_end = None

    _lip = row.get("Lipidation_sites") if hasattr(row, "get") else None
    if _chain_start is not None and _chain_end is not None and _lip is not None and str(_lip).strip() and str(_lip).strip().lower() != "nan":
        _lip_ranges = _extract_ranges_from_text(_lip)
        
        TERMINAL_THRESHOLD = 5
        
        start_overlaps = [(a,b) for (a,b) in _lip_ranges if a <= _chain_start + TERMINAL_THRESHOLD]
        end_overlaps   = [(a,b) for (a,b) in _lip_ranges if b >= _chain_end - TERMINAL_THRESHOLD]
        
        adjustment_made = False
        temp_start = _chain_start
        temp_end = _chain_end
        
        if start_overlaps:
            max_b = max(b for (_,b) in start_overlaps)
            temp_start = max_b + 5 # new insertion site 5aa downstream of N-term
            adjustment_made = True
        
        if end_overlaps:
            min_a = min(a for (a,_) in end_overlaps)
            temp_end = min_a - 5 # new insertion site 5aa upstream of C-term
            adjustment_made = True
        
        # Only set adjusted values if actually different from original
        if adjustment_made and temp_start <= temp_end:
            if temp_start != _chain_start or temp_end != _chain_end:
                adjusted_start = temp_start
                adjusted_end = temp_end

    # Now use adjusted coordinates for genomic calculations if they exist
    if adjusted_start is not None and adjusted_end is not None:
        aa_start_for_coords = adjusted_start
        aa_end_for_coords = adjusted_end
    else:
        aa_start_for_coords = aa_start
        aa_end_for_coords = aa_end

    # --- Override C-term to Ensembl length if UniProt shorter and no C features ---
    effective_aa_end = aa_end_for_coords
    override_note = None
    if mode == "chain" and aa_len and aa_end_for_coords < aa_len:
        c_features_present = any_nonempty(row.get("Propeptides")) or any_nonempty(row.get("Transit_peptides"))
        if not c_features_present:
            effective_aa_end = aa_len
            override_note = f"override_c_end_to_ensembl_len({aa_len})"

    # Reasons (side-specific) using the effective end
    n_reason, c_reason = build_chain_reason_pair(row, mode == "chain", int(aa_start_for_coords) if aa_len else None, int(effective_aa_end) if aa_len else None, int(aa_len) if aa_len else None)
    print(f"[debug] mode={mode} aa_len={aa_len} chain=({original_chain_start},{original_chain_end}) adjusted=({adjusted_start},{adjusted_end}) coords=({aa_start_for_coords},{aa_end_for_coords}) eff_end={effective_aa_end} for tx={tx_id}; reasons N='{n_reason}' C='{c_reason}' override={override_note}")

    # AA -> CDS nt span (using coordinates for genomic mapping)
    def aa_to_cds_nt_span(aidx: int) -> Tuple[int,int]:
        s = (aidx - 1) * 3 + 1
        e = s + 2
        return s, e

    n_cds_span = aa_to_cds_nt_span(int(aa_start_for_coords)) if cds_len_nt else None
    c_cds_span = aa_to_cds_nt_span(int(effective_aa_end)) if cds_len_nt else None

    # --- cDNA genomic endpoints (first and last coding nucleotides) ---
    cdna_start_chrom = cdna_start_pos = None
    cdna_end_chrom = cdna_end_pos = None
    try:
        if cds_len_nt and tx_id:
            s_span = aa_to_cds_nt_span(1)
            e_span = aa_to_cds_nt_span(int(aa_len))
            s_maps = map_cds_to_genome(session, tx_id, s_span[0], s_span[1])
            e_maps = map_cds_to_genome(session, tx_id, e_span[0], e_span[1])
            if s_maps:
                s_mb0 = s_maps[0]
                s_strand = int(s_mb0.get("strand", 0) or 0)
                s_chrom = s_mb0.get("seq_region_name")
                s_start_gen = min(int(m.get("start")) for m in s_maps if "start" in m)
                s_end_gen   = max(int(m.get("end"))   for m in s_maps if "end" in m)
                cdna_start_chrom = s_chrom
                cdna_start_pos   = s_start_gen if s_strand == 1 else s_end_gen
            if e_maps:
                e_mb0 = e_maps[0]
                e_strand = int(e_mb0.get("strand", 0) or 0)
                e_chrom = e_mb0.get("seq_region_name")
                e_start_gen = min(int(m.get("start")) for m in e_maps if "start" in m)
                e_end_gen   = max(int(m.get("end"))   for m in e_maps if "end" in m)
                cdna_end_chrom = e_chrom
                cdna_end_pos   = e_end_gen if e_strand == 1 else e_start_gen

    except Exception:
        cdna_start_chrom = cdna_start_pos = cdna_end_chrom = cdna_end_pos = None

    chrom_N = pos_N = strand_symbol = None
    chrom_C = pos_C = None

    if n_cds_span:
        print(f"[debug] N aa={aa_start_for_coords} cds_span={n_cds_span}")
        n_maps = map_cds_to_genome(session, tx_id, n_cds_span[0], n_cds_span[1])
        if n_maps:
            chrom_N, pos_N, strand_symbol = insertion_from_codon_blocks(n_maps, 'N')
            print(f"[debug] N mapped -> {chrom_N}:{pos_N} strand={strand_symbol}")
        else:
            print("[warn] N cds map failed")

    if c_cds_span:
        print(f"[debug] C aa={effective_aa_end} cds_span={c_cds_span}")
        c_maps = map_cds_to_genome(session, tx_id, c_cds_span[0], c_cds_span[1])
        if c_maps:
            chrom_C, pos_C, strand_symbol2 = insertion_from_codon_blocks(c_maps, 'C')
            print(f"[debug] C mapped -> {chrom_C}:{pos_C} strand={strand_symbol2}")
            if strand_symbol is None:
                strand_symbol = strand_symbol2
        else:
            print("[warn] C cds map failed")

    warn = None
    if chrom_N is None or chrom_C is None:
        warn = "cds_map_failed; used cdna fallback"
        cdna_start = int(tr["start"])
        cdna_end   = int(tr["end"])
        if mode == "chain":
            boundary_cdna_n = cdna_start + (int(aa_start_for_coords) - 1) * 3
            boundary_cdna_c = cdna_start + (int(effective_aa_end) * 3)
        else:
            boundary_cdna_n = cdna_start
            boundary_cdna_c = cdna_end + 1
        bfr_n = map_cdna_to_genomic(session, tx_id, boundary_cdna_n - 1 if boundary_cdna_n > 1 else None)
        bfr_c = map_cdna_to_genomic(session, tx_id, boundary_cdna_c - 1 if boundary_cdna_c > 1 else None)
        if not chrom_N and bfr_n:
            chrom_N, pos_N, strand3 = bfr_n
            strand_symbol = "+" if strand3 == 1 else "-"
        if not chrom_C and bfr_c:
            chrom_C, pos_C, strand3 = bfr_c
            strand_symbol = "+" if strand3 == 1 else "-"

    # Build base record with shared fields
    base_record = {
        "gene": gene_name,
        "WBGene": wbgene,
        "input_transcript": t_in,
        "uniprot": uniprot,
        "protein_len": row.get("protein_len"),
        "mode": mode,
        "chain_aa_start": int(original_chain_start) if aa_len else None,
        "Chain_adjusted_start": adjusted_start, 
        "Chain_adjusted_end": adjusted_end,   
        "Chain_adjusted_regions": f"{adjusted_start}_{adjusted_end}" if (adjusted_start is not None and adjusted_end is not None) else None, 
        "strand": strand_symbol,
        "cdna_start_chrom": cdna_start_chrom, 
        "cdna_start_pos": cdna_start_pos,
        "cdna_end_chrom": cdna_end_chrom, 
        "cdna_end_pos": cdna_end_pos,
        "source": "cds->genome" if not warn else "mixed",
        "warning": warn
    }
    
    # N-terminal specific fields
    out = []
    out.append({
        **base_record,
        "site_type": "Nterm",
        "chain_reason": n_reason,
        "chain_aa_end": int(original_chain_end) if aa_len else None,
        "chrom": chrom_N, 
        "pos": pos_N,
        "bases_from_cdna_start": (n_cds_span[0]) if n_cds_span else None,
        "bases_from_cdna_end": (cds_len_nt - n_cds_span[0]) if (cds_len_nt and n_cds_span) else None,
        "boundary_cdna": n_cds_span[0] if n_cds_span else None, 
        "before_cdna": (n_cds_span[0] - 1) if (n_cds_span and n_cds_span[0] > 1) else None,
        "boundary_chrom": chrom_N, 
        "boundary_pos": pos_N, 
        "boundary_strand": strand_symbol,
    })
    
    # C-terminal specific fields
    out.append({
        **base_record,
        "site_type": "Cterm",
        "chain_reason": c_reason,
        "chain_aa_end": int(original_chain_end) if aa_len else None,
        "chrom": chrom_C, 
        "pos": pos_C,
        "bases_from_cdna_start": (c_cds_span[1]) if c_cds_span else None,
        "bases_from_cdna_end": (cds_len_nt - c_cds_span[1]) if (cds_len_nt and c_cds_span) else None,
        "boundary_cdna": c_cds_span[1] if c_cds_span else None, 
        "before_cdna": (c_cds_span[1] - 1) if (c_cds_span and c_cds_span[1] > 1) else None,
        "boundary_chrom": chrom_C, 
        "boundary_pos": pos_C, 
        "boundary_strand": strand_symbol,
    })

    # --- UNIFIED enrichment pass: warnings, reasons, features ---
    _enrich_site_records(out, row)
    
    return out


def _enrich_site_records(out: List[Dict], row: pd.Series) -> None:
    """
    Single unified pass to enrich site records with:
    - Unknown boundary warnings
    - CAAX/PTS/ER motif warnings  
    - Lipidation and CAAX mentions in reasons
    - Features summary
    - Side-specific reason cleanup
    - Warning deduplication
    """
    # 1) Collect warnings for CAAX/PTS/ER motifs
    warn_msgs = []
    
    motif_warnings = [
        ("PTS2", "Predicted PTS2 motif - consider avoiding N-terminal tag"),
        ("PTS1", "Predicted PTS1 motif - consider avoiding C-terminal tag"),
        ("ER_KDEL_motif", "Predicted ER KDEL motif (C-terminal retention) - consider avoiding C-terminal tag"),
        ("ER_diLys_motif", "Predicted ER di-Lys motif (KKXX/KXKXX) - consider avoiding C-terminal tag"),
        ("CAAX", "Predicted CaaX prenylation site - consider avoiding C-terminal tag"),
    ]
    for col, msg in motif_warnings:
        val = row.get(col) if hasattr(row, "get") else None
        if val is not None and str(val).strip() and str(val).strip().lower() != "nan":
            if msg not in warn_msgs:
                warn_msgs.append(msg)

    # 2) Collect unknown boundary warnings
    coord_cols = ["Transit_peptides", "Chain_regions", "SP_regions", "Propeptides", "Lipidation_sites", "Motifs"]
    for col in coord_cols:
        val = row.get(col) if hasattr(row, "get") else None
        details = _detect_unknown_terminals(val)
        for detail in details:
            if "N-terminal" in detail:
                msg = f"Unknown N-terminal boundary in {col}: {val}"
            elif "C-terminal" in detail:
                msg = f"Unknown C-terminal boundary in {col}: {val}"
            else:
                msg = f"Unknown coordinate in {col}: {val}"
            if msg not in warn_msgs:
                warn_msgs.append(msg)

    # 3) Build features summary string
    feat_cols = ["SP_regions", "Propeptides", "Transit_peptides", "iMET_events",
                  "Lipidation_sites", "Motifs", "Chain_regions", "PTS1", "PTS2",
                  "ER_KDEL_motif", "ER_diLys_motif", "CAAX"]
    pairs = []
    for col in feat_cols:
        try:
            val = row.get(col) if hasattr(row, "get") else None
        except Exception:
            val = None
        pairs.append(f"{col}={val}")
    features_str = "features: " + "; ".join(pairs)

    # 4) Enrich each record
    for rec in out:
        site_type = str(rec.get("site_type") or "").lower()
        
        # Add warnings
        if warn_msgs:
            prev = rec.get("warning")
            msg = " | ".join(warn_msgs)
            rec["warning"] = _merge_warning(prev, msg)
        
        # Set features summary
        if not rec.get("features_summary"):
            rec["features_summary"] = features_str
        
        # Enrich chain_reason with lipidation and CAAX info
        side_reason = rec.get("chain_reason") or ""
        
        # Add lipidation mention (only once, if not already present)
        if "lipidation site" not in side_reason:
            lip = row.get("Lipidation_sites") if hasattr(row, "get") else None
            if lip is not None and str(lip).strip() and str(lip).strip().lower() != "nan":
                lip_ranges = _extract_ranges_from_text(lip)
                protein_len = rec.get("protein_len")
                try:
                    protein_len = int(protein_len) if protein_len == protein_len and protein_len is not None else None
                except Exception:
                    protein_len = None
                
                lip_note = "lipidation site present"
                if protein_len and lip_ranges:
                    # CAAX is typically at L-3; treat anything within last 5 aa as C-terminal
                    if any((a <= protein_len and b >= protein_len-5) or (a >= protein_len-5) for a,b in lip_ranges):
                        lip_note = "C-terminal lipidation site"
                    elif any(a <= 5 for a,b in lip_ranges):
                        lip_note = "N-terminal lipidation site"
                side_reason = (side_reason + (" | " if side_reason else "") + lip_note)

        # Add CAAX motif text (only once, if not already present)
        if "CAAX motif" not in side_reason and "CaaX motif" not in side_reason:
            caax = row.get("CAAX") if hasattr(row, "get") else None
            if caax is not None and str(caax).strip() and str(caax).strip().lower() != "nan":
                side_reason = (side_reason + (" | " if side_reason else "") + f"CAAX motif {str(caax).strip()}")
        
        # Side-specific cleanup: remove inappropriate mentions
        tokens = [t.strip() for t in side_reason.split("|") if t.strip()]
        if site_type.startswith("n"):
            # Remove C-terminal specific mentions from N-term
            tokens = [t for t in tokens if not t.lower().startswith("c-terminal lipidation site") 
                     and not t.lower().startswith("caax motif")]
        elif site_type.startswith("c"):
            # Remove CAAX from C-term reason (it's already in warnings)
            tokens = [t for t in tokens if not t.lower().startswith("caax motif")]
        
        # Deduplicate tokens while preserving order
        seen = set()
        deduped = []
        for t in tokens:
            if t not in seen:
                seen.add(t)
                deduped.append(t)
        
        rec["chain_reason"] = " | ".join(deduped)
        
        # Final warning deduplication
        _dedupe_warning_field(rec)

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

def _existing_done_pairs(out_path: str) -> set[tuple]:
    """Return {(WBGene, canonical_transcript)} pairs that already have both Nterm & Cterm."""
    done = set()
    try:
        if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
            df = pd.read_csv(out_path, usecols=["WBGene","input_transcript","site_type"], dtype=str, engine="python")
            df = df.dropna(subset=["WBGene","input_transcript","site_type"])
            df["WBGene"] = df["WBGene"].str.strip().str.lstrip("\ufeff")
            df["input_transcript"] = df["input_transcript"].str.strip()
            grp = df.groupby(["WBGene","input_transcript"])["site_type"].apply(lambda s: set((x or "").lower() for x in s))
            for (wb, tx), st in grp.items():
                if {"nterm","cterm"}.issubset(st):
                    done.add((wb, tx))
    except Exception:
        pass
    return done

# ===============================================================================
# Main
# ===============================================================================
def main():
    t0 = time.perf_counter()
    ap = argparse.ArgumentParser(description="Get chain-aware N/C-terminal insertion sites")
    ap.add_argument("--input", required=True, help="Input CSV")
    ap.add_argument("--output", required=True, help="Output CSV path")
    ap.add_argument("--flush-every", type=int, default=10, help="Write results to CSV every N genes (default: 10)")
    ap.add_argument("--no-resume", action="store_true", help="Do not resume; overwrite output")
    ap.add_argument("--max-rows", type=int, default=None, help="Limit rows (debugging)")
    args = ap.parse_args()

    df = pd.read_csv(args.input)
    if args.max_rows:
        df = df.head(args.max_rows)

    # Prepare output for incremental write/resume
    cols_for_write = ['gene', 'WBGene', 'input_transcript', 'uniprot', 'protein_len', 
                    'site_type', 'mode', 'chain_reason', 'chain_aa_start', 'chain_aa_end', 
                    'Chain_adjusted_start', 'Chain_adjusted_end', 'Chain_adjusted_regions',
                    'chrom', 'pos', 'strand', 'cdna_start_pos', 'bases_from_cdna_start', 
                    'cdna_end_pos', 'bases_from_cdna_end', 'features_summary', 'source', 'warning']
    out_path = args.output
    # initialize/overwrite per --no-resume
    if args.no_resume and os.path.exists(out_path):
        open(out_path, "w").close()
    # create header if file is empty/missing
    if not os.path.exists(out_path) or os.path.getsize(out_path) == 0:
        pd.DataFrame(columns=cols_for_write).to_csv(out_path, index=False)
    # determine already-completed (WBGene, input_transcript)
    done_pairs = set() if args.no_resume else _existing_done_pairs(out_path)
    if done_pairs and not args.no_resume:
        print(f"[resume] Detected {len(done_pairs)} completed gene(s); will skip those.", file=sys.stderr)

    print("[start] Building HTTP session...", file=sys.stderr)
    session = build_http_session()
    print("[start] Reading gene list...", file=sys.stderr)

    rows_out: List[Dict] = []
    buffer: List[Dict] = []
    flush_every = max(1, int(args.flush_every))
    total = len(df)
    
    print(f"[start] Processing {total} genes...", file=sys.stderr)
    
    for i, row in df.iterrows():
        wb = row.get('WBGene', '')
        tx = row.get('canonical_transcript', '')
        gene_name = get_gene_name(row)
        
        print(f"\n=== [{i+1}/{total}] {wb} ===", file=sys.stderr)
        print(f"  [gene] name={gene_name} transcript={tx}", file=sys.stderr)
        
        key = (str(wb).strip(), str(tx).strip())
        if key[0] and key[1] and key in done_pairs:
            print(f"  [skip] already have both N/C sites", file=sys.stderr)
            continue
        
        try:
            rows = compute_sites_for_row(session, row)
            rows_out.extend(rows)
            buffer.extend(rows)
            
            # Print summary of what was computed
            nterm = [r for r in rows if (r.get('site_type') or '').lower() == 'nterm']
            cterm = [r for r in rows if (r.get('site_type') or '').lower() == 'cterm']
            
            if nterm:
                n = nterm[0]
                print(f"  [Nterm] mode={n.get('mode')} chrom={n.get('chrom')}:{n.get('pos')} chain_aa={n.get('chain_aa_start')}-{n.get('chain_aa_end')}", file=sys.stderr)
                if n.get('warning'):
                    print(f"    [warn] {n.get('warning')}", file=sys.stderr)
            
            if cterm:
                c = cterm[0]
                print(f"  [Cterm] mode={c.get('mode')} chrom={c.get('chrom')}:{c.get('pos')} chain_aa={c.get('chain_aa_start')}-{c.get('chain_aa_end')}", file=sys.stderr)
                if c.get('warning'):
                    print(f"    [warn] {c.get('warning')}", file=sys.stderr)
            
            print(f"  [done] produced {len(rows)} site records", file=sys.stderr)
            
        except Exception as e:
            print(f"  [error] {type(e).__name__}: {e}", file=sys.stderr)
            for st in ("Nterm","Cterm"):
                buffer.append({
                    "gene": gene_name,
                    "WBGene": wb,
                    "input_transcript": tx,
                    "uniprot": row.get("uniprot"),
                    "protein_len": row.get("protein_len"),
                    "site_type": st,
                    "mode": "exception",
                    "chain_reason": "",
                    "chain_aa_start": None,
                    "chain_aa_end": None,
                    "chrom": None, "pos": None, "strand": None,
                    "cdna_start_pos": None,
                    "bases_from_cdna_start": None,
                    "cdna_end_pos": None,
                    "bases_from_cdna_end": None,
                    "features_summary": None,
                    "source": "exception",
                    "warning": f"error: {type(e).__name__}: {e}"
                })
        
        # mark as done if both sites present
        try:
            sts = {(r.get("site_type") or "").lower() for r in buffer[-2:]}
            if {"nterm","cterm"}.issubset(sts) and key[0] and key[1]:
                done_pairs.add(key)
        except Exception:
            pass
        
        if ((i+1) % flush_every) == 0:
            print(f"\n[flush] Processed {i+1}/{total} genes, writing {len(buffer)} records to disk...", file=sys.stderr)
            _write_rows(out_path, buffer, cols_for_write)
            buffer.clear()

    # final flush
    if buffer:
        print(f"\n[flush] Final flush: writing {len(buffer)} records...", file=sys.stderr)
        _write_rows(out_path, buffer, cols_for_write)
        buffer.clear()

    print(f"\n[done] Wrote {len(rows_out)} site records to {args.output}", file=sys.stderr)
    # timing
    dt = time.perf_counter() - t0
    rate = (total / dt) if dt > 0 else float('inf')
    print(f"[timing] Elapsed {dt:.2f}s | {total} genes | {rate:.2f} genes/s", file=sys.stderr)

if __name__ == "__main__":
    main()