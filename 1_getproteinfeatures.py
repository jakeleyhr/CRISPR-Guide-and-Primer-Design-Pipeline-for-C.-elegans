#!/usr/bin/env python3
"""
Protein Feature Annotation Tool for C. elegans Protein-Coding Genes

Retrieves protein features for WormBase genes by integrating Ensembl transcript data
with and UniProt protein features data.

Input:
  - CSV file containing WBGene IDs (column auto-detected: 'WBGene', 'WBGeneID', 
    'gene', or any column containing WBGene identifiers)

Output:
  - CSV file with one row per WBGene containing:
      • Gene & transcript identifiers (WBGene, canonical_transcript, UniProt ID)
      • Match type (canonical/non-canonical) and protein length
      • Sequence motifs: CAAX, PTS1, PTS2, ER retention signals (KDEL, KKxx)
      • UniProt features: signal peptides, transmembrane domains, lipidation sites,
        propeptides, transit peptides, initiator methionine events, chain processing
      • Compact feature summary tag for quick filtering

Transcript Selection Strategy:
  1. Canonical match: Uses Ensembl canonical transcript with corresponding UniProt entry
  2. Non-canonical match: If canonical lacks UniProt, finds longest UniProt sequence
     for the gene and matches it to the corresponding transcript
  3. If no Uniprot entry exists for the gene, outputs N/A with message

Features:
  - If interrupted, automatically resumes from where it left off by skipping genes already 
     present in the output file
  - Sequence-based validation to ensure Ensembl-UniProt entries match

Required arguments:
  --input PATH      Input CSV with WBIDs
  --output PATH     Output CSV path with added protein features
  
Optional arguments:
  --gene-col COLUMN     Specify WBGene column name (default: auto-detect)
  --flush-every N       Write results every N genes (default: 10)
  --no-resume           Overwrite existing output instead of resuming

Example:
  python 1_getproteinfeatures.py \
    --input genes_list.csv \
    --output 1.protein_features.csv \
    --flush-every 100
"""

import argparse, os, re, sys, time, requests
import pandas as pd
from typing import Dict, List, Optional, Tuple
from requests.adapters import HTTPAdapter, Retry

ENSEMBL_REST = "https://rest.ensembl.org"
SPECIES = "caenorhabditis_elegans"  # WBcel235
UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb/"
ALIPHATIC = set("AVILPM")  # for CAAX rule

# ===============================================================================
# UniProt Features of Interest
# ===============================================================================

KEEP = {
    "Signal": "SP",
    "Transit peptide": "TP",
    "Transmembrane": "TM",
    "Lipidation": "LIP",
    "Motif": "MOTIF",
    "Chain": "CHAIN",
    "Initiator methionine": "iMET",
    "Modified residue": "MOD",
    "Cross-link": "XLINK",
    "Glycosylation": "GLYCO",
    "Propeptide": "PRO",
    "Peptide": "PEP",
}

_FIELDNAMES = ["WBGene","canonical_transcript","uniprot","match_type","protein_len","CAAX","PTS1","PTS2","ER_KDEL_motif","ER_diLys_motif","iMET_events","SP_regions","Propeptides","Transit_peptides","Lipidation_sites","Motifs","TM_regions","Chain_regions","features","notes"]

# ===============================================================================
# HTTP Helpers
# ===============================================================================

def make_session() -> requests.Session:
    s = requests.Session()
    retries = Retry(
        total=6,
        backoff_factor=0.5,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET", "POST"),
    )
    s.mount("https://", HTTPAdapter(max_retries=retries))
    s.headers.update({
        "User-Agent": "celegans-canonical-feature-screener/1.1",
        "Accept": "application/json",
    })
    return s

def get_json(session: requests.Session, path: str, params: Dict = None) -> Optional[dict]:
    url = f"{ENSEMBL_REST}{path}"
    r = session.get(url, params=params or {}, timeout=45)
    if r.status_code == 200:
        return r.json()
    try:
        print(f"[http] {path} -> {r.status_code}: {r.json()}", file=sys.stderr)
    except Exception:
        print(f"[http] {path} -> {r.status_code}", file=sys.stderr)
    return None

def get_text(session: requests.Session, path: str, params: Dict = None) -> Optional[str]:
    url = f"{ENSEMBL_REST}{path}"
    r = session.get(url, params=params or {}, headers={"Accept":"text/plain"}, timeout=45)
    if r.status_code == 200:
        return r.text.strip()
    print(f"[http] {path} -> {r.status_code}", file=sys.stderr)
    return None

def resolve_gene_id(session: requests.Session, token: str) -> Optional[str]:
    token = (token or "").strip()
    if not token:
        return None
    # direct
    j = get_json(session, f"/lookup/id/{token}", {"expand": 0})
    if j and j.get("object_type") == "Gene":
        print(f"  [resolve] {token} -> {j['id']} (direct)", file=sys.stderr)
        return j["id"]
    # symbol
    xs = get_json(session, f"/xrefs/symbol/{SPECIES}/{token}")
    if isinstance(xs, list):
        for x in xs:
            if x.get("object_type") == "Gene":
                print(f"  [resolve] {token} -> {x['id']} (symbol)", file=sys.stderr)
                return x.get("id")
    # name
    xs = get_json(session, f"/xrefs/name/{SPECIES}/{token}")
    if isinstance(xs, list):
        for x in xs:
            if x.get("object_type") == "Gene":
                print(f"  [resolve] {token} -> {x['id']} (name)", file=sys.stderr)
                return x.get("id")
    print(f"  [resolve][warn] could not resolve: {token}", file=sys.stderr)
    return None

def get_all_transcripts_with_translations(session: requests.Session, gene_id: str) -> List[Tuple[str, str, int]]:
    """
    Get all translated transcripts for a gene.
    Returns list of (transcript_id, translation_id, aa_len) tuples.
    """
    gj = get_json(session, f"/lookup/id/{gene_id}", {"expand": 1})
    if not gj:
        print("  [lookup][warn] /lookup failed", file=sys.stderr)
        return []
    
    txs = gj.get("Transcript", []) or []
    print(f"  [lookup] transcripts={len(txs)}", file=sys.stderr)

    pcs: List[Tuple[str, str, int]] = []
    for tx in txs:
        tr = tx.get("Translation")
        if isinstance(tr, dict) and tr.get("id"):
            aa_len = tr.get("length")
            if aa_len is None and "start" in tr and "end" in tr:
                aa_len = (int(tr["end"]) - int(tr["start"]) + 1) // 3
            if aa_len and int(aa_len) > 0:
                pcs.append((tx["id"], tr["id"], int(aa_len)))

    print(f"  [lookup] translated_transcripts={len(pcs)}", file=sys.stderr)
    return pcs

def pick_canonical_transcript_and_translation(session: requests.Session, gene_id: str) -> Optional[Tuple[str,str,int,Optional[str],str]]:
    """
    Return (transcript_id, translation_id, aa_len, uniprot_id, match_type) for best transcript.
    
    Uses EXACT same logic as original code for transcript selection:
    1. Try canonical transcript, get UniProt for it
    2. If no canonical, use longest CDS
    
    NEW BEHAVIOR: If the chosen transcript has NO UniProt match, find longest UniProt 
    for the gene and use its corresponding transcript instead (marked as "non-canonical").
    
    match_type is either "canonical" or "non-canonical"
    """
    gj = get_json(session, f"/lookup/id/{gene_id}", {"expand": 1})
    if not gj:
        print("  [lookup][warn] /lookup failed", file=sys.stderr)
        return None
    
    txs = gj.get("Transcript", []) or []
    print(f"  [lookup] transcripts={len(txs)}", file=sys.stderr)

    pcs: List[Tuple[str, str, int]] = []
    for tx in txs:
        tr = tx.get("Translation")
        if isinstance(tr, dict) and tr.get("id"):
            aa_len = tr.get("length")
            if aa_len is None and "start" in tr and "end" in tr:
                aa_len = (int(tr["end"]) - int(tr["start"]) + 1) // 3
            if aa_len and int(aa_len) > 0:
                pcs.append((tx["id"], tr["id"], int(aa_len)))

    print(f"  [lookup] translated_transcripts={len(pcs)}", file=sys.stderr)
    if not pcs:
        return None

    # Find canonical transcript
    can_tx = gj.get("canonical_transcript")
    chosen_tx = None
    if can_tx:
        for tid, trid, aa in pcs:
            if tid == can_tx:
                chosen_tx = (tid, trid, aa)
                print(f"  [canonical] transcript={tid} protein={trid} aa={aa}", file=sys.stderr)
                break
    
    if not chosen_tx:
        pcs.sort(key=lambda x: x[2], reverse=True)
        chosen_tx = pcs[0]
        print(f"  [fallback] longest_cds transcript={chosen_tx[0]} protein={chosen_tx[1]} aa={chosen_tx[2]}", file=sys.stderr)
    
    # Now get UniProt ID for the chosen transcript
    tid, trid, aa = chosen_tx
    uniprot_id = _get_matching_uniprot_for_transcript(session, gene_id, trid, tid)
    
    print(f"  [uniprot_match] transcript={tid} -> uniprot={uniprot_id}", file=sys.stderr)
    
    # If UniProt found, return as canonical. Otherwise, try non-canonical strategy.
    if uniprot_id:
        # Original behavior: canonical transcript has UniProt
        print(f"  [canonical] using canonical transcript with UniProt", file=sys.stderr)
        return (tid, trid, aa, uniprot_id, "canonical")
    
    # If no UniProt for canonical transcript - find longest UniProt and its transcript
    print(f"  [non-canonical] canonical transcript has no UniProt, searching for longest UniProt...", file=sys.stderr)
    
    # Get all UniProt IDs for this gene from xrefs
    all_uniprot_ids = _get_all_uniprot_ids_for_gene(session, gene_id)
    
    if not all_uniprot_ids:
        print(f"  [non-canonical] no UniProt IDs found for gene at all", file=sys.stderr)
        # Return the canonical transcript with N/A for UniProt and informative match_type
        return (tid, trid, aa, "N/A", "No UniProt match found, using canonical transcript")
    
    # Find the longest UniProt sequence
    longest_uniprot = None
    longest_length = 0
    longest_seq = None
    
    for uniprot_id in all_uniprot_ids:
        try:
            # Fetch sequence
            resp = session.get(f"{UNIPROT_BASE}{uniprot_id}.fasta", timeout=20)
            if resp.status_code == 200:
                lines = resp.text.strip().split('\n')
                seq = ''.join(lines[1:])
                seq_len = len(seq)
                
                if seq_len > longest_length:
                    longest_length = seq_len
                    longest_uniprot = uniprot_id
                    longest_seq = seq
                    print(f"  [non-canonical] found {uniprot_id} with length {seq_len}", file=sys.stderr)
        except Exception as e:
            print(f"  [non-canonical][warn] failed to fetch {uniprot_id}: {e}", file=sys.stderr)
    
    if not longest_uniprot:
        print(f"  [non-canonical] could not fetch any UniProt sequences", file=sys.stderr)
        # Return canonical transcript with N/A for UniProt and informative match_type
        return (tid, trid, aa, "N/A", "No UniProt match found, using canonical transcript")
    
    print(f"  [non-canonical] longest UniProt: {longest_uniprot} ({longest_length} aa)", file=sys.stderr)
    
    # Now find which transcript matches this UniProt sequence
    for tx_id, tr_id, tx_aa in pcs:
        ensembl_seq = get_text(session, f"/sequence/id/{tr_id}", {"type": "protein"})
        if ensembl_seq and ensembl_seq.strip() == longest_seq.strip():
            print(f"  [non-canonical] MATCH! transcript={tx_id} matches {longest_uniprot}", file=sys.stderr)
            return (tx_id, tr_id, tx_aa, longest_uniprot, "non-canonical")
    
    # No exact transcript match found - use longest UniProt with its length, keep original transcript
    print(f"  [non-canonical] no matching transcript found, using canonical transcript={tid} with longest UniProt={longest_uniprot}", file=sys.stderr)
    return (tid, trid, longest_length, longest_uniprot, "non-canonical")

def _get_all_uniprot_ids_for_gene(session: requests.Session, gene_id: str) -> List[str]:
    """
    Get all UniProt IDs associated with a gene from Ensembl xrefs.
    """
    uniprot_ids = []
    
    # Get gene-level xrefs
    gene_xrefs = get_json(session, f"/xrefs/id/{gene_id}")
    if isinstance(gene_xrefs, list):
        uniprot_dbnames = ["Uniprot_gn", "Uniprot/SWISSPROT", "Uniprot/SPTREMBL", 
                          "UniProtKB/Swiss-Prot", "UniProtKB/TrEMBL", "Uniprot"]
        for xref in gene_xrefs:
            dbname = xref.get("dbname", "")
            if dbname in uniprot_dbnames:
                acc = xref.get("primary_id")
                if acc and acc not in uniprot_ids:
                    uniprot_ids.append(acc)
    
    print(f"    [uniprot_ids] found {len(uniprot_ids)} UniProt IDs for gene: {', '.join(uniprot_ids)}", file=sys.stderr)
    return uniprot_ids

def _query_uniprot_id_mapping(session: requests.Session, wbgene_id: str, gene_name: str, 
                              transcript_id: str, ensembl_protein_id: str, protein_sequence: str) -> List[str]:
    """
    Query UniProt's ID mapping service and sequence search to find UniProt IDs.
    Tries multiple approaches:
    1. WormBase Gene ID (WBGene)
    2. Gene name
    3. Transcript ID
    4. Ensembl Protein ID
    5. Direct sequence search
    """
    all_uniprot_ids = []
    
    # Try different ID types
    id_attempts = [
        (wbgene_id, "WormBase", "WBGene IDs"),
        (gene_name, "Gene_Name", "Gene_Name"),
        (transcript_id, "WormBase_Transcript", "WormBase_Transcript"),
        (ensembl_protein_id, "Ensembl_Protein", "Ensembl_Protein"),
    ]
    
    for identifier, from_db, db_label in id_attempts:
        if not identifier:
            continue
            
        try:
            print(f"    [id_mapping] trying {db_label}: {identifier}", file=sys.stderr)
            
            # Submit ID mapping job
            url = "https://rest.uniprot.org/idmapping/run"
            data = {
                "ids": identifier,
                "from": from_db,
                "to": "UniProtKB"
            }
            
            response = session.post(url, data=data, timeout=30)
            if response.status_code != 200:
                print(f"    [id_mapping] {db_label} submission failed: {response.status_code}", file=sys.stderr)
                continue
            
            job_id = response.json().get("jobId")
            if not job_id:
                print(f"    [id_mapping] {db_label} no job ID returned", file=sys.stderr)
                continue
            
            # Poll for results
            results_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
            max_attempts = 8
            for attempt in range(max_attempts):
                time.sleep(0.5)  # Wait 0.5 seconds between polls
                status_resp = session.get(results_url, timeout=30)
                if status_resp.status_code == 200:
                    status_data = status_resp.json()
                    if "results" in status_data or "failedIds" in status_data:
                        # Job complete
                        results = status_data.get("results", [])
                        for result in results:
                            to_entry = result.get("to", {})
                            if isinstance(to_entry, dict):
                                uniprot_id = to_entry.get("primaryAccession")
                                if uniprot_id:
                                    all_uniprot_ids.append(uniprot_id)
                                    print(f"    [id_mapping] {db_label} found: {uniprot_id}", file=sys.stderr)
                        break
            
        except Exception as e:
            print(f"    [id_mapping] {db_label} error: {e}", file=sys.stderr)
            continue
    
    # If found IDs via ID mapping, return them
    if all_uniprot_ids:
        # Remove duplicates while preserving order
        return list(dict.fromkeys(all_uniprot_ids))
    
    # Last resort: sequence-based search
    if protein_sequence:
        try:
            print(f"    [id_mapping] trying sequence search (length={len(protein_sequence)} aa)", file=sys.stderr)
            
            # Use UniProt's advanced search with sequence
            # Search for C. elegans proteins with similar length
            search_url = "https://rest.uniprot.org/uniprotkb/search"
            params = {
                "query": f"(organism_id:6239) AND (length:[{len(protein_sequence)-5} TO {len(protein_sequence)+5}])",
                "format": "json",
                "size": 100  # Get up to 100 candidates
            }
            
            search_resp = session.get(search_url, params=params, timeout=30)
            if search_resp.status_code == 200:
                results = search_resp.json().get("results", [])
                print(f"    [id_mapping] sequence search found {len(results)} length-matched candidates", file=sys.stderr)
                
                # Return all candidates for sequence validation
                candidates = []
                for entry in results:
                    acc = entry.get("primaryAccession")
                    if acc:
                        candidates.append(acc)
                
                if candidates:
                    print(f"    [id_mapping] returning {len(candidates)} candidates for sequence validation", file=sys.stderr)
                    return candidates[:20]  # Limit to 20 candidates for performance
                    
        except Exception as e:
            print(f"    [id_mapping] sequence search error: {e}", file=sys.stderr)
    
    print(f"    [id_mapping] no mappings found via any method", file=sys.stderr)
    return []

def _get_matching_uniprot_for_transcript(session: requests.Session, gene_id: str, translation_id: str, transcript_id: str = None) -> Optional[str]:
    """
    Get UniProt ID for a transcript by:
    1. First trying Ensembl xrefs (gene-level and transcript-level)
    2. If not found, using UniProt's ID mapping service with the Ensembl protein ID
    3. Validating by comparing sequences
    """
    
    # APPROACH 1: Try Ensembl xrefs first
    all_xrefs = []
    
    # Get gene-level xrefs
    gene_xrefs = get_json(session, f"/xrefs/id/{gene_id}")
    if isinstance(gene_xrefs, list):
        all_xrefs.extend(gene_xrefs)
        print(f"    [uniprot_match] found {len(gene_xrefs)} gene-level xrefs", file=sys.stderr)
    
    # Get transcript-level xrefs (if transcript_id provided)
    if transcript_id:
        tx_xrefs = get_json(session, f"/xrefs/id/{transcript_id}")
        if isinstance(tx_xrefs, list):
            all_xrefs.extend(tx_xrefs)
            print(f"    [uniprot_match] found {len(tx_xrefs)} transcript-level xrefs", file=sys.stderr)
    
    # Extract all UniProt IDs from xrefs (from any UniProt database)
    uniprot_ids = []
    uniprot_dbnames = ["Uniprot_gn", "Uniprot/SWISSPROT", "Uniprot/SPTREMBL", 
                       "UniProtKB/Swiss-Prot", "UniProtKB/TrEMBL", "Uniprot"]
    for xref in all_xrefs:
        dbname = xref.get("dbname", "")
        if dbname in uniprot_dbnames:
            acc = xref.get("primary_id")
            if acc:
                uniprot_ids.append(acc)
    
    # Remove duplicates while preserving order
    uniprot_ids = list(dict.fromkeys(uniprot_ids))
    
    # APPROACH 2: If no UniProt found in xrefs, try UniProt's ID mapping service
    if not uniprot_ids:
        print(f"    [uniprot_match] no UniProt entries in Ensembl xrefs, trying UniProt ID mapping...", file=sys.stderr)
        
        # Get the protein sequence first (needed for ID mapping)
        ensembl_seq = get_text(session, f"/sequence/id/{translation_id}", {"type": "protein"})
        if not ensembl_seq:
            print(f"    [uniprot_match] could not fetch Ensembl sequence for {translation_id}", file=sys.stderr)
            return None
        
        # Extract gene name from xrefs (often in wormbase_gseqname or display_id)
        gene_name = None
        for xref in all_xrefs:
            if xref.get("dbname") == "wormbase_gseqname":
                gene_name = xref.get("display_id")
                break
        if not gene_name:
            # Try to extract from transcript ID (e.g., C26D10.6a.1 -> C26D10.6)
            if transcript_id and "." in transcript_id:
                gene_name = transcript_id.rsplit(".", 1)[0].rstrip("ab")  # Remove isoform suffix
        
        print(f"    [uniprot_match] attempting ID mapping with: WBGene={gene_id}, gene_name={gene_name}, transcript={transcript_id}", file=sys.stderr)
        
        uniprot_ids = _query_uniprot_id_mapping(session, gene_id, gene_name, transcript_id, translation_id, ensembl_seq)
        if uniprot_ids:
            print(f"    [uniprot_match] found {len(uniprot_ids)} UniProt ID(s) via ID mapping: {', '.join(uniprot_ids)}", file=sys.stderr)
        else:
            print(f"    [uniprot_match] no UniProt mapping found", file=sys.stderr)
            return None
    else:
        print(f"    [uniprot_match] found {len(uniprot_ids)} UniProt candidates from xrefs: {', '.join(uniprot_ids)}", file=sys.stderr)
        # Get the Ensembl sequence for validation
        ensembl_seq = get_text(session, f"/sequence/id/{translation_id}", {"type": "protein"})
        if not ensembl_seq:
            print(f"    [uniprot_match] could not fetch Ensembl sequence for {translation_id}", file=sys.stderr)
            return None
    
    # At this point we have ensembl_seq and uniprot_ids ready for validation
    
    # Test all candidates (even if just one) - must match sequence
    print(f"    [uniprot_match] testing sequences to find match...", file=sys.stderr)
    for uniprot_id in uniprot_ids:
        try:
            # Fetch UniProt sequence
            resp = session.get(f"{UNIPROT_BASE}{uniprot_id}.fasta", timeout=20)
            if resp.status_code == 200:
                lines = resp.text.strip().split('\n')
                uniprot_seq = ''.join(lines[1:])  # Skip header
                
                if uniprot_seq.strip() == ensembl_seq.strip(): # Compares the amino acid sequences
                    print(f"    [uniprot_match] {uniprot_id}: MATCH (both {len(ensembl_seq)} aa)", file=sys.stderr)
                    return uniprot_id
                else:
                    print(f"    [uniprot_match] {uniprot_id}: mismatch (UniProt:{len(uniprot_seq)} vs Ensembl:{len(ensembl_seq)})", file=sys.stderr)
        except Exception as e:
            print(f"    [uniprot_match] error fetching {uniprot_id}: {e}", file=sys.stderr)
    
    # No match found
    print(f"    [uniprot_match] no sequence match found", file=sys.stderr)
    return None

def _rng(loc):
    """Return a start_end string even when end is unknown; include non-EXACT modifiers."""
    try:
        start = loc.get("start") or {}
        end = loc.get("end") or {}

        def fmt(p):
            # p like {"value": 1, "modifier": "EXACT"} or {"value": None, "modifier": "UNKNOWN"}
            if not p or p.get("value") is None:
                val = "?"
                mod = p.get("modifier")
            else:
                try:
                    val = str(int(p.get("value")))
                except Exception:
                    val = str(p.get("value"))
                mod = p.get("modifier")
            return f"{val}({mod})" if mod and mod != "EXACT" else val

        s = fmt(start)
        e = fmt(end)
        return f"{s}_{e}"
    except Exception:
        return ""

def fetch_uniprot_features_rich(uniprot_id: Optional[str], session: requests.Session) -> Dict[str, List[str]]:
    out = {k: [] for k in KEEP.values()}
    if not uniprot_id or uniprot_id == "N/A":
        return out
    try:
        r = session.get(f"{UNIPROT_BASE}{uniprot_id}.json", timeout=20)
        r.raise_for_status()
        data = r.json()
    except Exception as e:
        print(f"  [uniprot][warn] feature fetch failed for {uniprot_id}: {e}", file=sys.stderr)
        return out

    feats = data.get("features", [])
    print(f"  [uniprot] features fetched: {len(feats)} items", file=sys.stderr)
    for feat in feats:
        tag = KEEP.get(feat.get("type"))
        if not tag:
            continue
        loc = feat.get("location", {})
        rng = _rng(loc)
        desc = feat.get("description", "")
        if rng:
            s = f"{rng}" + (f" ({desc})" if desc else "")
        else:
            pos = (loc.get("position") or {}).get("value")
            s = (f"{pos}" if pos else "").strip()
            if desc:
                s = f"{s} ({desc})" if s else desc
        if s:
            out[tag].append(s)
    return out


# ===============================================================================
# Sequence Motifs
# ===============================================================================

def has_caax(seq: str) -> Tuple[bool, str]:
    """Return (has_caax, motif) where motif is last 4 AA if True else ''."""
    if not seq or len(seq) < 4:
        return False, ""
    motif = seq[-4:].upper()
    ok = (motif[0] == 'C') and ((motif[1] in ALIPHATIC) or (motif[2] in ALIPHATIC)) and motif[3].isalpha()
    return ok, (motif if ok else "")

def detect_pts1(seq: str):
    """
    C-terminal PTS1: classic [S/A/C/G/N/T][K/R/H][L/I/M].
    Returns (ok: bool, tripeptide: str).
    """
    if not seq or len(seq) < 3:
        return False, ""
    tri = seq[-3:].upper()
    if re.match(r'^[SACGNT][KRH][LIM]$', tri):
        return True, tri
    return False, ""

def detect_pts2(seq: str):
    """
    N-terminal PTS2 within first 40 aa: [R/K][L/V/I/Q]X{5}[H/Q][L/A/F].
    Returns (ok: bool, motif: str, start_index_1based: int)
    """
    if not seq:
        return False, "", -1
    s = seq[:40].upper()
    for i in range(len(s) - 8):
        window = s[i:i+9]
        if re.match(r'^[RK][LVIQ][A-Z]{5}[HQ][LAF]$', window):
            return True, window, i + 1
    return False, "", -1

def detect_er_kdel(seq: str):
    """
    Luminal ER-retention motif at C-terminus: (K|H|R)DEL at absolute end.
    Returns (ok: bool, motif: str)
    """
    if not seq or len(seq) < 4:
        return False, ""
    tail4 = seq[-4:].upper()
    if re.match(r'^[KHR]DEL$', tail4):
        return True, tail4
    return False, ""

def detect_er_di_lysine(seq: str):
    """
    Cytosolic ER retrieval motifs at C-terminus: KKxx or KxKxx at end.
    Returns (ok: bool, motif: str)
    """
    if not seq or len(seq) < 4:
        return False, ""
    up = seq.upper()
    tail4 = up[-4:]
    tail5 = up[-5:] if len(up) >= 5 else ""
    if re.match(r'^KK..$', tail4):
        return True, tail4
    if tail5 and re.match(r'^K.K..$', tail5):
        return True, tail5
    return False, ""

def autodetect_gene_col(df: pd.DataFrame) -> str:
    for cand in ("WBGene", "WBGeneID", "WBID", "wb_gene", "wb_gene_id"):
        if cand in df.columns:
            return cand
    for c in df.columns:
        if df[c].astype(str).str.contains(r"WBGene\d{8}", regex=True, na=False).any():
            return c
    raise ValueError("Could not find a WBGene ID column. Add --gene-col or rename a column to 'WBGeneID'.")

def _existing_done_set(out_path: str) -> set:
    try:
        if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
            df = pd.read_csv(out_path, usecols=["WBGene"])
            return set(df["WBGene"].astype(str).tolist())
    except Exception as e:
        print(f"[resume][warn] could not read existing output '{out_path}': {e}", file=sys.stderr)
    return set()

def _write_rows(out_path: str, rows: list):
    if not rows:
        return
    df = pd.DataFrame(rows)
    # ensure columns and order
    for col in _FIELDNAMES:
        if col not in df.columns:
            df[col] = "" if col not in ["protein_len"] else None
    df = df[_FIELDNAMES]
    header = not (os.path.exists(out_path) and os.path.getsize(out_path) > 0)
    df.to_csv(out_path, mode="a", header=header, index=False)
    print(f"[write] appended {len(df)} row(s) -> {out_path}", file=sys.stderr)

def _dedupe_file(out_path: str):
    """Keep only the last row per WBGene (append-safe compaction)."""
    try:
        df = pd.read_csv(out_path, dtype={"WBGene": str}, engine="python", on_bad_lines="skip")
        if "WBGene" not in df.columns:
            return
        # normalize
        df["WBGene"] = df["WBGene"].astype(str).str.strip().str.lstrip("\ufeff")
        df = df.dropna(subset=["WBGene"])
        df = df[df["WBGene"] != ""]
        df = df.drop_duplicates(subset=["WBGene"], keep="last")
        df.to_csv(out_path, index=False)
        print(f"[dedupe] compacted to {len(df)} unique WBGene rows", file=sys.stderr)
    except Exception as e:
        print(f"[dedupe][warn] could not dedupe: {e}", file=sys.stderr)

# ===============================================================================
# Main
# ===============================================================================
def main():
    t0 = time.perf_counter()  # Start timer

    ap = argparse.ArgumentParser(description="Get canonical isoform protein features per WBGene (Ensembl canonical transcript/protein)")
    ap.add_argument("--input", required=True, help="CSV with WBGene IDs (column auto-detected)")
    ap.add_argument("--output", required=True, help="Output CSV file path")
    ap.add_argument("--gene-col", default=None, help="Column name with WBGene IDs (optional)")
    ap.add_argument("--flush-every", type=int, default=10, help="Write results to CSV every N genes (default: 10)")
    ap.add_argument("--no-resume", action="store_true", help="Do not resume; overwrite output")
    args = ap.parse_args()

    out_path = args.output
    # Ensure .csv extension
    if not out_path.endswith('.csv'):
        out_path = f"{out_path}.csv"
    
    # Ensure output directory exists
    parent = os.path.dirname(os.path.abspath(out_path))
    if parent and not os.path.exists(parent):
        os.makedirs(parent, exist_ok=True)
    
    # Initialize output file with header if missing
    if not os.path.exists(out_path):
        pd.DataFrame(columns=_FIELDNAMES).to_csv(out_path, index=False)
        print(f"[init] created new output file with header: {os.path.abspath(out_path)}", file=sys.stderr)
        done = set()
    elif args.no_resume:
        print(f"[resume] --no-resume set: overwriting existing {out_path}", file=sys.stderr)
        open(out_path, "w").close()
        pd.DataFrame(columns=_FIELDNAMES).to_csv(out_path, index=False)
        done = set()
    else:
        done = _existing_done_set(out_path)
        if done:
            print(f"[resume] detected {len(done)} completed gene(s); will skip those.", file=sys.stderr)

    _dedupe_file(out_path)  # one-time deduplication on startup

    print("[start] reading gene list …", file=sys.stderr)
    df_in = pd.read_csv(args.input)
    gene_col = args.gene_col or autodetect_gene_col(df_in)
    wbgenes: List[str] = []
    for v in df_in[gene_col].astype(str):
        m = re.search(r"(WBGene\d{8})", v)
        if m: wbgenes.append(m.group(1))
        else: print(f"[warn] skipping value without WBGene: {v}", file=sys.stderr)
    if not wbgenes:
        print("[error] no valid WBGene IDs found", file=sys.stderr)
        sys.exit(2)
    print(f"[start] total genes: {len(wbgenes)}", file=sys.stderr)

    sess = make_session()
    buffer = []
    flush_every = max(1, int(args.flush_every))
    for i, wbg in enumerate(wbgenes, 1):
        if wbg in done:
            print(f"=== [skip] already have {wbg} ===", file=sys.stderr)
            continue
        print(f"\n=== [{i}/{len(wbgenes)}] {wbg} ===", file=sys.stderr)
        egid = resolve_gene_id(sess, wbg) or wbg
        result = pick_canonical_transcript_and_translation(sess, egid)
        if not result:
            print("  [warn] no translated transcript found", file=sys.stderr)
            buffer.append({
                "WBGene": wbg,
                "canonical_transcript": None,
                "uniprot": "N/A",
                "match_type": "No translated transcript found",
                "protein_len": None,
                "CAAX": "",
                "PTS1": "",
                "PTS2": "",
                "ER_KDEL_motif": "",
                "ER_diLys_motif": "",
                "TM_regions": "",
                "Motifs": "",
                "Lipidation_sites": "",
                "iMET_events": "",
                "SP_regions": "",
                "Propeptides": "",
                "Transit_peptides": "",
                "Chain_regions": "",
                "features": "NONE",
                "notes": "No translated transcript found",
            })
            if len(buffer) >= flush_every:
                done.update([r["WBGene"] for r in buffer])
                _write_rows(out_path, buffer)
                buffer.clear()
                _dedupe_file(out_path)
            continue

        tx_id, tr_id, aa_len, up, match_type = result
        print(f"  [chosen] transcript={tx_id} protein={tr_id} aa={aa_len} uniprot={up} type={match_type}", file=sys.stderr)

        # protein sequence (for CAAX)
        seq = get_text(sess, f"/sequence/id/{tr_id}", {"type": "protein"}) or ""
        print(f"  [seq] got_protein_seq={bool(seq)} len={len(seq)}", file=sys.stderr)
        has_c, caax = has_caax(seq)
        # Added: PTS/ER motifs from primary sequence
        pts1_ok, pts1_tri = detect_pts1(seq)
        pts2_ok, pts2_nonapeptide, pts2_pos = detect_pts2(seq)
        kdel_ok, kdel_tail = detect_er_kdel(seq)
        kk_ok, kk_tail = detect_er_di_lysine(seq)
        if pts1_ok: print(f"  [PTS1] motif={pts1_tri}", file=sys.stderr)
        if pts2_ok: print(f"  [PTS2] {pts2_nonapeptide} @ {pts2_pos}", file=sys.stderr)
        if kdel_ok: print(f"  [KDEL] tail={kdel_tail}", file=sys.stderr)
        if kk_ok:   print(f"  [KKxx] tail={kk_tail}", file=sys.stderr)
        if has_c:   print(f"  [CAAX] motif={caax}", file=sys.stderr)

        print(f"  [uniprot] accession={up} | match_type={match_type}", file=sys.stderr)

        # features from UniProt
        feats = fetch_uniprot_features_rich(up, sess)
        sp_regions = feats.get("SP", [])
        tm_regions = [re.sub(r"\s*\(.*?\)\s*$", "", x) for x in feats.get("TM", [])]  # strip "(Helical)"
        lip = feats.get("LIP", [])
        motif = feats.get("MOTIF", [])
        chain = feats.get("CHAIN", [])
        imet = feats.get("iMET", [])
        propeptide = feats.get("PRO", [])
        transit = feats.get("TP", [])

        # lipidation summary
        gpi  = any(re.search(r"\bGPI\b|glycosylphosphatidylinositol", x, re.I) for x in lip)
        myr  = any(re.search(r"myrist", x, re.I) for x in lip)
        palm = any(re.search(r"palmit", x, re.I) for x in lip)
        
        # processing events
        met_removed = any(re.search(r"removed", x, re.I) for x in imet)
        has_chain_processing = len(chain) > 0
        
        print(f"  [features] SP={len(sp_regions)} TM={len(tm_regions)} LIP={len(lip)} CHAIN={len(chain)} iMET={len(imet)}, PRO={len(propeptide)}, TP={len(transit)}", file=sys.stderr)
        print(f"  [features] GPI={gpi}, N-myristoyl={myr}, Palmitoyl={palm}, Met-removed={met_removed}, Chain-processing={has_chain_processing}", file=sys.stderr)

        tags = []
        if sp_regions: tags.append("SP")
        if tm_regions: tags.append(f"TMx{len(tm_regions)}")
        if gpi:        tags.append("GPI")
        if myr:        tags.append("N-myristoyl")
        if palm:       tags.append("Palmitoyl")
        if has_c:      tags.append("CAAX")
        if met_removed: tags.append("Met-removed")
        if propeptide: tags.append("Propeptide")
        if transit: tags.append("Transit Peptide")
        if has_chain_processing: tags.append("Chain-processed")
        feat_tag = ";".join(tags) if tags else "NONE"
        print(f"  [summary] tag={feat_tag}", file=sys.stderr)

        buffer.append({
            "WBGene": wbg,
            "canonical_transcript": tx_id,
            "uniprot": up,
            "match_type": match_type,
            "protein_len": aa_len,
            "CAAX": caax,
            "PTS1": (pts1_tri if pts1_ok else ""),
            "PTS2": (f"{pts2_nonapeptide}@{pts2_pos}" if pts2_ok else ""),
            "ER_KDEL_motif": (kdel_tail if kdel_ok else ""),
            "ER_diLys_motif": (kk_tail if kk_ok else ""),
            "TM_regions": ",".join(tm_regions),
            "Motifs": ",".join(motif),
            "Lipidation_sites": ",".join(lip),
            "Transit_peptides": ",".join(transit),
            "iMET_events": ",".join(imet),
            "SP_regions": ",".join(sp_regions),
            "Propeptides": ",".join(propeptide),
            "Chain_regions": ",".join(chain),
            "features": feat_tag,
            "notes": "",
        })

        if len(buffer) >= flush_every:
            done.update([r["WBGene"] for r in buffer])
            _write_rows(out_path, buffer)
            buffer.clear()
            _dedupe_file(out_path)
    
    # Final flush
    if buffer:
        _write_rows(out_path, buffer)
        done.update([r["WBGene"] for r in buffer])
        buffer.clear()
    
    print("[done]", file=sys.stderr)

    elapsed = time.perf_counter() - t0
    rate = (len(wbgenes) / elapsed) if elapsed > 0 else float('inf')
    print(f"[timing] Elapsed {elapsed:.2f}s | {len(wbgenes)} genes | {rate:.2f} genes/s")

if __name__ == "__main__":
    main()