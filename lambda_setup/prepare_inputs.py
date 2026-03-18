#!/usr/bin/env python3
"""
Converts test_sequences.csv → Boltz YAML input files (one per target).
Also writes chain_map.json used by convert_submission.py.

Handles:
  - Pure RNA (single or multi-chain / homo-multimers)
  - RNA-protein complexes: protein chains → 'protein' type in YAML
  - RNA-protein-DNA complexes: DNA chains → 'dna' type in YAML
  - Ligands from SMILES column

The submission only needs C1' coordinates for RNA residues.
Protein and DNA chains are included in the Boltz YAML for better
context-aware RNA structure prediction but excluded from chain_map.json.
"""

import csv
import json
import re
import sys
from pathlib import Path

import yaml  # PyYAML


# --------------------------------------------------------------------------- #
# Mol-type detection
# --------------------------------------------------------------------------- #

PROTEIN_EXCLUSIVE_CHARS = set("MDEFHIKLPQRSVWY")


def classify_sequence(header: str, original_seq: str) -> str:
    """
    Returns 'rna', 'dna', or 'protein'.

    Priority:
      1. 'DNA' in FASTA header → dna
      2. Protein-exclusive amino acid letters present → protein
      3. Has T but no U in sequence → dna (original T = thymine)
      4. Otherwise → rna
    """
    seq_upper = original_seq.upper()
    header_upper = header.upper()

    if "DNA" in header_upper:
        return "dna"

    if any(c in seq_upper for c in PROTEIN_EXCLUSIVE_CHARS):
        return "protein"

    has_T = "T" in seq_upper
    has_U = "U" in seq_upper

    if has_T and not has_U:
        return "dna"

    return "rna"


# --------------------------------------------------------------------------- #
# FASTA parsing
# --------------------------------------------------------------------------- #

def parse_all_sequences(all_seqs_str: str) -> list[tuple[str, str, int, str]]:
    """
    Parse the `all_sequences` FASTA field.

    Returns: [(header, sequence, n_copies, mol_type), ...]
        mol_type: 'rna' | 'dna' | 'protein'
        sequence: for RNA → T replaced with U; for DNA/protein → original case
        n_copies: number of auth chain labels in header
    """
    entries: list[tuple[str, str, int, str]] = []
    current_header: str | None = None
    current_seq_parts: list[str] = []

    def flush():
        if current_header is None or not current_seq_parts:
            return
        raw_seq = "".join(current_seq_parts).upper()
        mol_type = classify_sequence(current_header, raw_seq)

        if mol_type == "rna":
            seq = raw_seq.replace("T", "U")
        else:
            seq = raw_seq  # keep DNA with T; protein as-is

        auth_labels = re.findall(r"\[auth\s+(\S+?)\]", current_header)
        n_copies = len(auth_labels) if auth_labels else 1
        entries.append((current_header, seq, n_copies, mol_type))

    for line in (all_seqs_str or "").split("\n"):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            current_header = line
            current_seq_parts = []
        else:
            current_seq_parts.append(line)

    flush()
    return entries


def expand_entries(entries: list[tuple[str, str, int, str]]) -> tuple[list, list, list]:
    """
    Expand FASTA entries to individual chains with Boltz chain IDs.

    RNA chains come first (A, B, C, …), then protein, then DNA.
    This ordering ensures the Boltz output chains A, B, … map directly to
    the RNA residues in chain_map.json.

    Returns:
        rna_chains:     [(boltz_chain_id, seq), …]
        protein_chains: [(boltz_chain_id, seq), …]
        dna_chains:     [(boltz_chain_id, seq), …]
    """
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"

    rna_entries = [(h, s, n) for h, s, n, t in entries if t == "rna"]
    prot_entries = [(h, s, n) for h, s, n, t in entries if t == "protein"]
    dna_entries = [(h, s, n) for h, s, n, t in entries if t == "dna"]

    rna_chains, protein_chains, dna_chains = [], [], []
    idx = 0

    for _, seq, n in rna_entries:
        for _ in range(n):
            rna_chains.append((alphabet[idx], seq))
            idx += 1

    for _, seq, n in prot_entries:
        for _ in range(n):
            protein_chains.append((alphabet[idx], seq))
            idx += 1

    for _, seq, n in dna_entries:
        for _ in range(n):
            dna_chains.append((alphabet[idx], seq))
            idx += 1

    return rna_chains, protein_chains, dna_chains


# --------------------------------------------------------------------------- #
# Ligand parsing
# --------------------------------------------------------------------------- #

def parse_ligands(ligand_ids_str: str, ligand_smiles_str: str) -> list[dict]:
    if not ligand_ids_str or not ligand_smiles_str:
        return []
    ids = [x.strip() for x in ligand_ids_str.split(";")]
    smiles_list = [x.strip() for x in ligand_smiles_str.split(";")]
    return [
        {"id": lid, "smiles": smi}
        for lid, smi in zip(ids, smiles_list)
        if smi and smi not in (".", "N/A", "")
    ]


# --------------------------------------------------------------------------- #
# YAML builder
# --------------------------------------------------------------------------- #

def build_boltz_yaml(rna_chains: list[tuple[str, str]],
                     protein_chains: list[tuple[str, str]],
                     dna_chains: list[tuple[str, str]],
                     ligands: list[dict]) -> dict:
    sequences = []
    for cid, seq in rna_chains:
        sequences.append({"rna": {"id": cid, "sequence": seq}})
    for cid, seq in protein_chains:
        sequences.append({"protein": {"id": cid, "sequence": seq}})
    for cid, seq in dna_chains:
        sequences.append({"dna": {"id": cid, "sequence": seq}})

    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
    lig_start = len(rna_chains) + len(protein_chains) + len(dna_chains)
    for i, lig in enumerate(ligands):
        sequences.append({"ligand": {"id": alphabet[lig_start + i], "smiles": lig["smiles"]}})

    return {"version": 1, "sequences": sequences}


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main():
    base_dir = Path(__file__).parent
    test_csv = base_dir / "test_sequences.csv"
    inputs_dir = base_dir / "inputs"
    inputs_dir.mkdir(exist_ok=True)

    if not test_csv.exists():
        print(f"ERROR: {test_csv} not found. Upload it first.", file=sys.stderr)
        sys.exit(1)

    with open(test_csv) as f:
        rows = list(csv.DictReader(f))

    print(f"Found {len(rows)} test targets.")
    chain_map: dict[str, list[dict]] = {}

    for row in rows:
        target_id = row["target_id"]
        full_rna_seq = row["sequence"].upper().replace("T", "U")

        entries = parse_all_sequences(row.get("all_sequences", ""))
        rna_chains, protein_chains, dna_chains = expand_entries(entries)

        # Validate: concatenated RNA chains must match the sequence column
        concat_rna = "".join(seq for _, seq in rna_chains)
        if concat_rna != full_rna_seq:
            print(f"  [WARN] {target_id}: RNA concat mismatch "
                  f"({len(concat_rna)} vs {len(full_rna_seq)} nt) — "
                  f"falling back to full sequence as single RNA chain.")
            rna_chains = [("A", full_rna_seq)]
            protein_chains = []
            dna_chains = []

        ligands = parse_ligands(
            row.get("ligand_ids", ""),
            row.get("ligand_SMILES", ""),
        )

        boltz_input = build_boltz_yaml(rna_chains, protein_chains, dna_chains, ligands)
        yaml_path = inputs_dir / f"{target_id}.yaml"
        with open(yaml_path, "w") as f:
            yaml.dump(boltz_input, f, default_flow_style=False, allow_unicode=True)

        # chain_map = RNA chains only (used by convert_submission.py)
        chain_map[target_id] = [
            {"chain_id": cid, "seq_len": len(seq)} for cid, seq in rna_chains
        ]

        print(
            f"  {target_id}: "
            f"{len(rna_chains)}×RNA  "
            f"{len(protein_chains)}×prot  "
            f"{len(dna_chains)}×DNA  "
            f"{len(full_rna_seq)}nt  "
            f"{len(ligands)} ligands"
        )

    chain_map_path = base_dir / "chain_map.json"
    with open(chain_map_path, "w") as f:
        json.dump(chain_map, f, indent=2)

    print(f"\nchain_map.json → {chain_map_path}")
    print(f"YAML inputs    → {inputs_dir}/")
    print("\nNext: bash run_boltz.sh")


if __name__ == "__main__":
    main()
