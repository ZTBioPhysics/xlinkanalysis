"""
Input/Output functions for cross-link data.

This module handles reading and writing crosslink data files,
including CSV files and PDB structures.
"""

import pandas as pd
from pathlib import Path
from typing import Optional, List
from Bio.PDB import PDBParser


def load_crosslinks(filepath: str, has_header: bool = False) -> pd.DataFrame:
    """
    Load crosslink data from a CSV file.

    Args:
        filepath: Path to the CSV file.
        has_header: If True, first row is treated as column headers.
                   If False, columns are named Residue1, Residue2, Spectral Count.

    Returns:
        DataFrame with columns: Residue1, Residue2, Spectral Count

    Example:
        >>> df = load_crosslinks("crosslinks.csv")
        >>> df.head()
           Residue1  Residue2  Spectral Count
        0       100       150              12
        1       200       250               8
    """
    if has_header:
        df = pd.read_csv(filepath)
    else:
        df = pd.read_csv(filepath, header=None)
        df.columns = ['Residue1', 'Residue2', 'Spectral Count']

    return df


def load_annotated_crosslinks(filepath: str) -> pd.DataFrame:
    """
    Load a previously annotated crosslink CSV file (with headers).

    Args:
        filepath: Path to the annotated CSV file.

    Returns:
        DataFrame with all columns from the file.
    """
    return pd.read_csv(filepath)


def save_crosslinks(df: pd.DataFrame, filepath: str, index: bool = False) -> None:
    """
    Save crosslink data to a CSV file.

    Args:
        df: DataFrame to save.
        filepath: Output file path.
        index: Whether to include the DataFrame index.
    """
    df.to_csv(filepath, index=index)
    print(f"Saved {len(df)} crosslinks to {filepath}")


def get_lysines_from_pdb(pdb_file: str) -> dict:
    """
    Extract lysine residue positions from a PDB file.

    Args:
        pdb_file: Path to the PDB file.

    Returns:
        Dictionary mapping residue numbers to 'LYS'.

    Example:
        >>> lysines = get_lysines_from_pdb("structure.pdb")
        >>> 100 in lysines
        True
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    lysine_residues = {}

    for model in structure:
        for chain in model:
            for residue in chain.get_residues():
                if residue.get_resname() == 'LYS':
                    res_num = residue.get_full_id()[3][1]
                    lysine_residues[res_num] = 'LYS'

    return lysine_residues


def get_lysines_from_sequence(sequence: str) -> set:
    """
    Find all lysine positions in a protein sequence.

    Args:
        sequence: Protein sequence (single-letter amino acid codes).

    Returns:
        Set of 1-based residue positions that are lysines.

    Example:
        >>> lysines = get_lysines_from_sequence("MKAKDEF")
        >>> lysines
        {2, 4}
    """
    return {i + 1 for i, aa in enumerate(sequence) if aa == 'K'}


def load_sequence_from_fasta(filepath: str, entry_index: int = 0) -> str:
    """
    Load a protein sequence from a FASTA file.

    Args:
        filepath: Path to the FASTA file.
        entry_index: Which sequence to return if file has multiple entries.

    Returns:
        Protein sequence as a string.
    """
    sequences = []
    current_seq = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    sequences.append(''.join(current_seq))
                    current_seq = []
            else:
                current_seq.append(line)

        if current_seq:
            sequences.append(''.join(current_seq))

    if entry_index >= len(sequences):
        raise IndexError(f"Entry index {entry_index} out of range. "
                        f"File contains {len(sequences)} sequences.")

    return sequences[entry_index]
