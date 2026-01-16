"""
Distance calculation functions for cross-link analysis.

This module provides functions for calculating C-alpha distances from
PDB structures and sequence distances between residue pairs.
"""

import numpy as np
import pandas as pd
from typing import Dict, Set, Optional


def get_ca_coordinates(pdb_file: str, residue_numbers: Set[int]) -> Dict[int, np.ndarray]:
    """
    Extract C-alpha coordinates from a PDB file for specified residues.

    Args:
        pdb_file: Path to the PDB file.
        residue_numbers: Set of residue numbers to extract coordinates for.

    Returns:
        Dictionary mapping residue numbers to their CA coordinates (x, y, z).

    Example:
        >>> coords = get_ca_coordinates("structure.pdb", {100, 150, 200})
        >>> coords[100]
        array([10.5, 20.3, 15.2])
    """
    coordinates = {}

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                residue_number = int(line[22:26].strip())
                if residue_number in residue_numbers:
                    atom_name = line[12:16].strip()
                    if atom_name == 'CA':
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        coordinates[residue_number] = np.array([x, y, z])

    return coordinates


def get_atom_coordinates(pdb_file: str, residue_numbers: Set[int],
                         atom_type: str = 'CA') -> Dict[int, np.ndarray]:
    """
    Extract coordinates for a specific atom type from a PDB file.

    Args:
        pdb_file: Path to the PDB file.
        residue_numbers: Set of residue numbers to extract coordinates for.
        atom_type: Atom name to extract (e.g., 'CA', 'CB', 'NZ').

    Returns:
        Dictionary mapping residue numbers to coordinates.
    """
    coordinates = {}

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                residue_number = int(line[22:26].strip())
                if residue_number in residue_numbers:
                    atom_name = line[12:16].strip()
                    if atom_name == atom_type:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        coordinates[residue_number] = np.array([x, y, z])

    return coordinates


def euclidean_distance(coord1: np.ndarray, coord2: np.ndarray) -> float:
    """
    Calculate Euclidean distance between two 3D points.

    Args:
        coord1: First coordinate (x, y, z).
        coord2: Second coordinate (x, y, z).

    Returns:
        Distance in the same units as input coordinates (typically Angstroms).
    """
    return float(np.linalg.norm(coord1 - coord2))


def calculate_ca_distances(pdb_file: str, df: pd.DataFrame,
                           precision: int = 1) -> pd.DataFrame:
    """
    Calculate C-alpha distances for all crosslinked residue pairs.

    Adds a 'CA Distance' column to the DataFrame.

    Args:
        pdb_file: Path to the PDB file.
        df: DataFrame with 'Residue1' and 'Residue2' columns.
        precision: Number of decimal places for rounding distances.

    Returns:
        DataFrame with added 'CA Distance' column.

    Example:
        >>> df = load_crosslinks("crosslinks.csv")
        >>> df = calculate_ca_distances("structure.pdb", df)
        >>> df['CA Distance'].head()
        0    15.2
        1    22.8
        2    18.5
    """
    df = df.copy()

    # Get all unique residues
    unique_residues = set(df['Residue1']).union(set(df['Residue2']))

    # Extract CA coordinates
    ca_coords = get_ca_coordinates(pdb_file, unique_residues)

    # Calculate distances
    distances = []
    for _, row in df.iterrows():
        res1, res2 = row['Residue1'], row['Residue2']

        if res1 in ca_coords and res2 in ca_coords:
            dist = euclidean_distance(ca_coords[res1], ca_coords[res2])
            distances.append(round(dist, precision))
        else:
            distances.append(None)

    df['CA Distance'] = distances
    return df


def add_sequence_distance(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add sequence distance (residue separation) to the DataFrame.

    Sequence distance is the absolute difference between residue numbers.

    Args:
        df: DataFrame with 'Residue1' and 'Residue2' columns.

    Returns:
        DataFrame with added 'Sequence Distance' column.

    Example:
        >>> df['Residue1'] = [100, 200]
        >>> df['Residue2'] = [150, 210]
        >>> df = add_sequence_distance(df)
        >>> df['Sequence Distance'].tolist()
        [50, 10]
    """
    df = df.copy()

    if 'Residue1' not in df.columns or 'Residue2' not in df.columns:
        raise ValueError("DataFrame must include 'Residue1' and 'Residue2' columns.")

    df['Sequence Distance'] = (df['Residue1'] - df['Residue2']).abs()
    return df


def calculate_correlation(df: pd.DataFrame,
                          col1: str = 'CA Distance',
                          col2: str = 'Spectral Count') -> float:
    """
    Calculate Pearson correlation coefficient between two columns.

    Args:
        df: DataFrame containing the columns.
        col1: First column name.
        col2: Second column name.

    Returns:
        Pearson correlation coefficient.

    Raises:
        ValueError: If required columns are missing.
    """
    if col1 not in df.columns or col2 not in df.columns:
        raise ValueError(f"Required columns ('{col1}' and '{col2}') not found.")

    return df[col1].corr(df[col2])
