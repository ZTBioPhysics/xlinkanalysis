"""
Domain annotation and analysis functions.

This module handles mapping residues to protein domains and
computing domain-level statistics for crosslinks.
"""

import pandas as pd
from typing import Dict, List, Optional, Tuple, Set
from .config import load_domains, get_domain_ranges


def find_domain(residue: int, domain_ranges: Dict[str, List[List[int]]]) -> str:
    """
    Determine which domain a residue belongs to.

    Args:
        residue: Residue number (1-based).
        domain_ranges: Dictionary mapping domain names to their ranges.

    Returns:
        Domain name, or "Unknown" if residue is not in any defined domain.

    Example:
        >>> ranges = {'NTD': [[1, 100]], 'CTD': [[101, 200]]}
        >>> find_domain(50, ranges)
        'NTD'
        >>> find_domain(150, ranges)
        'CTD'
    """
    for domain, ranges in domain_ranges.items():
        for start, end in ranges:
            if start <= residue <= end:
                return domain
    return "Unknown"


def annotate_domains(df: pd.DataFrame,
                     domains_config: Optional[str] = None,
                     domain_ranges: Optional[Dict] = None) -> pd.DataFrame:
    """
    Add domain association annotations to crosslink data.

    Adds a 'Domain Association' column describing which domains each
    crosslink connects (e.g., "NTD to CTD").

    Args:
        df: DataFrame with 'Residue1' and 'Residue2' columns.
        domains_config: Path to domains YAML config file.
        domain_ranges: Alternatively, provide domain ranges directly.
                      Must provide either domains_config or domain_ranges.

    Returns:
        DataFrame with added 'Domain Association' column.

    Example:
        >>> df = annotate_domains(df, "config/domains.yaml")
        >>> df['Domain Association'].head()
        0    NTD to NTD
        1    NTD to CTD
        2    CTD to CTD
    """
    df = df.copy()

    # Get domain ranges from config or direct input
    if domain_ranges is None:
        if domains_config is None:
            raise ValueError("Must provide either domains_config or domain_ranges")
        domains = load_domains(domains_config)
        domain_ranges = get_domain_ranges(domains)

    # Find domain for each residue
    df['_Domain1'] = df['Residue1'].apply(lambda x: find_domain(x, domain_ranges))
    df['_Domain2'] = df['Residue2'].apply(lambda x: find_domain(x, domain_ranges))
    df['Domain Association'] = df['_Domain1'] + " to " + df['_Domain2']

    # Clean up temporary columns
    df.drop(['_Domain1', '_Domain2'], axis=1, inplace=True)

    return df


def get_domain_residues(domain_ranges: Dict[str, List[List[int]]]) -> Dict[str, Set[int]]:
    """
    Get all residue numbers belonging to each domain.

    Args:
        domain_ranges: Dictionary mapping domain names to their ranges.

    Returns:
        Dictionary mapping domain names to sets of residue numbers.
    """
    domain_residues = {}
    for domain, ranges in domain_ranges.items():
        residues = set()
        for start, end in ranges:
            residues.update(range(start, end + 1))
        domain_residues[domain] = residues
    return domain_residues


def domain_crosslink_stats(df: pd.DataFrame,
                           sequence: str,
                           domains_config: Optional[str] = None,
                           domain_ranges: Optional[Dict] = None) -> pd.DataFrame:
    """
    Calculate crosslink statistics for each protein domain.

    Args:
        df: DataFrame with crosslink data (Residue1, Residue2).
        sequence: Protein sequence (single-letter codes).
        domains_config: Path to domains YAML config file.
        domain_ranges: Alternatively, provide domain ranges directly.

    Returns:
        DataFrame with columns:
        - Domain
        - # of Residues
        - # of Lysines
        - Crosslinked Lysines
        - Unique Crosslinks
        - % Lysines Crosslinked
        - Crosslinks per Lysine

    Example:
        >>> stats = domain_crosslink_stats(df, sequence, "config/domains.yaml")
        >>> print(stats)
    """
    # Get domain ranges
    if domain_ranges is None:
        if domains_config is None:
            raise ValueError("Must provide either domains_config or domain_ranges")
        domains = load_domains(domains_config)
        domain_ranges = get_domain_ranges(domains)

    # Find all lysines in sequence
    lysine_residues = {i + 1 for i, aa in enumerate(sequence) if aa == "K"}

    results = []

    for domain, ranges in domain_ranges.items():
        # Get all residues in this domain
        domain_residues = set()
        for start, end in ranges:
            domain_residues.update(range(start, end + 1))

        # Lysines in this domain
        lysines_in_domain = domain_residues.intersection(lysine_residues)

        # Crosslinks involving this domain
        domain_crosslinks = df[
            df["Residue1"].isin(domain_residues) |
            df["Residue2"].isin(domain_residues)
        ]

        # Unique residue pairs
        unique_pairs = {
            tuple(sorted([row["Residue1"], row["Residue2"]]))
            for _, row in domain_crosslinks.iterrows()
        }

        # Lysines that participate in crosslinks
        lysines_in_crosslinks = {
            residue
            for pair in unique_pairs
            for residue in pair
            if residue in lysines_in_domain
        }

        # Calculate statistics
        total_residues = len(domain_residues)
        total_lysines = len(lysines_in_domain)
        crosslinked_lysines = len(lysines_in_crosslinks)
        unique_xlinks = len(unique_pairs)

        pct_lysines = (crosslinked_lysines / total_lysines * 100) if total_lysines > 0 else 0
        xlinks_per_lys = (unique_xlinks / crosslinked_lysines) if crosslinked_lysines > 0 else 0

        results.append({
            "Domain": domain,
            "# of Residues": total_residues,
            "# of Lysines": total_lysines,
            "Crosslinked Lysines": crosslinked_lysines,
            "Unique Crosslinks": unique_xlinks,
            "% Lysines Crosslinked": round(pct_lysines, 1),
            "Crosslinks per Lysine": round(xlinks_per_lys, 2)
        })

    return pd.DataFrame(results)


def summarize_by_domain_interaction(df: pd.DataFrame) -> pd.DataFrame:
    """
    Summarize crosslink statistics grouped by domain interaction type.

    Args:
        df: DataFrame with 'Domain Association', 'CA Distance', 'Spectral Count' columns.

    Returns:
        DataFrame with columns:
        - Domain Interaction
        - # of Unique Crosslinks
        - % C-alpha Distance <= 20 A
        - % C-alpha Distance <= 50 A
        - % Spectral Count >= 5
        - % Spectral Count >= 10
    """
    df = df.copy()

    # Add canonical pair for deduplication
    df["Canonical Pair"] = df.apply(
        lambda x: tuple(sorted([x["Residue1"], x["Residue2"]])), axis=1
    )

    # Group and aggregate
    summary_df = df.groupby("Domain Association").agg(
        **{
            "# of Unique Crosslinks": ("Canonical Pair", "nunique"),
            "% CA Distance <= 20": ("CA Distance", lambda x: (x <= 20).mean() * 100),
            "% CA Distance <= 50": ("CA Distance", lambda x: (x <= 50).mean() * 100),
            "% Spectral Count >= 5": ("Spectral Count", lambda x: (x >= 5).mean() * 100),
            "% Spectral Count >= 10": ("Spectral Count", lambda x: (x >= 10).mean() * 100),
        }
    )

    summary_df = summary_df.reset_index().rename(
        columns={"Domain Association": "Domain Interaction"}
    )

    # Round percentages
    for col in summary_df.columns:
        if col.startswith('%'):
            summary_df[col] = summary_df[col].round(1)

    return summary_df


def is_inter_domain(domain_association: str) -> bool:
    """
    Check if a domain association represents an inter-domain crosslink.

    Args:
        domain_association: String like "NTD to CTD" or "NTD to NTD".

    Returns:
        True if crosslink is between different domains.
    """
    domains = domain_association.split(" to ")
    return len(set(domains)) > 1

