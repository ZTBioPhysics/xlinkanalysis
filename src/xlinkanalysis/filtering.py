"""
Filtering functions for crosslink data.

This module provides various filtering operations for crosslink DataFrames,
including filtering by spectral count, distance, residue, and domain.
"""

import pandas as pd
from typing import Optional, List


def by_spectral_count(df: pd.DataFrame, min_count: int) -> pd.DataFrame:
    """
    Filter crosslinks by minimum spectral count.

    Args:
        df: DataFrame with 'Spectral Count' column.
        min_count: Minimum spectral count threshold (inclusive).

    Returns:
        Filtered DataFrame.

    Example:
        >>> df_filtered = by_spectral_count(df, min_count=5)
    """
    if 'Spectral Count' not in df.columns:
        raise ValueError("DataFrame must include a 'Spectral Count' column.")

    return df[df['Spectral Count'] >= min_count].copy()


def by_ca_distance(df: pd.DataFrame,
                   threshold: float,
                   comparison: str = 'less') -> pd.DataFrame:
    """
    Filter crosslinks by C-alpha distance.

    Args:
        df: DataFrame with 'CA Distance' column.
        threshold: Distance threshold in Angstroms.
        comparison: 'less' for <= threshold, 'greater' for >= threshold.

    Returns:
        Filtered DataFrame.

    Example:
        >>> df_close = by_ca_distance(df, 20, comparison='less')
        >>> df_far = by_ca_distance(df, 50, comparison='greater')
    """
    if 'CA Distance' not in df.columns:
        raise ValueError("DataFrame must include a 'CA Distance' column.")

    if comparison == 'less':
        return df[df['CA Distance'] <= threshold].copy()
    elif comparison == 'greater':
        return df[df['CA Distance'] >= threshold].copy()
    else:
        raise ValueError("comparison must be 'less' or 'greater'.")


def by_sequence_distance(df: pd.DataFrame,
                         threshold: int,
                         comparison: str = 'less') -> pd.DataFrame:
    """
    Filter crosslinks by sequence distance (residue separation).

    Args:
        df: DataFrame with 'Sequence Distance' column.
        threshold: Sequence distance threshold.
        comparison: 'less' for <= threshold, 'greater' for >= threshold.

    Returns:
        Filtered DataFrame.
    """
    if 'Sequence Distance' not in df.columns:
        raise ValueError("DataFrame must include a 'Sequence Distance' column.")

    if comparison == 'less':
        return df[df['Sequence Distance'] <= threshold].copy()
    elif comparison == 'greater':
        return df[df['Sequence Distance'] >= threshold].copy()
    else:
        raise ValueError("comparison must be 'less' or 'greater'.")


def by_residue(df: pd.DataFrame, residue: int) -> pd.DataFrame:
    """
    Filter to crosslinks involving a specific residue.

    Returns unique crosslinks (deduplicates pairs like 100-150 and 150-100).

    Args:
        df: DataFrame with 'Residue1' and 'Residue2' columns.
        residue: Residue number to filter for.

    Returns:
        Filtered DataFrame with unique crosslinks involving the residue.

    Example:
        >>> df_k100 = by_residue(df, 100)
    """
    if 'Residue1' not in df.columns or 'Residue2' not in df.columns:
        raise ValueError("DataFrame must include 'Residue1' and 'Residue2' columns.")

    df = df.copy()

    # Create canonical pair for deduplication
    df['_Canonical'] = df.apply(
        lambda x: tuple(sorted([x['Residue1'], x['Residue2']])), axis=1
    )

    # Deduplicate
    df_unique = df.drop_duplicates(subset='_Canonical')

    # Filter for residue
    filtered = df_unique[
        (df_unique['Residue1'] == residue) | (df_unique['Residue2'] == residue)
    ]

    # Clean up
    return filtered.drop(columns=['_Canonical']).copy()


def by_domain(df: pd.DataFrame,
              domain: Optional[str] = None,
              exclude_intra: bool = False,
              exclude_domains: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Filter crosslinks by domain association.

    Args:
        df: DataFrame with 'Domain Association' column.
        domain: Specific domain or domain pair to include.
               Examples: "NTD", "NTD to CTD"
               If None, includes all domains (subject to other filters).
        exclude_intra: If True, exclude intra-domain crosslinks.
        exclude_domains: List of domains to exclude from results.

    Returns:
        Filtered DataFrame.

    Examples:
        >>> # Get only inter-domain crosslinks
        >>> df_inter = by_domain(df, exclude_intra=True)

        >>> # Get crosslinks involving NTD
        >>> df_ntd = by_domain(df, domain="NTD")

        >>> # Get NTD to CTD crosslinks specifically
        >>> df_pair = by_domain(df, domain="NTD to CTD")
    """
    if "Domain Association" not in df.columns:
        raise ValueError("DataFrame must include a 'Domain Association' column.")

    def passes_filter(association: str) -> bool:
        domains = association.split(" to ")

        # Check intra-domain exclusion
        if exclude_intra and len(set(domains)) == 1:
            return False

        # Check excluded domains
        if exclude_domains:
            if any(d in exclude_domains for d in domains):
                return False

        # Check specific domain filter
        if domain:
            if " to " in domain:
                # Looking for specific pair (bidirectional)
                parts = domain.split(" to ")
                inverse = f"{parts[1]} to {parts[0]}"
                if association != domain and association != inverse:
                    return False
            else:
                # Looking for any crosslink involving this domain
                if domain not in domains:
                    return False

        return True

    return df[df['Domain Association'].apply(passes_filter)].copy()


def satisfied_crosslinks(df: pd.DataFrame,
                         max_distance: float = 20.0,
                         min_spectral_count: int = 1) -> pd.DataFrame:
    """
    Filter to "satisfied" crosslinks within expected distance constraints.

    Args:
        df: DataFrame with 'CA Distance' and 'Spectral Count' columns.
        max_distance: Maximum C-alpha distance to consider satisfied.
        min_spectral_count: Minimum spectral count.

    Returns:
        Filtered DataFrame of satisfied crosslinks.
    """
    return df[
        (df['CA Distance'] <= max_distance) &
        (df['Spectral Count'] >= min_spectral_count)
    ].copy()


def violated_crosslinks(df: pd.DataFrame,
                        min_distance: float = 20.0,
                        min_spectral_count: int = 5) -> pd.DataFrame:
    """
    Filter to "violated" crosslinks exceeding expected distance constraints.

    These are high-confidence crosslinks (high spectral count) that exceed
    the crosslinker's maximum span, potentially indicating conformational
    changes or dynamics.

    Args:
        df: DataFrame with 'CA Distance' and 'Spectral Count' columns.
        min_distance: Minimum C-alpha distance to consider violated.
        min_spectral_count: Minimum spectral count to ensure high confidence.

    Returns:
        Filtered DataFrame of violated crosslinks.
    """
    return df[
        (df['CA Distance'] > min_distance) &
        (df['Spectral Count'] >= min_spectral_count)
    ].copy()

