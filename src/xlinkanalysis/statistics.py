"""
Statistical analysis functions for crosslink data.

This module provides functions for computing statistics, comparing datasets,
and generating summary reports.
"""

import pandas as pd
import numpy as np
from typing import Tuple, Optional, List


def column_stats(df: pd.DataFrame, column: str) -> Tuple[float, float]:
    """
    Calculate mean and standard deviation for a numeric column.

    Args:
        df: DataFrame containing the column.
        column: Column name to analyze.

    Returns:
        Tuple of (mean, standard deviation).

    Raises:
        ValueError: If column doesn't exist or contains non-numeric data.
    """
    if column not in df.columns:
        raise ValueError(f"Column '{column}' not found in DataFrame.")

    try:
        mean_val = df[column].mean()
        std_val = df[column].std()
        return mean_val, std_val
    except TypeError:
        raise ValueError(f"Column '{column}' contains non-numeric data.")


def percentage_within_threshold(df: pd.DataFrame,
                                column: str,
                                threshold: float,
                                comparison: str = 'less') -> float:
    """
    Calculate percentage of values within a threshold.

    Args:
        df: DataFrame containing the column.
        column: Column name to analyze.
        threshold: Threshold value.
        comparison: 'less' for <=, 'greater' for >=.

    Returns:
        Percentage of values meeting the criteria.

    Example:
        >>> pct = percentage_within_threshold(df, 'CA Distance', 26, 'less')
        >>> print(f"{pct:.1f}% of crosslinks are within 26 Angstroms")
    """
    if column not in df.columns:
        raise ValueError(f"Column '{column}' not found in DataFrame.")

    total = len(df)
    if total == 0:
        return 0.0

    if comparison == 'less':
        count = (df[column] <= threshold).sum()
    elif comparison == 'greater':
        count = (df[column] >= threshold).sum()
    else:
        raise ValueError("comparison must be 'less' or 'greater'.")

    return (count / total) * 100


def satisfaction_rate(df: pd.DataFrame,
                      distance_threshold: float = 26.0,
                      spectral_count_min: int = 1) -> float:
    """
    Calculate the percentage of crosslinks satisfied by the structure.

    A crosslink is "satisfied" if its C-alpha distance is within the
    expected crosslinker span.

    Args:
        df: DataFrame with 'CA Distance' and 'Spectral Count' columns.
        distance_threshold: Maximum allowed C-alpha distance.
        spectral_count_min: Minimum spectral count to consider.

    Returns:
        Satisfaction rate as a percentage.
    """
    filtered = df[df['Spectral Count'] >= spectral_count_min]

    if len(filtered) == 0:
        return 0.0

    satisfied = (filtered['CA Distance'] <= distance_threshold).sum()
    return (satisfied / len(filtered)) * 100


def find_common_crosslinks(df1: pd.DataFrame, df2: pd.DataFrame) -> pd.DataFrame:
    """
    Find crosslinks present in both DataFrames.

    Spectral counts are averaged from both datasets.

    Args:
        df1: First DataFrame with crosslink data.
        df2: Second DataFrame with crosslink data.

    Returns:
        DataFrame containing common crosslinks with averaged spectral counts.
    """
    required = {'Residue1', 'Residue2', 'Spectral Count'}
    if not required.issubset(df1.columns) or not required.issubset(df2.columns):
        raise ValueError("DataFrames missing required columns.")

    df1 = df1.copy()
    df2 = df2.copy()

    # Create canonical pairs
    df1['_Canonical'] = df1.apply(
        lambda x: tuple(sorted([x['Residue1'], x['Residue2']])), axis=1
    )
    df2['_Canonical'] = df2.apply(
        lambda x: tuple(sorted([x['Residue1'], x['Residue2']])), axis=1
    )

    # Merge on canonical pairs
    merged = pd.merge(df1, df2, on='_Canonical', suffixes=('_1', '_2'))

    # Average spectral counts
    merged['Spectral Count'] = (
        merged['Spectral Count_1'] + merged['Spectral Count_2']
    ) / 2

    # Select and rename columns
    result = merged[['Residue1_1', 'Residue2_1', 'Spectral Count']].copy()
    result.columns = ['Residue1', 'Residue2', 'Spectral Count']

    # Add other columns from df1 if present
    if 'CA Distance' in df1.columns:
        ca_map = dict(zip(df1['_Canonical'], df1['CA Distance']))
        result['CA Distance'] = merged['_Canonical'].map(ca_map)

    if 'Sequence Distance' in df1.columns:
        seq_map = dict(zip(df1['_Canonical'], df1['Sequence Distance']))
        result['Sequence Distance'] = merged['_Canonical'].map(seq_map)

    if 'Domain Association' in df1.columns:
        dom_map = dict(zip(df1['_Canonical'], df1['Domain Association']))
        result['Domain Association'] = merged['_Canonical'].map(dom_map)

    return result


def find_unique_crosslinks(df1: pd.DataFrame,
                           df2: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Find crosslinks unique to each DataFrame.

    Args:
        df1: First DataFrame with crosslink data.
        df2: Second DataFrame with crosslink data.

    Returns:
        Tuple of (crosslinks unique to df1, crosslinks unique to df2).
    """
    # Create canonical pairs
    pairs1 = df1.apply(
        lambda x: tuple(sorted([x['Residue1'], x['Residue2']])), axis=1
    )
    pairs2 = df2.apply(
        lambda x: tuple(sorted([x['Residue1'], x['Residue2']])), axis=1
    )

    set1 = set(pairs1)
    set2 = set(pairs2)

    unique_to_1 = set1 - set2
    unique_to_2 = set2 - set1

    df1_unique = df1[pairs1.isin(unique_to_1)].copy()
    df2_unique = df2[pairs2.isin(unique_to_2)].copy()

    return df1_unique, df2_unique


def merge_crosslink_datasets(*dataframes: pd.DataFrame,
                             aggregate_spectral: str = 'sum') -> pd.DataFrame:
    """
    Merge multiple crosslink DataFrames, combining duplicate entries.

    Args:
        *dataframes: Variable number of DataFrames to merge.
        aggregate_spectral: How to aggregate spectral counts for duplicates.
                           Options: 'sum', 'mean', 'max'.

    Returns:
        Merged DataFrame with unique crosslinks.
    """
    if not dataframes:
        raise ValueError("At least one DataFrame required.")

    # Add canonical pairs to each
    dfs_with_canonical = []
    for df in dataframes:
        df = df.copy()
        df['_Canonical'] = df.apply(
            lambda x: tuple(sorted([x['Residue1'], x['Residue2']])), axis=1
        )
        dfs_with_canonical.append(df)

    # Concatenate all
    combined = pd.concat(dfs_with_canonical, ignore_index=True)

    # Aggregation function
    agg_funcs = {
        'Residue1': 'first',
        'Residue2': 'first',
        'Spectral Count': aggregate_spectral,
    }

    # Add optional columns
    if 'CA Distance' in combined.columns:
        agg_funcs['CA Distance'] = 'first'
    if 'Sequence Distance' in combined.columns:
        agg_funcs['Sequence Distance'] = 'first'
    if 'Domain Association' in combined.columns:
        agg_funcs['Domain Association'] = 'first'

    # Group and aggregate
    result = combined.groupby('_Canonical').agg(agg_funcs).reset_index()
    result.drop(columns=['_Canonical'], inplace=True)

    return result


def summary_report(df: pd.DataFrame,
                   spectral_thresholds: List[int] = [1, 5, 10],
                   distance_thresholds: List[float] = [20, 26, 50]) -> pd.DataFrame:
    """
    Generate a summary report of crosslink statistics.

    Args:
        df: DataFrame with crosslink data.
        spectral_thresholds: Spectral count thresholds to report.
        distance_thresholds: C-alpha distance thresholds to report.

    Returns:
        DataFrame with summary statistics.
    """
    results = []

    for sc_thresh in spectral_thresholds:
        filtered = df[df['Spectral Count'] >= sc_thresh]
        n_crosslinks = len(filtered)

        row = {
            'Spectral Count >=': sc_thresh,
            'N Crosslinks': n_crosslinks,
        }

        if 'CA Distance' in df.columns and n_crosslinks > 0:
            for dist in distance_thresholds:
                pct = (filtered['CA Distance'] <= dist).mean() * 100
                row[f'% CA <= {dist}A'] = round(pct, 1)

        results.append(row)

    return pd.DataFrame(results)
