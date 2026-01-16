"""
Visualization functions for crosslink analysis.

This module provides plotting functions for visualizing crosslink data,
including distance distributions, scatter plots, and cumulative curves.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, List, Tuple, Dict


def set_publication_style(font_size: int = 14) -> None:
    """
    Set matplotlib parameters for publication-quality figures.

    Args:
        font_size: Base font size for labels.
    """
    plt.rcParams.update({
        'font.size': font_size,
        'axes.labelsize': font_size + 2,
        'axes.titlesize': font_size + 2,
        'xtick.labelsize': font_size,
        'ytick.labelsize': font_size,
        'legend.fontsize': font_size,
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
    })


def plot_ca_vs_spectral(df: pd.DataFrame,
                        ca_threshold: Optional[float] = 26,
                        spectral_threshold: Optional[int] = 5,
                        figsize: Tuple[int, int] = (10, 6),
                        title: str = 'C-alpha Distance vs Spectral Count',
                        save_path: Optional[str] = None) -> plt.Figure:
    """
    Scatter plot of C-alpha distance vs spectral count.

    Args:
        df: DataFrame with 'CA Distance' and 'Spectral Count' columns.
        ca_threshold: Draw horizontal line at this distance (or None).
        spectral_threshold: Draw vertical line at this count (or None).
        figsize: Figure size (width, height).
        title: Plot title.
        save_path: If provided, save figure to this path.

    Returns:
        Matplotlib Figure object.
    """
    fig, ax = plt.subplots(figsize=figsize)

    ax.scatter(df['CA Distance'], df['Spectral Count'], alpha=0.7)
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)

    if ca_threshold:
        ax.axvline(x=ca_threshold, color='r', linestyle='--', linewidth=2,
                   label=f'CA = {ca_threshold} A')
    if spectral_threshold:
        ax.axhline(y=spectral_threshold, color='r', linestyle='--', linewidth=2,
                   label=f'SC = {spectral_threshold}')

    ax.set_xlabel('C-alpha Distance (A)')
    ax.set_ylabel('Spectral Count')
    ax.set_title(title)

    if ca_threshold or spectral_threshold:
        ax.legend()

    if save_path:
        fig.savefig(save_path)

    return fig


def plot_sequence_vs_ca(df: pd.DataFrame,
                        highlight_df: Optional[pd.DataFrame] = None,
                        highlight_label: str = 'Filtered',
                        ca_threshold: Optional[float] = 26,
                        figsize: Tuple[int, int] = (10, 6),
                        title: str = 'Sequence Distance vs C-alpha Distance',
                        save_path: Optional[str] = None) -> plt.Figure:
    """
    Scatter plot of sequence distance vs C-alpha distance.

    Args:
        df: DataFrame with 'Sequence Distance' and 'CA Distance' columns.
        highlight_df: Optional second DataFrame to highlight.
        highlight_label: Label for highlighted points.
        ca_threshold: Draw horizontal line at this distance.
        figsize: Figure size.
        title: Plot title.
        save_path: If provided, save figure to this path.

    Returns:
        Matplotlib Figure object.
    """
    fig, ax = plt.subplots(figsize=figsize)

    ax.scatter(df['Sequence Distance'], df['CA Distance'],
               alpha=0.7, label='All Data')

    if highlight_df is not None:
        ax.scatter(highlight_df['Sequence Distance'], highlight_df['CA Distance'],
                   alpha=0.7, color='r', label=highlight_label)

    ax.grid(True, alpha=0.3)

    if ca_threshold:
        ax.axhline(y=ca_threshold, color='red', linestyle='--', linewidth=2)

    ax.set_xlabel('Sequence Distance (residues)')
    ax.set_ylabel('C-alpha Distance (A)')
    ax.set_title(title)
    ax.legend()

    if save_path:
        fig.savefig(save_path)

    return fig


def plot_cumulative_distribution(df: pd.DataFrame,
                                 column: str = 'CA Distance',
                                 label: Optional[str] = None,
                                 ax: Optional[plt.Axes] = None,
                                 color: Optional[str] = None) -> plt.Axes:
    """
    Plot cumulative distribution of a column.

    Can be called multiple times on the same axes to compare datasets.

    Args:
        df: DataFrame containing the column.
        column: Column to plot distribution for.
        label: Legend label for this dataset.
        ax: Existing axes to plot on (creates new if None).
        color: Line color.

    Returns:
        Matplotlib Axes object.

    Example:
        >>> fig, ax = plt.subplots()
        >>> plot_cumulative_distribution(df_all, label="All")
        >>> plot_cumulative_distribution(df_filtered, label="Filtered", ax=ax)
        >>> plt.show()
    """
    if column not in df.columns:
        raise ValueError(f"Column '{column}' not found in DataFrame.")

    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))

    # Calculate cumulative distribution
    values = df[column].dropna().sort_values()
    cumulative = np.arange(1, len(values) + 1) / len(values) * 100

    plot_kwargs = {'drawstyle': 'steps-post'}
    if label:
        plot_kwargs['label'] = label
    if color:
        plot_kwargs['color'] = color

    ax.plot(values, cumulative, **plot_kwargs)
    ax.grid(True, alpha=0.3)

    if label:
        ax.legend()

    return ax


def plot_ca_cumulative_comparison(dataframes: List[pd.DataFrame],
                                  labels: List[str],
                                  threshold_line: Optional[float] = 26,
                                  figsize: Tuple[int, int] = (10, 6),
                                  title: str = 'Cumulative Distribution of C-alpha Distances',
                                  save_path: Optional[str] = None) -> plt.Figure:
    """
    Compare cumulative C-alpha distance distributions across multiple datasets.

    Args:
        dataframes: List of DataFrames to compare.
        labels: Labels for each DataFrame (include n in label if desired).
        threshold_line: Draw vertical line at this distance.
        figsize: Figure size.
        title: Plot title.
        save_path: If provided, save figure to this path.

    Returns:
        Matplotlib Figure object.

    Example:
        >>> dfs = [df_all, df_sc5, df_sc10]
        >>> labels = [f"All (n={len(df_all)})",
        ...           f"SC>=5 (n={len(df_sc5)})",
        ...           f"SC>=10 (n={len(df_sc10)})"]
        >>> plot_ca_cumulative_comparison(dfs, labels)
    """
    fig, ax = plt.subplots(figsize=figsize)

    for df, label in zip(dataframes, labels):
        plot_cumulative_distribution(df, 'CA Distance', label=label, ax=ax)

    if threshold_line:
        ax.axvline(x=threshold_line, color='r', linestyle='--',
                   label=f'{threshold_line} A threshold')

    ax.set_xlabel('C-alpha Distance (A)')
    ax.set_ylabel('Cumulative Percentage (%)')
    ax.set_title(title)
    ax.legend()

    if save_path:
        fig.savefig(save_path)

    return fig


def plot_domain_heatmap(df: pd.DataFrame,
                        value_column: str = 'count',
                        figsize: Tuple[int, int] = (10, 8),
                        cmap: str = 'Blues',
                        title: str = 'Inter-domain Crosslinks',
                        save_path: Optional[str] = None) -> plt.Figure:
    """
    Create a heatmap of crosslinks between domains.

    Args:
        df: DataFrame with 'Domain Association' column, or a pre-pivoted matrix.
        value_column: Column to use for heatmap values ('count' creates counts).
        figsize: Figure size.
        cmap: Colormap name.
        title: Plot title.
        save_path: If provided, save figure to this path.

    Returns:
        Matplotlib Figure object.
    """
    if 'Domain Association' in df.columns:
        # Create count matrix from associations
        df = df.copy()
        parts = df['Domain Association'].str.split(' to ', expand=True)
        df['Domain1'] = parts[0]
        df['Domain2'] = parts[1]

        # Create pivot table
        matrix = pd.crosstab(df['Domain1'], df['Domain2'])

        # Make symmetric (add transpose where not already counted)
        for i in matrix.index:
            for j in matrix.columns:
                if i != j:
                    total = matrix.loc[i, j] + matrix.loc[j, i] if j in matrix.index and i in matrix.columns else matrix.loc[i, j]
                    matrix.loc[i, j] = total
                    if j in matrix.index and i in matrix.columns:
                        matrix.loc[j, i] = total
    else:
        matrix = df

    fig, ax = plt.subplots(figsize=figsize)

    im = ax.imshow(matrix, cmap=cmap, aspect='auto')

    # Labels
    ax.set_xticks(range(len(matrix.columns)))
    ax.set_yticks(range(len(matrix.index)))
    ax.set_xticklabels(matrix.columns, rotation=45, ha='right')
    ax.set_yticklabels(matrix.index)

    # Colorbar
    plt.colorbar(im, ax=ax, label='Number of Crosslinks')

    # Annotate cells
    for i in range(len(matrix.index)):
        for j in range(len(matrix.columns)):
            val = matrix.iloc[i, j]
            if val > 0:
                ax.text(j, i, str(int(val)), ha='center', va='center',
                       color='white' if val > matrix.values.max()/2 else 'black')

    ax.set_title(title)

    if save_path:
        fig.savefig(save_path)

    return fig


def save_figure(fig: plt.Figure,
                path: str,
                formats: List[str] = ['png', 'pdf'],
                dpi: int = 300) -> None:
    """
    Save a figure in multiple formats.

    Args:
        fig: Matplotlib Figure object.
        path: Base path without extension.
        formats: List of formats to save (e.g., ['png', 'pdf', 'svg']).
        dpi: Resolution for raster formats.
    """
    from pathlib import Path
    base = Path(path)

    for fmt in formats:
        output = base.with_suffix(f'.{fmt}')
        fig.savefig(output, dpi=dpi, bbox_inches='tight')
        print(f"Saved: {output}")
