#!/usr/bin/env python3
"""
Cross-linking Mass Spectrometry Analysis Pipeline

This script provides a command-line interface for analyzing XL-MS data.
It can be run directly or used as a template for custom analysis scripts.

Usage:
    python run_analysis.py --input crosslinks.csv --pdb structure.pdb --output results/

Example:
    python run_analysis.py \
        --input /path/to/crosslinks.csv \
        --pdb /path/to/structure.pdb \
        --sequence /path/to/sequence.fasta \
        --domains config/domains.yaml \
        --config config/analysis_config.yaml \
        --output results/
"""

import argparse
import sys
from pathlib import Path

# Add src to path if running from scripts directory
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from xlinkanalysis import io, distances, domains, filtering, statistics, visualization, chimera
from xlinkanalysis.config import load_analysis_config, load_domains


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Analyze cross-linking mass spectrometry data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    # Required arguments
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Path to crosslink CSV file (Residue1, Residue2, Spectral Count)"
    )
    parser.add_argument(
        "--pdb", "-p",
        required=True,
        help="Path to PDB structure file"
    )

    # Optional arguments
    parser.add_argument(
        "--sequence", "-s",
        help="Path to FASTA file with protein sequence (for domain stats)"
    )
    parser.add_argument(
        "--domains", "-d",
        default="config/domains.yaml",
        help="Path to domains YAML config (default: config/domains.yaml)"
    )
    parser.add_argument(
        "--config", "-c",
        default="config/analysis_config.yaml",
        help="Path to analysis config YAML (default: config/analysis_config.yaml)"
    )
    parser.add_argument(
        "--output", "-o",
        default="results/",
        help="Output directory for results (default: results/)"
    )
    parser.add_argument(
        "--has-header",
        action="store_true",
        help="Input CSV has header row"
    )
    parser.add_argument(
        "--spectral-min",
        type=int,
        default=None,
        help="Override minimum spectral count threshold"
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip generating plots"
    )
    parser.add_argument(
        "--chimera-script",
        action="store_true",
        help="Generate ChimeraX visualization script"
    )

    return parser.parse_args()


def main():
    """Main analysis pipeline."""
    args = parse_args()

    # Setup paths
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "figures").mkdir(exist_ok=True)
    (output_dir / "tables").mkdir(exist_ok=True)
    (output_dir / "chimera").mkdir(exist_ok=True)

    # Load configuration
    print(f"Loading configuration from {args.config}...")
    try:
        config = load_analysis_config(args.config)
    except FileNotFoundError:
        print(f"Warning: Config file not found, using defaults")
        config = {'filtering': {'spectral_count_min': 5, 'ca_distance_max': 26}}

    # Get thresholds
    spectral_min = args.spectral_min or config.get('filtering', {}).get('spectral_count_min', 5)

    # Load crosslink data
    print(f"\nLoading crosslinks from {args.input}...")
    df = io.load_crosslinks(args.input, has_header=args.has_header)
    print(f"  Loaded {len(df)} crosslinks")

    # Calculate C-alpha distances
    print(f"\nCalculating C-alpha distances from {args.pdb}...")
    df = distances.calculate_ca_distances(args.pdb, df)

    # Add sequence distance
    df = distances.add_sequence_distance(df)

    # Annotate domains
    print(f"\nAnnotating domains from {args.domains}...")
    try:
        df = domains.annotate_domains(df, args.domains)
    except FileNotFoundError:
        print(f"Warning: Domains config not found, skipping domain annotation")

    # Save annotated data
    annotated_path = output_dir / "tables" / "crosslinks_annotated.csv"
    io.save_crosslinks(df, annotated_path)

    # Generate summary statistics
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)

    print(f"\nTotal crosslinks: {len(df)}")

    if 'CA Distance' in df.columns:
        ca_mean, ca_std = statistics.column_stats(df, 'CA Distance')
        print(f"C-alpha distance: {ca_mean:.1f} +/- {ca_std:.1f} A")

        pct_26 = statistics.percentage_within_threshold(df, 'CA Distance', 26)
        pct_50 = statistics.percentage_within_threshold(df, 'CA Distance', 50)
        print(f"  <= 26 A: {pct_26:.1f}%")
        print(f"  <= 50 A: {pct_50:.1f}%")

    # Filter by spectral count
    print(f"\nFiltering by spectral count >= {spectral_min}...")
    df_filtered = filtering.by_spectral_count(df, spectral_min)
    print(f"  {len(df_filtered)} crosslinks pass filter")

    if 'CA Distance' in df_filtered.columns and len(df_filtered) > 0:
        sat_rate = statistics.satisfaction_rate(df_filtered, distance_threshold=26)
        print(f"  Satisfaction rate (26 A): {sat_rate:.1f}%")

    # Save filtered data
    filtered_path = output_dir / "tables" / f"crosslinks_sc{spectral_min}.csv"
    io.save_crosslinks(df_filtered, filtered_path)

    # Domain statistics
    if 'Domain Association' in df.columns and args.sequence:
        print("\nDomain-level statistics:")
        sequence = io.load_sequence_from_fasta(args.sequence)
        domain_stats = domains.domain_crosslink_stats(df_filtered, sequence, args.domains)
        print(domain_stats.to_string(index=False))

        domain_stats_path = output_dir / "tables" / "domain_statistics.csv"
        io.save_crosslinks(domain_stats, domain_stats_path)

        # Domain interaction summary
        interaction_summary = domains.summarize_by_domain_interaction(df_filtered)
        print("\nDomain interactions:")
        print(interaction_summary.to_string(index=False))

        interaction_path = output_dir / "tables" / "domain_interactions.csv"
        io.save_crosslinks(interaction_summary, interaction_path)

    # Summary report
    report = statistics.summary_report(df)
    print("\nSummary by spectral count threshold:")
    print(report.to_string(index=False))
    report_path = output_dir / "tables" / "summary_report.csv"
    io.save_crosslinks(report, report_path)

    # Generate plots
    if not args.no_plots and 'CA Distance' in df.columns:
        print("\nGenerating plots...")
        visualization.set_publication_style()

        # CA vs Spectral Count
        fig = visualization.plot_ca_vs_spectral(df, ca_threshold=26, spectral_threshold=spectral_min)
        visualization.save_figure(fig, output_dir / "figures" / "ca_vs_spectral")

        # Sequence vs CA distance
        fig = visualization.plot_sequence_vs_ca(df, highlight_df=df_filtered,
                                                highlight_label=f"SC >= {spectral_min}")
        visualization.save_figure(fig, output_dir / "figures" / "sequence_vs_ca")

        # Cumulative distribution comparison
        df_sc5 = filtering.by_spectral_count(df, 5)
        df_sc10 = filtering.by_spectral_count(df, 10)
        dfs = [df, df_sc5, df_sc10]
        labels = [f"All (n={len(df)})",
                  f"SC>=5 (n={len(df_sc5)})",
                  f"SC>=10 (n={len(df_sc10)})"]
        fig = visualization.plot_ca_cumulative_comparison(dfs, labels)
        visualization.save_figure(fig, output_dir / "figures" / "cumulative_distribution")

        print("  Plots saved to results/figures/")

    # Generate ChimeraX script
    if args.chimera_script:
        print("\nGenerating ChimeraX visualization script...")
        script = chimera.generate_visualization_script(
            df_filtered,
            args.pdb,
            color_by_satisfaction=True,
            threshold=26
        )
        script_path = output_dir / "chimera" / "visualize_crosslinks.cxc"
        chimera.save_chimera_script(script, script_path)

    print("\n" + "="*60)
    print("Analysis complete!")
    print(f"Results saved to: {output_dir}")
    print("="*60)


if __name__ == "__main__":
    main()
