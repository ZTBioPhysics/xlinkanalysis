#!/usr/bin/env python3
"""
Example Workflow: Clean XL-MS Analysis

This script demonstrates the recommended workflow for analyzing crosslink data.
Instead of creating many intermediate variables (df, filtered_df1, filtered_df2...),
we use a staged approach with meaningful names.

STAGE 1: Create annotated dataset (run once, save the result)
STAGE 2: Load annotated data and answer specific questions

The key insight: Once you have the annotated CSV with all distances and domains,
you never need to recalculate them. Just load and filter as needed.
"""

# =============================================================================
# SETUP
# =============================================================================

from xlinkanalysis import io, distances, domains, filtering, statistics, visualization, chimera
import matplotlib.pyplot as plt

# Paths - EDIT THESE FOR YOUR DATA
RAW_CROSSLINKS = "/path/to/your/crosslinks.csv"  # Your crosslink CSV file
PDB_FILE = "/path/to/your/structure.pdb"          # Your PDB structure file
DOMAINS_CONFIG = "config/domains.yaml"            # Relative path to domains config
OUTPUT_DIR = "results"                            # Relative path to output directory

# =============================================================================
# STAGE 1: CREATE ANNOTATED DATASET (Run this once)
# =============================================================================
# This is the expensive step - calculating distances from PDB.
# Do it once, save the result, then use that file for all subsequent analysis.

def create_annotated_dataset():
    """
    Load raw crosslinks, calculate all distances, annotate domains, and save.

    Run this ONCE when you get new data. After that, just load the annotated file.
    """
    print("=" * 60)
    print("STAGE 1: Creating annotated dataset")
    print("=" * 60)

    # Load raw data
    print(f"\n1. Loading raw crosslinks from:\n   {RAW_CROSSLINKS}")
    df = io.load_crosslinks(RAW_CROSSLINKS, has_header=False)
    print(f"   Loaded {len(df)} crosslinks")

    # Calculate C-alpha distances (this is the slow step - reads PDB)
    print(f"\n2. Calculating C-alpha distances from PDB...")
    df = distances.calculate_ca_distances(PDB_FILE, df)

    # Add sequence distance
    print("3. Adding sequence distances...")
    df = distances.add_sequence_distance(df)

    # Annotate domains
    print("4. Annotating domains...")
    df = domains.annotate_domains(df, DOMAINS_CONFIG)

    # Save the annotated dataset - THIS IS YOUR NEW STARTING POINT
    output_file = f"{OUTPUT_DIR}/tables/HDL_crosslinks_annotated.csv"
    io.save_crosslinks(df, output_file)

    print(f"\n✓ Annotated dataset saved to:\n  {output_file}")
    print("\nFrom now on, start your analysis by loading this file!")

    return df


# =============================================================================
# STAGE 2: ANALYSIS (Load annotated data and ask questions)
# =============================================================================
# This is where you do your actual analysis. Since everything is pre-calculated,
# filtering is instant. No need to track multiple DataFrames - just filter on demand.

def analyze_crosslinks(df):
    """
    Example analyses you can run on the annotated dataset.

    KEY CONCEPT: Don't create lots of intermediate variables.
    Instead, filter directly when you need something specific.
    """
    print("\n" + "=" * 60)
    print("STAGE 2: Analysis")
    print("=" * 60)

    # -------------------------------------------------------------------------
    # BASIC STATISTICS (on all data)
    # -------------------------------------------------------------------------
    print("\n--- BASIC STATISTICS ---")
    print(f"Total crosslinks: {len(df)}")

    ca_mean, ca_std = statistics.column_stats(df, 'CA Distance')
    print(f"C-alpha distance: {ca_mean:.1f} ± {ca_std:.1f} Å")

    # What percentage are within crosslinker span?
    pct_26 = statistics.percentage_within_threshold(df, 'CA Distance', 26)
    pct_50 = statistics.percentage_within_threshold(df, 'CA Distance', 50)
    print(f"Within 26 Å: {pct_26:.1f}%")
    print(f"Within 50 Å: {pct_50:.1f}%")

    # -------------------------------------------------------------------------
    # QUESTION 1: What are the high-confidence crosslinks?
    # -------------------------------------------------------------------------
    print("\n--- HIGH-CONFIDENCE CROSSLINKS (Spectral Count ≥ 5) ---")

    # Filter and analyze in one go - no intermediate variable needed
    high_conf = filtering.by_spectral_count(df, min_count=5)
    print(f"Count: {len(high_conf)}")

    sat_rate = statistics.satisfaction_rate(high_conf, distance_threshold=26)
    print(f"Satisfaction rate (≤26 Å): {sat_rate:.1f}%")

    # Save this subset if you need it
    io.save_crosslinks(high_conf, f"{OUTPUT_DIR}/tables/HDL_high_confidence_sc5.csv")

    # -------------------------------------------------------------------------
    # QUESTION 2: What are the inter-domain crosslinks?
    # -------------------------------------------------------------------------
    print("\n--- INTER-DOMAIN CROSSLINKS ---")

    # Chain filters together for complex queries
    inter_domain = filtering.by_domain(
        filtering.by_spectral_count(df, min_count=5),  # Start with high-confidence
        exclude_intra=True  # Only inter-domain
    )
    print(f"High-confidence inter-domain crosslinks: {len(inter_domain)}")

    # -------------------------------------------------------------------------
    # QUESTION 3: Which crosslinks are "violated" (too far apart)?
    # -------------------------------------------------------------------------
    print("\n--- VIOLATED CROSSLINKS (High SC but >26 Å) ---")

    violated = filtering.violated_crosslinks(df, min_distance=26, min_spectral_count=5)
    print(f"Count: {len(violated)}")

    if len(violated) > 0:
        print("\nThese might indicate conformational dynamics:")
        print(violated[['Residue1', 'Residue2', 'CA Distance', 'Spectral Count', 'Domain Association']].head(10).to_string(index=False))

    # -------------------------------------------------------------------------
    # QUESTION 4: What's happening in a specific domain?
    # -------------------------------------------------------------------------
    print("\n--- VWF-DOMAIN CROSSLINKS ---")

    vwf_links = filtering.by_domain(df, domain="VWF-domain")
    vwf_high_conf = filtering.by_spectral_count(vwf_links, min_count=5)
    print(f"Total involving VWF-domain: {len(vwf_links)}")
    print(f"High-confidence: {len(vwf_high_conf)}")

    # -------------------------------------------------------------------------
    # QUESTION 5: What crosslinks involve a specific residue?
    # -------------------------------------------------------------------------
    print("\n--- CROSSLINKS INVOLVING RESIDUE 1500 ---")

    res_1500 = filtering.by_residue(df, 1500)
    print(f"Crosslinks involving residue 1500: {len(res_1500)}")
    if len(res_1500) > 0:
        print(res_1500.to_string(index=False))

    # -------------------------------------------------------------------------
    # DOMAIN-LEVEL SUMMARY
    # -------------------------------------------------------------------------
    print("\n--- DOMAIN INTERACTION SUMMARY ---")

    # Use high-confidence data for domain analysis
    high_conf = filtering.by_spectral_count(df, min_count=5)
    interaction_summary = domains.summarize_by_domain_interaction(high_conf)
    print(interaction_summary.to_string(index=False))

    return df


def create_plots(df):
    """Generate publication-quality plots."""
    print("\n" + "=" * 60)
    print("GENERATING PLOTS")
    print("=" * 60)

    visualization.set_publication_style()

    # Plot 1: CA Distance vs Spectral Count
    print("\n1. CA Distance vs Spectral Count...")
    fig = visualization.plot_ca_vs_spectral(df, ca_threshold=26, spectral_threshold=5)
    visualization.save_figure(fig, f"{OUTPUT_DIR}/figures/ca_vs_spectral")
    plt.close()

    # Plot 2: Cumulative distribution comparison
    print("2. Cumulative distribution...")
    df_sc5 = filtering.by_spectral_count(df, 5)
    df_sc10 = filtering.by_spectral_count(df, 10)

    fig = visualization.plot_ca_cumulative_comparison(
        [df, df_sc5, df_sc10],
        [f"All (n={len(df)})", f"SC≥5 (n={len(df_sc5)})", f"SC≥10 (n={len(df_sc10)})"],
        threshold_line=26
    )
    visualization.save_figure(fig, f"{OUTPUT_DIR}/figures/cumulative_distribution")
    plt.close()

    # Plot 3: Sequence vs CA distance
    print("3. Sequence distance vs CA distance...")
    fig = visualization.plot_sequence_vs_ca(df, highlight_df=df_sc5, highlight_label="SC ≥ 5")
    visualization.save_figure(fig, f"{OUTPUT_DIR}/figures/sequence_vs_ca")
    plt.close()

    print(f"\n✓ Plots saved to: {OUTPUT_DIR}/figures/")


def generate_chimera_script(df):
    """Generate ChimeraX visualization script for high-confidence crosslinks."""
    print("\n" + "=" * 60)
    print("GENERATING CHIMERAX SCRIPT")
    print("=" * 60)

    high_conf = filtering.by_spectral_count(df, min_count=5)

    script = chimera.generate_visualization_script(
        high_conf,
        PDB_FILE,
        color_by_satisfaction=True,
        threshold=26
    )

    script_path = f"{OUTPUT_DIR}/chimera/visualize_HDL_crosslinks.cxc"
    chimera.save_chimera_script(script, script_path)

    print(f"\nTo visualize in ChimeraX:")
    print(f"  open {script_path}")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    # OPTION 1: Create annotated dataset (run this first time)
    df = create_annotated_dataset()

    # OPTION 2: Load existing annotated dataset (use this after first run)
    # df = io.load_annotated_crosslinks(f"{OUTPUT_DIR}/tables/HDL_crosslinks_annotated.csv")

    # Run analysis
    analyze_crosslinks(df)

    # Generate plots
    create_plots(df)

    # Generate ChimeraX script
    generate_chimera_script(df)

    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE!")
    print("=" * 60)
