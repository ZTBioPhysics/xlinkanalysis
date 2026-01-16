"""
xlinkanalysis - Cross-linking Mass Spectrometry Analysis Package

A Python package for analyzing cross-linking mass spectrometry (XL-MS) data,
calculating distances from protein structures, and generating visualizations.

Example usage:
    from xlinkanalysis import io, distances, domains, filtering

    # Load crosslink data
    df = io.load_crosslinks("crosslinks.csv")

    # Calculate C-alpha distances from PDB
    df = distances.calculate_ca_distances("structure.pdb", df)

    # Annotate with domain information
    df = domains.annotate_domains(df, "config/domains.yaml")

    # Filter by spectral count
    df_filtered = filtering.by_spectral_count(df, min_count=5)
"""

__version__ = "1.0.0"
__author__ = "Berndsen Lab"

from . import io
from . import distances
from . import domains
from . import filtering
from . import statistics
from . import visualization
from . import chimera
from . import config
