# xlinkanalysis

A Python package for analyzing cross-linking mass spectrometry (XL-MS) data, calculating distances from protein structures, and generating visualizations.

## Features

- **Distance Calculations**: Calculate C-alpha distances from PDB structures for crosslinked residue pairs
- **Domain Annotation**: Map residues to protein domains and analyze domain-level crosslinking patterns
- **Flexible Filtering**: Filter crosslinks by spectral count, distance, residue, or domain association
- **Statistical Analysis**: Generate summary statistics and compare datasets
- **Publication-Ready Plots**: Create cumulative distributions, scatter plots, and heatmaps
- **ChimeraX Integration**: Generate visualization scripts for UCSF ChimeraX

## Installation

### From source (recommended for development)

```bash
git clone https://github.com/berndsenlab/xlinkanalysis.git
cd xlinkanalysis
pip install -e .
```

### Dependencies only

```bash
pip install -r requirements.txt
```

## Quick Start

### Command Line

```bash
# Basic analysis
python scripts/run_analysis.py \
    --input data/raw/crosslinks.csv \
    --pdb data/external/structure.pdb \
    --output results/

# With all options
python scripts/run_analysis.py \
    --input data/raw/crosslinks.csv \
    --pdb data/external/structure.pdb \
    --sequence data/external/sequence.fasta \
    --domains config/domains.yaml \
    --config config/analysis_config.yaml \
    --output results/ \
    --spectral-min 5 \
    --chimera-script
```

### Python API

```python
from xlinkanalysis import io, distances, domains, filtering, visualization

# Load data
df = io.load_crosslinks("crosslinks.csv")

# Calculate distances from structure
df = distances.calculate_ca_distances("structure.pdb", df)
df = distances.add_sequence_distance(df)

# Annotate with domain information
df = domains.annotate_domains(df, "config/domains.yaml")

# Filter high-confidence crosslinks
df_filtered = filtering.by_spectral_count(df, min_count=5)
df_satisfied = filtering.satisfied_crosslinks(df_filtered, max_distance=26)

# Generate plots
fig = visualization.plot_ca_vs_spectral(df, ca_threshold=26)
fig.savefig("ca_vs_spectral.png")

# Save results
io.save_crosslinks(df_filtered, "filtered_crosslinks.csv")
```

## Project Structure

```
xlinkanalysis/
├── config/                 # Configuration files
│   ├── domains.yaml       # Protein domain definitions
│   └── analysis_config.yaml
├── src/xlinkanalysis/     # Main package
│   ├── io.py              # Data import/export
│   ├── distances.py       # Distance calculations
│   ├── domains.py         # Domain annotation
│   ├── filtering.py       # Data filtering
│   ├── statistics.py      # Statistical analysis
│   ├── visualization.py   # Plotting functions
│   └── chimera.py         # ChimeraX integration
├── scripts/               # Command-line tools
│   ├── run_analysis.py    # Main analysis pipeline
│   └── extract_residues.py
├── data/                  # Data directories
│   ├── raw/              # Original MS output
│   ├── processed/        # Intermediate files
│   └── external/         # Reference files (PDB, FASTA)
├── results/              # Output directories
│   ├── figures/
│   ├── tables/
│   └── chimera/
└── tests/                # Unit tests
```

## Configuration

### Domain Definitions (`config/domains.yaml`)

Define protein domains for annotation:

```yaml
domains:
  NTD:
    ranges: [[1, 994]]
    color: "#1f77b4"
    description: "N-terminal domain"

  beta-belt:
    ranges: [[995, 2378]]
    color: "#ff7f0e"
```

### Analysis Parameters (`config/analysis_config.yaml`)

Configure thresholds and plotting options:

```yaml
filtering:
  spectral_count_min: 5
  ca_distance_max: 26

plotting:
  figure_size: [10, 6]
  dpi: 300
```

## Input File Format

The basic input format is a CSV with three columns (no header):

```
Residue1,Residue2,Spectral Count
100,150,12
200,250,8
```

Or with header:
```
Residue1,Residue2,Spectral Count
100,150,12
200,250,8
```

## ChimeraX Visualization

Generate and use ChimeraX scripts:

```python
from xlinkanalysis import chimera

# Generate distance commands
commands = chimera.format_distance_commands(df, model=1)
print(commands)

# Generate full visualization script
script = chimera.generate_visualization_script(df, "structure.pdb")
chimera.save_chimera_script(script, "visualize.cxc")
```

In ChimeraX:
```
open visualize.cxc
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

MIT License - see LICENSE file for details.

## Citation

If you use this software in your research, please cite:

[Citation information will be added upon publication]

## Contact

Berndsen Lab - University of Missouri
