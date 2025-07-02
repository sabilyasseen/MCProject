# Monte Carlo Visualization Suite

A modular Python toolkit for analyzing and visualizing Monte Carlo simulation results.

## Overview

The original large notebook has been broken down into focused, testable modules:

- **color_utils.py** - Color palette and utility functions
- **filename_factory.py** - Filename and label generation  
- **cluster_analysis.py** - Cluster size analysis and plotting
- **energy_analysis.py** - Energy convergence analysis
- **grid_visualization.py** - Static grid configuration plotting
- **visualization_drivers.py** - Main orchestration functions

## Quick Start

1. Install dependencies:
```bash
pip install -r requirements.txt
```

2. Test the system:
```bash
python test_visualization_suite.py
```

3. Basic usage:
```python
from visualization_drivers import simple_batch

simple_batch(
    "BJ=-0.44_Kawasaki_0_tol0.00_scp0.50",
    "BJ=-0.44_DR_1_tol0.50_scp0.50",
    plot_energy=True,
    plot_clusters=True,
    plot_grids=True
)
```

## Features

- **Energy Analysis**: Convergence plots, PDFs, statistical tables
- **Cluster Analysis**: Dynamic detection, boundary/interior tracking
- **Grid Visualization**: Species mapping, cluster overlays
- **Modular Design**: Clean separation, independent testing

## Data Format

CSV files should follow this naming convention:
```
BJ={value}_{algorithm}_{type}_{tolerance}_{scp}.csv
```

Required columns include Energy, BoundarySize_Cluster{N}, InteriorSize_Cluster{N}.

## Main Functions

- `simple_batch()` - Process pre-configured filenames
- `run_visualizations()` - Complex parameter combinations  
- `batch_cluster_analysis()` - Cluster analysis only
- `plot_grid_from_file()` - Single grid visualization

## Testing

Run the test suite to verify all components work correctly:
```bash
python test_visualization_suite.py
```