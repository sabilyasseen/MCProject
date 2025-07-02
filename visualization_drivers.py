"""
Main driver functions for Monte Carlo visualization suite.

This module provides high-level orchestration functions that coordinate
all visualization components including energy analysis, cluster analysis,
and grid visualization.
"""

import sys
import math
from pathlib import Path
from typing import List, Optional
import pandas as pd

# Import our modular components
from filename_factory import build_filenames
from energy_analysis import compare_runs
from cluster_analysis import analyze_cluster_sizes
from grid_visualization import plot_grid, parse_grid_config_string
from color_utils import get_distinct_colors

def run_visualizations(*, 
                      BJ: str,
                      include_kawasaki: bool = True,
                      base_variants: Optional[List[str]] = None,
                      signifiers: Optional[List[List[str]]] = None,
                      cluster_types: Optional[List[str]] = None,
                      plot_energy: bool = True,
                      plot_grids: bool = False,
                      plot_moves: bool = False,
                      plot_clusters: bool = True):
    """
    Complex combinatoric workflow for generating filenames and visualizations.
    
    Args:
        BJ: Beta*J parameter for filename generation
        include_kawasaki: Whether to include Kawasaki algorithm files
        base_variants: Base algorithm variants (e.g., ["NewMC", "DR"])
        signifiers: Parameter combinations (e.g., [["tol0.01"], ["scp0.50"]])
        cluster_types: Cluster type variants (e.g., ["1", "2"])
        plot_energy: Whether to plot energy and convergence analysis
        plot_grids: Whether to plot static grid configurations
        plot_moves: Whether to plot move statistics (not yet implemented)
        plot_clusters: Whether to plot cluster size analysis
    """
    base_variants = base_variants or [""]
    signifiers = signifiers or [[""]]
    cluster_types = cluster_types or ["1"]

    csv_files, tbl, col, leg, leg_col, expand = build_filenames(
        BJ=BJ, include_kawasaki=include_kawasaki,
        base_variants=base_variants, signifiers=signifiers, cluster_types=cluster_types)

    if plot_energy:
        compare_runs(csv_files=csv_files, table_labels=tbl, curve_colours=col,
                    legend_labels=leg, legend_colours=leg_col)
    
    if plot_clusters:
        print("\n--- Running Cluster Analysis for Complex Batch ---", file=sys.stderr)
        analyze_cluster_sizes(csv_files=csv_files, table_labels=tbl, curve_colours=col)
    
    if plot_grids:
        print("\n--- Running Static Grid Plots for Complex Batch ---", file=sys.stderr)
        _plot_static_grids(csv_files, tbl)

    if plot_moves:
        print("\n--- Move Analysis not yet implemented ---", file=sys.stderr)


def simple_batch(*m_batch_filenames, 
                plot_energy: bool = True, 
                plot_grids: bool = False, 
                plot_clusters: bool = True):
    """
    Streamlined workflow for generating visualizations from a list of pre-configured filenames.
    Uses a simpler, distinct color scheme for plots.

    Args:
        m_batch_filenames: Base filenames (e.g., ["BJ=1.0_MC_1_tol0.01_scp0.05"])
        plot_energy: Whether to plot energy and convergence data
        plot_grids: Whether to plot static grid configurations
        plot_clusters: Whether to plot cluster size analysis
    """
    csv_files = [f"{bf}.csv" for bf in m_batch_filenames]
    table_labels = list(m_batch_filenames)  # Use base filename as label
    curve_colours = get_distinct_colors(len(m_batch_filenames))
    legend_labels = list(m_batch_filenames)  # Use base filename as label
    legend_colours = curve_colours

    if plot_energy:
        print("\n--- Running Energy Plots for Simple Batch ---", file=sys.stderr)
        compare_runs(csv_files=csv_files, table_labels=table_labels, curve_colours=curve_colours,
                    legend_labels=legend_labels, legend_colours=legend_colours,
                    title="Simple Batch Energy & Convergence Analysis")
    
    if plot_clusters:
        print("\n--- Running Cluster Analysis for Simple Batch ---", file=sys.stderr)
        analyze_cluster_sizes(csv_files=csv_files, table_labels=table_labels, curve_colours=curve_colours)
    
    if plot_grids:
        print("\n--- Running Static Grid Plots for Simple Batch ---", file=sys.stderr)
        _plot_static_grids(csv_files, table_labels, list(m_batch_filenames))


def _plot_static_grids(csv_files: List[str], table_labels: List[str], base_filenames: Optional[List[str]] = None):
    """
    Internal function to plot static grids for a batch of files.
    
    Args:
        csv_files: List of CSV file paths
        table_labels: Labels for each file
        base_filenames: Optional base filenames for grid config files
    """
    for i, (csv_file, lab) in enumerate(zip(csv_files, table_labels)):
        main_csv_file = Path(csv_file)
        
        # Try to find grid config file
        if base_filenames:
            grid_cfg_file = Path(f"{base_filenames[i]}_GridConfig.csv")
        else:
            grid_cfg_file = main_csv_file.with_name(main_csv_file.stem + "_GridConfig.csv")
        
        df_grid_to_plot = pd.DataFrame()
        
        if grid_cfg_file.exists():
            print(f"Debug: Plotting static grid from {grid_cfg_file}", file=sys.stderr)
            try:
                df_grid_to_plot = pd.read_csv(grid_cfg_file)
                plot_grid(df_grid_to_plot, title=f"{lab} (Static Grid)", 
                         save_path=grid_cfg_file.with_suffix(".png"))
            except Exception as e:
                print(f"Error loading or plotting static grid from {grid_cfg_file}: {e}", file=sys.stderr)
                
        elif main_csv_file.exists():
            print(f"Debug: Plotting static grid from main CSV {main_csv_file} (fallback)", file=sys.stderr)
            try:
                _plot_grid_from_main_csv(main_csv_file, lab)
            except Exception as e:
                print(f"Error loading or plotting static grid from main CSV {main_csv_file}: {e}", file=sys.stderr)
        else:
            print(f"Grid config file and main CSV not found for static plot: {grid_cfg_file}", file=sys.stderr)


def _plot_grid_from_main_csv(main_csv_file: Path, label: str):
    """
    Plot grid from the main CSV file by parsing GridConfig columns.
    
    Args:
        main_csv_file: Path to the main CSV file
        label: Label for the plot
    """
    dtype_mapping = {
        "PossiblyChanged": str, "ProposedAtoms": str, "AcceptedAtoms": str,
        "AcceptedPossiblyChanged": str, "GridConfigStarting": str, "GridConfigPotential": str
    }
    
    df_sim_last_row = pd.read_csv(main_csv_file, dtype=dtype_mapping).tail(1)
    if df_sim_last_row.empty:
        print(f"Warning: Main CSV file {main_csv_file.name} is empty for static plot.", file=sys.stderr)
        return

    grid_config_str = None
    N_for_static = 0

    # Determine N from the first available config string in the last row
    for _, row in df_sim_last_row.iterrows():  # This loop will run once for the tail(1) row
        potential_str = str(row.get("GridConfigPotential", "")).strip('"').strip("'")
        starting_str = str(row.get("GridConfigStarting", "")).strip('"').strip("'")

        if potential_str and potential_str not in ['nan', 'NaN', '']:
            num_cells = len(potential_str.split(';'))
            N_for_static = int(math.sqrt(num_cells))
            if N_for_static * N_for_static != num_cells:
                print(f"Warning: Potential GridConfig string implies non-square grid or malformed data. "
                     f"Num cells: {num_cells}. Trying to proceed with N={N_for_static}", file=sys.stderr)
            grid_config_str = potential_str
            title_suffix = "(from Final Potential GridConfig in Main CSV)"
            break
        elif starting_str and starting_str not in ['nan', 'NaN', '']:
            num_cells = len(starting_str.split(';'))
            N_for_static = int(math.sqrt(num_cells))
            if N_for_static * N_for_static != num_cells:
                print(f"Warning: Starting GridConfig string implies non-square grid or malformed data. "
                     f"Num cells: {num_cells}. Trying to proceed with N={N_for_static}", file=sys.stderr)
            grid_config_str = starting_str
            title_suffix = "(from Final Starting GridConfig in Main CSV)"
            break
    
    if grid_config_str and N_for_static > 0:
        df_grid_to_plot = parse_grid_config_string(grid_config_str, N_for_static)
        if not df_grid_to_plot.empty:
            plot_grid(df_grid_to_plot, title=f"{label} {title_suffix}",
                     save_path=main_csv_file.with_suffix(".static_grid.png"),
                     grid_size=N_for_static)
        else:
            print(f"Warning: GridConfig column empty or malformed in {main_csv_file.name} for static plot.", 
                  file=sys.stderr)
    else:
        print(f"Warning: Neither 'GridConfigPotential' nor 'GridConfigStarting' found/valid in "
              f"{main_csv_file.name} for static plot.", file=sys.stderr)


def demo_simple_batch():
    """
    Demonstration function showing how to use simple_batch with example files.
    """
    example_filenames = [
        "BJ=-0.44_Kawasaki_0_tol0.00_scp0.50",
        "BJ=-0.44_DR_1_tol0.50_scp0.50"
    ]

    print("--- Starting Simple Batch Visualization Demo ---", file=sys.stderr)
    simple_batch(
        *example_filenames,  # Use * to unpack the list to positional arguments
        plot_energy=True,    # Set to True to generate energy plots
        plot_grids=True,     # Set to True to generate static grid plots
        plot_clusters=True   # Set to True to generate cluster size analysis
    )
    print("--- Simple Batch Visualization Demo Complete ---", file=sys.stderr)


def demo_complex_batch():
    """
    Demonstration function showing how to use run_visualizations with complex parameters.
    """
    print("--- Starting Complex Batch Visualization Demo ---", file=sys.stderr)
    run_visualizations(
        BJ="-0.44",
        include_kawasaki=True,
        base_variants=["NewMC", "DR"],
        signifiers=[["tol0.01", "tol0.05"], ["scp0.00", "scp0.10"]],
        cluster_types=["1", "2"],
        plot_energy=True,
        plot_grids=False,
        plot_clusters=True
    )
    print("--- Complex Batch Visualization Demo Complete ---", file=sys.stderr)


if __name__ == '__main__':
    # Run demo when executed directly
    print("🎯 MONTE CARLO VISUALIZATION SUITE")
    print("=" * 50)
    print("Running demonstration with example data...")
    demo_simple_batch()