"""
Grid visualization module for Monte Carlo visualization suite.

This module provides functions for plotting static grid configurations
with cluster information, species visualization, and statistical overlays.
"""

import sys
import re
import math
from pathlib import Path
from typing import Optional
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from color_utils import get_distinct_colors

def parse_grid_config_string(grid_config_str: str, N: int) -> pd.DataFrame:
    """
    Parse a grid configuration string into a DataFrame.
    
    Args:
        grid_config_str: String containing grid configuration data
        N: Grid size (N x N)
        
    Returns:
        DataFrame with columns: x, y, species, cluster_id, cell_type
    """
    if not grid_config_str or grid_config_str in ['nan', 'NaN', '']:
        return pd.DataFrame()
    
    try:
        # Remove quotes and split by semicolons
        config_str = str(grid_config_str).strip('"').strip("'")
        if not config_str:
            return pd.DataFrame()
            
        cells = config_str.split(';')
        
        data = []
        for i, cell_str in enumerate(cells):
            if not cell_str or cell_str.strip() == '':
                continue
                
            # Calculate x, y coordinates from index
            x = i % N
            y = i // N
            
            # Parse cell information
            # Expected format: something like "species:cluster_id:cell_type" or similar
            # We'll try to extract numbers and categorize them
            numbers = re.findall(r'-?\d+', cell_str)
            
            if len(numbers) >= 2:
                species = int(numbers[0])
                cluster_id = int(numbers[1])
            elif len(numbers) == 1:
                species = int(numbers[0])
                cluster_id = 0  # Default cluster
            else:
                species = 0
                cluster_id = 0
            
            # Determine cell type based on string content
            if 'BOUNDARY' in cell_str.upper():
                cell_type = 'BOUNDARY'
            elif 'INTERIOR' in cell_str.upper():
                cell_type = 'INTERIOR'
            else:
                cell_type = 'ANY'
            
            data.append({
                'x': x,
                'y': y,
                'species': species,
                'cluster_id': cluster_id,
                'cell_type': cell_type
            })
        
        return pd.DataFrame(data)
        
    except Exception as e:
        print(f"Error parsing grid config string: {e}", file=sys.stderr)
        return pd.DataFrame()

def plot_grid(cfg_df: pd.DataFrame, *, 
              title: Optional[str] = None, 
              figsize=(8, 8),
              ax=None, 
              save_path: Optional[Path] = None, 
              grid_size: Optional[int] = None):
    """
    Plot the grid configuration from GridConfig CSV data.
    
    Args:
        cfg_df: DataFrame with grid configuration data
        title: Plot title
        figsize: Figure size tuple
        ax: Matplotlib axis (if None, creates new figure)
        save_path: Path to save the plot
        grid_size: Override grid size detection
    """
    # Validate input data
    if cfg_df.empty:
        print("Warning: Empty DataFrame provided to plot_grid", file=sys.stderr)
        fig, ax_empty = plt.subplots(figsize=figsize)
        ax_empty.text(0.5, 0.5, "Empty Grid (No Config)", ha='center', va='center', 
                     transform=ax_empty.transAxes, fontsize=16)
        ax_empty.set(xlim=(0, 1), ylim=(0, 1), xticks=[], yticks=[])
        if save_path:
            fig.savefig(save_path, dpi=250, bbox_inches="tight")
            plt.close(fig)
        else:
            plt.show()
        return
    
    # Validate required columns
    required_cols = ["x", "y", "species", "cluster_id", "cell_type"]
    missing_cols = [col for col in required_cols if col not in cfg_df.columns]
    if missing_cols:
        print(f"Error: Missing required columns in GridConfig: {missing_cols}", file=sys.stderr)
        print(f"Available columns: {list(cfg_df.columns)}", file=sys.stderr)
        return
    
    # Determine grid size
    if grid_size is None:
        N = int(max(cfg_df["x"].max(), cfg_df["y"].max()) + 1) if not cfg_df.empty else 0
    else:
        N = grid_size
    
    print(f"Debug: Grid size determined as {N}x{N} ({len(cfg_df)} cells)", file=sys.stderr)
    
    # Create figure
    create_new_fig = ax is None
    if create_new_fig:
        fig, ax_local = plt.subplots(1, 1, figsize=figsize)
    else:
        ax_local = ax
        fig = ax.figure
        
    ax_local.set(xlim=(-0.5, N-0.5), ylim=(-0.5, N-0.5), aspect="equal", xticks=[], yticks=[])
    ax_local.set_title(title or "Grid Configuration", fontweight="bold", fontsize=14)
    
    # Use a copy to avoid modifying original data
    df_display = cfg_df.copy()
    
    # Get unique cluster IDs for color mapping
    unique_clusters = sorted(df_display["cluster_id"].unique())
    cluster_colors = get_distinct_colors(len(unique_clusters))
    cluster_color_map = {cluster_id: cluster_colors[i] for i, cluster_id in enumerate(unique_clusters)}
    
    # Draw grid cells
    for _, r in df_display.iterrows():
        # Determine cell color based on species
        if r["species"] == 0:
            base_color = "#d62728"  # Red for species 0
        elif r["species"] == 1:  
            base_color = "#1f77b4"  # Blue for species 1
        else:
            base_color = "#808080"  # Gray for unknown species
        
        # Determine transparency based on cell type
        if r["cell_type"] == "INTERIOR":
            alpha = 0.4
        elif r["cell_type"] == "BOUNDARY":
            alpha = 0.8
        else:  # "ANY" or other
            alpha = 0.6
        
        # Create rectangle for the cell
        rect = plt.Rectangle((r["x"]-0.4, r["y"]-0.4), 0.8, 0.8,
                           facecolor=base_color, alpha=alpha, 
                           edgecolor='black', linewidth=0.5)
        ax_local.add_patch(rect)
        
        # Add cluster ID text in the center of the cell
        txt = ax_local.text(r["x"], r["y"], str(r["cluster_id"]),
                          ha="center", va="center", fontsize=8, 
                          color='white', weight='bold')
        ax_local.add_artist(txt)
    
    # Add legend
    legend_elements = [
        Patch(facecolor="#d62728", alpha=0.8, label='Species 0 (Boundary)'),
        Patch(facecolor="#d62728", alpha=0.4, label='Species 0 (Interior)'),
        Patch(facecolor="#1f77b4", alpha=0.8, label='Species 1 (Boundary)'),
        Patch(facecolor="#1f77b4", alpha=0.4, label='Species 1 (Interior)')
    ]
    ax_local.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=10)
    
    # Add grid lines for better visualization
    for i in range(N+1):
        ax_local.axhline(i-0.5, color='gray', linewidth=0.3, alpha=0.5)
        ax_local.axvline(i-0.5, color='gray', linewidth=0.3, alpha=0.5)
    
    # Add cluster statistics as text
    cluster_stats = []
    for cluster_id in unique_clusters:
        cluster_cells = df_display[df_display["cluster_id"] == cluster_id]
        boundary_count = len(cluster_cells[cluster_cells["cell_type"] == "BOUNDARY"])
        interior_count = len(cluster_cells[cluster_cells["cell_type"] == "INTERIOR"])
        total_count = len(cluster_cells)
        cluster_stats.append(f"Cluster {cluster_id}: B={boundary_count}, I={interior_count}, Total={total_count}")
    
    stats_text = "\n".join(cluster_stats)
    ax_local.text(1.02, 0.5, stats_text, transform=ax_local.transAxes, 
                 fontsize=9, verticalalignment='center',
                 bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
    
    if create_new_fig:
        plt.tight_layout()
    
    # Handle saving and displaying  
    if save_path:
        fig.savefig(save_path, dpi=250, bbox_inches="tight")
        print(f"Grid plot saved to: {save_path}", file=sys.stderr)
    
    # Always display the plot visually in addition to saving (only if we created the figure)
    if create_new_fig:
        plt.show()


def plot_grid_from_file(csv_file_path: str, *, 
                       title: Optional[str] = None,
                       save_path: Optional[str] = None):
    """
    Convenience function to plot grid directly from a CSV file.
    
    Args:
        csv_file_path: Path to the CSV file containing grid configuration
        title: Title for the plot
        save_path: Path to save the plot
    """
    try:
        df = pd.read_csv(csv_file_path)
        save_path_obj = Path(save_path) if save_path else None
        plot_grid(df, title=title or f"Grid from {Path(csv_file_path).name}", 
                 save_path=save_path_obj)
    except Exception as e:
        print(f"Error plotting grid from file {csv_file_path}: {e}", file=sys.stderr)