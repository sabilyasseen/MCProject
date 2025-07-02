"""
Cluster analysis module for Monte Carlo visualization suite.

This module provides functions for analyzing cluster sizes (boundary and interior)
over simulation iterations, with support for dynamic cluster detection and visualization.
"""

import sys
from pathlib import Path
from typing import List, Optional
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from color_utils import get_distinct_colors

def analyze_cluster_sizes(csv_files: List[str], *, 
                         table_labels: List[str], 
                         curve_colours: List[str], 
                         legend_labels: Optional[List[str]] = None, 
                         legend_colours: Optional[List[str]] = None):
    """
    Analyze and plot cluster boundary and interior sizes as a function of iterations.
    Creates individual charts for each file with dynamic cluster detection.
    
    Args:
        csv_files: List of CSV file paths
        table_labels: Labels for each CSV file
        curve_colours: Colors for each CSV file
        legend_labels: Optional legend labels
        legend_colours: Optional legend colors
    """
    
    print(f"\n--- Starting Cluster Size Analysis for {len(csv_files)} files ---", file=sys.stderr)
    
    plt.rcParams.update({'font.size': 10})
    
    files_processed = 0
    
    # Process each file individually
    for csv_file, table_label, curve_color in zip(csv_files, table_labels, curve_colours):
        file_path = Path(csv_file)
        if not file_path.is_file():
            print(f"[WARNING] Skipping {csv_file}: File not found", file=sys.stderr)
            continue
            
        try:
            # Read CSV and detect cluster columns
            df = pd.read_csv(file_path)
            
            # Find all cluster boundary and interior columns
            boundary_cols = [col for col in df.columns if col.startswith('BoundarySize_Cluster')]
            interior_cols = [col for col in df.columns if col.startswith('InteriorSize_Cluster')]
            
            # Extract cluster numbers and sort them
            cluster_numbers = []
            for col in boundary_cols:
                cluster_num = int(col.replace('BoundarySize_Cluster', ''))
                cluster_numbers.append(cluster_num)
            cluster_numbers = sorted(list(set(cluster_numbers)))
            
            num_clusters = len(cluster_numbers)
            if num_clusters == 0:
                print(f"[WARNING] No cluster data found in {csv_file}", file=sys.stderr)
                continue
                
            print(f"[INFO] Processing {csv_file}: Found {num_clusters} clusters", file=sys.stderr)
            
            # Create figure for this file
            fig = plt.figure(figsize=(16, 10))
            fig.suptitle(f'Cluster Size Analysis: {table_label}', fontsize=16, fontweight='bold')
            
            # Calculate grid layout
            if num_clusters <= 2:
                rows, cols = 1, num_clusters
            elif num_clusters <= 4:
                rows, cols = 2, 2
            elif num_clusters <= 6:
                rows, cols = 2, 3
            elif num_clusters <= 9:
                rows, cols = 3, 3
            else:
                rows, cols = 4, int(np.ceil(num_clusters / 4))
            
            # Get iterations (x-axis)
            iterations = np.arange(len(df))
            
            # Generate distinct colors for clusters
            cluster_colors = get_distinct_colors(num_clusters)
            
            # Plot each cluster
            for i, cluster_num in enumerate(cluster_numbers):
                ax = fig.add_subplot(rows, cols, i + 1)
                
                boundary_col = f'BoundarySize_Cluster{cluster_num}'
                interior_col = f'InteriorSize_Cluster{cluster_num}'
                
                if boundary_col in df.columns and interior_col in df.columns:
                    # Get data, handling NaN values
                    boundary_data = df[boundary_col].fillna(0)
                    interior_data = df[interior_col].fillna(0)
                    
                    # Plot boundary and interior sizes
                    ax.plot(iterations, boundary_data, 
                           label='Boundary', 
                           color=cluster_colors[i], 
                           linewidth=2, 
                           alpha=0.8)
                    
                    ax.plot(iterations, interior_data, 
                           label='Interior', 
                           color=cluster_colors[i], 
                           linewidth=2, 
                           linestyle='--', 
                           alpha=0.8)
                    
                    # Calculate total cluster size
                    total_size = boundary_data + interior_data
                    ax.plot(iterations, total_size, 
                           label='Total', 
                           color='gray', 
                           linewidth=1, 
                           alpha=0.6)
                    
                    # Formatting
                    ax.set_title(f'Cluster {cluster_num}', fontweight='bold')
                    ax.set_xlabel('Iteration')
                    ax.set_ylabel('Number of Cells')
                    ax.grid(True, alpha=0.3)
                    ax.legend(fontsize=8)
                    
                    # Add statistics text
                    final_boundary = boundary_data.iloc[-1] if len(boundary_data) > 0 else 0
                    final_interior = interior_data.iloc[-1] if len(interior_data) > 0 else 0
                    final_total = final_boundary + final_interior
                    
                    stats_text = f'Final: B={int(final_boundary)}, I={int(final_interior)}, T={int(final_total)}'
                    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
                           fontsize=8, verticalalignment='top',
                           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                else:
                    ax.text(0.5, 0.5, f'No data for\nCluster {cluster_num}', 
                           transform=ax.transAxes, ha='center', va='center')
                    ax.set_title(f'Cluster {cluster_num} (No Data)', fontweight='bold')
            
            # Remove empty subplots
            for i in range(len(cluster_numbers), rows * cols):
                fig.delaxes(fig.add_subplot(rows, cols, i + 1))
            
            plt.tight_layout()
            plt.show()
            
            # Print summary statistics
            print(f"\n📈 Summary for {table_label}:")
            print(f"   Total clusters: {num_clusters}")
            for cluster_num in cluster_numbers:
                boundary_col = f'BoundarySize_Cluster{cluster_num}'
                interior_col = f'InteriorSize_Cluster{cluster_num}'
                
                if boundary_col in df.columns and interior_col in df.columns:
                    boundary_data = df[boundary_col].fillna(0)
                    interior_data = df[interior_col].fillna(0)
                    
                    final_boundary = boundary_data.iloc[-1] if len(boundary_data) > 0 else 0
                    final_interior = interior_data.iloc[-1] if len(interior_data) > 0 else 0
                    
                    avg_boundary = boundary_data.mean()
                    avg_interior = interior_data.mean()
                    
                    print(f"   Cluster {cluster_num}: Final B={int(final_boundary)}, I={int(final_interior)} | "
                          f"Avg B={avg_boundary:.1f}, I={avg_interior:.1f}")
            
            files_processed += 1
            
        except Exception as e:
            print(f"[ERROR] Error processing {csv_file}: {e}", file=sys.stderr)
            continue
    
    print(f"\n--- Cluster Size Analysis Complete ---", file=sys.stderr)


def batch_cluster_analysis(*base_filenames):
    """
    Convenience function to analyze cluster sizes for multiple files.
    
    Args:
        base_filenames: Base names of CSV files (without .csv extension)
    """
    csv_files = [f"{base_name}.csv" for base_name in base_filenames]
    table_labels = list(base_filenames)
    
    # Generate colors for the files
    num_files = len(csv_files)
    curve_colours = get_distinct_colors(num_files)
    
    analyze_cluster_sizes(
        csv_files,
        table_labels=table_labels,
        curve_colours=curve_colours
    )