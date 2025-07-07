# MOVE ANALYSIS FIXED - CORRECTED VERSION
from __future__ import annotations
from pathlib import Path
from typing import List, Tuple, Dict, Any, Set, Optional
import math, colorsys, itertools, random, re, sys, time, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import ListedColormap, to_rgb, to_hex
from matplotlib.patches import Patch, Rectangle
from mpl_toolkits.mplot3d import Axes3D

BASE_COLOURS = {
    "NewMC": "#ffcc99", "DR": "#ffcccc", "NewMC_Expand": "#ccffcc",
    "NewMC_Restricted": "#ccccff", "NewMC_Biased": "#ffffcc",
    "NewMC_Restricted_Expand": "#ffccff", "NewMC_Biased_Expand": "#ccffff",
    "Kawasaki": "#ffcccc", "Kawasaki_Black": "#000000",
}

def get_distinct_colors(n):
    colors = ["#000000"]
    if n > 1:
        for i in range(n - 1):
            hue = (i / (n - 1)) * 0.8 + 0.1
            saturation = 0.8
            value = 0.75
            rgb = colorsys.hsv_to_rgb(hue, saturation, value)
            colors.append(to_hex(rgb))
    return colors

def moves_analysis_driver(csv_files, table_labels, curve_colours):
    """
    FIXED VERSION: Corrected acceptance rate normalization and scatter plot visualization
    """
    for csv_file, table_label, color in zip(csv_files, table_labels, curve_colours):
        if not Path(csv_file).is_file(): 
            print(f"File not found: {csv_file}")
            continue
        
        try:
            df = pd.read_csv(csv_file)
            required_cols = ["Total_Acceptances", "DeltaE", "NumAtomsSwapped"]
            if not all(col in df.columns for col in required_cols): 
                print(f"Missing required columns in {csv_file}")
                continue
            
            # Calculate acceptance status for each move
            df['IsAccepted'] = (df['Total_Acceptances'].diff().fillna(0) > 0).astype(int)
            
            # Calculate normalized acceptance rates by move type buckets
            # Group by number of atoms swapped to create buckets
            bucket_stats = df.groupby('NumAtomsSwapped').agg({
                'IsAccepted': ['sum', 'count'],
                'DeltaE': ['mean', 'std']
            }).round(4)
            
            # Flatten column names
            bucket_stats.columns = ['_'.join(col).strip() for col in bucket_stats.columns.values]
            bucket_stats['AcceptanceRate'] = (bucket_stats['IsAccepted_sum'] / bucket_stats['IsAccepted_count']).fillna(0)
            
            # Create figure with scatter plot as main focus
            fig, axes = plt.subplots(2, 2, figsize=(16, 12))
            fig.suptitle(f'Move Analysis: {table_label}', fontsize=16, fontweight='bold')
            
            # Main scatter plot: Number of atoms vs DeltaE with acceptance status
            accepted_mask = df['IsAccepted'] == 1
            rejected_mask = df['IsAccepted'] == 0
            
            # Use low DPI markers for better performance
            axes[0, 0].scatter(df[accepted_mask]['NumAtomsSwapped'], df[accepted_mask]['DeltaE'], 
                             c='green', alpha=0.6, s=8, label='Accepted', rasterized=True)
            axes[0, 0].scatter(df[rejected_mask]['NumAtomsSwapped'], df[rejected_mask]['DeltaE'], 
                             c='red', alpha=0.6, s=8, label='Rejected', rasterized=True)
            
            # Set appropriate axis limits with minimum ranges but allow expansion for data
            min_atoms = df['NumAtomsSwapped'].min()
            max_atoms = df['NumAtomsSwapped'].max()
            atom_range = max_atoms - min_atoms
            x_min = max(0, min_atoms - 0.1 * atom_range) if atom_range > 0 else 0
            x_max = max_atoms + 0.1 * atom_range if atom_range > 0 else max(10, max_atoms)
            
            min_deltaE = df['DeltaE'].min()
            max_deltaE = df['DeltaE'].max()
            deltaE_range = max_deltaE - min_deltaE
            y_min = min_deltaE - 0.1 * deltaE_range if deltaE_range > 0 else min_deltaE - 1
            y_max = max_deltaE + 0.1 * deltaE_range if deltaE_range > 0 else max_deltaE + 1
            
            axes[0, 0].set_xlim(x_min, x_max)
            axes[0, 0].set_ylim(y_min, y_max)
            axes[0, 0].set_xlabel('Number of Atoms Swapped')
            axes[0, 0].set_ylabel('ΔE (Energy Change)')
            axes[0, 0].set_title('Move Acceptance by Atoms Swapped vs Energy Change')
            axes[0, 0].legend()
            axes[0, 0].grid(True, alpha=0.3)
            
            # Normalized acceptance rates by move type buckets
            if len(bucket_stats) > 0:
                bucket_indices = bucket_stats.index
                acceptance_rates = bucket_stats['AcceptanceRate']
                
                axes[0, 1].bar(bucket_indices, acceptance_rates, color=color, alpha=0.7)
                axes[0, 1].set_xlabel('Number of Atoms Swapped (Move Type)')
                axes[0, 1].set_ylabel('Normalized Acceptance Rate')
                axes[0, 1].set_title('Acceptance Rate by Move Type Bucket')
                axes[0, 1].set_ylim(0, 1)
                axes[0, 1].grid(True, alpha=0.3)
                
                # Add count labels on bars
                for idx, rate in zip(bucket_indices, acceptance_rates):
                    count = bucket_stats.loc[idx, 'IsAccepted_count']
                    axes[0, 1].text(idx, rate + 0.02, f'n={count}', ha='center', va='bottom', fontsize=8)
            
            # Energy change distribution
            axes[1, 0].hist(df['DeltaE'], bins=50, alpha=0.7, color=color, edgecolor='black')
            axes[1, 0].set_xlabel('ΔE (Energy Change)')
            axes[1, 0].set_ylabel('Frequency')
            axes[1, 0].set_title('Energy Change Distribution')
            axes[1, 0].grid(True, alpha=0.3)
            
            # Running acceptance rate (properly normalized)
            iterations = np.arange(len(df))
            running_acceptance = df['Total_Acceptances'] / (iterations + 1)
            axes[1, 1].plot(iterations, running_acceptance, color=color, linewidth=2)
            axes[1, 1].set_xlabel('Iteration')
            axes[1, 1].set_ylabel('Cumulative Acceptance Rate')
            axes[1, 1].set_title('Running Acceptance Rate')
            axes[1, 1].set_ylim(0, 1)
            axes[1, 1].grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.show()
            
            # Print summary statistics
            print(f"\n=== Move Analysis Summary for {table_label} ===")
            print(f"Total moves: {len(df)}")
            print(f"Total accepted: {df['IsAccepted'].sum()}")
            print(f"Overall acceptance rate: {df['IsAccepted'].mean():.4f}")
            print(f"Move type buckets (by atoms swapped):")
            for idx, row in bucket_stats.iterrows():
                print(f"  {idx} atoms: {row['AcceptanceRate']:.4f} rate ({row['IsAccepted_sum']:.0f}/{row['IsAccepted_count']:.0f})")
            
        except Exception as e:
            print(f"Error processing {csv_file}: {e}")
            continue

def simple_batch_processor(*filenames, run_moves=True):
    """Simple batch processor for copy-paste filenames"""
    if not filenames: return
    
    csv_files = [f"{base_name}.csv" for base_name in filenames]
    table_labels = list(filenames)
    curve_colours = get_distinct_colors(len(csv_files))
    
    existing_files = [f for f in csv_files if Path(f).is_file()]
    existing_labels = [table_labels[i] for i, f in enumerate(csv_files) if Path(f).is_file()]
    existing_colours = [curve_colours[i] for i, f in enumerate(csv_files) if Path(f).is_file()]
    
    if not existing_files: 
        print("No existing files found!")
        return
    
    if run_moves:
        moves_analysis_driver(existing_files, existing_labels, existing_colours)

# USAGE EXAMPLE:
if __name__ == "__main__":
    # Example usage - replace with your actual filenames
    filenames_to_analyze = [
        "BJ=-0.44_Kawasaki_0_tol0.00_scp0.50",
        "BJ=-0.44_DR_1_tol0.50_scp0.50",
        "BJ=-0.44_DR_2_tol0.50_scp0.50"
    ]
    
    simple_batch_processor(*filenames_to_analyze, run_moves=True)