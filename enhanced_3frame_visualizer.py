#!/usr/bin/env python3
"""
Ultra-Optimized 3-Frame Ising Model Animation System
=======================================================

This system creates animations with 3 frames per Monte Carlo iteration:
1. Plain Grid: Current state of the lattice
2. Yellow Outline: Highlighting possible move positions  
3. Accept/Reject Overlay: Green for accepted moves, Orange for rejected moves
"""

from __future__ import annotations
from pathlib import Path
from typing import List, Tuple, Dict, Any, Set, Optional
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle
from matplotlib.colors import ListedColormap
import sys
import time
import os
from dataclasses import dataclass

# Set matplotlib backend for better performance
plt.switch_backend('Agg')

@dataclass
class IterationData:
    """Container for all data needed for one iteration's 3-frame sequence"""
    grid_starting: np.ndarray
    grid_potential: Optional[np.ndarray]
    possibly_changed: List[Tuple[int, int]]
    proposed_atoms: List[Tuple[int, int]]
    accepted_atoms: List[Tuple[int, int]]
    is_accepted: bool
    iteration_num: int

def parse_atom_indices_to_coords(indices_raw: any, grid_size: int) -> List[Tuple[int, int]]:
    """Parse atom indices string to (x, y) coordinates with error handling"""
    if pd.isna(indices_raw) or str(indices_raw).strip() in ["", '""', "''"]:
        return []
    
    s = str(indices_raw).strip('"\\\' ')
    if not s:
        return []
    
    try:
        # Handle both comma and semicolon separators
        separators = [',', ';']
        for sep in separators:
            if sep in s:
                parts = s.split(sep)
                break
        else:
            parts = [s]
            
        indices = []
        for x in parts:
            x = x.strip()
            if x:
                try:
                    idx = int(x)
                    if 0 <= idx < grid_size * grid_size:
                        indices.append(idx)
                except ValueError:
                    continue
        
        return [(idx % grid_size, idx // grid_size) for idx in indices]
    except Exception:
        return []

def parse_grid_config_to_species_array(config_raw: any, N: int) -> Optional[np.ndarray]:
    """Parse grid configuration string to species array"""
    if pd.isna(config_raw) or str(config_raw).strip() in ["", '""', "''"]:
        return None
    
    s = str(config_raw).strip('"\\\' ')
    if not s:
        return None
        
    parts = s.split(';')
    if len(parts) != N * N:
        return None
        
    species = np.zeros((N, N), dtype=np.int8)
    try:
        for i, cell_data_str in enumerate(parts):
            if not cell_data_str.strip():
                continue
            cell_parts = cell_data_str.split('|')
            if len(cell_parts) == 0:
                continue
            species_val = int(cell_parts[0])
            if species_val not in [0, 1]:
                species_val = 0  # Default to 0 for invalid values
            y, x = i // N, i % N
            species[y, x] = species_val
    except (ValueError, IndexError):
        return None
    return species

def load_and_process_csv_data(csv_file: Path, start_iter: int, end_iter: int) -> List[IterationData]:
    """Load CSV data and process into IterationData objects"""
    try:
        # First, determine grid size
        df_sample = pd.read_csv(csv_file, dtype=str, nrows=1)
        if df_sample.empty or 'GridConfigStarting' not in df_sample.columns:
            return []
            
        config_str = df_sample["GridConfigStarting"].iloc[0]
        num_cells = len(str(config_str).strip('"\\\' ').split(';'))
        N = int(math.sqrt(num_cells))
        
        if N * N != num_cells:
            return []
            
        # Load the iteration range
        nrows_to_read = end_iter - start_iter + 1
        if nrows_to_read <= 0:
            return []
            
        df = pd.read_csv(
            csv_file,
            dtype=str,
            skiprows=range(1, start_iter + 1),
            nrows=nrows_to_read
        )
        
        if df.empty:
            return []
            
        # Required columns
        required_cols = ['GridConfigStarting', 'GridConfigPotential', 'PossiblyChanged', 
                        'ProposedAtoms', 'AcceptedAtoms']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            return []
            
        # Process each iteration
        iterations = []
        for idx, row in df.iterrows():
            # Parse grid states
            grid_starting = parse_grid_config_to_species_array(row['GridConfigStarting'], N)
            grid_potential = parse_grid_config_to_species_array(row['GridConfigPotential'], N)
            
            if grid_starting is None:
                continue
                
            # Parse atom coordinates
            possibly_changed = parse_atom_indices_to_coords(row['PossiblyChanged'], N)
            proposed_atoms = parse_atom_indices_to_coords(row['ProposedAtoms'], N)
            accepted_atoms = parse_atom_indices_to_coords(row['AcceptedAtoms'], N)
            
            # Determine if move was accepted
            is_accepted = len(accepted_atoms) > 0
            
            iteration_data = IterationData(
                grid_starting=grid_starting,
                grid_potential=grid_potential,
                possibly_changed=possibly_changed,
                proposed_atoms=proposed_atoms,
                accepted_atoms=accepted_atoms,
                is_accepted=is_accepted,
                iteration_num=start_iter + len(iterations)
            )
            iterations.append(iteration_data)
            
        return iterations
        
    except Exception as e:
        return []

def create_ultra_optimized_3frame_animation(iterations: List[IterationData], 
                                          output_filename: str, 
                                          fps: int = 15,
                                          dpi: int = 100):
    """Create 3-frame animation with maximum optimization"""
    if not iterations:
        return
        
    N = iterations[0].grid_starting.shape[0]
    total_frames = len(iterations) * 3  # 3 frames per iteration
    
    # Setup figure with optimized parameters
    fig, ax = plt.subplots(figsize=(10, 10), dpi=dpi)
    ax.set_xlim(-0.5, N - 0.5)
    ax.set_ylim(-0.5, N - 0.5)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Ising Model Monte Carlo Simulation', fontsize=16, pad=20, weight='bold')
    
    # Invert y-axis to match matrix indexing
    ax.invert_yaxis()
    plt.tight_layout()
    
    # Color scheme - vibrant colors for better visibility
    species_cmap = ListedColormap(['#FF6B6B', '#4ECDC4'])  # Coral Red=0, Teal=1
    
    # Pre-create all rectangle patches for maximum efficiency
    grid_patches = []
    for y in range(N):
        for x in range(N):
            rect = Rectangle((x - 0.5, y - 0.5), 1, 1, 
                           facecolor=species_cmap(0), 
                           edgecolor='black', 
                           linewidth=0.5)
            ax.add_patch(rect)
            grid_patches.append(rect)
    
    # Overlay patches for highlighting (created on demand)
    overlay_patches = []
    status_text = ax.text(0.5, 1.08, '', transform=ax.transAxes, 
                         ha='center', va='bottom', fontsize=14, weight='bold',
                         bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))
    
    def clear_overlays():
        """Remove all overlay patches"""
        for patch in overlay_patches:
            patch.remove()
        overlay_patches.clear()
    
    def update_grid_colors(grid_state: np.ndarray):
        """Update grid colors efficiently"""
        flat_grid = grid_state.flatten()
        for i, patch in enumerate(grid_patches):
            patch.set_facecolor(species_cmap(flat_grid[i]))
    
    def add_outline_overlay(coords: List[Tuple[int, int]], color: str, linewidth: int = 6):
        """Add colored outline overlay"""
        for x, y in coords:
            if 0 <= x < N and 0 <= y < N:
                rect = Rectangle((x - 0.5, y - 0.5), 1, 1,
                               facecolor='none',
                               edgecolor=color,
                               linewidth=linewidth,
                               zorder=10)
                ax.add_patch(rect)
                overlay_patches.append(rect)
    
    def add_fill_overlay(coords: List[Tuple[int, int]], color: str):
        """Add colored fill overlay (no transparency)"""
        for x, y in coords:
            if 0 <= x < N and 0 <= y < N:
                rect = Rectangle((x - 0.5, y - 0.5), 1, 1,
                               facecolor=color,
                               edgecolor='black',
                               linewidth=2,
                               zorder=5)
                ax.add_patch(rect)
                overlay_patches.append(rect)
    
    def animate(frame_idx: int):
        """Animation function for 3-frame sequence"""
        iteration_idx = frame_idx // 3
        frame_type = frame_idx % 3
        
        if iteration_idx >= len(iterations):
            return grid_patches + overlay_patches + [status_text]
            
        iteration = iterations[iteration_idx]
        
        # Clear previous overlays
        clear_overlays()
        
        if frame_type == 0:
            # Frame 1: Plain grid (starting state)
            update_grid_colors(iteration.grid_starting)
            status_text.set_text(f"Iteration {iteration.iteration_num:04d} - Starting State")
            status_text.set_color('black')
            
        elif frame_type == 1:
            # Frame 2: Grid with yellow outline showing possible moves
            update_grid_colors(iteration.grid_starting)
            if iteration.possibly_changed:
                add_outline_overlay(iteration.possibly_changed, '#FFD700', linewidth=8)  # Gold/Yellow
            status_text.set_text(f"Iteration {iteration.iteration_num:04d} - Possible Move Highlighted")
            status_text.set_color('#FF8C00')  # Dark Orange
            
        elif frame_type == 2:
            # Frame 3: Result with accept/reject overlay
            if iteration.is_accepted and iteration.grid_potential is not None:
                # Show final state and green overlay for accepted moves
                update_grid_colors(iteration.grid_potential)
                if iteration.accepted_atoms:
                    add_fill_overlay(iteration.accepted_atoms, '#32CD32')  # LimeGreen
                status_text.set_text(f"Iteration {iteration.iteration_num:04d} - [ACCEPTED]")
                status_text.set_color('#228B22')  # Forest Green
            else:
                # Show starting state and orange overlay for rejected moves
                update_grid_colors(iteration.grid_starting)
                if iteration.proposed_atoms:
                    add_fill_overlay(iteration.proposed_atoms, '#FF4500')  # OrangeRed
                status_text.set_text(f"Iteration {iteration.iteration_num:04d} - [REJECTED]")
                status_text.set_color('#DC143C')  # Crimson
        
        return grid_patches + overlay_patches + [status_text]
    
    # Create animation with optimal settings
    try:
        # Check if ffmpeg is available
        if 'ffmpeg' not in animation.writers.avail:
            writer = animation.PillowWriter(fps=fps)
            output_filename = output_filename.replace('.mp4', '.gif')
        else:
            # Use the fastest possible ffmpeg settings
            Writer = animation.writers['ffmpeg']
            writer = Writer(
                fps=fps,
                metadata=dict(
                    title='3-Frame Ising Model Animation',
                    artist='Ultra-Optimized Visualization System',
                    comment=f'{len(iterations)} iterations, {total_frames} frames'
                ),
                bitrate=-1,  # Variable bitrate for best quality/speed balance
                extra_args=[
                    '-preset', 'ultrafast',       # Fastest encoding preset
                    '-crf', '18',                # High quality (lower is better)
                    '-pix_fmt', 'yuv420p',       # Compatibility format
                    '-movflags', '+faststart',   # Fast web streaming
                    '-tune', 'animation',        # Optimize for animation content
                    '-threads', '0',             # Use all available CPU threads
                    '-loglevel', 'error'         # Reduce ffmpeg output
                ]
            )
        
        # Create animation
        anim = animation.FuncAnimation(
            fig, animate, frames=total_frames, 
            interval=1000//fps, blit=False, repeat=True
        )
        
        # Save animation
        start_time = time.time()
        anim.save(output_filename, writer=writer, dpi=dpi)
        save_time = time.time() - start_time
        
        if os.path.exists(output_filename):
            file_size = os.path.getsize(output_filename) / 1024 / 1024
            return {
                'success': True,
                'filename': output_filename,
                'size_mb': file_size,
                'render_time': save_time,
                'total_frames': total_frames,
                'duration': total_frames/fps
            }
        else:
            return {'success': False, 'error': 'Failed to save animation file'}
            
    except Exception as e:
        return {'success': False, 'error': str(e)}
    finally:
        plt.close(fig)

def batch_process_3frame_animations(
    *base_filenames: str,
    start_iteration: int = 0,
    end_iteration: int = 50,
    fps: int = 15,
    dpi: int = 100,
    output_dir: str = "outputs"
):
    """Process multiple CSV files into 3-frame animations"""
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    results = []
    
    for base_name in base_filenames:
        csv_file = Path(f"{base_name}.csv")
        if not csv_file.exists():
            results.append({'filename': base_name, 'success': False, 'error': 'File not found'})
            continue
            
        # Load and process data
        iterations = load_and_process_csv_data(csv_file, start_iteration, end_iteration)
        if not iterations:
            results.append({'filename': base_name, 'success': False, 'error': 'No valid data'})
            continue
        
        # Create animation
        output_filename = output_path / f"{base_name}_3frame_iter_{start_iteration}_{end_iteration}.mp4"
        result = create_ultra_optimized_3frame_animation(
            iterations, str(output_filename), fps=fps, dpi=dpi
        )
        result['filename'] = base_name
        results.append(result)
    
    return results