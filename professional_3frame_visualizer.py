#!/usr/bin/env python3
"""
Professional 3-Frame Ising Model Animation System
=================================================

This system creates high-performance animations with a structured 3-frame sequence
per Monte Carlo iteration for comprehensive move visualization.

Frame Structure:
1. Base State: Current lattice configuration
2. Proposed Move: Highlighted positions showing potential changes (yellow outline)  
3. Result State: Final outcome with color-coded acceptance status
   - Green overlay: Accepted moves
   - Orange overlay: Rejected moves

Technical Features:
- Maximum rendering efficiency with optimized FFmpeg settings
- No transparency overlays for faster processing
- Efficient memory management with on-demand patch creation
- Professional error handling and status reporting
- Comprehensive CSV data validation

Encoding Specifications:
- H.264 with ultrafast preset for maximum speed
- CRF 18 for high visual quality
- Multi-threaded processing with animation tuning
- YUV420P pixel format for universal compatibility
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

# Configure matplotlib for optimal performance
plt.switch_backend('Agg')
plt.rcParams.update({'font.size': 12})

@dataclass
class IterationData:
    """Container for iteration data required for 3-frame sequence generation"""
    grid_starting: np.ndarray
    grid_potential: Optional[np.ndarray]
    possibly_changed: List[Tuple[int, int]]
    proposed_atoms: List[Tuple[int, int]]
    accepted_atoms: List[Tuple[int, int]]
    is_accepted: bool
    iteration_num: int

def parse_atom_indices_to_coords(indices_raw: any, grid_size: int) -> List[Tuple[int, int]]:
    """Convert atom indices string to coordinate pairs with robust error handling"""
    if pd.isna(indices_raw) or str(indices_raw).strip() in ["", '""', "''"]:
        return []
    
    s = str(indices_raw).strip('"\\\' ')
    if not s:
        return []
    
    try:
        # Handle multiple separator types
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
    """Parse grid configuration string to species array with validation"""
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
                species_val = 0  # Default fallback for invalid values
            y, x = i // N, i % N
            species[y, x] = species_val
    except (ValueError, IndexError):
        return None
    return species

def load_and_process_csv_data(csv_file: Path, start_iter: int, end_iter: int) -> List[IterationData]:
    """Load and process CSV data into structured iteration objects"""
    print(f"Loading simulation data from {csv_file.name}")
    print(f"Processing iterations {start_iter} to {end_iter}")
    
    try:
        # Determine grid dimensions from first row
        df_sample = pd.read_csv(csv_file, dtype=str, nrows=1)
        if df_sample.empty or 'GridConfigStarting' not in df_sample.columns:
            print(f"Error: Invalid CSV structure in {csv_file.name}")
            return []
            
        config_str = df_sample["GridConfigStarting"].iloc[0]
        num_cells = len(str(config_str).strip('"\\\' ').split(';'))
        N = int(math.sqrt(num_cells))
        
        if N * N != num_cells:
            print(f"Error: Grid is not square. Detected {num_cells} cells")
            return []
            
        print(f"Detected grid dimensions: {N}x{N}")
            
        # Load specified iteration range
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
            print(f"Warning: No data found in specified iteration range")
            return []
            
        # Validate required columns
        required_cols = ['GridConfigStarting', 'GridConfigPotential', 'PossiblyChanged', 
                        'ProposedAtoms', 'AcceptedAtoms']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"Error: Missing required columns: {missing_cols}")
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
            
            # Determine acceptance status
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
            
        print(f"Successfully processed {len(iterations)} iterations")
        return iterations
        
    except Exception as e:
        print(f"Error during CSV processing: {e}")
        return []

def create_professional_3frame_animation(
    iterations: List[IterationData], 
    output_filename: str, 
    fps: int = 12,
    dpi: int = 120
):
    """Create optimized 3-frame animation with professional quality settings"""
    if not iterations:
        print("Error: No iterations available for animation")
        return
        
    N = iterations[0].grid_starting.shape[0]
    total_frames = len(iterations) * 3  # 3 frames per iteration
    
    print(f"Initializing animation system")
    print(f"Grid size: {N}x{N}")
    print(f"Iterations: {len(iterations)}")
    print(f"Total frames: {total_frames}")
    
    # Setup figure with professional parameters
    fig, ax = plt.subplots(figsize=(10, 10), dpi=dpi)
    ax.set_xlim(-0.5, N - 0.5)
    ax.set_ylim(-0.5, N - 0.5)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Monte Carlo Simulation Progress', fontsize=16, pad=20, weight='bold')
    
    # Invert y-axis for proper matrix visualization
    ax.invert_yaxis()
    plt.tight_layout()
    
    # Professional color scheme
    species_cmap = ListedColormap(['#E74C3C', '#3498DB'])  # Red=0, Blue=1
    
    # Pre-create grid patches for optimal performance
    grid_patches = []
    for y in range(N):
        for x in range(N):
            rect = Rectangle((x - 0.5, y - 0.5), 1, 1, 
                           facecolor=species_cmap(0), 
                           edgecolor='black', 
                           linewidth=0.5)
            ax.add_patch(rect)
            grid_patches.append(rect)
    
    # Dynamic overlay management
    overlay_patches = []
    status_text = ax.text(0.5, 1.08, '', transform=ax.transAxes, 
                         ha='center', va='bottom', fontsize=14, weight='bold',
                         bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.9))
    
    def clear_overlays():
        """Remove all overlay patches for clean frame transitions"""
        for patch in overlay_patches:
            patch.remove()
        overlay_patches.clear()
    
    def update_grid_colors(grid_state: np.ndarray):
        """Update grid colors with optimized batch processing"""
        flat_grid = grid_state.flatten()
        for i, patch in enumerate(grid_patches):
            patch.set_facecolor(species_cmap(flat_grid[i]))
    
    def add_outline_overlay(coords: List[Tuple[int, int]], color: str, linewidth: int = 6):
        """Add colored outline overlay for highlighting"""
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
        """Add solid fill overlay for result indication"""
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
        """Animation function implementing 3-frame sequence logic"""
        iteration_idx = frame_idx // 3
        frame_type = frame_idx % 3
        
        if iteration_idx >= len(iterations):
            return grid_patches + overlay_patches + [status_text]
            
        iteration = iterations[iteration_idx]
        
        # Clear previous overlays
        clear_overlays()
        
        if frame_type == 0:
            # Frame 1: Base state
            update_grid_colors(iteration.grid_starting)
            status_text.set_text(f"Iteration {iteration.iteration_num:04d} - Initial State")
            status_text.set_color('black')
            
        elif frame_type == 1:
            # Frame 2: Proposed move highlighting
            update_grid_colors(iteration.grid_starting)
            if iteration.possibly_changed:
                add_outline_overlay(iteration.possibly_changed, '#FFD700', linewidth=8)
            status_text.set_text(f"Iteration {iteration.iteration_num:04d} - Proposed Move")
            status_text.set_color('#FF8C00')
            
        elif frame_type == 2:
            # Frame 3: Result state with acceptance indication
            if iteration.is_accepted and iteration.grid_potential is not None:
                # Show accepted state
                update_grid_colors(iteration.grid_potential)
                if iteration.accepted_atoms:
                    add_fill_overlay(iteration.accepted_atoms, '#2ECC71')
                status_text.set_text(f"Iteration {iteration.iteration_num:04d} - ACCEPTED")
                status_text.set_color('#27AE60')
            else:
                # Show rejected state
                update_grid_colors(iteration.grid_starting)
                if iteration.proposed_atoms:
                    add_fill_overlay(iteration.proposed_atoms, '#E67E22')
                status_text.set_text(f"Iteration {iteration.iteration_num:04d} - REJECTED")
                status_text.set_color('#C0392B')
        
        # Progress reporting
        progress = (frame_idx + 1) / total_frames * 100
        if frame_idx % max(1, total_frames // 20) == 0:
            print(f"Animation progress: {progress:.1f}%")
        
        return grid_patches + overlay_patches + [status_text]
    
    # Create animation with professional encoding settings
    try:
        # Check FFmpeg availability
        if 'ffmpeg' not in animation.writers.avail:
            print("Warning: FFmpeg not available, using fallback writer")
            writer = animation.PillowWriter(fps=fps)
            output_filename = output_filename.replace('.mp4', '.gif')
        else:
            # Professional FFmpeg configuration
            Writer = animation.writers['ffmpeg']
            writer = Writer(
                fps=fps,
                metadata=dict(
                    title='Professional 3-Frame Ising Model Animation',
                    artist='Monte Carlo Visualization System',
                    comment=f'{len(iterations)} iterations processed'
                ),
                bitrate=-1,
                extra_args=[
                    '-preset', 'ultrafast',
                    '-crf', '18',
                    '-pix_fmt', 'yuv420p',
                    '-movflags', '+faststart',
                    '-tune', 'animation',
                    '-threads', '0',
                    '-loglevel', 'error'
                ]
            )
        
        # Generate animation
        print(f"Configuring animation at {fps} FPS")
        anim = animation.FuncAnimation(
            fig, animate, frames=total_frames, 
            interval=1000//fps, blit=False, repeat=True
        )
        
        # Save with performance monitoring
        print(f"Rendering animation to {output_filename}")
        start_time = time.time()
        anim.save(output_filename, writer=writer, dpi=dpi)
        render_time = time.time() - start_time
        
        if os.path.exists(output_filename):
            file_size = os.path.getsize(output_filename) / 1024 / 1024
            print(f"Animation completed successfully")
            print(f"Output file: {output_filename}")
            print(f"File size: {file_size:.1f} MB")
            print(f"Render time: {render_time:.1f} seconds")
            print(f"Animation duration: {total_frames/fps:.1f} seconds")
        else:
            print("Error: Animation file was not created")
            
    except Exception as e:
        print(f"Error during animation creation: {e}")
        print("Ensure FFmpeg is properly installed and accessible")
        import traceback
        traceback.print_exc()
    finally:
        plt.close(fig)

def batch_process_professional_animations(
    *base_filenames: str,
    start_iteration: int = 0,
    end_iteration: int = 50,
    fps: int = 12,
    dpi: int = 120,
    output_dir: str = "outputs"
):
    """Process multiple simulation files with professional batch processing"""
    print("Starting Professional Animation Batch Processing")
    print(f"Configuration: Iterations {start_iteration}-{end_iteration}, {fps} FPS, {dpi} DPI")
    print(f"Output directory: {output_dir}")
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    processing_results = []
    
    for base_name in base_filenames:
        print(f"\nProcessing simulation: {base_name}")
        
        csv_file = Path(f"{base_name}.csv")
        if not csv_file.exists():
            error_msg = f"CSV file not found: {csv_file}"
            print(f"Error: {error_msg}")
            processing_results.append(f"FAILED - {base_name}: File not found")
            continue
            
        # Load and validate data
        iterations = load_and_process_csv_data(csv_file, start_iteration, end_iteration)
        if not iterations:
            error_msg = f"No valid iterations found in {base_name}"
            print(f"Error: {error_msg}")
            processing_results.append(f"FAILED - {base_name}: No valid data")
            continue
        
        # Generate animation
        output_filename = output_path / f"{base_name}_professional_3frame_{start_iteration}_{end_iteration}.mp4"
        create_professional_3frame_animation(
            iterations, str(output_filename), fps=fps, dpi=dpi
        )
        processing_results.append(f"SUCCESS - {base_name}: Animation created")
    
    print("\nBatch Processing Complete")
    print("Processing Summary:")
    for result in processing_results:
        print(f"  {result}")

def main():
    """Main execution function for demonstration"""
    print("Professional 3-Frame Ising Model Animation System")
    print("=" * 60)
    
    # Example configuration
    test_files = [
        "BJ=-0.44_Kawasaki_0_tol0.00_scp0.50",
        "BJ=-0.44_DR_1_tol0.50_scp0.50"
    ]
    
    # Professional settings
    config = {
        'start_iteration': 0,
        'end_iteration': 30,
        'fps': 10,
        'dpi': 120,
        'output_dir': "outputs"
    }
    
    print("Configuration parameters:")
    for key, value in config.items():
        print(f"  {key}: {value}")
    
    batch_process_professional_animations(
        *test_files,
        **config
    )

if __name__ == '__main__':
    main()