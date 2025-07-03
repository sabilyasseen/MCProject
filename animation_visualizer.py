import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle
import ast
import math

def parse_atom_list(atom_str):
    """Parse comma-separated atom indices from CSV string format"""
    if pd.isna(atom_str) or atom_str == '':
        return []
    try:
        # Remove quotes if present
        atom_str = atom_str.strip('"')
        if atom_str == '':
            return []
        return [int(x.strip()) for x in atom_str.split(',') if x.strip()]
    except:
        return []

def parse_grid_config(config_str, grid_size):
    """Parse grid configuration string into species and cluster arrays"""
    if pd.isna(config_str) or config_str == '':
        return np.zeros((grid_size, grid_size)), np.zeros((grid_size, grid_size))
    
    try:
        # Remove quotes if present
        config_str = config_str.strip('"')
        cells = config_str.split(';')
        
        species_grid = np.zeros((grid_size, grid_size))
        cluster_grid = np.zeros((grid_size, grid_size))
        
        for i, cell_data in enumerate(cells):
            if i >= grid_size * grid_size:
                break
            parts = cell_data.split('|')
            if len(parts) >= 2:
                row = i // grid_size
                col = i % grid_size
                species_grid[row, col] = int(parts[0])
                cluster_grid[row, col] = int(parts[1])
        
        return species_grid, cluster_grid
    except Exception as e:
        print(f"Warning: Failed to parse grid config: {e}")
        return np.zeros((grid_size, grid_size)), np.zeros((grid_size, grid_size))

def atoms_to_coordinates(atom_indices, grid_size):
    """Convert flat atom indices to (row, col) coordinates"""
    coords = []
    for atom in atom_indices:
        if atom >= 0 and atom < grid_size * grid_size:
            row = atom // grid_size
            col = atom % grid_size
            coords.append((row, col))
    return coords

def create_frame_data(iteration_data, grid_size, frame_type):
    """Create frame data for a specific frame type"""
    if frame_type == 'starting':
        species_grid, cluster_grid = parse_grid_config(iteration_data['GridConfigStarting'], grid_size)
        highlight_atoms = parse_atom_list(iteration_data['PossiblyChanged'])
        highlight_color = 'yellow'
        title = f"Iteration {iteration_data.name}: Starting State (Possibly Changed)"
    elif frame_type == 'potential':
        species_grid, cluster_grid = parse_grid_config(iteration_data['GridConfigPotential'], grid_size)
        highlight_atoms = parse_atom_list(iteration_data['ProposedAtoms'])
        highlight_color = 'orange'
        title = f"Iteration {iteration_data.name}: Potential State (Proposed)"
    else:  # final
        # Use potential grid config for final state
        species_grid, cluster_grid = parse_grid_config(iteration_data['GridConfigPotential'], grid_size)
        accepted_atoms = parse_atom_list(iteration_data['AcceptedAtoms'])
        proposed_atoms = parse_atom_list(iteration_data['ProposedAtoms'])
        
        # Rejected atoms are proposed atoms that weren't accepted
        rejected_atoms = [atom for atom in proposed_atoms if atom not in accepted_atoms]
        
        return species_grid, cluster_grid, accepted_atoms, rejected_atoms, f"Iteration {iteration_data.name}: Final State"
    
    return species_grid, cluster_grid, highlight_atoms, highlight_color, title

def create_animation_frame(ax, species_grid, cluster_grid, grid_size, iteration_num, 
                         highlight_atoms=None, highlight_color='yellow', 
                         accepted_atoms=None, rejected_atoms=None, title=""):
    """Create a single animation frame"""
    ax.clear()
    
    # Create base visualization using species
    cmap = plt.cm.RdYlBu_r
    im = ax.imshow(species_grid, cmap=cmap, vmin=0, vmax=1, interpolation='nearest')
    
    # Add highlights for different atom types
    if highlight_atoms:
        coords = atoms_to_coordinates(highlight_atoms, grid_size)
        for row, col in coords:
            rect = Rectangle((col-0.4, row-0.4), 0.8, 0.8, 
                           linewidth=2, edgecolor=highlight_color, facecolor='none')
            ax.add_patch(rect)
    
    if accepted_atoms:
        coords = atoms_to_coordinates(accepted_atoms, grid_size)
        for row, col in coords:
            rect = Rectangle((col-0.4, row-0.4), 0.8, 0.8, 
                           linewidth=3, edgecolor='green', facecolor='none')
            ax.add_patch(rect)
    
    if rejected_atoms:
        coords = atoms_to_coordinates(rejected_atoms, grid_size)
        for row, col in coords:
            rect = Rectangle((col-0.4, row-0.4), 0.8, 0.8, 
                           linewidth=3, edgecolor='red', facecolor='none')
            ax.add_patch(rect)
    
    ax.set_title(title, fontsize=12)
    ax.set_xlim(-0.5, grid_size-0.5)
    ax.set_ylim(-0.5, grid_size-0.5)
    ax.set_aspect('equal')
    
    # Remove ticks for cleaner look
    ax.set_xticks([])
    ax.set_yticks([])
    
    return im

def run_animation(base_filenames, batch_mode=True, input_dir='outputs', file_suffix='.csv'):
    """
    Create animations for simulation data showing atom state changes
    
    Args:
        base_filenames: List of base filenames to process
        batch_mode: Whether to run in batch mode
        input_dir: Directory containing input CSV files
        file_suffix: File suffix for input files
    """
    
    for base in base_filenames:
        try:
            # Construct full filename
            full_filename = os.path.join(input_dir, base + file_suffix)
            
            if not os.path.exists(full_filename):
                print(f"[ERROR] File not found: {full_filename}")
                continue
            
            print(f"[INFO] Processing {full_filename}")
            
            # Read CSV data
            try:
                df = pd.read_csv(full_filename)
            except Exception as e:
                print(f"[ERROR] Failed to read CSV file {full_filename}: {e}")
                continue
            
            if df.empty:
                print(f"[WARNING] Empty CSV file: {full_filename}")
                continue
            
            # Determine grid size from first valid grid configuration
            grid_size = None
            for idx, row in df.iterrows():
                config_str = row.get('GridConfigStarting', '')
                if pd.notna(config_str) and config_str != '':
                    config_str = config_str.strip('"')
                    cells = config_str.split(';')
                    total_cells = len(cells)
                    grid_size = int(math.sqrt(total_cells))
                    break
            
            if grid_size is None:
                print(f"[ERROR] Could not determine grid size from {full_filename}")
                continue
            
            print(f"[INFO] Detected grid size: {grid_size}x{grid_size}")
            
            # Limit iterations for performance (can be adjusted)
            max_iterations = min(len(df), 100)
            df_subset = df.head(max_iterations)
            
            # Create figure and animation
            fig, ax = plt.subplots(figsize=(10, 10))
            
            frames_data = []
            
            # Prepare frame data (3 frames per iteration)
            for idx, row in df_subset.iterrows():
                # Frame 1: Starting state with possibly changed atoms
                try:
                    species_grid, cluster_grid, highlight_atoms, highlight_color, title = create_frame_data(row, grid_size, 'starting')
                    frames_data.append({
                        'species_grid': species_grid,
                        'cluster_grid': cluster_grid,
                        'highlight_atoms': highlight_atoms,
                        'highlight_color': highlight_color,
                        'title': title,
                        'type': 'starting'
                    })
                except Exception as e:
                    print(f"[WARNING] Failed to create starting frame for iteration {idx}: {e}")
                
                # Frame 2: Potential state with proposed atoms
                try:
                    species_grid, cluster_grid, highlight_atoms, highlight_color, title = create_frame_data(row, grid_size, 'potential')
                    frames_data.append({
                        'species_grid': species_grid,
                        'cluster_grid': cluster_grid,
                        'highlight_atoms': highlight_atoms,
                        'highlight_color': highlight_color,
                        'title': title,
                        'type': 'potential'
                    })
                except Exception as e:
                    print(f"[WARNING] Failed to create potential frame for iteration {idx}: {e}")
                
                # Frame 3: Final state with accepted (green) and rejected (red) atoms
                try:
                    species_grid, cluster_grid, accepted_atoms, rejected_atoms, title = create_frame_data(row, grid_size, 'final')
                    frames_data.append({
                        'species_grid': species_grid,
                        'cluster_grid': cluster_grid,
                        'accepted_atoms': accepted_atoms,
                        'rejected_atoms': rejected_atoms,
                        'title': title,
                        'type': 'final'
                    })
                except Exception as e:
                    print(f"[WARNING] Failed to create final frame for iteration {idx}: {e}")
            
            if not frames_data:
                print(f"[ERROR] No valid frames created for {full_filename}")
                continue
            
            def animate(frame_idx):
                frame_data = frames_data[frame_idx]
                
                if frame_data['type'] == 'final':
                    return create_animation_frame(
                        ax, frame_data['species_grid'], frame_data['cluster_grid'], 
                        grid_size, frame_idx // 3,
                        accepted_atoms=frame_data.get('accepted_atoms'),
                        rejected_atoms=frame_data.get('rejected_atoms'),
                        title=frame_data['title']
                    )
                else:
                    return create_animation_frame(
                        ax, frame_data['species_grid'], frame_data['cluster_grid'],
                        grid_size, frame_idx // 3,
                        highlight_atoms=frame_data.get('highlight_atoms'),
                        highlight_color=frame_data.get('highlight_color'),
                        title=frame_data['title']
                    )
            
            # Create animation
            anim = animation.FuncAnimation(
                fig, animate, frames=len(frames_data), 
                interval=1000, repeat=True, blit=False
            )
            
            # Save animation as MP4
            output_filename = os.path.join(input_dir, f"{base}_animation.mp4")
            
            try:
                Writer = animation.writers['ffmpeg']
                writer = Writer(fps=1, metadata=dict(artist='MC Simulation'), bitrate=1800)
                anim.save(output_filename, writer=writer)
                print(f"[SUCCESS] Animation saved as: {output_filename}")
            except Exception as e:
                print(f"[ERROR] Failed to save animation as MP4: {e}")
                # Try saving as GIF as fallback
                try:
                    gif_filename = os.path.join(input_dir, f"{base}_animation.gif")
                    anim.save(gif_filename, writer='pillow', fps=1)
                    print(f"[SUCCESS] Animation saved as GIF: {gif_filename}")
                except Exception as e2:
                    print(f"[ERROR] Failed to save animation as GIF: {e2}")
            
            plt.close(fig)
            
        except Exception as e:
            print(f"[ERROR] Failed to process {base}: {e}")
            continue
    
    print("[INFO] Animation generation complete")

if __name__ == "__main__":
    # Example usage
    if len(sys.argv) > 1:
        base_filenames = sys.argv[1:]
    else:
        # Default test case
        base_filenames = ["test_simulation"]
    
    run_animation(base_filenames)