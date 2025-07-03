import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def run_energy_analysis(base_filenames, batch_mode=True, input_dir='outputs', file_suffix='.csv'):
    """
    Enhanced energy analysis with multiple plot types:
    - Regular energy (linear and log x-axis) with all files overlayed
    - Rolling average energy (linear and log x-axis) with all files overlayed
    - Change in rolling average at each iteration (log scale)
    - Change in rolling average slope over 1000 iterations (log scale)
    """
    
    # Prepare data containers
    all_data = {}
    colors = plt.cm.tab10(np.linspace(0, 1, len(base_filenames)))
    
    # Load all data files
    for idx, base in enumerate(base_filenames):
        file_path = os.path.join(input_dir, base + file_suffix)
        if not os.path.exists(file_path):
            # Try without input_dir (direct filename)
            file_path = base + file_suffix
            if not os.path.exists(file_path):
                print(f"[ERROR] File not found: {base + file_suffix}", file=sys.stderr)
                continue
        
        try:
            df = pd.read_csv(file_path)
            if 'total_energy' in df.columns:
                energy_col = 'total_energy'
            elif 'Energy' in df.columns:
                energy_col = 'Energy'
            elif 'energy' in df.columns:
                energy_col = 'energy'
            else:
                print(f"[ERROR] No energy column found in {file_path}", file=sys.stderr)
                continue
                
            energy_data = df[energy_col].ffill().bfill()
            all_data[base] = {
                'energy': energy_data,
                'iterations': np.arange(len(energy_data)),
                'color': colors[idx]
            }
            
        except Exception as e:
            print(f"[ERROR] Failed to process {file_path}: {e}", file=sys.stderr)
    
    if not all_data:
        print("[ERROR] No valid data files found", file=sys.stderr)
        return
    
    # Create the comprehensive plot layout
    fig = plt.figure(figsize=(20, 16))
    fig.suptitle('Comprehensive Energy Analysis', fontsize=16, fontweight='bold')
    
    # Plot 1: Regular Energy - Linear X-axis
    ax1 = fig.add_subplot(3, 2, 1)
    for base, data in all_data.items():
        ax1.plot(data['iterations'], data['energy'], 
                label=base, color=data['color'], linewidth=1.5, alpha=0.8)
    ax1.set_xlabel('Iteration')
    ax1.set_ylabel('Energy')
    ax1.set_title('Energy vs Iteration (Linear Scale)')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Regular Energy - Log X-axis
    ax2 = fig.add_subplot(3, 2, 2)
    for base, data in all_data.items():
        ax2.semilogx(data['iterations'][1:], data['energy'][1:], 
                    label=base, color=data['color'], linewidth=1.5, alpha=0.8)
    ax2.set_xlabel('Iteration (Log Scale)')
    ax2.set_ylabel('Energy')
    ax2.set_title('Energy vs Iteration (Log X-axis)')
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Rolling Average Energy - Linear X-axis
    ax3 = fig.add_subplot(3, 2, 3)
    for base, data in all_data.items():
        rolling_avg = data['energy'].rolling(window=min(100, len(data['energy'])//10), min_periods=1).mean()
        ax3.plot(data['iterations'], rolling_avg, 
                label=base, color=data['color'], linewidth=2, alpha=0.9)
    ax3.set_xlabel('Iteration')
    ax3.set_ylabel('Rolling Average Energy')
    ax3.set_title('Rolling Average Energy (Linear Scale)')
    ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Rolling Average Energy - Log X-axis
    ax4 = fig.add_subplot(3, 2, 4)
    for base, data in all_data.items():
        rolling_avg = data['energy'].rolling(window=min(100, len(data['energy'])//10), min_periods=1).mean()
        ax4.semilogx(data['iterations'][1:], rolling_avg[1:], 
                    label=base, color=data['color'], linewidth=2, alpha=0.9)
    ax4.set_xlabel('Iteration (Log Scale)')
    ax4.set_ylabel('Rolling Average Energy')
    ax4.set_title('Rolling Average Energy (Log X-axis)')
    ax4.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax4.grid(True, alpha=0.3)
    
    # Plot 5: Change in Rolling Average at Each Iteration (Log Scale)
    ax5 = fig.add_subplot(3, 2, 5)
    for base, data in all_data.items():
        rolling_avg = data['energy'].rolling(window=min(100, len(data['energy'])//10), min_periods=1).mean()
        change = np.abs(np.diff(rolling_avg))
        change[change == 0] = 1e-10  # Avoid log(0)
        ax5.semilogy(data['iterations'][1:], change, 
                    label=base, color=data['color'], linewidth=1.5, alpha=0.8)
    ax5.set_xlabel('Iteration')
    ax5.set_ylabel('|Change in Rolling Average| (Log Scale)')
    ax5.set_title('Change in Rolling Average per Iteration')
    ax5.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax5.grid(True, alpha=0.3)
    
    # Plot 6: Change in Rolling Average over 1000 iterations 'Slope' (Log Scale)
    ax6 = fig.add_subplot(3, 2, 6)
    for base, data in all_data.items():
        rolling_avg = data['energy'].rolling(window=min(100, len(data['energy'])//10), min_periods=1).mean()
        window_1000 = min(1000, len(rolling_avg)//5)
        if window_1000 < 10:
            window_1000 = min(10, len(rolling_avg)//2)
        
        slopes = []
        iterations_slope = []
        
        for i in range(window_1000, len(rolling_avg)):
            y_vals = rolling_avg[i-window_1000:i]
            x_vals = np.arange(len(y_vals))
            if len(y_vals) > 1:
                slope = np.polyfit(x_vals, y_vals, 1)[0]
                slopes.append(abs(slope))
                iterations_slope.append(i)
        
        if slopes:
            slopes = np.array(slopes)
            slopes[slopes == 0] = 1e-10  # Avoid log(0)
            ax6.semilogy(iterations_slope, slopes, 
                        label=base, color=data['color'], linewidth=1.5, alpha=0.8)
    
    ax6.set_xlabel('Iteration')
    ax6.set_ylabel('|Slope of Rolling Average| (Log Scale)')
    ax6.set_title(f'Rolling Average Slope over {window_1000} Iterations')
    ax6.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()

# For compatibility with different interfaces
def energy_analysis_driver(csv_files, table_labels, curve_colours, title="Energy Analysis"):
    """
    Driver function that matches the interface expected by visualizer_revamped
    """
    base_filenames = [Path(f).stem for f in csv_files]
    run_energy_analysis(base_filenames, batch_mode=True, input_dir='', file_suffix='.csv')