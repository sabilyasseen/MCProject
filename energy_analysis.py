import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def _assign_colors(num, user_colors=None):
    """Assign colors, always making the first one black."""
    if user_colors is not None and len(user_colors) >= num:
        colors = list(user_colors)
        colors[0] = 'black'
        return colors
    colors = list(plt.cm.tab10(np.linspace(0, 1, max(num, 2))))
    colors[0] = 'black'
    return colors

def _truncate_label(label, maxlen=30):
    """Truncate a label (such as a filename) for table display."""
    if len(label) <= maxlen:
        return label
    else:
        # Show start and end, hide middle
        return label[:maxlen//2-2] + "..." + label[-maxlen//2+1:]

def _build_summary_table(all_data, rolling_window=100, label_truncate_len=30):
    """
    Build a summary table DataFrame with final energy, final rolling average, and allow for easy extension.
    Truncates long labels for display.
    """
    summary = []
    for label, data in all_data.items():
        energy = data['energy']
        rolling_avg = energy.rolling(window=min(rolling_window, max(1, len(energy)//10)), min_periods=1).mean()
        row = {
            "Label": _truncate_label(str(label), maxlen=label_truncate_len),
            "Full Label": str(label),
            "Final Energy": energy.iloc[-1] if len(energy) > 0 else np.nan,
            "Final Rolling Avg": rolling_avg.iloc[-1] if len(rolling_avg) > 0 else np.nan,
        }
        summary.append(row)
    df = pd.DataFrame(summary)
    df = df.set_index("Label")
    return df

def _save_summary_table_figure(summary_df, output_filename, font_size=8):
    """
    Save the summary table as its own figure, using a simple, flexible, small-text implementation.
    """
    summary_df = summary_df.copy()
    if "Full Label" in summary_df.columns:
        summary_df_display = summary_df.drop(columns=["Full Label"])
    else:
        summary_df_display = summary_df

    n_rows, n_cols = summary_df_display.shape

    # Simple flexible sizing: width/height per cell, but small font
    cell_width = 1.5  # inches per column
    cell_height = 0.4  # inches per row
    fig_width = max(4, n_cols * cell_width + 1.5)
    fig_height = max(2, n_rows * cell_height + 1.2)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.axis('off')

    # Prepare cell text, rounding numbers but keeping strings as is
    cell_text = []
    for i in range(n_rows):
        row = []
        for j, col in enumerate(summary_df_display.columns):
            val = summary_df_display.iloc[i, j]
            if isinstance(val, float) or isinstance(val, np.floating):
                row.append(np.round(val, 6))
            else:
                row.append(str(val))
        cell_text.append(row)

    table = ax.table(
        cellText=cell_text,
        rowLabels=summary_df_display.index,
        colLabels=summary_df_display.columns,
        loc='center',
        cellLoc='center',
        rowLoc='center'
    )
    table.auto_set_font_size(False)
    table.set_fontsize(font_size)
    table.scale(1, 1)
    plt.tight_layout()
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    fig.savefig(output_filename, bbox_inches='tight', dpi=200)
    print(f"Energy analysis table saved as: {output_filename}")
    plt.close(fig)

def _run_energy_analysis_core(all_data, output_filename, title="Comprehensive Energy Analysis", restrict_histogram_to=None):
    """
    Core plotting logic for energy analysis (no table in this figure).
    If restrict_histogram_to is not None, it should be a dict with the same keys as all_data,
    and values are the shortened energy arrays to use for the histogram.
    """
    fig = plt.figure(figsize=(22, 18))
    fig.suptitle(title, fontsize=22, fontweight='bold')
    plt.subplots_adjust(hspace=0.28, wspace=0.18, top=0.90, bottom=0.08)

    # Plot 1: Energy vs Iteration
    ax1 = fig.add_subplot(4, 2, 1)
    for idx, (base, data) in enumerate(all_data.items()):
        ax1.plot(data['iterations'], data['energy'],
                label=_truncate_label(str(base)), color=data['color'], linewidth=3.0, alpha=0.9)
    ax1.set_xlabel('Iteration')
    ax1.set_ylabel('Energy')
    ax1.set_title('Energy vs Iteration', fontsize=20)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper right', fontsize=12, frameon=True)

    # Plot 2: Energy vs Iteration (Logarithmic X-axis)
    ax2 = fig.add_subplot(4, 2, 2)
    for idx, (base, data) in enumerate(all_data.items()):
        ax2.semilogx(data['iterations'][1:], data['energy'][1:],
                    color=data['color'], linewidth=3.0, alpha=0.9)
    ax2.set_xlabel('Iteration (Log Scale)')
    ax2.set_ylabel('Energy')
    ax2.set_title('Energy vs Iteration (Logarithmic X-axis)', fontsize=20)
    ax2.grid(True, alpha=0.3)

    # Plot 3: Average Energy vs Iteration
    ax3 = fig.add_subplot(4, 2, 3)
    for idx, (base, data) in enumerate(all_data.items()):
        rolling_avg = data['energy'].rolling(window=min(100, max(1, len(data['energy'])//10)), min_periods=1).mean()
        ax3.plot(data['iterations'], rolling_avg,
                color=data['color'], linewidth=3.0, alpha=0.95)
        # Add weighted rolling average if weights are available
        if 'weights' in data:
            weighted_rolling_avg = (data['energy'] * data['weights']).rolling(window=min(100, max(1, len(data['energy'])//10)), min_periods=1).sum() / data['weights'].rolling(window=min(100, max(1, len(data['energy'])//10)), min_periods=1).sum()
            ax3.plot(data['iterations'], weighted_rolling_avg,
                    color=data['color'], linewidth=3.0, alpha=0.95, linestyle='--')
    ax3.set_xlabel('Iteration')
    ax3.set_ylabel('Average Energy')
    ax3.set_title('Average Energy vs Iteration (Solid: Regular, Dashed: Weighted)', fontsize=20)
    ax3.grid(True, alpha=0.3)

    # Plot 4: Average Energy vs Iteration (Logarithmic X-axis)
    ax4 = fig.add_subplot(4, 2, 4)
    for idx, (base, data) in enumerate(all_data.items()):
        rolling_avg = data['energy'].rolling(window=min(100, max(1, len(data['energy'])//10)), min_periods=1).mean()
        ax4.semilogx(data['iterations'][1:], rolling_avg[1:],
                    color=data['color'], linewidth=3.0, alpha=0.95)
        # Add weighted rolling average if weights are available
        if 'weights' in data:
            weighted_rolling_avg = (data['energy'] * data['weights']).rolling(window=min(100, max(1, len(data['energy'])//10)), min_periods=1).sum() / data['weights'].rolling(window=min(100, max(1, len(data['energy'])//10)), min_periods=1).sum()
            ax4.semilogx(data['iterations'][1:], weighted_rolling_avg[1:],
                        color=data['color'], linewidth=3.0, alpha=0.95, linestyle='--')
    ax4.set_xlabel('Iteration (Log Scale)')
    ax4.set_ylabel('Average Energy')
    ax4.set_title('Average Energy vs Iteration (Logarithmic X-axis)', fontsize=20)
    ax4.grid(True, alpha=0.3)

    # Plot 5: Change in Average Energy at Each Iteration (Logarithmic Y-axis) -- DOTS
    ax5 = fig.add_subplot(4, 2, 5)
    for idx, (base, data) in enumerate(all_data.items()):
        rolling_avg = data['energy'].rolling(window=min(100, max(1, len(data['energy'])//10)), min_periods=1).mean()
        change = np.abs(np.diff(rolling_avg))
        change[change == 0] = 1e-10  # Avoid log(0)
        # Use dots instead of lines
        ax5.semilogy(data['iterations'][1:], change,
                    color=data['color'], marker='o', linestyle='None', markersize=3, alpha=0.8, label=_truncate_label(str(base)))
    ax5.set_xlabel('Iteration')
    ax5.set_ylabel('|Change in Average Energy| (Log Scale)')
    ax5.set_title('Change in Average Energy per Iteration (Logarithmic Y-axis)', fontsize=20)
    ax5.grid(True, alpha=0.3)
    ax5.legend(loc='best', fontsize=10, frameon=True)

    # Plot 6: Change in Average Energy over 1000 iterations 'Slope' (Logarithmic Y-axis) -- DOTS
    ax6 = fig.add_subplot(4, 2, 6)
    window_1000 = 1000
    for idx, (base, data) in enumerate(all_data.items()):
        rolling_avg = data['energy'].rolling(window=min(100, max(1, len(data['energy'])//10)), min_periods=1).mean()
        actual_window = min(window_1000, max(1, len(rolling_avg)//5))
        if actual_window < 10:
            actual_window = min(10, max(1, len(rolling_avg)//2))

        slopes = []
        iterations_slope = []

        for i in range(actual_window, len(rolling_avg)):
            y_vals = rolling_avg[i-actual_window:i]
            x_vals = np.arange(len(y_vals))
            if len(y_vals) > 1:
                slope = np.polyfit(x_vals, y_vals, 1)[0]
                slopes.append(abs(slope))
                iterations_slope.append(i)

        if slopes:
            slopes = np.array(slopes)
            slopes[slopes == 0] = 1e-10  # Avoid log(0)
            # Use dots instead of lines
            ax6.semilogy(iterations_slope, slopes,
                        color=data['color'], marker='o', linestyle='None', markersize=3, alpha=0.8, label=_truncate_label(str(base)))

    ax6.set_xlabel('Iteration')
    ax6.set_ylabel('|Slope of Average Energy| (Log Scale)')
    ax6.set_title(f'Average Energy Slope over ~{window_1000} Iterations (Logarithmic Y-axis)', fontsize=20)
    ax6.grid(True, alpha=0.3)
    ax6.legend(loc='best', fontsize=10, frameon=True)

    # Plot 7: Histogram of all energy values (make it span the entire width of the figure)
    ax_hist = fig.add_subplot(4, 1, 4)
    all_energies = []
    all_weights = []
    if restrict_histogram_to is None:
        for base, data in all_data.items():
            all_energies.append(np.asarray(data['energy']))
            if 'weights' in data:
                all_weights.append(np.asarray(data['weights']))
            else:
                all_weights.append(np.ones(len(data['energy'])))
    else:
        for base in all_data.keys():
            all_energies.append(np.asarray(restrict_histogram_to[base]))
            if 'weights' in all_data[base]:
                # Use weights corresponding to the restricted energy values
                weights_restricted = np.asarray(all_data[base]['weights'])[:len(restrict_histogram_to[base])]
                all_weights.append(weights_restricted)
            else:
                all_weights.append(np.ones(len(restrict_histogram_to[base])))
    
    if all_energies:
        all_energies_flat = np.concatenate(all_energies)
        all_weights_flat = np.concatenate(all_weights)
        min_energy = np.nanmin(all_energies_flat)
        max_energy = np.nanmax(all_energies_flat)
        start = 5 + 50 * np.floor((min_energy - 5) / 50)
        end = 5 + 50 * np.ceil((max_energy - 5) / 50)
        bins = np.arange(start, end + 50, 50)
    else:
        bins = 50  # fallback

    if restrict_histogram_to is None:
        for idx, (base, data) in enumerate(all_data.items()):
            weights = data.get('weights', np.ones(len(data['energy'])))
            ax_hist.hist(data['energy'], bins=bins, weights=weights, alpha=0.6, color=data['color'], 
                        label=_truncate_label(str(base)), histtype='stepfilled', linewidth=2, density=True)
    else:
        for idx, (base, data) in enumerate(all_data.items()):
            energy_values = restrict_histogram_to[base]
            if 'weights' in data:
                weights = np.asarray(data['weights'])[:len(energy_values)]
            else:
                weights = np.ones(len(energy_values))
            ax_hist.hist(energy_values, bins=bins, weights=weights, alpha=0.6, color=data['color'], 
                        label=_truncate_label(str(base)), histtype='stepfilled', linewidth=2, density=True)
    ax_hist.set_xlabel('Energy Value')
    ax_hist.set_ylabel('Probability Density')
    ax_hist.set_title('Probability Density of Energy Values (Weighted)', fontsize=20)
    ax_hist.grid(True, alpha=0.3)
    ax_hist.legend(loc='best', fontsize=10, frameon=True)

    output_path = os.path.join('outputs', output_filename)
    fig.savefig(output_path, bbox_inches='tight')
    print(f"Energy analysis plot saved as: {output_path}")
    plt.close(fig)

def run_energy_analysis(base_filenames, batch_mode=True, input_dir='outputs', file_suffix='.csv'):
    """
    Enhanced energy analysis with multiple plot types:
    - Energy (linear and log x-axis) with all files overlayed
    - Average energy (linear and log x-axis) with all files overlayed
    - Change in average energy at each iteration (log scale)
    - Change in average energy slope over 1000 iterations (log scale)
    - Summary table of final energies (saved as a separate figure)
    - Histogram of all energy values
    - Also, recursively repeat the analysis with the first 20% of moves removed
    """
    print(f"[EnergyAnalysis][INFO] Running energy analysis for {base_filenames}")

    all_data = {}
    colors = _assign_colors(len(base_filenames))

    for idx, base in enumerate(base_filenames):
        file_path = os.path.join(input_dir, base + file_suffix)
        if not os.path.exists(file_path):
            file_path = base + file_suffix
            if not os.path.exists(file_path):
                print(f"[EnergyAnalysis][ERROR] File not found: {base + file_suffix}", file=sys.stderr)
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
                print(f"[EnergyAnalysis][ERROR] No energy column found in {file_path}", file=sys.stderr)
                continue

            energy_data = df[energy_col].ffill().bfill()
            
            # Calculate weights from Norm column if available
            weights = None
            if 'Norm' in df.columns:
                norm_data = df['Norm'].ffill().bfill()
                # Calculate weights: 1 + (1-Norm)/Norm
                weights = 1 + (1 - norm_data) / norm_data
                weights = weights.fillna(1.0)  # Handle division by zero
            
            all_data[base] = {
                'energy': energy_data,
                'iterations': np.arange(len(energy_data)),
                'color': colors[idx]
            }
            
            if weights is not None:
                all_data[base]['weights'] = weights

        except Exception as e:
            print(f"[EnergyAnalysis][ERROR] Failed to process {file_path}: {e}", file=sys.stderr)

    if not all_data:
        print("[EnergyAnalysis][ERROR] No valid data files found", file=sys.stderr)
        return

    base_concat = ''.join(str(b) for b in base_filenames)
    output_filename = base_concat + '_energy_analysis.png'
    table_filename = base_concat + '_energy_analysis_table.png'

    _run_energy_analysis_core(all_data, output_filename, title="Comprehensive Energy Analysis")

    summary_df = _build_summary_table(all_data, label_truncate_len=30)
    table_output_path = os.path.join('outputs', table_filename)
    _save_summary_table_figure(summary_df, table_output_path)

    shortened_all_data = {}
    shortened_hist_data = {}
    for base, data in all_data.items():
        n = len(data['energy'])
        start_idx = int(np.floor(0.2 * n))
        if n - start_idx < 2:
            continue
        shortened_energy = data['energy'].iloc[start_idx:].reset_index(drop=True)
        shortened_iterations = np.arange(len(shortened_energy))
        shortened_all_data[base] = {
            'energy': shortened_energy,
            'iterations': shortened_iterations,
            'color': data['color']
        }
        if 'weights' in data:
            shortened_weights = data['weights'].iloc[start_idx:].reset_index(drop=True)
            shortened_all_data[base]['weights'] = shortened_weights
        shortened_hist_data[base] = shortened_energy

    if shortened_all_data:
        shortened_output_filename = output_filename.replace('.png', '_shortened.png')
        shortened_table_filename = table_filename.replace('.png', '_shortened.png')
        _run_energy_analysis_core(shortened_all_data, shortened_output_filename,
                                 title="Energy Analysis (After Removing First 20% of Moves)",
                                 restrict_histogram_to=shortened_hist_data)
        shortened_summary_df = _build_summary_table(shortened_all_data, label_truncate_len=30)
        shortened_table_output_path = os.path.join('outputs', shortened_table_filename)
        _save_summary_table_figure(shortened_summary_df, shortened_table_output_path)
    else:
        print("[EnergyAnalysis][WARN] Not enough data to perform shortened analysis.")

def energy_analysis_driver(csv_files, table_labels, curve_colours, title="Energy Analysis"):
    """
    Driver function that matches the interface expected by visualizer_revamped
    This function directly processes the csv_files as provided by the notebook
    """
    all_data = {}
    colors = _assign_colors(len(csv_files), user_colors=curve_colours)

    for idx, (csv_file, table_label) in enumerate(zip(csv_files, table_labels)):
        if not Path(csv_file).is_file():
            print(f"[EnergyAnalysis][ERROR] File not found: {csv_file}", file=sys.stderr)
            continue

        try:
            df = pd.read_csv(csv_file)
            if 'energy' in df.columns:
                energy_col = 'energy'
            else:
                print(f"[EnergyAnalysis][ERROR] No energy column found in {csv_file}", file=sys.stderr)
                continue

            energy_data = df[energy_col].ffill().bfill()
            
            # Calculate weights from Norm column if available
            weights = None
            if 'Norm' in df.columns:
                norm_data = df['Norm'].ffill().bfill()
                # Calculate weights: 1 + (1-Norm)/Norm
                weights = 1 + (1 - norm_data) / norm_data
                weights = weights.fillna(1.0)  # Handle division by zero
            
            all_data[table_label] = {
                'energy': energy_data,
                'iterations': np.arange(len(energy_data)),
                'color': colors[idx]
            }
            
            if weights is not None:
                all_data[table_label]['weights'] = weights

        except Exception as e:
            print(f"[EnergyAnalysis][ERROR] Failed to process {csv_file}: {e}", file=sys.stderr)
            continue

    if not all_data:
        print("[EnergyAnalysis][ERROR] No valid data files found", file=sys.stderr)
        return

    base_concat = ''.join(str(f) for f in csv_files)
    output_filename = base_concat + '_energy_analysis.png'
    table_filename = base_concat + '_energy_analysis_table.png'

    _run_energy_analysis_core(all_data, output_filename, title=title)

    summary_df = _build_summary_table(all_data, label_truncate_len=30)
    table_output_path = os.path.join('outputs', table_filename)
    _save_summary_table_figure(summary_df, table_output_path)

    shortened_all_data = {}
    shortened_hist_data = {}
    for label, data in all_data.items():
        n = len(data['energy'])
        start_idx = int(np.floor(0.2 * n))
        if n - start_idx < 2:
            continue
        shortened_energy = data['energy'].iloc[start_idx:].reset_index(drop=True)
        shortened_iterations = np.arange(len(shortened_energy))
        shortened_all_data[label] = {
            'energy': shortened_energy,
            'iterations': shortened_iterations,
            'color': data['color']
        }
        if 'weights' in data:
            shortened_weights = data['weights'].iloc[start_idx:].reset_index(drop=True)
            shortened_all_data[label]['weights'] = shortened_weights
        shortened_hist_data[label] = shortened_energy

    if shortened_all_data:
        shortened_output_filename = output_filename.replace('.png', '_shortened.png')
        shortened_table_filename = table_filename.replace('.png', '_shortened.png')
        _run_energy_analysis_core(shortened_all_data, shortened_output_filename,
                                 title=title + " (After Removing First 20% of Moves)",
                                 restrict_histogram_to=shortened_hist_data)
        shortened_summary_df = _build_summary_table(shortened_all_data, label_truncate_len=30)
        shortened_table_output_path = os.path.join('outputs', shortened_table_filename)
        _save_summary_table_figure(shortened_summary_df, shortened_table_output_path)
    else:
        print("[EnergyAnalysis][WARN] Not enough data to perform shortened analysis.")