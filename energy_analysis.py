"""
Energy analysis module for Monte Carlo visualization suite.

This module provides functions for analyzing energy convergence, plotting energy trajectories,
and generating convergence tables and PDFs from simulation data.
"""

import sys
import math
import itertools
from pathlib import Path
from typing import List, Optional, Dict, Any
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def compare_runs(csv_files: List[str], *, 
                table_labels: List[str], 
                curve_colours: List[str],
                legend_labels: Optional[List[str]], 
                legend_colours: Optional[List[str]],
                title: str = "Energy & Convergence Analysis"):
    """
    Compare multiple simulation runs through energy plots, slopes, and convergence analysis.
    
    Args:
        csv_files: List of CSV file paths
        table_labels: Labels for each run
        curve_colours: Colors for plotting curves
        legend_labels: Labels for legend
        legend_colours: Colors for legend
        title: Title for the analysis plots
    """
    
    plt.rcParams.update({'font.size': 13})

    # --- Figure 1: Original Energy Plots ---
    fig1 = plt.figure(figsize=(20, 15))
    gs1 = fig1.add_gridspec(3, 3, height_ratios=[1, 1, 1], hspace=0.4)

    ax_raw = fig1.add_subplot(gs1[0, :])
    ax_E, ax_dE, ax_simple = fig1.add_subplot(gs1[1, 0]), fig1.add_subplot(gs1[1, 1]), fig1.add_subplot(gs1[1, 2])
    ax_slope = fig1.add_subplot(gs1[2, :])

    fig1.suptitle(title + " - Energy Plots (Full Data)", fontsize=22, fontweight="bold", y=0.95)

    lengths = [len(pd.read_csv(f)) for f in csv_files if Path(f).is_file()]
    if not lengths: 
        return
    
    W = max(50, min(1000, min(lengths)//10))
    thr = 1e-4
    streak = max(3, min(10, W//100))
    max_steps_full = 0  # Max steps for full data plots

    T_full_data = {k: [] for k in ("run", "Efin", "conv_iter", "conv_E",
                                   "norm_Efin", "norm_conv_iter", "norm_conv_E", "acceptance_rate")}
    all_filtered_energies = {}  # Dictionary to store energies and weights for each run

    # Loop for plotting full data
    for fp, lbl_tbl, col, lbl_leg in zip(csv_files, table_labels, curve_colours,
                                        itertools.cycle(legend_labels or [None])):
        p = Path(fp)
        if not p.is_file(): 
            continue
        try:
            df = pd.read_csv(p)
            if "Energy" not in df.columns: 
                continue

            steps = np.arange(len(df))
            max_steps_full = max(max_steps_full, len(df))
            E = df["Energy"].ffill().bfill()
            ax_raw.plot(steps, E, lw=1, color=col, alpha=0.5, label=lbl_leg)

            cum = E.expanding().mean()
            ax_E.plot(steps, cum, lw=2, color=col, alpha=0.7)
            
            # Simple energy plot with log x-axis
            ax_simple.plot(steps, E, lw=1, color=col, alpha=0.7, label=lbl_leg)

            # Calculate acceptance rate
            if "Total_Acceptances" in df.columns and "Move_Attempts" in df.columns:
                acceptance_rate = (df["Total_Acceptances"].iloc[-1] * 100.0 / df["Move_Attempts"].iloc[-1]) if df["Move_Attempts"].iloc[-1] > 0 else 0.0
            else:
                acceptance_rate = float('nan')

            wn = None
            weights = None
            if "Norm" in df.columns:
                N_norm = df["Norm"].replace([np.inf, -np.inf, 0], np.nan).ffill().bfill()
                eff = (1-N_norm)*0/N_norm + 1  # This line seems to set eff to 1 always if N_norm is a number
                wn = (E*eff).cumsum()/eff.cumsum()
                weights = eff  # Store weights for histogram
                ax_E.plot(steps, wn, lw=2, ls="--", color=col, alpha=0.7)

            stride = max(1, len(steps)//1000)
            ax_dE.scatter(steps[1::stride], np.abs(np.diff(cum))[::stride],
                         s=10, alpha=.7, color=col)
            
            conv = nconv = None
            if len(cum) >= 2*W:
                idx = np.arange(0, len(cum)-W, W)
                sl = np.array([(cum[i+W]-cum[i])/W for i in idx])
                ax_slope.scatter(idx+W/2, np.abs(sl), s=24, alpha=.8, color=col)
                conv = next((int(idx[i]) for i in range(len(sl)-streak+1)
                            if np.all(np.abs(sl[i:i+streak]) < thr)), None)
            
            if wn is not None and len(wn) >= 2*W:
                idxN = np.arange(0, len(wn)-W, W)
                slN = np.array([(wn[i+W]-wn[i])/W for i in idxN])
                ax_slope.scatter(idxN+W/2, np.abs(slN),
                               marker="x", s=26, alpha=.8, color=col)
                nconv = next((int(idxN[i]) for i in range(len(slN)-streak+1)
                             if np.all(np.abs(slN[i:i+streak]) < thr)), None)

            # Store filtered energies with their color and weights
            start_index_for_80_percent = int(0.2 * len(df))
            all_filtered_energies[lbl_tbl] = {
                'energies': df["Energy"].iloc[start_index_for_80_percent:].tolist(),
                'color': col,
                'weights': weights.iloc[start_index_for_80_percent:].tolist() if weights is not None else None,
                'steps': steps[start_index_for_80_percent:].tolist(),
                'cum': cum.iloc[start_index_for_80_percent:].tolist(),
                'wn': wn.iloc[start_index_for_80_percent:].tolist() if wn is not None else None
            }

            # Populate table data for full run
            T_full_data["run"].append(lbl_tbl)
            T_full_data["Efin"].append(f"{cum.iloc[-1]:.4g}")
            T_full_data["conv_iter"].append(conv if conv is not None else "¬conv")
            T_full_data["conv_E"].append(f"{cum[conv]:.4g}" if conv is not None else "N/A")
            T_full_data["acceptance_rate"].append(f"{acceptance_rate:.2f}%" if not np.isnan(acceptance_rate) else "N/A")
            
            if wn is not None:
                T_full_data["norm_Efin"].append(f"{wn.iloc[-1]:.4g}")
                T_full_data["norm_conv_iter"].append(nconv if nconv is not None else "N/A")
                T_full_data["norm_conv_E"].append(f"{wn[nconv]:.4g}" if nconv is not None else "N/A")
            else:
                T_full_data["norm_Efin"].append("N/A")
                T_full_data["norm_conv_iter"].append("N/A")
                T_full_data["norm_conv_E"].append("N/A")

        except Exception as e: 
            print(f"[warn] {p.name}: {e}", file=sys.stderr)

    # Cosmetics for Figure 1 axes (Full Data)
    ax_raw.set(xlabel="step", ylabel="E")
    ax_raw.grid(ls=":")
    if legend_labels and any(legend_labels):
        ax_raw.legend(fontsize=12, handleheight=8, handlelength=8)
        for legline in ax_raw.get_legend().get_lines():
            legline.set_linewidth(6.0)

    ax_E.set(xscale="log", xlabel="step", ylabel="⟨E⟩")
    ax_E.grid(ls=":")

    ax_dE.set(xscale="log", yscale="log", xlim=(1, max_steps_full),
              xlabel="step", ylabel="|Δ⟨E⟩|")
    ax_dE.grid(ls=":")

    ax_simple.set(xscale="log", xlim=(1, max_steps_full), xlabel="step", ylabel="E")
    ax_simple.grid(ls=":")
    if legend_labels and any(legend_labels):
        ax_simple.legend(fontsize=10, loc='best')

    ax_slope.set(xscale="log", yscale="log",
                 xlabel="step", ylabel=f"|slope| (W={W})")
    ax_slope.grid(ls=":")

    plt.tight_layout()
    fig1.savefig("energy_plots_full_data.png", dpi=300, bbox_inches="tight",
                 pad_inches=0.5, transparent=False)
    plt.show()

    # --- Figure 2: Original Convergence Table ---
    _create_convergence_table(T_full_data, title + " - Convergence Data (Full Data)", 
                             "convergence_table_full_data.png")

    # --- Figure 3: Filtered Data Plots ---
    if all_filtered_energies:
        _create_filtered_plots(all_filtered_energies, title, W, thr, streak, legend_labels)

    # --- Figure 5: PDF Plot ---
    if all_filtered_energies:
        _create_pdf_plot(all_filtered_energies, title)


def _create_convergence_table(table_data: Dict[str, List], title: str, filename: str):
    """Create and display a convergence table."""
    num_rows = len(table_data["run"])
    row_height = 0.5
    fig_height = max(8.5, 2 + (num_rows * row_height))

    fig = plt.figure(figsize=(20, fig_height))
    ax_tbl = fig.add_subplot(111)
    ax_tbl.axis("off")

    fig.suptitle(title, fontsize=22, fontweight="bold", y=0.95)

    cell_text = list(zip(table_data["run"], table_data["Efin"], table_data["conv_iter"], table_data["conv_E"],
                        table_data["norm_Efin"], table_data["norm_conv_iter"], table_data["norm_conv_E"],
                        table_data["acceptance_rate"]))

    col_labels = ["run", "final E", "conv iter", "conv E",
                  "Norm Efin", "Norm conv iter", "Norm conv E",
                  "Accept Rate"]

    cell_scale = min(2.0, 20/num_rows)

    tbl = ax_tbl.table(cellText=cell_text,
                       colLabels=col_labels,
                       loc="center",
                       cellLoc="center")

    font_size = min(12, 160/num_rows)
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(font_size)
    tbl.scale(cell_scale, cell_scale * 1.2)

    for cell in tbl._cells.values():
        cell.set_edgecolor('black')
        cell.set_facecolor('white')
        cell.set_alpha(0.8)

    plt.tight_layout()
    fig.savefig(filename, dpi=300, bbox_inches="tight", pad_inches=0.5, transparent=False)
    plt.show()


def _create_filtered_plots(all_filtered_energies: Dict, title: str, W: int, thr: float, streak: int, legend_labels: Optional[List[str]]):
    """Create filtered data plots."""
    fig3 = plt.figure(figsize=(20, 17))
    gs3 = fig3.add_gridspec(3, 3, height_ratios=[1, 1, 1], hspace=0.4)

    ax_raw_f = fig3.add_subplot(gs3[0, :])
    ax_E_f, ax_dE_f, ax_simple_f = fig3.add_subplot(gs3[1, 0]), fig3.add_subplot(gs3[1, 1]), fig3.add_subplot(gs3[1, 2])
    ax_slope_f = fig3.add_subplot(gs3[2, :])

    fig3.suptitle(title + " - Energy Plots (Feed-in Removed)", fontsize=22, fontweight="bold", y=0.95)

    max_steps_filtered = 0
    T_filtered_data = {k: [] for k in ("run", "Efin", "conv_iter", "conv_E",
                                      "norm_Efin", "norm_conv_iter", "norm_conv_E",
                                      "acceptance_rate")}

    for run_name, run_data in all_filtered_energies.items():
        steps = run_data['steps']
        max_steps_filtered = max(max_steps_filtered, len(steps))
        
        # Plot raw energies
        ax_raw_f.plot(steps, run_data['energies'], lw=1, color=run_data['color'], alpha=0.5, label=run_name)
        
        # Plot cumulative average
        ax_E_f.plot(steps, run_data['cum'], lw=2, color=run_data['color'], alpha=0.7)
        
        if run_data['wn'] is not None:
            ax_E_f.plot(steps, run_data['wn'], lw=2, ls="--", color=run_data['color'], alpha=0.7)

        # Plot simple energy with log x-axis (filtered data)
        ax_simple_f.plot(steps, run_data['energies'], lw=1, color=run_data['color'], alpha=0.7, label=run_name)

        # Plot energy differences
        stride = max(1, len(steps)//1000)
        cum_array = np.array(run_data['cum'])
        ax_dE_f.scatter(steps[1::stride], np.abs(np.diff(cum_array))[::stride],
                       s=10, alpha=.7, color=run_data['color'])

        # Plot slopes
        if len(cum_array) >= 2*W:
            idx = np.arange(0, len(cum_array)-W, W)
            sl = np.array([(cum_array[i+W]-cum_array[i])/W for i in idx])
            ax_slope_f.scatter(idx+W/2, np.abs(sl), s=24, alpha=.8, color=run_data['color'])
            conv = next((int(idx[i]) for i in range(len(sl)-streak+1)
                        if np.all(np.abs(sl[i:i+streak]) < thr)), None)
        else:
            conv = None

        # Populate filtered data table
        T_filtered_data["run"].append(run_name)
        T_filtered_data["Efin"].append(f"{run_data['cum'][-1]:.4g}")
        T_filtered_data["conv_iter"].append(conv if conv is not None else "¬conv")
        T_filtered_data["conv_E"].append(f"{run_data['cum'][conv]:.4g}" if conv is not None else "N/A")
        
        if run_data['wn'] is not None:
            T_filtered_data["norm_Efin"].append(f"{run_data['wn'][-1]:.4g}")
            T_filtered_data["norm_conv_iter"].append("N/A")  # Simplified for filtered data
            T_filtered_data["norm_conv_E"].append("N/A")    # Simplified for filtered data
        else:
            T_filtered_data["norm_Efin"].append("N/A")
            T_filtered_data["norm_conv_iter"].append("N/A")
            T_filtered_data["norm_conv_E"].append("N/A")
        
        # Calculate acceptance rate for filtered data (simplified)
        T_filtered_data["acceptance_rate"].append("N/A")

    # Cosmetics for Figure 3 axes (Filtered Data)
    ax_raw_f.set(xlabel="step", ylabel="E")
    ax_raw_f.grid(ls=":")
    if any(legend_labels or []):
        ax_raw_f.legend(fontsize=12, handleheight=8, handlelength=8)
        for legline in ax_raw_f.get_legend().get_lines():
            legline.set_linewidth(6.0)

    ax_E_f.set(xscale="log", xlabel="step", ylabel="⟨E⟩")
    ax_E_f.grid(ls=":")

    # Fix x-axis alignment - use proper range for filtered data
    min_step_filtered = min([min(run_data['steps']) for run_data in all_filtered_energies.values()])
    ax_dE_f.set(xscale="log", yscale="log", xlim=(min_step_filtered, max_steps_filtered),
                xlabel="step", ylabel="|Δ⟨E⟩|")
    ax_dE_f.grid(ls=":")

    ax_simple_f.set(xscale="log", xlim=(min_step_filtered, max_steps_filtered), xlabel="step", ylabel="E")
    ax_simple_f.grid(ls=":")
    if any(legend_labels or []):
        ax_simple_f.legend(fontsize=10, loc='best')

    ax_slope_f.set(xscale="log", yscale="log",
                   xlabel="step", ylabel=f"|slope| (W={W})")
    ax_slope_f.grid(ls=":")

    plt.tight_layout()
    fig3.savefig("energy_plots_filtered_data.png", dpi=300, bbox_inches="tight",
                 pad_inches=0.5, transparent=False)
    plt.show()

    # Create filtered data convergence table
    _create_convergence_table(T_filtered_data, title + " - Convergence Data (Feed-in Removed)", 
                             "convergence_table_filtered_data.png")


def _create_pdf_plot(all_filtered_energies: Dict, title: str):
    """Create PDF plot for filtered energy data."""
    fig5 = plt.figure(figsize=(10, 7))
    ax_pdf = fig5.add_subplot(111)
    fig5.suptitle(title + " - PDF of Raw Energy (Feed-in Removed)", fontsize=18, fontweight="bold")

    # Create bins from min to max energy
    all_energies = [e for run in all_filtered_energies.values() for e in run['energies']]
    max_e = max(all_energies)
    min_e = min(all_energies)
    bin_start = math.floor(min_e/500) * 500 - 5
    bin_end = math.floor(max_e/500) * 500 - 5
    bins = np.arange(bin_start, bin_end + 500, 500)

    # Plot weighted histogram for each run
    for run_name, run_data in all_filtered_energies.items():
        weights = run_data['weights'] if run_data['weights'] is not None else None
        ax_pdf.hist(run_data['energies'], bins=bins, density=True, alpha=0.5,
                   weights=weights, color=run_data['color'], 
                   label=run_name, edgecolor='black')

    ax_pdf.set_xlabel("Energy")
    ax_pdf.set_ylabel("Probability Density")
    ax_pdf.grid(ls=":")
    ax_pdf.legend(fontsize=8)
    plt.tight_layout()
    fig5.savefig("energy_pdf_feed_in_removed.png", dpi=300, bbox_inches="tight",
                 pad_inches=0.5, transparent=False)
    plt.show()