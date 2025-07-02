import sys
import math
import itertools
from pathlib import Path
from typing import List, Optional, Dict, Any
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def compare_runs(csv_files: List[str], *, table_labels: List[str], curve_colours: List[str], legend_labels: Optional[List[str]], legend_colours: Optional[List[str]], title: str = "Energy & Convergence Analysis"):
    plt.rcParams.update({'font.size': 13})
    fig1 = plt.figure(figsize=(20, 15))
    gs1 = fig1.add_gridspec(3, 3, height_ratios=[1, 1, 1], hspace=0.4)
    
    ax_raw = fig1.add_subplot(gs1[0, :])
    ax_E, ax_dE, ax_simple = fig1.add_subplot(gs1[1, 0]), fig1.add_subplot(gs1[1, 1]), fig1.add_subplot(gs1[1, 2])
    ax_slope = fig1.add_subplot(gs1[2, :])
    
    fig1.suptitle(title + " - Energy Plots", fontsize=22, fontweight="bold", y=0.95)
    
    max_steps_full = 0
    
    for fp, lbl_tbl, col, lbl_leg in zip(csv_files, table_labels, curve_colours, itertools.cycle(legend_labels or [None])):
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
            
            ax_simple.plot(steps, E, lw=1, color=col, alpha=0.7, label=lbl_leg)
            
            stride = max(1, len(steps)//1000)
            ax_dE.scatter(steps[1::stride], np.abs(np.diff(cum))[::stride], s=10, alpha=.7, color=col)
            
        except Exception as e: 
            continue

    ax_raw.set(xlabel="step", ylabel="E")
    ax_raw.grid(ls=":")
    if legend_labels and any(legend_labels):
        ax_raw.legend(fontsize=12)

    ax_E.set(xscale="log", xlabel="step", ylabel="⟨E⟩")
    ax_E.grid(ls=":")

    ax_dE.set(xscale="log", yscale="log", xlim=(1, max_steps_full), xlabel="step", ylabel="|Δ⟨E⟩|")
    ax_dE.grid(ls=":")

    ax_simple.set(xscale="log", xlim=(1, max_steps_full), xlabel="step", ylabel="E")
    ax_simple.grid(ls=":")
    if legend_labels and any(legend_labels):
        ax_simple.legend(fontsize=10, loc='best')

    ax_slope.set(xscale="log", yscale="log", xlabel="step", ylabel="|slope|")
    ax_slope.grid(ls=":")

    plt.tight_layout()
    plt.show()
