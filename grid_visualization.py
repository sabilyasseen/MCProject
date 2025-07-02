import sys
import re
import math
from pathlib import Path
from typing import Optional
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle

def parse_grid_config_string(grid_config_str: str, N: int) -> pd.DataFrame:
    if not grid_config_str or grid_config_str in ['nan', 'NaN', '']:
        return pd.DataFrame()
    
    try:
        config_str = str(grid_config_str).strip('"').strip("'")
        if not config_str:
            return pd.DataFrame()
            
        cells = config_str.split(';')
        
        data = []
        for i, cell_str in enumerate(cells):
            if not cell_str or cell_str.strip() == '':
                continue
                
            x = i % N
            y = i // N
            
            numbers = re.findall(r'-?\d+', cell_str)
            
            if len(numbers) >= 2:
                species = int(numbers[0])
                cluster_id = int(numbers[1])
            elif len(numbers) == 1:
                species = int(numbers[0])
                cluster_id = 0
            else:
                species = 0
                cluster_id = 0
            
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
        return pd.DataFrame()

def plot_grid(cfg_df: pd.DataFrame, *, title: Optional[str] = None, figsize=(8, 8), ax=None, save_path: Optional[Path] = None, grid_size: Optional[int] = None):
    if cfg_df.empty:
        return
    
    required_cols = ["x", "y", "species", "cluster_id", "cell_type"]
    missing_cols = [col for col in required_cols if col not in cfg_df.columns]
    if missing_cols:
        return
    
    if grid_size is None:
        N = int(max(cfg_df["x"].max(), cfg_df["y"].max()) + 1) if not cfg_df.empty else 0
    else:
        N = grid_size
    
    create_new_fig = ax is None
    if create_new_fig:
        fig, ax_local = plt.subplots(1, 1, figsize=figsize)
    else:
        ax_local = ax
        fig = ax.figure
        
    ax_local.set(xlim=(-0.5, N-0.5), ylim=(-0.5, N-0.5), aspect="equal", xticks=[], yticks=[])
    ax_local.set_title(title or "Grid Configuration", fontweight="bold", fontsize=14)
    
    df_display = cfg_df.copy()
    
    for _, r in df_display.iterrows():
        if r["species"] == 0:
            base_color = "#d62728"
        elif r["species"] == 1:  
            base_color = "#1f77b4"
        else:
            base_color = "#808080"
        
        if r["cell_type"] == "INTERIOR":
            alpha = 0.4
        elif r["cell_type"] == "BOUNDARY":
            alpha = 0.8
        else:
            alpha = 0.6
        
        rect = plt.Rectangle((r["x"]-0.4, r["y"]-0.4), 0.8, 0.8, facecolor=base_color, alpha=alpha, edgecolor='black', linewidth=0.5)
        ax_local.add_patch(rect)
        
        txt = ax_local.text(r["x"], r["y"], str(r["cluster_id"]), ha="center", va="center", fontsize=8, color='white', weight='bold')
    
    if create_new_fig:
        plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=250, bbox_inches="tight")
    
    if create_new_fig:
        plt.show()

def plot_grid_from_file(csv_file_path: str, *, title: Optional[str] = None, save_path: Optional[str] = None):
    try:
        df = pd.read_csv(csv_file_path)
        save_path_obj = Path(save_path) if save_path else None
        plot_grid(df, title=title or f"Grid from {Path(csv_file_path).name}", save_path=save_path_obj)
    except Exception as e:
        pass
