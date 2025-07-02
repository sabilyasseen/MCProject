"""
Filename and label factory for Monte Carlo visualization suite.

This module provides functions for generating consistent filenames, labels, and colors
for different algorithm configurations and parameter combinations.
"""

import sys
from typing import List, Tuple
from color_utils import _anchor_colour, _tint_red, _tint_green, BASE_COLOURS

def build_filenames(*, BJ: str, include_kawasaki: bool,
                   base_variants: List[str], 
                   signifiers: List[List[str]], 
                   cluster_types: List[str]) -> Tuple[List[str], List[str], List[str], List[str], List[str], List[str]]:
    """
    Build comprehensive lists of filenames, labels, and colors for visualization.
    
    Args:
        BJ: Beta*J parameter value as string
        include_kawasaki: Whether to include Kawasaki algorithm files
        base_variants: Base algorithm variants (e.g., ["NewMC", "DR"])
        signifiers: Parameter combinations (e.g., [["tol0.01"], ["scp0.50"]])  
        cluster_types: Cluster type variants (e.g., ["1", "2"])
        
    Returns:
        Tuple containing:
        - csv_files: List of CSV filenames
        - table_labels: List of table labels
        - curve_colours: List of colors for curves
        - legend_labels: List of legend labels
        - legend_colours: List of colors for legend
        - expand_files: List of all expanded filenames
    """
    csv, tbl, col, leg, leg_col, expand = [], [], [], [], [], []

    if include_kawasaki:
        # Kawasaki specific filename format
        fn = f"BJ={BJ}_Kawasaki_0_tol0.00_scp0.00.csv"
        print(f"Trying file: {fn}", file=sys.stderr)
        csv.append(fn)
        tbl.append("Kawasaki")
        col.append(BASE_COLOURS["Kawasaki_Black"])
        leg.append("Kawasaki")
        leg_col.append(BASE_COLOURS["Kawasaki_Black"])

    # Iterate through base variants, cluster types, and signifier combinations
    for stem in base_variants:
        for ct in cluster_types:
            anchor = _anchor_colour(stem)
            name = f"{stem}_{ct}"
            # Assuming signifiers[0] is for 'tol' and signifiers[1] is for 'scp'
            for j, signifier_1 in enumerate(signifiers[0]):  # e.g., tol values
                for k, signifier_2 in enumerate(signifiers[1]):  # e.g., scp values
                    fn_signifier = f"BJ={BJ}_{name}_{signifier_1}_{signifier_2}.csv"
                    csv.append(fn_signifier)
                    tbl.append(f"{name}_{signifier_1}_{signifier_2}")

                    # Calculate tint factors based on position in both loops
                    # Using much smaller tint ranges for subtler variation
                    tint1 = 0.8 + (0.4 * j/len(signifiers[0]))  # First tint varies from 0.8 to 1.2
                    tint2 = 0.8 + (0.4 * k/len(signifiers[1]))  # Second tint varies from 0.8 to 1.2
                    # Apply both tints sequentially
                    tinted_color = _tint_red(_tint_green(anchor, tint1), tint2)
                    col.append(tinted_color)

                    # Only add to legend for specific iterations (e.g., first and last of the inner loop)
                    # This helps prevent an overly cluttered legend
                    if k == 0 or k == len(signifiers[1])-1:
                        leg.append(f"{name}_{signifier_1}_{signifier_2}")
                        leg_col.append(tinted_color)
                    expand.append(fn_signifier)
    
    return csv, tbl, col, leg, leg_col, expand