import sys
from typing import List, Tuple
from color_utils import _anchor_colour, _tint_red, _tint_green, BASE_COLOURS

def build_filenames(*, BJ: str, include_kawasaki: bool, base_variants: List[str], signifiers: List[List[str]], cluster_types: List[str]) -> Tuple[List[str], List[str], List[str], List[str], List[str], List[str]]:
    csv, tbl, col, leg, leg_col, expand = [], [], [], [], [], []

    if include_kawasaki:
        fn = f"BJ={BJ}_Kawasaki_0_tol0.00_scp0.00.csv"
        csv.append(fn)
        tbl.append("Kawasaki")
        col.append(BASE_COLOURS["Kawasaki_Black"])
        leg.append("Kawasaki")
        leg_col.append(BASE_COLOURS["Kawasaki_Black"])

    for stem in base_variants:
        for ct in cluster_types:
            anchor = _anchor_colour(stem)
            name = f"{stem}_{ct}"
            for j, signifier_1 in enumerate(signifiers[0]):
                for k, signifier_2 in enumerate(signifiers[1]):
                    fn_signifier = f"BJ={BJ}_{name}_{signifier_1}_{signifier_2}.csv"
                    csv.append(fn_signifier)
                    tbl.append(f"{name}_{signifier_1}_{signifier_2}")
                    tint1 = 0.8 + (0.4 * j/len(signifiers[0]))
                    tint2 = 0.8 + (0.4 * k/len(signifiers[1]))
                    tinted_color = _tint_red(_tint_green(anchor, tint1), tint2)
                    col.append(tinted_color)
                    if k == 0 or k == len(signifiers[1])-1:
                        leg.append(f"{name}_{signifier_1}_{signifier_2}")
                        leg_col.append(tinted_color)
                    expand.append(fn_signifier)
    
    return csv, tbl, col, leg, leg_col, expand
