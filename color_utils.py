import colorsys
from matplotlib.colors import to_rgb, to_hex

BASE_COLOURS = {
    "NewMC": "#ffcc99", "DR": "#ffcccc", "NewMC_Expand": "#ccffcc",
    "NewMC_Restricted": "#ccccff", "NewMC_Biased": "#ffffcc",
    "NewMC_Restricted_Expand": "#ffccff", "NewMC_Biased_Expand": "#ccffff",
    "Kawasaki": "#ffcccc", "Kawasaki_Black": "#000000",
}

def get_distinct_colors(n):
    colors = ["#000000"]
    if n > 1:
        for i in range(n - 1):
            hue = (i / (n - 1)) * 0.8 + 0.1
            saturation = 0.8
            value = 0.75
            rgb = colorsys.hsv_to_rgb(hue, saturation, value)
            colors.append(to_hex(rgb))
    return colors

def _anchor_colour(stem: str) -> str:
    for k in sorted(BASE_COLOURS, key=len, reverse=True):
        if stem.startswith(k): 
            return BASE_COLOURS[k]
    return "#808080"

def _tint_red(hex_colour: str, factor: float) -> str:
    r, g, b = to_rgb(hex_colour)
    r = min(1.0, r * factor)
    return to_hex((r, g, b))

def _tint_green(hex_colour: str, factor: float) -> str:
    r, g, b = to_rgb(hex_colour)
    g = min(1.0, g * factor)
    return to_hex((r, g, b))
