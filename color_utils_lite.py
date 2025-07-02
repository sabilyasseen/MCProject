"""
Lightweight color utilities for Monte Carlo visualization suite.

This version works without matplotlib dependencies for testing.
"""

import colorsys

# Base color scheme for different algorithm types
BASE_COLOURS = {
    "NewMC":              "#ffcc99",  # Light orange
    "DR":                 "#ffcccc",  # Light red
    "NewMC_Expand":       "#ccffcc",  # Light green
    "NewMC_Restricted":   "#ccccff",  # Light blue
    "NewMC_Biased":       "#ffffcc",  # Light yellow
    "NewMC_Restricted_Expand": "#ffccff", # Light magenta
    "NewMC_Biased_Expand":    "#ccffff", # Light cyan
    "Kawasaki":           "#ffcccc",  # Light red
    "Kawasaki_Black":     "#000000",  # Dark black
}

def _anchor_colour(stem: str) -> str:
    """
    Find the base color for a given algorithm stem name.
    
    Args:
        stem: Algorithm name stem (e.g., "NewMC", "DR", "Kawasaki")
        
    Returns:
        Hex color string for the base color
    """
    for k in sorted(BASE_COLOURS, key=len, reverse=True):
        if stem.startswith(k): 
            return BASE_COLOURS[k]
    return "#808080"  # Default gray

def _simple_to_rgb(hex_color: str):
    """Simple hex to RGB conversion without matplotlib."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16)/255.0 for i in (0, 2, 4))

def _simple_to_hex(rgb):
    """Simple RGB to hex conversion without matplotlib."""
    r, g, b = rgb
    return f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}"

def _tint_red(hex_colour: str, factor: float) -> str:
    """Apply red tint to a color."""
    r, g, b = _simple_to_rgb(hex_colour)
    r = min(1.0, r * factor)
    return _simple_to_hex((r, g, b))

def _tint_green(hex_colour: str, factor: float) -> str:
    """Apply green tint to a color."""
    r, g, b = _simple_to_rgb(hex_colour)
    g = min(1.0, g * factor)
    return _simple_to_hex((r, g, b))

def _tint_blue(hex_colour: str, factor: float) -> str:
    """Apply blue tint to a color."""
    r, g, b = _simple_to_rgb(hex_colour)
    b = min(1.0, b * factor)
    return _simple_to_hex((r, g, b))

def get_distinct_colors(n):
    """
    Generate n distinct colors for visualization.
    
    Args:
        n: Number of distinct colors needed
        
    Returns:
        List of hex color strings, with black as the first color
    """
    colors = ["#000000"]  # Black for the first color
    if n > 1:
        # Use HSV to get perceptually distinct colors
        # Start hue from a non-zero value to avoid very dark colors close to black
        for i in range(n - 1):
            hue = (i / (n - 1)) * 0.8 + 0.1  # Spread hues between 0.1 and 0.9
            # Ensure good saturation and value for visibility
            saturation = 0.8  # Fixed saturation
            value = 0.75      # Fixed value
            rgb = colorsys.hsv_to_rgb(hue, saturation, value)
            colors.append(_simple_to_hex(rgb))
    return colors