#!/usr/bin/env python3
"""
Demo script for the Monte Carlo Visualization Suite.

This script demonstrates how to use the modular visualization components.
It can be run in a Jupyter notebook by copying cells, or as a standalone script.
"""

# =============================================================================
# CELL 1: Import and Setup
# =============================================================================

print("=" * 60)
print("🎯 MONTE CARLO VISUALIZATION SUITE DEMO")
print("=" * 60)

# Check if we're in a Jupyter environment
try:
    from IPython import get_ipython
    if get_ipython() is not None:
        print("📘 Running in Jupyter environment")
        # Enable inline plotting for Jupyter
        get_ipython().run_line_magic('matplotlib', 'inline')
    else:
        print("🐍 Running as standalone Python script")
except ImportError:
    print("🐍 Running as standalone Python script")

# Test imports with graceful fallback
print("\n📦 Testing imports...")

try:
    # Test our modular components
    from color_utils import get_distinct_colors, BASE_COLOURS
    print("  ✅ color_utils imported")
    
    from filename_factory import build_filenames
    print("  ✅ filename_factory imported")
    
    # Only import the heavy modules if dependencies are available
    dependencies_available = True
    try:
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        print("  ✅ numpy, pandas, matplotlib available")
        
        from cluster_analysis import batch_cluster_analysis
        from energy_analysis import compare_runs
        from grid_visualization import plot_grid_from_file
        from visualization_drivers import simple_batch, run_visualizations
        print("  ✅ All visualization modules imported")
        
    except ImportError as e:
        dependencies_available = False
        print(f"  ⚠️  Heavy dependencies not available: {e}")
        print("  💡 Install with: pip install -r requirements.txt")
        
except ImportError as e:
    print(f"  ❌ Module import failed: {e}")
    print("  💡 Ensure all .py files are in the same directory")
    exit(1)

# =============================================================================
# CELL 2: Color Utilities Demo
# =============================================================================

print("\n" + "=" * 60)
print("🎨 COLOR UTILITIES DEMO")
print("=" * 60)

# Demonstrate color generation
print("Generating distinct color palettes:")
for n in [3, 5, 8]:
    colors = get_distinct_colors(n)
    print(f"  {n} colors: {colors}")

# Show base algorithm colors
print(f"\nBase algorithm colors ({len(BASE_COLOURS)} total):")
for algo, color in BASE_COLOURS.items():
    print(f"  {algo:20} -> {color}")

if dependencies_available:
    print("\n📊 Creating color palette visualization...")
    try:
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        
        for i, n in enumerate([3, 5, 8]):
            colors = get_distinct_colors(n)
            ax = axes[i]
            
            # Create a simple bar chart showing the colors
            bars = ax.bar(range(n), [1] * n, color=colors)
            ax.set_title(f'{n} Distinct Colors')
            ax.set_xlabel('Color Index')
            ax.set_ylabel('Value')
            
            # Add color hex values as labels
            for j, (bar, color) in enumerate(zip(bars, colors)):
                ax.text(j, 0.5, color, rotation=90, ha='center', va='center', 
                       fontsize=8, color='white', weight='bold')
        
        plt.tight_layout()
        plt.show()
        print("  ✅ Color palette visualization created")
        
    except Exception as e:
        print(f"  ⚠️  Visualization failed: {e}")

# =============================================================================
# CELL 3: Filename Factory Demo
# =============================================================================

print("\n" + "=" * 60)
print("📁 FILENAME FACTORY DEMO")
print("=" * 60)

print("Generating filenames for parameter combinations...")

# Example parameter combinations
BJ = "-0.44"
base_variants = ["NewMC", "DR"]
signifiers = [["tol0.01", "tol0.05"], ["scp0.00", "scp0.50"]]
cluster_types = ["1", "2"]

csv_files, tbl, col, leg, leg_col, expand = build_filenames(
    BJ=BJ,
    include_kawasaki=True,
    base_variants=base_variants,
    signifiers=signifiers,
    cluster_types=cluster_types
)

print(f"Generated {len(csv_files)} filename combinations:")
for i, (filename, label, color) in enumerate(zip(csv_files[:10], tbl[:10], col[:10])):
    print(f"  {i+1:2d}. {filename}")
    print(f"      Label: {label}")
    print(f"      Color: {color}")
    if i >= 4:  # Show first 5 only
        print(f"      ... and {len(csv_files)-5} more")
        break

# =============================================================================
# CELL 4: File Discovery Demo
# =============================================================================

print("\n" + "=" * 60)
print("🔍 FILE DISCOVERY DEMO")
print("=" * 60)

from pathlib import Path

print("Scanning current directory for data files...")

# Look for CSV files
csv_files_found = list(Path(".").glob("*.csv"))
print(f"📊 Found {len(csv_files_found)} CSV files total")

# Look for Monte Carlo specific files
mc_files = list(Path(".").glob("BJ=*_*.csv"))
print(f"🎯 Found {len(mc_files)} Monte Carlo files")

# Look for grid config files
grid_files = list(Path(".").glob("*_GridConfig.csv"))
print(f"🔲 Found {len(grid_files)} grid configuration files")

if csv_files_found:
    print("\n📋 Sample files found:")
    for f in sorted(csv_files_found)[:5]:
        print(f"  - {f.name}")
    if len(csv_files_found) > 5:
        print(f"  ... and {len(csv_files_found)-5} more")
else:
    print("\n💡 No CSV files found. To test with data:")
    print("   1. Place your simulation output files in this directory")
    print("   2. Files should follow pattern: BJ=X_Algorithm_Y_tolZ_scpW.csv")

# =============================================================================
# CELL 5: Modular Architecture Demo
# =============================================================================

print("\n" + "=" * 60)
print("🏗️  MODULAR ARCHITECTURE DEMO")
print("=" * 60)

print("The visualization suite consists of these modules:")
print()

modules = [
    ("color_utils.py", "Color palette management and utility functions"),
    ("filename_factory.py", "Automatic filename and label generation"),
    ("cluster_analysis.py", "Cluster size tracking and boundary analysis"),
    ("energy_analysis.py", "Energy convergence and statistical analysis"),
    ("grid_visualization.py", "Static grid configuration plotting"),
    ("visualization_drivers.py", "High-level orchestration functions")
]

for module, description in modules:
    status = "✅" if Path(module).exists() else "❌"
    print(f"{status} {module:25} - {description}")

print("\n🔧 Key Benefits of Modular Design:")
print("  • Each component can be tested independently")
print("  • Easy to modify or extend individual features")
print("  • Cleaner error handling and debugging")
print("  • Reusable components for different analysis workflows")
print("  • Reduced memory usage compared to monolithic notebook")

# =============================================================================
# CELL 6: Usage Examples (syntax only if no dependencies)
# =============================================================================

print("\n" + "=" * 60)
print("📚 USAGE EXAMPLES")
print("=" * 60)

if dependencies_available:
    print("🎉 Full functionality available! Here are working examples:")
    
    print("\n1. Simple Batch Processing:")
    print("   from visualization_drivers import simple_batch")
    print("   simple_batch('file1', 'file2', plot_energy=True, plot_clusters=True)")
    
    print("\n2. Complex Parameter Combinations:")
    print("   from visualization_drivers import run_visualizations")
    print("   run_visualizations(BJ='-0.44', base_variants=['DR'], plot_energy=True)")
    
    print("\n3. Individual Component Usage:")
    print("   from cluster_analysis import batch_cluster_analysis")
    print("   batch_cluster_analysis('file1', 'file2')")
    
    # Try to run a simple example if files exist
    if mc_files:
        print(f"\n🎬 Demo with available files...")
        try:
            # Get first few files for demo
            demo_files = [f.stem for f in mc_files[:2]]
            print(f"   Demo files: {demo_files}")
            
            # This would run a real demo - commented out to avoid heavy processing
            # simple_batch(*demo_files, plot_energy=False, plot_clusters=True, plot_grids=False)
            print("   (Demo execution skipped to save processing time)")
            
        except Exception as e:
            print(f"   Demo failed: {e}")
    
else:
    print("⚠️  Dependencies not available. Here's what you could do:")
    print()
    print("1. Install dependencies:")
    print("   pip install -r requirements.txt")
    print()
    print("2. Then run these examples:")
    print("   from visualization_drivers import simple_batch")
    print("   simple_batch('your_file1', 'your_file2', plot_energy=True)")
    print()
    print("3. Or use the complex workflow:")
    print("   from visualization_drivers import run_visualizations")
    print("   run_visualizations(BJ='-0.44', base_variants=['DR', 'NewMC'])")

# =============================================================================
# CELL 7: Summary and Next Steps
# =============================================================================

print("\n" + "=" * 60)
print("🎯 SUMMARY AND NEXT STEPS")
print("=" * 60)

print("✅ Modular visualization suite successfully created!")
print()
print("📁 Files created:")
print("  • 6 core Python modules")
print("  • Test suite for validation")
print("  • Requirements file")
print("  • Documentation")
print()
print("🧪 Testing:")
print("  • Syntax test: python3 test_syntax_only.py")
print("  • Full test: python3 test_visualization_suite.py")
print()
print("🚀 Ready to use:")
if dependencies_available:
    print("  • All dependencies available")
    print("  • Ready for data analysis")
else:
    print("  • Install dependencies: pip install -r requirements.txt")
    print("  • Then ready for data analysis")

print("\n💡 Tips:")
print("  • Place your CSV files in the current directory")
print("  • Use simple_batch() for straightforward analysis")
print("  • Use run_visualizations() for complex parameter sweeps")
print("  • Each module can be used independently")

print("\n🎉 Happy analyzing!")

if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("📝 This demo completed successfully!")
    print("Copy the code blocks above into Jupyter notebook cells to run interactively.")
    print("=" * 60)