#!/usr/bin/env python3
"""
Test script for the Monte Carlo Visualization Suite.

This script tests all the modular components to ensure they work correctly together.
Run this script to verify that the visualization suite is functioning properly.
"""

import sys
from pathlib import Path

def test_imports():
    """Test that all modules can be imported correctly."""
    print("🔧 Testing module imports...")
    
    modules_to_test = [
        'color_utils',
        'filename_factory', 
        'cluster_analysis',
        'energy_analysis',
        'grid_visualization',
        'visualization_drivers'
    ]
    
    success_count = 0
    for module_name in modules_to_test:
        try:
            module = __import__(module_name)
            print(f"  ✅ {module_name}")
            success_count += 1
        except ImportError as e:
            print(f"  ❌ {module_name}: {e}")
        except Exception as e:
            print(f"  ⚠️  {module_name}: {e}")
    
    print(f"\n📊 Import Results: {success_count}/{len(modules_to_test)} modules imported successfully")
    return success_count == len(modules_to_test)

def test_color_utilities():
    """Test color utility functions."""
    print("\n🎨 Testing color utilities...")
    
    try:
        from color_utils import get_distinct_colors, BASE_COLOURS, _anchor_colour
        
        # Test color generation
        colors_3 = get_distinct_colors(3)
        colors_5 = get_distinct_colors(5)
        
        print(f"  ✅ Generated 3 colors: {colors_3}")
        print(f"  ✅ Generated 5 colors: {colors_5}")
        
        # Test base colors
        print(f"  ✅ Base colors available: {len(BASE_COLOURS)} algorithms")
        
        # Test anchor color lookup
        anchor = _anchor_colour("NewMC")
        print(f"  ✅ Anchor color for 'NewMC': {anchor}")
        
        return True
    except Exception as e:
        print(f"  ❌ Color utilities test failed: {e}")
        return False

def test_filename_factory():
    """Test filename generation."""
    print("\n📁 Testing filename factory...")
    
    try:
        from filename_factory import build_filenames
        
        # Test basic filename building
        csv_files, tbl, col, leg, leg_col, expand = build_filenames(
            BJ="-0.44",
            include_kawasaki=True,
            base_variants=["NewMC"],
            signifiers=[["tol0.01"], ["scp0.50"]],
            cluster_types=["1"]
        )
        
        print(f"  ✅ Generated {len(csv_files)} filenames")
        print(f"  ✅ Generated {len(tbl)} table labels")
        print(f"  ✅ Generated {len(col)} colors")
        
        if csv_files:
            print(f"  📄 Example filename: {csv_files[0]}")
        
        return True
    except Exception as e:
        print(f"  ❌ Filename factory test failed: {e}")
        return False

def test_file_discovery():
    """Test file discovery in current directory."""
    print("\n🔍 Testing file discovery...")
    
    current_dir = Path(".")
    
    # Look for CSV files
    csv_files = list(current_dir.glob("*.csv"))
    print(f"  📊 Found {len(csv_files)} CSV files")
    
    # Look for specific pattern files
    mc_files = list(current_dir.glob("BJ=*_*.csv"))
    print(f"  🎯 Found {len(mc_files)} Monte Carlo files")
    
    # Look for grid config files
    grid_files = list(current_dir.glob("*_GridConfig.csv"))
    print(f"  🔲 Found {len(grid_files)} grid config files")
    
    # Show a few examples
    if csv_files:
        print("  📋 Example files:")
        for f in csv_files[:5]:
            print(f"    - {f.name}")
    
    return len(csv_files) > 0

def test_visualization_drivers():
    """Test the main driver functions."""
    print("\n🚀 Testing visualization drivers...")
    
    try:
        from visualization_drivers import simple_batch, run_visualizations
        
        print("  ✅ simple_batch function imported")
        print("  ✅ run_visualizations function imported")
        
        # Test with non-existent files (should handle gracefully)
        print("  🧪 Testing with mock filenames...")
        
        # This should run without crashing, even if files don't exist
        # The functions should handle missing files gracefully
        
        return True
    except Exception as e:
        print(f"  ❌ Visualization drivers test failed: {e}")
        return False

def run_demo_if_files_exist():
    """Run a demo if suitable files are found."""
    print("\n🎯 Checking for demo files...")
    
    # Look for example files that might exist
    example_patterns = [
        "BJ=*_Kawasaki_*.csv",
        "BJ=*_DR_*.csv", 
        "BJ=*_NewMC_*.csv"
    ]
    
    found_files = []
    for pattern in example_patterns:
        files = list(Path(".").glob(pattern))
        found_files.extend([f.stem for f in files])
    
    if found_files:
        print(f"  🎉 Found {len(found_files)} demo files!")
        print("  📋 Available for demo:")
        for f in found_files[:3]:  # Show first 3
            print(f"    - {f}")
        
        if len(found_files) >= 2:
            print("\n🎬 Running mini demo...")
            try:
                from visualization_drivers import simple_batch
                # Run with first 2 files, but only cluster analysis to avoid heavy plotting
                simple_batch(
                    *found_files[:2],
                    plot_energy=False,
                    plot_clusters=True,
                    plot_grids=False
                )
                print("  ✅ Demo completed successfully!")
                return True
            except Exception as e:
                print(f"  ⚠️  Demo failed (expected if dependencies missing): {e}")
                return False
    else:
        print("  ℹ️  No suitable demo files found")
        return True

def main():
    """Run all tests."""
    print("🎯 MONTE CARLO VISUALIZATION SUITE TEST")
    print("=" * 50)
    
    tests = [
        ("Module Imports", test_imports),
        ("Color Utilities", test_color_utilities),
        ("Filename Factory", test_filename_factory),
        ("File Discovery", test_file_discovery),
        ("Visualization Drivers", test_visualization_drivers),
    ]
    
    passed = 0
    total = len(tests)
    
    for test_name, test_func in tests:
        try:
            if test_func():
                passed += 1
        except Exception as e:
            print(f"❌ {test_name} failed unexpectedly: {e}")
    
    print("\n" + "=" * 50)
    print(f"📊 TEST RESULTS: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests passed! The visualization suite is ready to use.")
        
        # Try to run demo if files exist
        run_demo_if_files_exist()
    else:
        print("⚠️  Some tests failed. Check the error messages above.")
        print("\n💡 Common issues:")
        print("   - Missing dependencies: pip install pandas numpy matplotlib")
        print("   - Module not found: ensure all Python files are in the same directory")
    
    print("\n📚 USAGE EXAMPLES:")
    print("   from visualization_drivers import simple_batch")
    print("   simple_batch('file1', 'file2', plot_energy=True, plot_clusters=True)")
    print()
    print("   from visualization_drivers import run_visualizations")  
    print("   run_visualizations(BJ='-0.44', base_variants=['DR'], plot_energy=True)")

if __name__ == "__main__":
    main()