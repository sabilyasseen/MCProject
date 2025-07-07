#!/usr/bin/env python3
"""
CORRECTED MOVE ANALYSIS SCRIPT
This script fixes the issues with acceptance rate normalization and creates the scatter plot visualization.

Run this script instead of the notebook to get the corrected analysis.
"""

import sys
import os
from pathlib import Path

# Add the script directory to the Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from Move_Analysis_Fixed import moves_analysis_driver, simple_batch_processor
    print("Successfully imported corrected move analysis functions")
except ImportError as e:
    print(f"Import error: {e}")
    print("Please ensure Move_Analysis_Fixed.py is in the same directory")
    sys.exit(1)

def main():
    print("=== CORRECTED MOVE ANALYSIS ===")
    print("This script provides the corrected move analysis with:")
    print("1. Fixed acceptance rate normalization by move type buckets")
    print("2. Scatter plot with green/red markers for accepted/rejected moves")
    print("3. Proper axis scaling with minimum ranges that expand for data")
    print("4. Low DPI markers for better performance")
    print("")
    
    # Example usage - replace with your actual filenames
    filenames_to_analyze = [
        "BJ=-0.44_Kawasaki_0_tol0.00_scp0.50",
        "BJ=-0.44_DR_1_tol0.50_scp0.50", 
        "BJ=-0.44_DR_2_tol0.50_scp0.50"
    ]
    
    print("Looking for files:")
    for filename in filenames_to_analyze:
        csv_file = f"{filename}.csv"
        if Path(csv_file).is_file():
            print(f"  ✓ Found: {csv_file}")
        else:
            print(f"  ✗ Not found: {csv_file}")
    print("")
    
    # Run the corrected analysis
    try:
        simple_batch_processor(*filenames_to_analyze, run_moves=True)
        print("\n=== ANALYSIS COMPLETE ===")
    except Exception as e:
        print(f"Error running analysis: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)