# How to Fix the Move Analysis in Visualizer_Revamped.ipynb

## Issue Summary
The original `moves_analysis_driver` function had several problems:
1. **Broken acceptance rate calculation** - not normalized by move type bucket sizes
2. **Too restrictive axis ranges** - didn't allow expansion for data
3. **Missing scatter plot** - no visualization of accepted vs rejected moves by atoms and energy

## Manual Fix Instructions

### Step 1: Open Visualizer_Revamped.ipynb
Find the `moves_analysis_driver` function in Cell 1 (around line 120-160)

### Step 2: Replace the entire function
Replace the existing `moves_analysis_driver` function with this corrected version:

```python
def moves_analysis_driver(csv_files, table_labels, curve_colours):
    """
    CORRECTED VERSION: Fixed acceptance rate normalization and scatter plot visualization
    """
    for csv_file, table_label, color in zip(csv_files, table_labels, curve_colours):
        if not Path(csv_file).is_file(): continue
        
        try:
            df = pd.read_csv(csv_file)
            required_cols = ["Total_Acceptances", "DeltaE", "NumAtomsSwapped"]
            if not all(col in df.columns for col in required_cols): continue
            
            # Calculate acceptance status for each move
            df['IsAccepted'] = (df['Total_Acceptances'].diff().fillna(0) > 0).astype(int)
            
            # Calculate normalized acceptance rates by move type buckets
            # Group by number of atoms swapped to create buckets
            bucket_stats = df.groupby('NumAtomsSwapped').agg({
                'IsAccepted': ['sum', 'count'],
                'DeltaE': ['mean', 'std']
            }).round(4)
            
            # Flatten column names
            bucket_stats.columns = ['_'.join(col).strip() for col in bucket_stats.columns.values]
            bucket_stats['AcceptanceRate'] = (bucket_stats['IsAccepted_sum'] / bucket_stats['IsAccepted_count']).fillna(0)
            
            # Create figure with scatter plot as main focus
            fig, axes = plt.subplots(2, 2, figsize=(16, 12))
            fig.suptitle(f'Move Analysis: {table_label}', fontsize=16, fontweight='bold')
            
            # Main scatter plot: Number of atoms vs DeltaE with acceptance status
            accepted_mask = df['IsAccepted'] == 1
            rejected_mask = df['IsAccepted'] == 0
            
            # Use low DPI markers for better performance
            axes[0, 0].scatter(df[accepted_mask]['NumAtomsSwapped'], df[accepted_mask]['DeltaE'], 
                             c='green', alpha=0.6, s=8, label='Accepted', rasterized=True)
            axes[0, 0].scatter(df[rejected_mask]['NumAtomsSwapped'], df[rejected_mask]['DeltaE'], 
                             c='red', alpha=0.6, s=8, label='Rejected', rasterized=True)
            
            # Set appropriate axis limits with minimum ranges but allow expansion for data
            min_atoms = df['NumAtomsSwapped'].min()
            max_atoms = df['NumAtomsSwapped'].max()
            atom_range = max_atoms - min_atoms
            x_min = max(0, min_atoms - 0.1 * atom_range) if atom_range > 0 else 0
            x_max = max_atoms + 0.1 * atom_range if atom_range > 0 else max(10, max_atoms)
            
            min_deltaE = df['DeltaE'].min()
            max_deltaE = df['DeltaE'].max()
            deltaE_range = max_deltaE - min_deltaE
            y_min = min_deltaE - 0.1 * deltaE_range if deltaE_range > 0 else min_deltaE - 1
            y_max = max_deltaE + 0.1 * deltaE_range if deltaE_range > 0 else max_deltaE + 1
            
            axes[0, 0].set_xlim(x_min, x_max)
            axes[0, 0].set_ylim(y_min, y_max)
            axes[0, 0].set_xlabel('Number of Atoms Swapped')
            axes[0, 0].set_ylabel('ΔE (Energy Change)')
            axes[0, 0].set_title('Move Acceptance by Atoms Swapped vs Energy Change')
            axes[0, 0].legend()
            axes[0, 0].grid(True, alpha=0.3)
            
            # Normalized acceptance rates by move type buckets
            if len(bucket_stats) > 0:
                bucket_indices = bucket_stats.index
                acceptance_rates = bucket_stats['AcceptanceRate']
                
                axes[0, 1].bar(bucket_indices, acceptance_rates, color=color, alpha=0.7)
                axes[0, 1].set_xlabel('Number of Atoms Swapped (Move Type)')
                axes[0, 1].set_ylabel('Normalized Acceptance Rate')
                axes[0, 1].set_title('Acceptance Rate by Move Type Bucket')
                axes[0, 1].set_ylim(0, 1)
                axes[0, 1].grid(True, alpha=0.3)
                
                # Add count labels on bars
                for idx, rate in zip(bucket_indices, acceptance_rates):
                    count = bucket_stats.loc[idx, 'IsAccepted_count']
                    axes[0, 1].text(idx, rate + 0.02, f'n={count}', ha='center', va='bottom', fontsize=8)
            
            # Energy change distribution
            axes[1, 0].hist(df['DeltaE'], bins=50, alpha=0.7, color=color, edgecolor='black')
            axes[1, 0].set_xlabel('ΔE (Energy Change)')
            axes[1, 0].set_ylabel('Frequency')
            axes[1, 0].set_title('Energy Change Distribution')
            axes[1, 0].grid(True, alpha=0.3)
            
            # Running acceptance rate (properly normalized)
            iterations = np.arange(len(df))
            running_acceptance = df['Total_Acceptances'] / (iterations + 1)
            axes[1, 1].plot(iterations, running_acceptance, color=color, linewidth=2)
            axes[1, 1].set_xlabel('Iteration')
            axes[1, 1].set_ylabel('Cumulative Acceptance Rate')
            axes[1, 1].set_title('Running Acceptance Rate')
            axes[1, 1].set_ylim(0, 1)
            axes[1, 1].grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.show()
            
            # Print summary statistics
            print(f"\n=== Move Analysis Summary for {table_label} ===")
            print(f"Total moves: {len(df)}")
            print(f"Total accepted: {df['IsAccepted'].sum()}")
            print(f"Overall acceptance rate: {df['IsAccepted'].mean():.4f}")
            print(f"Move type buckets (by atoms swapped):")
            for idx, row in bucket_stats.iterrows():
                print(f"  {idx} atoms: {row['AcceptanceRate']:.4f} rate ({row['IsAccepted_sum']:.0f}/{row['IsAccepted_count']:.0f})")
            
        except Exception as e:
            print(f"Error processing {csv_file}: {e}")
            continue
```

## Alternative: Use the Corrected Scripts

If you don't want to manually edit the notebook, you can use the corrected scripts I created:

### Option 1: Run the corrected analysis directly
```bash
python run_corrected_move_analysis.py
```

### Option 2: Import the corrected functions
```python
from Move_Analysis_Fixed import moves_analysis_driver, simple_batch_processor

# Use the corrected functions
filenames = ["your_file1", "your_file2"]
simple_batch_processor(*filenames, run_moves=True)
```

## Key Improvements in the Corrected Version

1. **Fixed Acceptance Rate Normalization**: 
   - Groups moves by number of atoms swapped (move type buckets)
   - Calculates acceptance rate for each bucket separately
   - Prevents artificially high rates caused by uneven move type distributions

2. **Proper Axis Scaling**:
   - Dynamically calculates appropriate ranges based on actual data
   - Sets minimum ranges but allows expansion to fit all data points
   - Uses 10% buffer around data range for better visualization

3. **Scatter Plot Visualization**:
   - Main plot shows number of atoms (x-axis) vs ΔE (y-axis)
   - Green markers for accepted moves, red for rejected
   - Low DPI markers (s=8, rasterized=True) for better performance

4. **Additional Analysis**:
   - Bar chart showing normalized acceptance rates by move type bucket
   - Energy change distribution histogram
   - Properly normalized running acceptance rate
   - Detailed summary statistics

## Files Created

- `Move_Analysis_Fixed.py` - Standalone corrected analysis functions
- `run_corrected_move_analysis.py` - Simple script to run the corrected analysis
- `NOTEBOOK_FIX_GUIDE.md` - This guide

Use whichever approach works best for your workflow!