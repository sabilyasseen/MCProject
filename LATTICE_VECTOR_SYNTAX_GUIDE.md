# Comprehensive Guide: Lattice Point Vector Saves and Grid Configuration Syntax

This guide documents all vector saves and grid configuration formats used in the Monte Carlo simulation for atomic species swapping.

## Overview

The simulation saves multiple types of data for each iteration:
- **Lattice point vectors**: Sets of atom indices involved in moves
- **Grid configuration snapshots**: Complete grid states at different move phases
- **Cluster size tracking**: Boundary and interior atom counts per cluster

---

## **LATTICE POINT VECTOR SAVES**

All lattice point vectors are saved as **comma-separated integers enclosed in quotes** in CSV format.

### **1. PossiblyChanged Vector**
- **Description**: All atoms whose local configurations were potentially affected by the move
- **Source**: `move.possibly_changed_vec`
- **Format**: `"atom1,atom2,atom3,..."`
- **When Populated**: Every iteration (both accepted and rejected moves)
- **Example**: `"45,46,78,79,123"`

### **2. ProposedAtoms Vector** 
- **Description**: All atoms that were proposed for swapping in the move (flattened from all tiers)
- **Source**: `move.primary_atoms` (flattened across all tiers)
- **Format**: `"atom1,atom2,atom3,..."`
- **When Populated**: Every iteration (both accepted and rejected moves)
- **Example**: `"45,78"` (for a single swap between atoms 45 and 78)

### **3. AcceptedAtoms Vector**
- **Description**: Primary atoms that were actually swapped (only for accepted moves)
- **Source**: `move.primary_atoms` (only if move accepted)
- **Format**: `"atom1,atom2,atom3,..."` (accepted moves) or `""` (rejected moves)
- **When Populated**: 
  - **Accepted moves**: Contains the same atoms as ProposedAtoms
  - **Rejected moves**: Empty string
- **Example**: 
  - Accepted: `"45,78"`
  - Rejected: `""`

### **4. AcceptedPossiblyChanged Vector** ⭐ *NEW*
- **Description**: Possibly changed atoms for accepted moves only
- **Source**: `move.possibly_changed_vec` (only if move accepted)
- **Format**: `"atom1,atom2,atom3,..."` (accepted moves) or `""` (rejected moves)
- **When Populated**:
  - **Accepted moves**: Contains the same atoms as PossiblyChanged
  - **Rejected moves**: Empty string
- **Example**:
  - Accepted: `"45,46,78,79,123"`
  - Rejected: `""`

---

## **GRID CONFIGURATION SAVES**

Grid configurations capture the complete state of the lattice at specific points during move processing.

### **Grid Configuration Format**
Each grid configuration is a string with the format:
```
"species|cluster_id|cell_type;species|cluster_id|cell_type;..."
```

**Per-atom format**: `species|cluster_id|cell_type`
- **species**: `0` or `1` (atomic species)
- **cluster_id**: Integer cluster ID (`-1` for unassigned)
- **cell_type**: `B` (boundary), `I` (interior), `A` (any/unassigned)
- **separator**: `;` between atoms

### **5. GridConfigStarting** ⭐ *NEW*
- **Description**: Complete grid state before any move processing begins
- **When Captured**: Before `deterministic_build_from_swap()` is called
- **Format**: Full grid configuration string (as described above)
- **Example**: `"0|2|B;1|2|B;0|1|I;1|1|I;..."`

### **6. GridConfigPotential** ⭐ *NEW*  
- **Description**: Complete grid state after move is applied but before accept/reject decision
- **When Captured**: After `deterministic_build_from_swap()` and move evaluation
- **Format**: Full grid configuration string or `""` for failed moves
- **Example**: 
  - Valid move: `"1|2|B;0|2|B;0|1|I;1|1|I;..."`
  - Failed move: `""`

---

## **CLUSTER SIZE TRACKING**

### **7. BoundarySize_ClusterN & InteriorSize_ClusterN**
- **Description**: Total boundary and interior atom counts per cluster
- **Format**: Integer values (one pair of columns per cluster)
- **Example**: For 3 clusters:
  ```
  BoundarySize_Cluster0,InteriorSize_Cluster0,BoundarySize_Cluster1,InteriorSize_Cluster1,BoundarySize_Cluster2,InteriorSize_Cluster2
  150,200,75,125,80,95
  ```

---

## **GRID CONFIGURATION FILE (GridConfig.csv)**

### **Separate Detailed Grid File**
- **Filename Pattern**: `{base_filename}_GridConfig.csv`
- **Format**: Per-atom detailed information
- **Columns**: `x,y,species,cluster_id,cell_type`
- **Example**:
  ```csv
  x,y,species,cluster_id,cell_type
  0,0,1,2,BOUNDARY
  1,0,0,2,BOUNDARY  
  2,0,1,1,INTERIOR
  ```

---

## **DATA RELATIONSHIP EXAMPLES**

### **Accepted Move Example**
```csv
PossiblyChanged,ProposedAtoms,AcceptedAtoms,AcceptedPossiblyChanged,GridConfigStarting,GridConfigPotential
"45,46,78,79","45,78","45,78","45,46,78,79","0|2|B;1|2|B;...","1|2|B;0|2|B;..."
```

### **Rejected Move Example**  
```csv
PossiblyChanged,ProposedAtoms,AcceptedAtoms,AcceptedPossiblyChanged,GridConfigStarting,GridConfigPotential
"45,46,78,79","45,78","","","0|2|B;1|2|B;...","1|2|B;0|2|B;..."
```

### **Failed Move Example**
```csv
PossiblyChanged,ProposedAtoms,AcceptedAtoms,AcceptedPossiblyChanged,GridConfigStarting,GridConfigPotential
"","","","","0|2|B;1|2|B;...","" 
```

---

## **KEY DISTINCTIONS**

| Vector Type | Always Populated | Only for Accepted | Purpose |
|-------------|-----------------|-------------------|---------|
| PossiblyChanged | ✅ | ❌ | Track all affected atoms |
| ProposedAtoms | ✅ | ❌ | Track attempted swaps |
| AcceptedAtoms | ❌ | ✅ | Track successful swaps |
| AcceptedPossiblyChanged | ❌ | ✅ | Track affected atoms in successful moves |

| Grid Config | When Captured | Purpose |
|-------------|---------------|---------|
| GridConfigStarting | Before move | Initial state |
| GridConfigPotential | After move application | Potential outcome |

---

## **USAGE NOTES**

1. **Move Tracking**: Use ProposedAtoms vs AcceptedAtoms to determine which moves were accepted
2. **Impact Analysis**: Use PossiblyChanged to analyze local configuration changes
3. **State Comparison**: Compare GridConfigStarting vs GridConfigPotential to see exact changes
4. **Success Analysis**: Use AcceptedPossiblyChanged to analyze only successful move impacts
5. **Lattice Position**: Atom indices correspond to `y * grid_size + x` positions

This comprehensive tracking system enables detailed analysis of move dynamics, acceptance patterns, and grid evolution throughout the simulation. 