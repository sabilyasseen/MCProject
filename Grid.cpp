#include "Grid.hpp"
#include <random>
#include <ctime>
#include <iostream>
#include <algorithm> // For std::sort
#include <omp.h>



// Constructor and Destructor
Grid::Grid(int N, bool minsize5_boundary_bool) {
    minsize5_boundary=minsize5_boundary_bool;
    Size = N;
    num_clusters = 0;
    array = new Cell[N * N];
    int num_species=2;
    total_boundary_atoms = new int[num_species];
    total_interior_atoms = new int[num_species];
    // Initialize all totals to 0
    for (int i = 0; i < num_species; i++) {
        total_boundary_atoms[i] = 0;
        total_interior_atoms[i] = 0;
    }
    // Initialize random number generator
    gen = std::mt19937(rd());
    species_dist = std::uniform_int_distribution<>(0, 1);
    
    // Step 1: Initialize a vector of all indices and set total number of cells
    int total_cells = N * N;
    std::vector<int> cell_indices(total_cells);
    for (int i = 0; i < total_cells; i++) {
        cell_indices[i] = i;
    }
    
    // Step 2: Shuffle the indices randomly
    std::shuffle(cell_indices.begin(), cell_indices.end(), gen);
    
    // Step 3: Assign first half to species 0, second half to species 1
    int half_cells = total_cells / 2;
    for (int i = 0; i < total_cells; i++) {
        int idx = cell_indices[i];
        int x = idx % N;
        int y = idx / N;
        
        // Assign species 0 to first half, species 1 to second half
        int species = (i < half_cells) ? 0 : 1;
        array[idx] = Cell(x, y, N, species);
    }
    
    // Step 4: Set up local configurations for all cells
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int current_idx = j * N + i;
            Cell& current_cell = array[current_idx];
            
            // Initialize local_config with center species
            current_cell.local_config.clear();
            current_cell.local_config.push_back(current_cell.cell_species);
            
            // Add up neighbor (j+1 with wraparound)
            int up_idx;
            if (j == N-1) {
                up_idx = 0 * N + i;  // Wrap to bottom (j=0)
            } else {
                up_idx = (j+1) * N + i;  // Normal case: one row up
            }
            current_cell.local_config.push_back(array[up_idx].cell_species);
            current_cell.neighbor_indexes.push_back(up_idx);
            
            // Add right neighbor
            int right_idx;
            if (i == N-1) {
                right_idx = j * N + 0;  // Wrap to left
            } else {
                right_idx = j * N + (i+1);
            }
            current_cell.local_config.push_back(array[right_idx].cell_species);
            current_cell.neighbor_indexes.push_back(right_idx);
            
            // Add down neighbor (j-1 with wraparound)
            int down_idx;
            if (j == 0) {
                down_idx = (N-1) * N + i;  // Wrap to top (j=N-1)
            } else {
                down_idx = (j-1) * N + i;  // Normal case: one row down
            }
            current_cell.local_config.push_back(array[down_idx].cell_species);
            current_cell.neighbor_indexes.push_back(down_idx);
            
            // Add left neighbor
            int left_idx;
            if (i == 0) {
                left_idx = j * N + (N-1);  // Wrap to right
            } else {
                left_idx = j * N + (i-1);
            }
            current_cell.local_config.push_back(array[left_idx].cell_species);
            current_cell.neighbor_indexes.push_back(left_idx);
            
            // Set neighbor indexes and pointers
        }
    }
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            int current_idx = j * N + i;
            Cell& current_cell = array[current_idx];
            current_cell.expansion_indexes = current_cell.neighbor_indexes;
        }
    }
    
    // Initialize lattice site filtering vectors
    // Create index 0 as the default case: all lattice sites in interior, empty boundary
    interior_lattice_sites.resize(1);
    boundary_lattice_sites.resize(1);
    
    // Add all lattice sites to interior set at index 0
    for (int i = 0; i < total_cells; i++) {
        interior_lattice_sites[0].insert(i);
    }
    
    // Leave boundary set at index 0 empty (initialized by default)
    // This ensures lookups in boundary_lattice_sites[0] return false for all atoms
}

Grid::~Grid() {
    delete[] array;
    delete[] total_boundary_atoms;
    delete[] total_interior_atoms;
}

// Grid operations
void Grid::add_cluster(int minsize) {
    // We have 2 species in this simulation (species 0 and 1)
    num_clusters=clusters.size();
    Cluster newcluster(num_clusters, minsize);
    clusters.push_back(newcluster);
    num_clusters=clusters.size();
}

void Grid::add_restricted_cluster(int minsize) {
    // Create a new restricted cluster with the same functionality as add_cluster
    // Use the current size of restricted_clusters as the ID
    int restricted_cluster_id = restricted_clusters.size();
    Cluster newcluster(restricted_cluster_id, minsize);
    restricted_clusters.push_back(newcluster);
}

    void Grid::add_atom(int atom, int cluster_index, CellType type) {
        int species = array[atom].cell_species; // Species is now directly 0 or 1, no adjustment needed
        
        if (type == CellType::BOUNDARY) {
            clusters[cluster_index].add_atom_boundary(atom, species);
            adjust_total_boundary(species, 1); // Increment boundary counter
        } else if (type == CellType::INTERIOR) {
            clusters[cluster_index].add_atom_interior(atom, species);
            adjust_total_interior(species, 1); // Increment interior counter
        }
        array[atom].cell_cluster_ID = cluster_index;
        array[atom].cell_type = type;
        array[atom].min_size = clusters[cluster_index].property.min_size;
    }

void Grid::remove_atom(int atom) {
    // Get current state from the cell
    int cluster_index = array[atom].cell_cluster_ID;
    CellType type = array[atom].cell_type;
    int species = array[atom].cell_species; // Species is now directly 0 or 1
    
    array[atom].min_size = 0;
    // Only proceed if the atom is actually assigned to a cluster
    if (cluster_index != -1) {
        if (type == CellType::BOUNDARY) {
            clusters[cluster_index].remove_atom_boundary(atom, species);
            adjust_total_boundary(species, -1); // Decrement boundary counter
        } else if (type == CellType::INTERIOR) {
            clusters[cluster_index].remove_atom_interior(atom, species);
            adjust_total_interior(species, -1); // Decrement interior counter
        }
        // Reset the atom's state
        array[atom].cell_cluster_ID = -1;
        array[atom].cell_type = CellType::CELL_ANY;
    }
}
void Grid::assign_initial_clusters() {
    // Step 1: Clear all existing cluster assignments and reset cluster sizes
    for (int i = 0; i < Size*Size; i++) {
        array[i].cell_cluster_ID = -1;
        array[i].cell_type = CellType::CELL_ANY;
        array[i].min_size = 0;
    }

    std::vector<int> changed_atoms_vec(Size*Size);
    for(int i = 0; i < Size*Size; i++){
        changed_atoms_vec[i]=i;
    }

    reclassify_cluster_status(changed_atoms_vec, 0);
    reclassify_boundary_status_minsize5(changed_atoms_vec, 0);
    reclassify_boundary_status_minsize2(changed_atoms_vec, 0);
    reclassify_boundary_status(changed_atoms_vec, 0);

    for(int i = 0; i < num_clusters; i++){
        clusters[i].requires_refresh=true;
    }

    for(int i = 0; i < Size*Size; i++){
        array[i].accept_cell();
    }
}


void Grid::apply_cluster_change(std::vector<int> atom_list, int filter_index) {
    for(int i = 0; i < num_clusters; i++){
        clusters[i].requires_refresh = true;
    }
    for(int atom : atom_list){
        // Filter atoms based on lattice site assignments if filter_index is valid
        if(filter_index >= 0) {
            bool atom_in_interior = false;
            bool atom_in_boundary = false;
            
            // Check if atom is in the interior lattice sites at this index
            if(filter_index < interior_lattice_sites.size()) {
                atom_in_interior = interior_lattice_sites[filter_index].find(atom) != interior_lattice_sites[filter_index].end();
            }
            
            // Check if atom is in the boundary lattice sites at this index
            if(filter_index < boundary_lattice_sites.size()) {
                atom_in_boundary = boundary_lattice_sites[filter_index].find(atom) != boundary_lattice_sites[filter_index].end();
            }
            
            // Skip processing this atom if it's not in either set for this filter index
            if(!atom_in_interior && !atom_in_boundary) {
                continue;
            }
        }
        
        array[atom].accept_cell();
    }
}

void Grid::undo_cluster_changes(std::vector<int> atom_list, int filter_index) {        
    for(int atom: atom_list) {
        bool atom_in_interior = false;
        bool atom_in_boundary = false;
        
        // Filter atoms based on lattice site assignments if filter_index is valid
        if(filter_index >= 0) {
            // Check if atom is in the interior lattice sites at this index
            if(filter_index < interior_lattice_sites.size()) {
                atom_in_interior = interior_lattice_sites[filter_index].find(atom) != interior_lattice_sites[filter_index].end();
            }
            
            // Check if atom is in the boundary lattice sites at this index
            if(filter_index < boundary_lattice_sites.size()) {
                atom_in_boundary = boundary_lattice_sites[filter_index].find(atom) != boundary_lattice_sites[filter_index].end();
            }
            
            // Skip processing this atom if it's not in either set for this filter index
            if(!atom_in_interior && !atom_in_boundary) {
                continue;
            }
        } else {
            // When no filtering, treat as if atom is in interior (for add_atom logic)
            atom_in_interior = true;
        }
        
        // Remove from current cluster if any
        remove_atom(atom);
        array[atom].revert_cell();
        
        // Only add back to previous cluster if atom is in interior set
        if(atom_in_interior) {
            add_atom(atom, array[atom].cell_cluster_ID, array[atom].cell_type);
        }
    }
}
void Grid::change_specie(int atom) {
    if(array[atom].cell_species == 0){
        array[atom].cell_species = 1;
    } else if(array[atom].cell_species == 1){
        array[atom].cell_species = 0;
    } 
    int new_specie = array[atom].cell_species;

    array[atom].local_config[0] = new_specie;

    // Get the cell's coordinates
    for(int neighbor: array[atom].neighbor_indexes){
        refresh_local_config(neighbor);
    }
}

void Grid::refresh_local_config(int atom){
    for(int i = 0; i < 4; i++){
        array[atom].local_config[i+1]=array[array[atom].neighbor_indexes[i]].cell_species;
    }
}

bool Grid::any_neighbors_id(int atom,int cluster_id){
    for(int neighbor: array[atom].neighbor_indexes){
        if(array[neighbor].cell_cluster_ID==cluster_id){
            return true;
        }
    }
    return false;
}
bool Grid::all_neighbors_id(int atom,int cluster_id){
    for(int neighbor: array[atom].neighbor_indexes){
        if(array[neighbor].cell_cluster_ID!=cluster_id){
            return false;
        }
    }
    return true;
}
bool Grid::any_neighbors_property(int atom,int cluster_id){
    for(int neighbor: array[atom].neighbor_indexes){
        if(clusters[cluster_id].satisfies_properties(array[neighbor].local_config)){
            return true;
        }
    }
    return false;
}
bool Grid::all_neighbors_property(int atom,int cluster_id){
    for(int neighbor: array[atom].neighbor_indexes){
        if(!clusters[cluster_id].satisfies_properties(array[neighbor].local_config)){
            return false;
        }
    }
    return true;
}

// -----------------------------------------------------------------------------
// New overload: assign_initial_clusters using a subset of lattice sites
// -----------------------------------------------------------------------------
void Grid::assign_initial_clusters(const std::vector<int>& allowed_indices) {
    // If the caller passes an empty vector, fall back to full-grid behaviour
    if (allowed_indices.empty()) {
        assign_initial_clusters();
        return;
    }

    // Step 1: Reset cell assignment state for the *whole* grid so that any
    // previous cluster configuration does not leak into the new one.
    // NOTE: this mirrors the behaviour of the original implementation.
    for (int i = 0; i < Size * Size; ++i) {
        array[i].cell_cluster_ID = -1;
        array[i].cell_type       = CellType::CELL_ANY;
        array[i].min_size        = 0;
    }

    // Reset global counters
    const int num_species = 2; // hard-coded in current simulation
    for (int s = 0; s < num_species; ++s) {
        total_boundary_atoms[s] = 0;
        total_interior_atoms[s] = 0;
    }

    // Step 2: Prepare the list of atoms that will be (re)classified – only the
    // allowed subset.
    std::vector<int> changed_atoms_vec = allowed_indices;

    // Step 3: Reclassify cluster & boundary status for the subset.
    reclassify_cluster_status(changed_atoms_vec, 0);
    reclassify_boundary_status_minsize5(changed_atoms_vec, 0);
    reclassify_boundary_status_minsize2(changed_atoms_vec, 0);
    reclassify_boundary_status(changed_atoms_vec, 0);

    // Step 4: Mark clusters as needing a refresh of their cached vectors.
    for (int i = 0; i < num_clusters; ++i) {
        clusters[i].requires_refresh = true;
    }

    // Step 5: Accept the tentative cell modifications for the allowed set.
    for (int idx : allowed_indices) {
        if (idx >= 0 && idx < Size * Size) {
            array[idx].accept_cell();
        }
    }
}

void Grid::adjust_total_boundary(int specie, int increment) {
    #pragma omp atomic
    total_boundary_atoms[specie] += increment;
}

void Grid::adjust_total_interior(int specie, int increment) {
    #pragma omp atomic 
    total_interior_atoms[specie] += increment;
}

bool Grid::reclassify_cluster_status(std::vector<int> atom_list, int filter_index){
    #pragma omp parallel for
    for(int atom: atom_list){
        // Filter atoms based on lattice site assignments if filter_index is valid
        if(filter_index >= 0) {
            bool atom_in_interior = false;
            bool atom_in_boundary = false;
            
            // Check if atom is in the interior lattice sites at this index
            if(filter_index < interior_lattice_sites.size()) {
                atom_in_interior = interior_lattice_sites[filter_index].find(atom) != interior_lattice_sites[filter_index].end();
            }
            
            // Check if atom is in the boundary lattice sites at this index
            if(filter_index < boundary_lattice_sites.size()) {
                atom_in_boundary = boundary_lattice_sites[filter_index].find(atom) != boundary_lattice_sites[filter_index].end();
            }
            
            // Skip processing this atom if it's not in either set for this filter index
            if(!atom_in_interior && !atom_in_boundary) {
                continue;
            }
        }
        
        remove_atom(atom);
        
        // ------------------------------------------------------------------
        // Check if this atom's local configuration violates any restricted pattern
        // (Check BEFORE positive cluster assignment)
        // ------------------------------------------------------------------
        const std::vector<int> &cfg = array[atom].local_config;
        for ( auto &rc : restricted_clusters) {
            if (rc.satisfies_properties(cfg)) {
                return false; // Violation detected - exit immediately
            }
        }
        
        for(int j = num_clusters-1; j >= 0; j--){
            if(array[atom].cell_cluster_ID == -1){
                // Get the minsize for this cluster
                
                // Check if atom itself satisfies the property
                if(clusters[j].satisfies_properties(array[atom].local_config)){
                    array[atom].cell_cluster_ID = j;
                    array[atom].cell_type = CellType::INTERIOR;
                    array[atom].min_size = clusters[j].property.min_size;
                    break;
                }
                if(clusters[j].property.min_size==8){
                    if(any_neighbors_property(atom,j)){
                        array[atom].cell_cluster_ID = j;
                        array[atom].cell_type = CellType::BOUNDARY;
                        array[atom].min_size = clusters[j].property.min_size;
                        break;
                    }
                }
                
                // For minsize 5, do additional neighbor check

            }
        }   
    }

    return true; // No violations detected
}
void Grid::reclassify_boundary_status_minsize5(std::vector<int> atom_list, int filter_index){
    if(!minsize5_boundary)return;
    #pragma omp parallel for
    for(int atom: atom_list){
        // Filter atoms based on lattice site assignments if filter_index is valid
        if(filter_index >= 0) {
            bool atom_in_interior = false;
            bool atom_in_boundary = false;
            
            // Check if atom is in the interior lattice sites at this index
            if(filter_index < interior_lattice_sites.size()) {
                atom_in_interior = interior_lattice_sites[filter_index].find(atom) != interior_lattice_sites[filter_index].end();
            }
            
            // Check if atom is in the boundary lattice sites at this index
            if(filter_index < boundary_lattice_sites.size()) {
                atom_in_boundary = boundary_lattice_sites[filter_index].find(atom) != boundary_lattice_sites[filter_index].end();
            }
            
            // Skip processing this atom if it's not in either set for this filter index
            if(!atom_in_interior && !atom_in_boundary) {
                continue;
            }
        }

        for(int neighbor: array[atom].neighbor_indexes){
            if(array[neighbor].min_size==5&& array[neighbor].cell_type==CellType::INTERIOR && (array[neighbor].cell_cluster_ID>array[atom].cell_cluster_ID)){
               array[atom].cell_cluster_ID = array[neighbor].cell_cluster_ID;
                array[atom].cell_type = CellType::BOUNDARY;
                array[atom].min_size = 5;
                break;
            }
        }
        
    }
}
void Grid::reclassify_boundary_status_minsize2(std::vector<int> atom_list, int filter_index){
    #pragma omp parallel for
    for(int atom: atom_list){
        // Filter atoms based on lattice site assignments if filter_index is valid
        if(filter_index >= 0) {
            bool atom_in_interior = false;
            bool atom_in_boundary = false;
            
            // Check if atom is in the interior lattice sites at this index
            if(filter_index < interior_lattice_sites.size()) {
                atom_in_interior = interior_lattice_sites[filter_index].find(atom) != interior_lattice_sites[filter_index].end();
            }
            
            // Check if atom is in the boundary lattice sites at this index
            if(filter_index < boundary_lattice_sites.size()) {
                atom_in_boundary = boundary_lattice_sites[filter_index].find(atom) != boundary_lattice_sites[filter_index].end();
            }
            
            // Skip processing this atom if it's not in either set for this filter index
            if(!atom_in_interior && !atom_in_boundary) {
                continue;
            }
        }
        
        // Add minsize2-specific logic here if needed
    }
}

void Grid::reclassify_boundary_status(std::vector<int> atom_list, int filter_index){
    #pragma omp parallel for
    for(int atom: atom_list){
        // Filter atoms based on lattice site assignments if filter_index is valid
        if(filter_index >= 0) {
            bool atom_in_interior = false;
            
            // Check if atom is in the interior lattice sites at this index
            if(filter_index < interior_lattice_sites.size()) {
                atom_in_interior = interior_lattice_sites[filter_index].find(atom) != interior_lattice_sites[filter_index].end();
            }
            
            // Skip processing this atom if it's not in the interior set for this filter index
            if(!atom_in_interior) {
                continue;
            }
        }
        
        if(all_neighbors_id(atom,array[atom].cell_cluster_ID) && array[atom].cell_type!=CellType::BOUNDARY){
            add_atom(atom, array[atom].cell_cluster_ID, CellType::INTERIOR);
        }   
        else{
            add_atom(atom, array[atom].cell_cluster_ID, CellType::BOUNDARY);
        }
    }
}