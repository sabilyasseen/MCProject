#include "Sim.hpp"
#include "globals.hpp"
#include "Move.hpp"
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <unordered_set>
#include <unordered_map>
#include <mutex>

// Custom hash function for std::vector<std::vector<int>>
struct VectorVectorHash {
    std::size_t operator()(const std::vector<std::vector<int>>& vv) const {
        std::size_t seed = vv.size();
        for (const auto& v : vv) {
            std::size_t h = 0;
            for (int i : v) {
                h ^= std::hash<int>{}(i) + 0x9e3779b9 + (h << 6) + (h >> 2);
            }
            seed ^= h + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

// Global map for tested cluster configurations and their energies
std::unordered_map<std::vector<std::vector<int>>, float, VectorVectorHash> tested_configurations;

// Function declarations
bool test_cluster_set_consistency(Grid& grid);
bool test_interior_cell_properties(Grid& grid);
bool test_boundary_cell_properties(Grid& grid);
bool test_interior_cell_neighbors(Grid& grid);
bool test_boundary_cell_neighbors(Grid& grid);
bool test_cluster_sizes(Grid& grid);
bool test_cluster_membership(Grid& grid);
bool test_neighbor_count(Grid& grid);
bool test_local_config_consistency(Grid& grid);
bool test_grid_total_consistency(Grid& grid);
void validate_all_cluster_properties(Grid& grid);
void validate_without_totals(Grid& grid);
std::string run_and_save(std::string selection_string, int iterations, int cluster_type, bool overload, float same_cluster_prob, float expand, float tolerance_override);
void display_cluster_composition(const Grid& grid);
std::string run_and_save_delayed_rejection(std::string selection_string, int iterations, int cluster_type, int max_num_swaps, bool overload, float same_cluster_prob, float expand, float tolerance_override);

// Add forward declaration for new helper
std::vector<std::vector<int>> build_cluster_tree_(const std::string& selection_string, 
                                    int min_clusters, int number_of_clusters,
                                    int min_properties_per_cluster, int max_properties_per_cluster,
                                    float cluster_bias, float same_cluster_prob, 
                                    int depth, int final_test_runs, bool use_minsize5_boundary,
                                    bool use_rotation_only = false);

// Function to prune and standardize a cluster conditions matrix
std::vector<std::vector<int>> prune(const std::vector<std::vector<int>>& input_matrix) {
    if (input_matrix.empty()) {
        return input_matrix;
    }
    
    std::vector<std::vector<int>> result = input_matrix;
    
    // Step 1: Process each vector - sort and remove duplicates
    for (auto& cluster : result) {
        if (!cluster.empty()) {
            // Sort the cluster
            std::sort(cluster.begin(), cluster.end());
            // Remove duplicates within this cluster
            cluster.erase(std::unique(cluster.begin(), cluster.end()), cluster.end());
        }
    }
    
    // Step 2: Remove duplicate rows across the entire matrix
    for (size_t i = 0; i < result.size(); ++i) {
        for (size_t j = i + 1; j < result.size(); ) {
            if (result[i] == result[j]) {
                result.erase(result.begin() + j);
            } else {
                ++j;
            }
        }
    }
    
    // Step 3: Sort clusters by their first element for consistent ordering
    std::sort(result.begin(), result.end(), [](const std::vector<int>& a, const std::vector<int>& b) {
        if (a.empty() && b.empty()) return false;
        if (a.empty()) return true;
        if (b.empty()) return false;
        return a[0] < b[0];
    });
    
    return result;
}

// Test individual cluster properties
bool test_cluster_set_consistency(Grid& grid) {
    int boundary_errors = 0;
    int interior_errors = 0;
    
    for (int i = 0; i < grid.Size * grid.Size; i++) {
        Cell& cell = grid.array[i];
        if (cell.cell_cluster_ID >= 0 && cell.cell_cluster_ID < grid.num_clusters) {
            Cluster& cluster = grid.clusters[cell.cell_cluster_ID];
            bool in_boundary = cluster.boundary_cells[cell.cell_species].find(i) != cluster.boundary_cells[cell.cell_species].end();
            bool in_interior = cluster.interior_cells[cell.cell_species].find(i) != cluster.interior_cells[cell.cell_species].end();
            
            if (cell.cell_type == CellType::BOUNDARY && !in_boundary) {
                boundary_errors++;
            }
            if (cell.cell_type == CellType::INTERIOR && !in_interior) {
                interior_errors++;
            }
        }
    }
    
    bool valid = (boundary_errors == 0 && interior_errors == 0);
    if (!valid) {
        if (boundary_errors > 0) {
            std::cout << "FAILED: Cluster set consistency - " << boundary_errors 
                      << " cells marked as boundary but not in their cluster's boundary set" << std::endl;
        }
        if (interior_errors > 0) {
            std::cout << "FAILED: Cluster set consistency - " << interior_errors 
                      << " cells marked as interior but not in their cluster's interior set" << std::endl;
        }
    }
    return valid;
}

bool test_interior_cell_properties(Grid& grid) {
    int property_errors_minsize1 = 0;
    int property_errors_minsize5 = 0;
    int total_interior_minsize1 = 0;
    int total_interior_minsize5 = 0;
    
    for (int i = 0; i < grid.Size * grid.Size; i++) {
        Cell& cell = grid.array[i];
        if (cell.cell_cluster_ID >= 0 && cell.cell_type == CellType::INTERIOR) {
            Cluster& cluster = grid.clusters[cell.cell_cluster_ID];
            int minsize = cluster.property.min_size;
            
            // Count total interior cells by minsize
            if (minsize == 1) {
                total_interior_minsize1++;
            } else if (minsize == 5) {
                total_interior_minsize5++;
            }
            
            // Check if cell satisfies its cluster's property
            if (!cluster.satisfies_properties(cell.local_config)) {
                // Track errors by minsize
                if (minsize == 1) {
                    property_errors_minsize1++;
                } else if (minsize == 5) {
                    if(!grid.any_neighbors_property(cell.index, cell.cell_cluster_ID)) {
                        property_errors_minsize5++;
                    }
                }
            }
        }
    }
    
    bool valid = (property_errors_minsize1 == 0 && property_errors_minsize5 == 0);
    if (!valid) {
        std::cout << "FAILED: Interior cell properties" << std::endl;
        if (property_errors_minsize1 > 0) {
            std::cout << "  Minsize 1 clusters: " << property_errors_minsize1 
                      << " errors out of " << total_interior_minsize1 << " interior cells" << std::endl;
        }
        if (property_errors_minsize5 > 0) {
            std::cout << "  Minsize 5 clusters: " << property_errors_minsize5 
                      << " errors out of " << total_interior_minsize5 << " interior cells" << std::endl;
        }
    }
    return valid;
}

bool test_boundary_cell_properties(Grid& grid) {
    int property_errors_minsize1 = 0;
    int property_errors_minsize5 = 0;
    int total_boundary_minsize1 = 0;
    int total_boundary_minsize5 = 0;
    
    for (int i = 0; i < grid.Size * grid.Size; i++) {
        Cell& cell = grid.array[i];
        if (cell.cell_cluster_ID >= 0 && cell.cell_type == CellType::BOUNDARY) {
            Cluster& cluster = grid.clusters[cell.cell_cluster_ID];
            int minsize = cluster.property.min_size;
            
            // Count total boundary cells by minsize
            if (minsize == 1) {
                total_boundary_minsize1++;
                
                // For minsize 1, boundary cells should still satisfy the cluster property
                if (!cluster.satisfies_properties(cell.local_config)) {
                    property_errors_minsize1++;
                }
            } else if (minsize == 5) {
                total_boundary_minsize5++;
                
                // For minsize 5, check if any neighbors satisfy the property
                bool any_neighbor_satisfies = false;
                for (int neighbor_idx : cell.neighbor_indexes) {
                    if (grid.array[neighbor_idx].cell_cluster_ID == cell.cell_cluster_ID &&
                        cluster.satisfies_properties(grid.array[neighbor_idx].local_config)) {
                        any_neighbor_satisfies = true;
                        break;
                    }
                }
                
                // Only count as error if no neighbors satisfy the property
                if (!any_neighbor_satisfies && !grid.clusters[cell.cell_cluster_ID].satisfies_properties(cell.local_config)) {
                    property_errors_minsize5++;
                }
            }
        }
    }
    
    // Only minsize 1 errors count as test failures
    bool valid = (property_errors_minsize1 == 0);
    if (!valid) {
        std::cout << "FAILED: Boundary cell properties" << std::endl;
        if (property_errors_minsize1 > 0) {
            std::cout << "  Minsize 1 clusters: " << property_errors_minsize1 
                      << " errors out of " << total_boundary_minsize1 << " boundary cells" << std::endl;
        }
    }
    
    // For minsize 5, report errors but don't fail the test
    if (property_errors_minsize5 > 0) {
        std::cout << "WARNING: " << property_errors_minsize5 
                  << " minsize 5 boundary cells have no neighbors or itself satisfying property (out of " 
                  << total_boundary_minsize5 << " total)" << std::endl;
    }
    
    return valid;
}

bool test_interior_cell_neighbors(Grid& grid) {
    int neighbor_errors_minsize1 = 0;
    int neighbor_errors_minsize5 = 0;
    int total_interior_minsize1 = 0;
    int total_interior_minsize5 = 0;
    
    for (int i = 0; i < grid.Size * grid.Size; i++) {
        Cell& cell = grid.array[i];
        if (cell.cell_cluster_ID >= 0 && cell.cell_type == CellType::INTERIOR) {
            int minsize = grid.clusters[cell.cell_cluster_ID].property.min_size;
            
            // Count total interior cells by minsize
            if (minsize == 1) {
                total_interior_minsize1++;
            } else if (minsize == 5) {
                total_interior_minsize5++;
            }
            
            // Check that all neighbors have the same cluster ID
            for (int neighbor_idx : cell.neighbor_indexes) {
                if (grid.array[neighbor_idx].cell_cluster_ID != cell.cell_cluster_ID) {
                    if (minsize == 1) {
                        neighbor_errors_minsize1++;
                    } else if (minsize == 5) {
                        neighbor_errors_minsize5++;
                    }
                    break; // Count each problematic cell once
                }
            }
        }
    }
    
    bool valid = (neighbor_errors_minsize1 == 0 && neighbor_errors_minsize5 == 0);
    if (!valid) {
        std::cout << "FAILED: Interior cell neighbors" << std::endl;
        if (neighbor_errors_minsize1 > 0) {
            std::cout << "  Minsize 1 clusters: " << neighbor_errors_minsize1 
                      << " errors out of " << total_interior_minsize1 << " interior cells" << std::endl;
        }
        if (neighbor_errors_minsize5 > 0) {
            std::cout << "  Minsize 5 clusters: " << neighbor_errors_minsize5
                      << " errors out of " << total_interior_minsize5 << " interior cells" << std::endl;
        }
    }
    return valid;
}

bool test_boundary_cell_neighbors(Grid& grid) {
    int neighbor_errors_minsize1 = 0;
    int neighbor_errors_minsize5 = 0;
    int total_boundary_minsize1 = 0;
    int total_boundary_minsize5 = 0;
    
    for (int i = 0; i < grid.Size * grid.Size; i++) {
        Cell& cell = grid.array[i];
        if (cell.cell_cluster_ID >= 0 && cell.cell_type == CellType::BOUNDARY) {
            int minsize = grid.clusters[cell.cell_cluster_ID].property.min_size;
            
            // Count total boundary cells by minsize
            if (minsize == 1) {
                total_boundary_minsize1++;
            } else if (minsize == 5) {
                total_boundary_minsize5++;
            }
            
            // Check that at least one neighbor has a different cluster ID
            bool has_different_cluster_neighbor = false;
            for (int neighbor_idx : cell.neighbor_indexes) {
                if (grid.array[neighbor_idx].cell_cluster_ID != cell.cell_cluster_ID) {
                    has_different_cluster_neighbor = true;
                    break;
                }
            }
            
            if (!has_different_cluster_neighbor) {
                if (minsize == 1) {
                    neighbor_errors_minsize1++;
                } else if (minsize == 5) {
                    if(grid.clusters[cell.cell_cluster_ID].satisfies_properties(cell.local_config)) {
                        neighbor_errors_minsize5++;
                    }
                }
            }
        }
    }
    
    bool valid = (neighbor_errors_minsize1 == 0 && neighbor_errors_minsize5 == 0);
    if (!valid) {
        std::cout << "FAILED: Boundary cell neighbors" << std::endl;
        if (neighbor_errors_minsize1 > 0) {
            std::cout << "  Minsize 1 clusters: " << neighbor_errors_minsize1 
                      << " cells with no neighbors from different clusters out of " << total_boundary_minsize1 
                      << " boundary cells" << std::endl;
        }
        if (neighbor_errors_minsize5 > 0) {
            std::cout << "  Minsize 5 clusters: " << neighbor_errors_minsize5 
                      << " cells with no neighbors from different clusters out of " << total_boundary_minsize5 
                      << " boundary cells" << std::endl;
        }
    }
    return valid;
}

bool test_cluster_sizes(Grid& grid) {
    int species_errors = 0;
    
    for (int i = 0; i < grid.num_clusters; i++) {
        // Count atoms of each species in boundary and interior sets for this cluster
        std::vector<int> boundary_species_counts(2, 0); // [species 0, species 1]
        std::vector<int> interior_species_counts(2, 0); // [species 0, species 1]
        
        for (int j = 0; j < grid.Size * grid.Size; j++) {
            if (grid.array[j].cell_cluster_ID == i) {
                int species = grid.array[j].cell_species;
                if (species < 0 || species > 1) {
                    std::cout << "Invalid species " << species << " found at atom " << j << std::endl;
                    continue;
                }
                
                if (grid.array[j].cell_type == CellType::BOUNDARY) {
                    boundary_species_counts[species]++;
                } else if (grid.array[j].cell_type == CellType::INTERIOR) {
                    interior_species_counts[species]++;
                }
            }
        }
        
        // Check sizes for each species in boundary sets
        for (int species = 0; species < 2; species++) {
            int actual_boundary_size = grid.clusters[i].boundary_cells[species].size();
            if (boundary_species_counts[species] != actual_boundary_size) {
                species_errors++;
                std::cout << "FAILED: Cluster " << i << ", Species " << species 
                          << " boundary size mismatch: counted " << boundary_species_counts[species]
                          << " but set has " << actual_boundary_size << std::endl;
            }
            
            int actual_interior_size = grid.clusters[i].interior_cells[species].size();
            if (interior_species_counts[species] != actual_interior_size) {
                species_errors++;
                std::cout << "FAILED: Cluster " << i << ", Species " << species 
                          << " interior size mismatch: counted " << interior_species_counts[species]
                          << " but set has " << actual_interior_size << std::endl;
            }
        }
    }
    
    bool valid = (species_errors == 0);
    if (!valid) {
        std::cout << "FAILED: " << species_errors << " cluster size mismatches found" << std::endl;
    }
    
    return valid;
}

bool test_cluster_membership(Grid& grid) {
    int unassigned_cells = 0;
    int duplicate_assignments = 0;
    int total_cells_in_sets = 0;
    int total_cells = grid.Size * grid.Size;
    
    // Create a map to track which atoms have been assigned to sets
    std::vector<int> cell_assignment_count(total_cells, 0);
    
    // Count cells in each cluster's sets, separated by species
    for (int i = 0; i < grid.num_clusters; i++) {
        Cluster& cluster = grid.clusters[i];
        
        // Count boundary and interior cells for each species
        for (int species = 0; species < 2; species++) {
            for (const int& atom : cluster.boundary_atoms_vector[species]) {
                if (atom >= 0 && atom < total_cells) {
                    cell_assignment_count[atom]++;
                    total_cells_in_sets++;
                } else {
                    std::cout << "FAILED: Invalid atom index " << atom 
                              << " in boundary_cells of cluster " << i 
                              << ", species " << species << std::endl;
                }
            }
            
            for (const int& atom : cluster.interior_atoms_vector[species]) {
                if (atom >= 0 && atom < total_cells) {
                    cell_assignment_count[atom]++;
                    total_cells_in_sets++;
                } else {
                    std::cout << "FAILED: Invalid atom index " << atom 
                              << " in interior_cells of cluster " << i 
                              << ", species " << species << std::endl;
                }
            }
        }
    }
    
    // Check how many cells are unassigned or have multiple assignments
    for (int i = 0; i < total_cells; i++) {
        if (grid.array[i].cell_cluster_ID == -1) {
            unassigned_cells++;
        }
        
        if (cell_assignment_count[i] == 0) {
            // Grid says this cell is not in any cluster's sets
            if (grid.array[i].cell_cluster_ID != -1) {
                std::cout << "FAILED: Cell " << i << " has cluster ID " << grid.array[i].cell_cluster_ID 
                          << " but is not in any cluster's sets" << std::endl;
            }
        } else if (cell_assignment_count[i] > 1) {
            // Cell appears in multiple sets
            duplicate_assignments++;
            std::cout << "FAILED: Cell " << i << " is assigned to " 
                      << cell_assignment_count[i] << " different sets" << std::endl;
        }
        
        // Verify that cell's cluster ID and type match its presence in sets
        if (grid.array[i].cell_cluster_ID != -1) {
            int cluster_id = grid.array[i].cell_cluster_ID;
            int species = grid.array[i].cell_species;
            CellType type = grid.array[i].cell_type;
            
            if (type == CellType::BOUNDARY) {
                if (grid.clusters[cluster_id].boundary_cells[species].find(i) == 
                    grid.clusters[cluster_id].boundary_cells[species].end()) {
                    std::cout << "FAILED: Cell " << i << " has cluster ID " << cluster_id 
                              << " and type BOUNDARY but is not in that cluster's boundary set" << std::endl;
                }
            } else if (type == CellType::INTERIOR) {
                if (grid.clusters[cluster_id].interior_cells[species].find(i) == 
                    grid.clusters[cluster_id].interior_cells[species].end()) {
                    std::cout << "FAILED: Cell " << i << " has cluster ID " << cluster_id 
                              << " and type INTERIOR but is not in that cluster's interior set" << std::endl;
                }
            }
        }
    }
    
    bool valid = (unassigned_cells == 0 && duplicate_assignments == 0);
    if (!valid) {
        if (unassigned_cells > 0) {
            std::cout << "FAILED: Cluster membership - " << unassigned_cells 
                      << " cells are not assigned to any cluster" << std::endl;
        }
        if (duplicate_assignments > 0) {
            std::cout << "FAILED: Cluster membership - " << duplicate_assignments 
                      << " cells are assigned to multiple cluster sets" << std::endl;
        }
        if (total_cells_in_sets != total_cells) {
            std::cout << "FAILED: Total cells in sets (" << total_cells_in_sets 
                      << ") doesn't match grid size (" << total_cells << ")" << std::endl;
        }
    }
    
    return valid;
}

// Test that each atom has exactly 4 neighbors
bool test_neighbor_count(Grid& grid) {
    int count_errors = 0;
    
    for (int i = 0; i < grid.Size * grid.Size; i++) {
        Cell& cell = grid.array[i];
        if (cell.neighbor_indexes.size() != 4) {
            count_errors++;
        }
    }
    
    bool valid = (count_errors == 0);
    if (!valid) {
        std::cout << "FAILED: Neighbor count - " << count_errors 
                  << " cells have an incorrect number of neighbors" << std::endl;
    }
    return valid;
}

// Test that each neighbor's species matches the corresponding local_config value
bool test_local_config_consistency(Grid& grid) {
    int size_errors = 0;
    int center_errors = 0;
    int neighbor_errors = 0;
    
    for (int i = 0; i < grid.Size * grid.Size; i++) {
        Cell& cell = grid.array[i];
        
        // Check that local_config has 5 elements: center + 4 neighbors
        if (cell.local_config.size() != 5) {
            size_errors++;
            continue;
        }
        
        // Check that center species matches
        if (cell.local_config[0] != cell.cell_species) {
            center_errors++;
        }
        
        // Check each neighbor's species matches the corresponding local_config
        bool cell_has_neighbor_error = false;
        for (int j = 0; j < 4; j++) {
            int neighbor_idx = cell.neighbor_indexes[j];
            if (grid.array[neighbor_idx].cell_species != cell.local_config[j + 1]) {
                cell_has_neighbor_error = true;
                break;
            }
        }
        if (cell_has_neighbor_error) {
            neighbor_errors++;
        }
    }
    
    bool valid = (size_errors == 0 && center_errors == 0 && neighbor_errors == 0);
    if (!valid) {
        if (size_errors > 0) {
            std::cout << "FAILED: Local config consistency - " << size_errors 
                      << " cells have incorrect local_config size" << std::endl;
        }
        if (center_errors > 0) {
            std::cout << "FAILED: Local config consistency - " << center_errors 
                      << " cells have mismatched center species in local_config" << std::endl;
        }
        if (neighbor_errors > 0) {
            std::cout << "FAILED: Local config consistency - " << neighbor_errors 
                      << " cells have neighbor species mismatches in local_config" << std::endl;
        }
    }
    return valid;
}

// Add this function declaration at the top with other test declarations
bool test_grid_total_consistency(Grid& grid);

// Add the new test function implementation
bool test_grid_total_consistency(Grid& grid) {
    //std::cout << "\n=== TESTING GRID TOTAL CONSISTENCY ===" << std::endl;
    
    // Initialize counters for manual count
    int actual_boundary_count[2] = {0, 0}; // For species 0 and 1
    int actual_interior_count[2] = {0, 0}; // For species 0 and 1
    
    // Count actual cells in the grid
    for (int i = 0; i < grid.Size * grid.Size; i++) {
        Cell& cell = grid.array[i];
        if (cell.cell_species < 0 || cell.cell_species > 1) {
            std::cout << "WARNING: Invalid species " << cell.cell_species << " at cell " << i << std::endl;
            continue;
        }
        
        if (cell.cell_type == CellType::BOUNDARY) {
            actual_boundary_count[cell.cell_species]++;
        } else if (cell.cell_type == CellType::INTERIOR) {
            actual_interior_count[cell.cell_species]++;
        }
    }
    
    // Get cluster-summed counts
    int cluster_boundary_count[2] = {0, 0};
    int cluster_interior_count[2] = {0, 0};
    
    for (int i = 0; i < grid.num_clusters; i++) {
        for (int species = 0; species < 2; species++) {
            cluster_boundary_count[species] += grid.clusters[i].boundary_cells[species].size();
            cluster_interior_count[species] += grid.clusters[i].interior_cells[species].size();
        }
    }
    
    // Calculate expected pair counts based on the actual counts
    int expected_bb_pairs = actual_boundary_count[0] * actual_boundary_count[1];
    int expected_bi_pairs = actual_boundary_count[0] * actual_interior_count[1] + 
                           actual_interior_count[0] * actual_boundary_count[1];
    int expected_ii_pairs = actual_interior_count[0] * actual_interior_count[1];
    
    // Check for errors
    bool has_errors = false;
    
    // Validate grid totals against actual counts
    for (int species = 0; species < 2; species++) {
        if (grid.total_boundary_atoms[species] != actual_boundary_count[species]) {
            std::cout << "ERROR: Grid total_boundary_atoms[" << species << "] = " 
                      << grid.total_boundary_atoms[species] 
                      << " doesn't match actual count " << actual_boundary_count[species] << std::endl;
            has_errors = true;
        }
        
        if (grid.total_interior_atoms[species] != actual_interior_count[species]) {
            std::cout << "ERROR: Grid total_interior_atoms[" << species << "] = " 
                      << grid.total_interior_atoms[species] 
                      << " doesn't match actual count " << actual_interior_count[species] << std::endl;
            has_errors = true;
        }
    }
    
    // Validate grid totals against cluster sums
    for (int species = 0; species < 2; species++) {
        if (grid.total_boundary_atoms[species] != cluster_boundary_count[species]) {
            std::cout << "ERROR: Grid total_boundary_atoms[" << species << "] = " 
                      << grid.total_boundary_atoms[species] 
                      << " doesn't match cluster sum " << cluster_boundary_count[species] << std::endl;
            has_errors = true;
        }
        
        if (grid.total_interior_atoms[species] != cluster_interior_count[species]) {
            std::cout << "ERROR: Grid total_interior_atoms[" << species << "] = " 
                      << grid.total_interior_atoms[species] 
                      << " doesn't match cluster sum " << cluster_interior_count[species] << std::endl;
            has_errors = true;
        }
    }
    
    // Only print summary if there are errors
    if (has_errors) {
        std::cout << "\nSummary of atom counts:" << std::endl;
        std::cout << "-----------------------------------------------------------------" << std::endl;
        std::cout << "Source        | Boundary Sp0 | Boundary Sp1 | Interior Sp0 | Interior Sp1" << std::endl;
        std::cout << "-----------------------------------------------------------------" << std::endl;
        std::cout << "Grid Totals   | " 
                  << std::setw(12) << grid.total_boundary_atoms[0] << " | " 
                  << std::setw(12) << grid.total_boundary_atoms[1] << " | "
                  << std::setw(12) << grid.total_interior_atoms[0] << " | "
                  << std::setw(12) << grid.total_interior_atoms[1] << std::endl;
        
        std::cout << "Actual Count  | " 
                  << std::setw(12) << actual_boundary_count[0] << " | " 
                  << std::setw(12) << actual_boundary_count[1] << " | "
                  << std::setw(12) << actual_interior_count[0] << " | "
                  << std::setw(12) << actual_interior_count[1] << std::endl;
        
        std::cout << "Cluster Sum   | " 
                  << std::setw(12) << cluster_boundary_count[0] << " | " 
                  << std::setw(12) << cluster_boundary_count[1] << " | "
                  << std::setw(12) << cluster_interior_count[0] << " | "
                  << std::setw(12) << cluster_interior_count[1] << std::endl;
        std::cout << "-----------------------------------------------------------------" << std::endl;
        
        // Output the expected pair counts
        std::cout << "\nExpected pair counts (based on actual atom counts):" << std::endl;
        std::cout << "  BB pairs: " << expected_bb_pairs << std::endl;
        std::cout << "  BI pairs: " << expected_bi_pairs << std::endl;
        std::cout << "  II pairs: " << expected_ii_pairs << std::endl;
    }
    
    return !has_errors;
}

// Function to validate all cluster properties
void validate_all_cluster_properties(Grid& grid) {
    bool all_valid = true;
    
    // Basic mechanics tests
    bool mechanics_valid = true;
    if (!test_neighbor_count(grid)) {
        std::cout << "FAILED: Basic mechanics - Incorrect neighbor count" << std::endl;
        mechanics_valid = false;
        all_valid = false;
    }
    if (!test_local_config_consistency(grid)) {
        std::cout << "FAILED: Basic mechanics - Local config inconsistency" << std::endl;
        mechanics_valid = false;
        all_valid = false;
    }
    
    if (!mechanics_valid) {
        std::cout << "Basic mechanics tests failed. Skipping cluster validation tests." << std::endl;
        return;
    }
    
    // New grid total consistency test
    if (!test_grid_total_consistency(grid)) {
        std::cout << "FAILED: Grid totals inconsistency detected" << std::endl;
        all_valid = false;
    }
    
    // Cluster validation tests
    if (!test_cluster_set_consistency(grid)) {
        all_valid = false;
    }
    if (!test_interior_cell_properties(grid)) {
        all_valid = false;
    }
    if (!test_boundary_cell_properties(grid)) {
        all_valid = false;
    }
    if (!test_interior_cell_neighbors(grid)) {
        all_valid = false;
    }
    if (!test_boundary_cell_neighbors(grid)) {
        all_valid = false;
    }
    if (!test_cluster_sizes(grid)) {
        all_valid = false;
    }
    if (!test_cluster_membership(grid)) {
        all_valid = false;
    }
    
    if (!all_valid) {
        std::cout << "Some validation tests failed" << std::endl;
    }
}

// Add this function to display cluster composition details
void display_cluster_composition(const Grid& grid) {
    std::cout << "\n=== CLUSTER COMPOSITION DETAILS ===" << std::endl;
    
    for (int i = 0; i < grid.num_clusters; i++) {
        std::cout << "Cluster " << i << " (min_size=" << grid.clusters[i].property.min_size << "):" << std::endl;
        
        // Display boundary cells counts for each species
        std::cout << "  Boundary cells: " << std::endl;
        int total_boundary = 0;
        for (int species = 0; species < 2; species++) {
            int boundary_count = grid.clusters[i].boundary_cells[species].size();
            total_boundary += boundary_count;
            std::cout << "    Species " << species << ": " << boundary_count << " cells" << std::endl;
        }
        std::cout << "    Total: " << total_boundary << " cells" << std::endl;
        
        // Display interior cells counts for each species
        std::cout << "  Interior cells: " << std::endl;
        int total_interior = 0;
        for (int species = 0; species < 2; species++) {
            int interior_count = grid.clusters[i].interior_cells[species].size();
            total_interior += interior_count;
            std::cout << "    Species " << species << ": " << interior_count << " cells" << std::endl;
        }
        std::cout << "    Total: " << total_interior << " cells" << std::endl;
        
        std::cout << "  Cluster total: " << (total_boundary + total_interior) << " cells" << std::endl;
    }
    
    // Display global totals
    int total_boundary_global = 0;
    int total_interior_global = 0;
    
    for (int species = 0; species < 2; species++) {
        total_boundary_global += grid.total_boundary_atoms[species];
        total_interior_global += grid.total_interior_atoms[species];
    }
    
    std::cout << "\nGlobal totals:" << std::endl;
    std::cout << "  Total boundary cells: " << total_boundary_global << std::endl;
    std::cout << "  Total interior cells: " << total_interior_global << std::endl;
    std::cout << "  Total cells: " << grid.Size * grid.Size << std::endl;
    std::cout << "=================================" << std::endl;
}

void validate_without_totals(Grid& grid) {
    bool all_valid = true;
    
    // Basic mechanics tests
    bool mechanics_valid = true;
    if (!test_neighbor_count(grid)) {
        std::cout << "FAILED: Basic mechanics - Incorrect neighbor count" << std::endl;
        mechanics_valid = false;
        all_valid = false;
    }
    if (!test_local_config_consistency(grid)) {
        std::cout << "FAILED: Basic mechanics - Local config inconsistency" << std::endl;
        mechanics_valid = false;
        all_valid = false;
    }
    
    if (!mechanics_valid) {
        std::cout << "Basic mechanics tests failed. Skipping cluster validation tests." << std::endl;
        return;
    }
    
    // Cluster validation tests
    if (!test_cluster_set_consistency(grid)) {
        all_valid = false;
    }
    if (!test_interior_cell_properties(grid)) {
        all_valid = false;
    }
    if (!test_boundary_cell_properties(grid)) {
        all_valid = false;
    }
    if (!test_interior_cell_neighbors(grid)) {
        all_valid = false;
    }
    if (!test_boundary_cell_neighbors(grid)) {
        all_valid = false;
    }
    if (!test_cluster_sizes(grid)) {
        all_valid = false;
    }
    if (!test_cluster_membership(grid)) {
        all_valid = false;
    }
    
    if (!all_valid) {
        std::cout << "Some validation tests failed" << std::endl;
    }
}

std::string run_and_save(std::string selection_string, int iterations, int cluster_type, bool overload, float same_cluster_prob, float expand, float tolerance_override) {
    // 1. Create a simulation object
    Sim sim;
    // 2. Initialize grid with cluster conditions and assign cells to clusters
    sim.conditions = SimConditions(SIM_CONDITIONS.at(selection_string));
    // Only override if overload is true
    if (overload) {
        sim.conditions.prob_same_cluster_base = same_cluster_prob;
        sim.conditions.prob_bb_expansion = expand;
        sim.conditions.prob_bi_expansion = expand;
        sim.conditions.prob_ii_expansion = expand;
    }
    // --- Tolerance override logic ---
    if (tolerance_override >= 0.0f) {
        sim.conditions.tolerance = tolerance_override;
    }
    // --- End tolerance override logic ---
    //std::cout << "[DEBUG] Initializing simulation parameters" << std::endl;
    sim.cluster_conditions(cluster_type);
    // --- Debug: print cluster compositions and weights ---
    std::cout << "[INFO] Cluster summary:" << std::endl;
    for (int ci = 0; ci < sim.grid.num_clusters; ++ci) {
        const Cluster &cl = sim.grid.clusters[ci];
        std::cout << "  Cluster " << ci
                  << " | min_size=" << cl.property.min_size
                  << " | weight="   << cl.cluster_weight << std::endl;
    }
    sim.initialize(sim.conditions);
    sim.calculate_initial_selection_probs_weighted(false);
    sim.system_energy = sim.compute_total_energy();
    //std::cout << "[DEBUG] Validating cluster properties" << std::endl;
    validate_all_cluster_properties(sim.grid);
    // Only run boundary cell property test at the start
    test_boundary_cell_properties(sim.grid);
    // 4. Run simulation
    //std::cout << "[DEBUG] Starting simulation loop" << std::endl;
    float total_delta_e = 0.0f;
    int total_accepted = 0;
    for (int i = 0; i < iterations; i++) {
        // Validation every 1000 iterations is now fully disabled
        // if (i % 1000 == 0) {
        //     test_grid_total_consistency(sim.grid);
        //     test_cluster_set_consistency(sim.grid);
        //     test_interior_cell_properties(sim.grid);
        //     // test_boundary_cell_properties(sim.grid); // Mute during loop
        //     test_interior_cell_neighbors(sim.grid);
        //     test_boundary_cell_neighbors(sim.grid);
        //     test_cluster_sizes(sim.grid);
        //     test_cluster_membership(sim.grid);
        // }
        sim.calculate_initial_selection_probs_weighted(false);
        Swap initial_swap = sim.select_initial_atoms_weighted(false);
        if (initial_swap.primary_atom == -1 || initial_swap.secondary_atom == -1) {
            std::cout << "  Primary atom: " << initial_swap.primary_atom << std::endl;
            std::cout << "  Secondary atom: " << initial_swap.secondary_atom << std::endl;
            std::cout << "  Swap type: " << (initial_swap.this_swap_type == SwapType::BB ? "BB" : initial_swap.this_swap_type == SwapType::BI ? "BI" : "II") << std::endl;
            std::cout << "  Same cluster: " << (initial_swap.same_cluster ? "YES" : "NO") << std::endl;
            std::cout << "[DEBUG] Exiting run_and_save due to invalid swap" << std::endl;
            test_grid_total_consistency(sim.grid);
            test_cluster_set_consistency(sim.grid);
            test_interior_cell_properties(sim.grid);
            // test_boundary_cell_properties(sim.grid); // Mute during error
            test_interior_cell_neighbors(sim.grid);
            test_boundary_cell_neighbors(sim.grid);
            test_cluster_sizes(sim.grid);
            test_cluster_membership(sim.grid);
            return "";
        }
        if (initial_swap.primary_atom < 0 || initial_swap.primary_atom >= sim.grid.Size * sim.grid.Size ||
            initial_swap.secondary_atom < 0 || initial_swap.secondary_atom >= sim.grid.Size * sim.grid.Size) {
            std::cout << "[DEBUG][ERROR] Atom indices out of bounds at iteration " << i << std::endl;
            std::cout << "  Primary atom: " << initial_swap.primary_atom << " (valid range: 0-" << (sim.grid.Size * sim.grid.Size - 1) << ")" << std::endl;
            std::cout << "  Secondary atom: " << initial_swap.secondary_atom << " (valid range: 0-" << (sim.grid.Size * sim.grid.Size - 1) << ")" << std::endl;
            std::cout << "[DEBUG] Exiting run_and_save due to out-of-bounds atom indices" << std::endl;
            test_grid_total_consistency(sim.grid);
            test_cluster_set_consistency(sim.grid);
            test_interior_cell_properties(sim.grid);
            // test_boundary_cell_properties(sim.grid); // Mute during error
            test_interior_cell_neighbors(sim.grid);
            test_boundary_cell_neighbors(sim.grid);
            test_cluster_sizes(sim.grid);
            test_cluster_membership(sim.grid);
            return "";
        }
        CellType primary_actual_type = sim.grid.array[initial_swap.primary_atom].cell_type;
        CellType secondary_actual_type = sim.grid.array[initial_swap.secondary_atom].cell_type;
        if (primary_actual_type != initial_swap.primary_atom_type || 
            secondary_actual_type != initial_swap.secondary_atom_type) {
            std::cout << "[DEBUG][ERROR] Atom type mismatch at iteration " << i << std::endl;
            std::cout << "  Primary atom " << initial_swap.primary_atom << ": expected " 
                      << (initial_swap.primary_atom_type == CellType::BOUNDARY ? "BOUNDARY" : "INTERIOR")
                      << ", actual " << (primary_actual_type == CellType::BOUNDARY ? "BOUNDARY" : "INTERIOR") << std::endl;
            std::cout << "  Secondary atom " << initial_swap.secondary_atom << ": expected " 
                      << (initial_swap.secondary_atom_type == CellType::BOUNDARY ? "BOUNDARY" : "INTERIOR")
                      << ", actual " << (secondary_actual_type == CellType::BOUNDARY ? "BOUNDARY" : "INTERIOR") << std::endl;
            std::cout << "[DEBUG] Exiting run_and_save due to atom type mismatch" << std::endl;
            test_grid_total_consistency(sim.grid);  
            test_cluster_set_consistency(sim.grid);
            test_interior_cell_properties(sim.grid);
            // test_boundary_cell_properties(sim.grid); // Mute during error
            test_interior_cell_neighbors(sim.grid);
            test_boundary_cell_neighbors(sim.grid);
            test_cluster_sizes(sim.grid);
            test_cluster_membership(sim.grid);
            return "";
        }
        sim.variables.backup();
        Move move(sim.conditions.J, sim.conditions.Beta, sim.conditions.tolerance);
        move.deterministic_build_from_swap(initial_swap, sim.grid, false);
        sim.calculate_initial_selection_probs_weighted(false);
        Swap reverse_swap = sim.calculate_reverse_swap(move);
        float acceptance_prob = move.get_final_acceptance_deterministic(sim.grid, reverse_swap, false);
        if (std::isnan(acceptance_prob) || acceptance_prob < 0.0f || acceptance_prob > 1.0f) {
            std::cout << "[DEBUG][ERROR] Invalid acceptance probability at iteration " << i << std::endl;
            std::cout << "  Acceptance probability: " << acceptance_prob << std::endl;
            std::cout << "  Forward probability: " << move.final_forward_probability << std::endl;
            std::cout << "  Reverse probability: " << move.final_reverse_probability << std::endl;
            std::cout << "  Delta E: " << move.DeltaE << std::endl;
            std::cout << "[DEBUG] Exiting run_and_save due to invalid acceptance probability" << std::endl;
            test_grid_total_consistency(sim.grid);
            test_cluster_set_consistency(sim.grid);
            test_interior_cell_properties(sim.grid);
            // test_boundary_cell_properties(sim.grid); // Mute during error
            test_interior_cell_neighbors(sim.grid);
            test_boundary_cell_neighbors(sim.grid);
            test_cluster_sizes(sim.grid);
            test_cluster_membership(sim.grid);
            return "";
        }
        std::uniform_real_distribution<float> dist(0.0f, 1.0f);
        float rand_num = dist(sim.gen);
        bool accepted = false;
        if (rand_num < acceptance_prob) {
            //std::cout << "[DEBUG] Move accepted at iteration " << i << std::endl;
            sim.accept_move(move);
            accepted = true;
            total_accepted++;
        } else {
            //std::cout << "[DEBUG] Move rejected at iteration " << i << std::endl;
            sim.undo_move(move);
        }
        if (iterations >= 5 && (i + 1) % (iterations / 5) == 0) {
            int percent_complete = ((i + 1) * 20) / (iterations / 5);
            float current_acceptance_rate = (total_accepted * 100.0f) / (i + 1);
            std::cout << "[PROGRESS] " << percent_complete << "% complete | Acceptance rate: "
                      << current_acceptance_rate << "% | Energy: " << sim.system_energy << std::endl;
        }
    }
    //std::cout << "[DEBUG] Simulation loop complete" << std::endl;
    std::cout << "\n====== SIMULATION COMPLETE ======" << std::endl;
    std::cout << "Final system energy: " << sim.system_energy << std::endl;
    std::cout << "Final validation:" << std::endl;
    validate_all_cluster_properties(sim.grid);
    // Only run boundary cell property test at the end
    test_boundary_cell_properties(sim.grid);
    float final_acceptance_rate = (total_accepted * 100.0f / iterations);
    std::cout << "\nFinal Statistics:" << std::endl;
    std::cout << "  Total iterations: " << iterations << std::endl;
    std::cout << "  Accepted moves: " << total_accepted << std::endl;
    std::cout << "  Acceptance rate: " << final_acceptance_rate << "%" << std::endl;
    std::string cluster_type_str = std::to_string(cluster_type);
    std::ostringstream tol_signifier;
    if (tolerance_override >= 0.0f) {
        tol_signifier << std::fixed << std::setprecision(2) << "_tol" << tolerance_override;
    }
    std::ostringstream scp_signifier;
    scp_signifier << std::fixed << std::setprecision(2) << "_scp" << same_cluster_prob;
    std::string results_file = "outputs/" + selection_string + "_" + cluster_type_str + tol_signifier.str() + scp_signifier.str() + ".csv";
    std::string config_file = "outputs/" + selection_string + "_" + cluster_type_str + tol_signifier.str() + scp_signifier.str() + "_GridConfig.csv";
    std::cout << "[DEBUG] Saving results to CSV: " << results_file << std::endl;
    sim.save_results_to_csv(results_file);
    sim.save_grid_configuration_to_csv(config_file);
    // Only return the results file (without .csv extension)
    size_t dot_pos = results_file.rfind(".csv");
    if (dot_pos != std::string::npos) {
        return results_file.substr(0, dot_pos);
    } else {
        return results_file;
    }
}

// Add a new function for delayed rejection iteration
std::string run_and_save_delayed_rejection(std::string selection_string, int iterations, int cluster_type, int max_num_swaps, bool overload, float same_cluster_prob, float expand, float tolerance_override) {
    std::cout << "[DEBUG] Starting run_and_save_delayed_rejection with " << iterations << " iterations, cluster_type=" << cluster_type << ", max_num_swaps=" << max_num_swaps << std::endl;
    Sim sim;
    sim.conditions = SimConditions(SIM_CONDITIONS.at(selection_string));
    if (overload) {
        sim.conditions.prob_same_cluster_base = same_cluster_prob;
        sim.conditions.prob_bb_expansion = expand;
        sim.conditions.prob_bi_expansion = expand;
        sim.conditions.prob_ii_expansion = expand;
    }
    if (tolerance_override >= 0.0f) {
        sim.conditions.tolerance = tolerance_override;
    }
   // std::cout << "[DEBUG] Setting up cluster conditions..." << std::endl;
   sim.conditions.cluster_bias=0.25f;
    sim.cluster_conditions(cluster_type);
    // --- Debug: print cluster compositions and weights ---
    //std::cout << "[INFO] Cluster summary:" << std::endl;
    for (int ci = 0; ci < sim.grid.num_clusters; ++ci) {
      //  const Cluster &cl = sim.grid.clusters[ci];
      //  std::cout << "  Cluster " << ci
      //            << " | min_size=" << cl.property.min_size
      //            << " | weight="   << cl.cluster_weight << std::endl;
    }
    std::cout << "[DEBUG] Initializing simulation..." << std::endl;
    sim.initialize(sim.conditions);
    sim.calculate_initial_selection_probs_weighted(false);
    sim.system_energy = sim.compute_total_energy();
    std::cout << "[DEBUG] Validating cluster properties..." << std::endl;
    validate_all_cluster_properties(sim.grid);
    test_boundary_cell_properties(sim.grid);
    float total_delta_e = 0.0f;
    int total_accepted = 0;
    for (int i = 0; i < iterations; i++) {
        bool debug_this_iter = (iterations >= 5 && (i + 1) % (iterations / 5) == 0); // debug only at 20% progress points
        debug_this_iter=false;
        bool accepted = sim.delayed_rejection_iteration(max_num_swaps, debug_this_iter);
        if (accepted) {
            total_accepted++;
        }
        if (iterations >= 5 && (i + 1) % (iterations / 5) == 0) {
            int percent_complete = ((i + 1) * 20) / (iterations / 5);
            float current_acceptance_rate = (total_accepted * 100.0f) / (i + 1);
            std::cout << "[PROGRESS] " << percent_complete << "% complete | Acceptance rate: "
                      << current_acceptance_rate << "% | Energy: " << sim.system_energy << std::endl;
        }
    }
    std::cout << "\n====== DELAYED REJECTION SIMULATION COMPLETE ======" << std::endl;
    std::cout << "Final system energy: " << sim.system_energy << std::endl;
    std::cout << "Final validation:" << std::endl;
    validate_all_cluster_properties(sim.grid);
    test_boundary_cell_properties(sim.grid);
    float final_acceptance_rate = (total_accepted * 100.0f / iterations);
    std::cout << "\nFinal Statistics:" << std::endl;
    std::cout << "  Total iterations: " << iterations << std::endl;
    std::cout << "  Accepted moves: " << total_accepted << std::endl;
    std::cout << "  Acceptance rate: " << final_acceptance_rate << "%" << std::endl;
    std::string cluster_type_str = std::to_string(cluster_type);
    std::ostringstream tol_signifier;
    if (tolerance_override >= 0.0f) {
        tol_signifier << std::fixed << std::setprecision(2) << "_tol" << tolerance_override;
    }
    std::ostringstream scp_signifier;
    scp_signifier << std::fixed << std::setprecision(2) << "_scp" << same_cluster_prob;
    std::string results_file = "outputs/" + selection_string  +"_"+ cluster_type_str + tol_signifier.str() + scp_signifier.str() + ".csv";
    std::string config_file = "outputs/" + selection_string  +"_"+ cluster_type_str + tol_signifier.str() + scp_signifier.str() + "_GridConfig.csv";
    std::cout << "[DEBUG] Saving delayed rejection results to CSV: " << results_file << std::endl;
    sim.save_results_to_csv(results_file);
    sim.save_grid_configuration_to_csv(config_file);
    size_t dot_pos = results_file.rfind(".csv");
    if (dot_pos != std::string::npos) {
        return results_file.substr(0, dot_pos);
    } else {
        return results_file;
    }
}

// ----------------- New helper implementation -----------------
std::vector<std::vector<int>> build_cluster_tree_(const std::string& selection_string, 
                                    int min_clusters, int number_of_clusters,
                                    int min_properties_per_cluster, int max_properties_per_cluster,
                                    float cluster_bias, float same_cluster_prob, 
                                    int depth, int final_test_runs, bool use_minsize5_boundary,
                                    bool use_rotation_only) {
    // Create a "same starting point" simulation with cluster condition 2
    Sim same_starting_point;
    same_starting_point.conditions = SimConditions(SIM_CONDITIONS.at(selection_string));
    same_starting_point.cluster_conditions(2);
    same_starting_point.initialize(same_starting_point.conditions);
    validate_all_cluster_properties(same_starting_point.grid);
    
    int properties_checked = 0;
    int branches_pruned = 0;
    float lowest = same_starting_point.system_energy;
    std::vector<std::vector<int>> lowest_cluster_matrix;
    bool lowest_used_catch_all = true;
    
    // Select pattern library based on use_rotation_only parameter
    int num_properties = use_rotation_only ? 
        static_cast<int>(GLOBAL_PATTERN_VECTORS_ROT.size()) : 
        static_cast<int>(GLOBAL_PATTERN_VECTORS.size());
    
    // Track which properties have been used to avoid repetition
    std::unordered_set<int> used_properties;
    
    // Helper function to check if a configuration is valid
    // A configuration is valid if it uses all available properties
    auto is_valid_configuration = [&](const std::vector<std::vector<int>>& config) -> bool {
        std::unordered_set<int> all_used;
        for (const auto& cluster : config) {
            for (int prop : cluster) {
                all_used.insert(prop);
            }
        }
        
        // Check if all properties are used
        for (int prop = 0; prop < num_properties; prop++) {
            if (all_used.find(prop) == all_used.end()) {
                return false; // Invalid - missing property
            }
        }
        return true; // Valid - all properties used
    };
    
    // Thread-safe tested configurations map with mutex
    std::mutex tested_configurations_mutex;
    const float PLACEHOLDER_VALUE = -999999.0f; // Placeholder while testing
    
    // Helper function to safely check and reserve a configuration for testing
    auto try_reserve_configuration = [&](const std::vector<std::vector<int>>& pruned_config) -> bool {
        std::lock_guard<std::mutex> lock(tested_configurations_mutex);
        auto it = tested_configurations.find(pruned_config);
        if (it != tested_configurations.end()) {
            return false; // Already tested or being tested
        }
        // Reserve this configuration with placeholder
        tested_configurations[pruned_config] = PLACEHOLDER_VALUE;
        return true; // Successfully reserved
    };
    
    // Helper function to safely update configuration result
    auto update_configuration_result = [&](const std::vector<std::vector<int>>& pruned_config, float energy) -> void {
        std::lock_guard<std::mutex> lock(tested_configurations_mutex);
        tested_configurations[pruned_config] = energy;
    };
    
    // Helper function to safely get configuration result
    auto get_configuration_result = [&](const std::vector<std::vector<int>>& pruned_config) -> std::pair<bool, float> {
        std::lock_guard<std::mutex> lock(tested_configurations_mutex);
        auto it = tested_configurations.find(pruned_config);
        if (it != tested_configurations.end() && it->second != PLACEHOLDER_VALUE) {
            return {true, it->second};
        }
        return {false, 0.0f};
    };
    
    // Thread-safe update of best configuration
    std::mutex best_config_mutex;
    auto update_best_if_better = [&](const std::vector<std::vector<int>>& config, float energy) -> void {
        std::lock_guard<std::mutex> lock(best_config_mutex);
        if (energy < lowest) {
            lowest = energy;
            lowest_cluster_matrix = config;
            lowest_used_catch_all = false; // No catch-all anymore
        }
    };

    // Updated test_configuration function with thread safety
    auto test_configuration = [&](const std::vector<std::vector<int>>& config, bool unused_param) -> void {
        // Use the configuration directly (no catch-all manipulation needed)
        std::vector<std::vector<int>> test_config = config;
        
        // Prune the configuration to get its standard form
        std::vector<std::vector<int>> pruned_config = prune(test_config);
        
        // Check if already tested and get result if available
        std::pair<bool, float> test_result = get_configuration_result(pruned_config);
        if (test_result.first) {
            update_best_if_better(config, test_result.second);
            return;
        }
        
        // Try to reserve this configuration for testing
        if (!try_reserve_configuration(pruned_config)) {
            return; // Another thread is testing this or it's already done
        }
        
        #pragma omp atomic
        properties_checked++;
        
        Sim sim;
        sim.conditions = SimConditions(SIM_CONDITIONS.at(selection_string));
        sim.conditions.minsize5_boundary_setting = use_minsize5_boundary;
        sim.conditions.prob_same_cluster_base = same_cluster_prob;
        sim.conditions.cluster_bias = cluster_bias;
        
        sim.cluster_conditions(2);
        sim.copy_grid(same_starting_point.grid);
        
        // Use the original test_config (not pruned) for simulation
        // Select appropriate cluster function and tier size based on use_rotation_only parameter
        if (use_rotation_only) {
            sim.cluster_conditions_from_indices_matrix_uninverted(test_config, 2);
        } else {
            sim.cluster_conditions_from_indices_matrix(test_config, 1);
        }
        sim.initialize(sim.conditions);
        
        float final_energy;
        
        // Check if variables.Norm is 0 - if so, skip iterations and store a high energy value
        if (sim.variables.Norm == 0.0) {
            final_energy = 1e6f;
        } else {
            // Run the specified number of iterations
            for(int iter = 0; iter < depth; iter++) {
                sim.iterate_improved();
            }
            final_energy = sim.system_energy;
        }
        
        // Update the result in the map
        update_configuration_result(pruned_config, final_energy);
        
        // Update best configuration if this is better
        update_best_if_better(config, final_energy);
    };
    
    // Debug flag for permutation generation
    bool debug_permutations = false; // Set to true to enable debug output
    
    // Recursive function to generate all permutations with unused properties tracking
    std::function<void(std::vector<std::vector<int>>&, std::vector<int>&)> generate_permutations = 
        [&](std::vector<std::vector<int>>& current_config, std::vector<int>& unused_properties) {
            if (debug_permutations && current_config.size() <= 2) {
                std::cout << "[DEBUG] Config size: " << current_config.size() << ", Unused properties: [";
                for (size_t i = 0; i < unused_properties.size(); i++) {
                    std::cout << unused_properties[i];
                    if (i < unused_properties.size() - 1) std::cout << ",";
                }
                std::cout << "]" << std::endl;
            }
            // Test current configuration if it meets minimum cluster requirements
            if (static_cast<int>(current_config.size()) >= min_clusters) {
                // Check if all clusters meet minimum properties requirement
                bool all_clusters_valid = true;
                for (const auto& cluster : current_config) {
                    if (static_cast<int>(cluster.size()) < min_properties_per_cluster) {
                        all_clusters_valid = false;
                        break;
                    }
                }
                
                if (all_clusters_valid) {
                    bool is_valid = is_valid_configuration(current_config);
                    
                    if (is_valid) {
                        // Configuration is valid, test it
                        test_configuration(current_config, false);
                    }
                    // Note: Invalid configurations are skipped
                }
            }
            
            // Check if we've reached maximum capacity (stop mutating)
            bool at_max_capacity = true;
            if (static_cast<int>(current_config.size()) < number_of_clusters) {
                at_max_capacity = false; // Can still add more vectors
            } else {
                // Check if any vector has less than max properties
                for (const auto& cluster : current_config) {
                    if (static_cast<int>(cluster.size()) < max_properties_per_cluster) {
                        at_max_capacity = false; // Can still add more properties
                        break;
                    }
                }
            }
            
            // Also check if we have unused properties available
            if (unused_properties.empty()) {
                at_max_capacity = true; // No properties left to use
                #pragma omp atomic
                branches_pruned++;
            }
            
            if (at_max_capacity) {
                return; // Stop mutating - at maximum capacity or no unused properties
            }
            
            // Check if this configuration (or any permutation of it) has already been processed
            std::vector<std::vector<int>> pruned_config = prune(current_config);
            
            // If this configuration has been processed, stop mutating this branch
            {
                std::lock_guard<std::mutex> lock(tested_configurations_mutex);
                if (tested_configurations.find(pruned_config) != tested_configurations.end()) {
                    return; // Stop mutating - already processed
                }
            }

            
            // Mutation 1: Add every possible cluster property to every existing vector 
            // that has less than max properties
            bool mutation1_applied = false;
            for (size_t cluster_idx = 0; cluster_idx < current_config.size(); cluster_idx++) {
                if (static_cast<int>(current_config[cluster_idx].size()) < max_properties_per_cluster) {
                    for (size_t prop_idx = 0; prop_idx < unused_properties.size(); prop_idx++) {
                        int prop = unused_properties[prop_idx];
                        // Create new configuration with this property added to this cluster
                        std::vector<std::vector<int>> new_config = current_config;
                        new_config[cluster_idx].push_back(prop);
                        
                        // Create new unused_properties vector without this property
                        std::vector<int> new_unused_properties = unused_properties;
                        new_unused_properties.erase(new_unused_properties.begin() + prop_idx);
                        
                        generate_permutations(new_config, new_unused_properties);
                        mutation1_applied = true;
                    }
                }
            }
            
            // Mutation 2: Add a new vector with each possible cluster property
            // (only if we haven't reached max number of vectors)
            bool mutation2_applied = false;
            if (static_cast<int>(current_config.size()) < number_of_clusters) {
                for (size_t prop_idx = 0; prop_idx < unused_properties.size(); prop_idx++) {
                    int prop = unused_properties[prop_idx];
                    // Create new configuration with a new single-property cluster
                    std::vector<std::vector<int>> new_config = current_config;
                    new_config.push_back({prop});
                    
                    // Create new unused_properties vector without this property
                    std::vector<int> new_unused_properties = unused_properties;
                    new_unused_properties.erase(new_unused_properties.begin() + prop_idx);
                    
                    generate_permutations(new_config, new_unused_properties);
                    mutation2_applied = true;
                }
            }
            
            // If neither mutation was applied, we're at a dead end - this prevents infinite recursion
            if (!mutation1_applied && !mutation2_applied) {
                return;
            }
        };
    
    // Generate all configurations in parallel
    std::vector<std::vector<std::vector<int>>> initial_configs;
    
    // First, generate all possible single-cluster configurations
    for (int cluster_size = min_properties_per_cluster; cluster_size <= max_properties_per_cluster; cluster_size++) {
        if (cluster_size > num_properties) break;
        
        // Generate all combinations of cluster_size properties
        std::function<void(std::vector<int>&, int, int)> generate_initial_combinations = 
            [&](std::vector<int>& current_combination, int start_idx, int remaining_props) {
                if (remaining_props == 0) {
                    initial_configs.push_back({current_combination});
                    return;
                }
                
                for (int i = start_idx; i <= num_properties - remaining_props; i++) {
                    current_combination.push_back(i);
                    generate_initial_combinations(current_combination, i + 1, remaining_props - 1);
                    current_combination.pop_back();
                }
            };
        
        std::vector<int> current_combination;
        generate_initial_combinations(current_combination, 0, cluster_size);
    }
    
    // Process all configurations in parallel
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < initial_configs.size(); i++) {
        std::vector<std::vector<int>> config = initial_configs[i];
        
        // Create initial unused properties vector excluding those already used in config
        std::unordered_set<int> used_properties;
        for (const auto& cluster : config) {
            for (int prop : cluster) {
                used_properties.insert(prop);
            }
        }
        
        std::vector<int> initial_unused_properties;
        // Add all regular properties (0 to num_properties-1)
        for (int prop = 0; prop < num_properties; prop++) {
            if (used_properties.find(prop) == used_properties.end()) {
                initial_unused_properties.push_back(prop);
            }
        }
        
        generate_permutations(config, initial_unused_properties);
    }
    
    std::cout << "[INFO] Total cluster configurations screened: " << properties_checked << std::endl;
    std::cout << "[INFO] Tested configurations: " << tested_configurations.size() << " unique configurations" << std::endl;
    std::cout << "[INFO] Configurations actually simulated: " << properties_checked << std::endl;
    std::cout << "[INFO] Branches pruned due to unused properties optimization: " << branches_pruned << std::endl;
    
    if(lowest_cluster_matrix.empty()) {
        std::cout << "[WARNING] No valid configurations found, skipping confirmation simulation." << std::endl;
        return std::vector<std::vector<int>>();
    }

    // Run confirmation simulation with the best configuration
    Sim best_sim;
    best_sim.conditions = SimConditions(SIM_CONDITIONS.at(selection_string));
    best_sim.conditions.minsize5_boundary_setting = use_minsize5_boundary;
    best_sim.conditions.prob_same_cluster_base = same_cluster_prob;
    best_sim.conditions.cluster_bias = cluster_bias;

    // Use the configuration directly (catch-all is already included if present)
    std::vector<std::vector<int>> final_config = lowest_cluster_matrix;

    // Select appropriate cluster function and tier size for confirmation simulation
    if (use_rotation_only) {
        best_sim.cluster_conditions_from_indices_matrix_uninverted(final_config, 2);
    } else {
        best_sim.cluster_conditions_from_indices_matrix(final_config, 1);
    }
    best_sim.initialize(best_sim.conditions);

    int confirm_iters = final_test_runs;
    int accepted_moves = 0;
    for (int t = 0; t < confirm_iters; ++t) {
        if (best_sim.iterate_improved()) ++accepted_moves;
    }

    // Calculate average energy over the run
    double avg_energy = 0.0;
    if(!best_sim.energy_history.empty()) {
        for(double e : best_sim.energy_history) avg_energy += e;
        avg_energy /= static_cast<double>(best_sim.energy_history.size());
    }

    std::cout << "[INFO] Confirmation sim completed: " << confirm_iters << " iterations, "
              << accepted_moves << " accepted. Final energy = "
              << best_sim.system_energy << " | Avg energy = " << avg_energy << std::endl;
    
    std::cout << "[INFO] Best cluster configuration:" << std::endl;
    
    // Show the cluster 0 (0-weight catch-all) and then all property clusters
    std::cout << "  Cluster 0: properties [0-weight catch-all]" << std::endl;
    
    for(size_t i = 0; i < lowest_cluster_matrix.size(); i++) {
        std::cout << "  Cluster " << (i + 1) << ": properties [";
        for(size_t j = 0; j < lowest_cluster_matrix[i].size(); j++) {
            std::cout << lowest_cluster_matrix[i][j];
            if(j < lowest_cluster_matrix[i].size() - 1) std::cout << ",";
        }
        std::cout << "]" << std::endl;
    }

    return lowest_cluster_matrix;
}

// Disable the original heavy main function for testing subset utilities
#if 0
int main() {
    initialize_energy();
    initialize_pattern_library();
    
    // Output cluster property sets from globals.cpp
    std::cout << "\n============ CLUSTER PROPERTY SETS OUTPUT ============" << std::endl;
    
    // Output GLOBAL_PATTERN_VECTORS (with inversions) - Index to Pattern mapping
    std::cout << "\n--- INDEX TO PATTERN MAPPING (with inversions) ---" << std::endl;
    std::cout << "Vector size: " << GLOBAL_PATTERN_VECTORS.size() << std::endl;
    for (size_t i = 0; i < GLOBAL_PATTERN_VECTORS.size(); ++i) {
        std::cout << "Index " << i << " -> Pattern [";
        for (size_t j = 0; j < GLOBAL_PATTERN_VECTORS[i].size(); ++j) {
            std::cout << GLOBAL_PATTERN_VECTORS[i][j];
            if (j < GLOBAL_PATTERN_VECTORS[i].size() - 1) std::cout << ",";
        }
        std::cout << "]" << std::endl;
    }
    
    // Initialize rotation-only pattern library
    initialize_pattern_library_rot();
    
    // Output GLOBAL_PATTERN_VECTORS_ROT (rotation-only) - Index to Pattern mapping
    std::cout << "\n--- INDEX TO PATTERN MAPPING (rotation-only) ---" << std::endl;
    std::cout << "Vector size: " << GLOBAL_PATTERN_VECTORS_ROT.size() << std::endl;
    for (size_t i = 0; i < GLOBAL_PATTERN_VECTORS_ROT.size(); ++i) {
        std::cout << "Index " << i << " -> Pattern [";
        for (size_t j = 0; j < GLOBAL_PATTERN_VECTORS_ROT[i].size(); ++j) {
            std::cout << GLOBAL_PATTERN_VECTORS_ROT[i][j];
            if (j < GLOBAL_PATTERN_VECTORS_ROT[i].size() - 1) std::cout << ",";
        }
        std::cout << "]" << std::endl;
    }
    
    std::cout << "\n============ END CLUSTER PROPERTY SETS OUTPUT ============" << std::endl;
    
    std::string BJ = "-0.44";
    
    // Run cluster tree optimization with standard patterns (includes inversions)
    // std::vector<std::vector<int>> optimal_cluster_config = build_cluster_tree_("BJ=-0.44_DR", 2, 2, 1, 1, 0.5f, 0.0f, 100, 10000, false, false);
    
    // Run cluster tree optimization with rotation-only patterns (no inversions)
    // std::vector<std::vector<int>> optimal_cluster_config_rot = build_cluster_tree_("BJ=-0.44_DR", 4, 4, 2, 2, 0.5f, 0.0f, 10, 10000, false, true);
    int iterations = 50;
    std::vector<std::string> all_filenames;
    std::cout << "\n============ RUNNING KAWASAKI SIMULATIONS ============" << std::endl;
    // Run Kawasaki simulation
    std::string Kawasaki_string = "BJ=" + BJ + "_Kawasaki";
    std::cout << "\nRunning Kawasaki simulation" << std::endl;
//  all_filenames.push_back(run_and_save(Kawasaki_string, iterations, 0, false, 0.0f, 0.0f, 0.0f));
    all_filenames.push_back(run_and_save_delayed_rejection(Kawasaki_string, iterations, 0, 1, true, 0.5f, 0.0f, 0.0f));
    std::cout << "\n============ RUNNING NEW MC SIMULATIONS ============" << std::endl;
    
    // Run with different NewMC variants
    std::vector<std::string> mc_variants = {
     "BJ=" + BJ + "_DR",
     "BJ=" + BJ + "_DR_Expand",
     "BJ=" + BJ + "_NewMC",
    //  "BJ=" + BJ + "_NewMC_Restricted", 
    //  "BJ=" + BJ + "_NewMC_Biased",
    //  "BJ=" + BJ + "_NewMC_Same",
    //   "BJ=" + BJ + "_NewMC_Expand",
    //   "BJ=" + BJ + "_NewMC_Restricted_Expand",
    //   "BJ=" + BJ + "_NewMC_Biased_Expand",
    //   "BJ=" + BJ + "_NewMC_Same_Expand"
    };
    // Run NewMC base case separately
    std::string base_variant = "BJ=" + BJ + "_NewMC";
    std::string base_variant_expand = "BJ=" + BJ + "_NewMC_Expand";
  //  std::cout << "\nRunning " << base_variant << " base case" << std::endl;
 //   std::vector<ClusterProperty> properties = build_cluster_tree_(BJ);

    //std::cout << "\nRunning Delayed Rejection simulation with max_num_swaps = ";
    int max_num_swaps = 1;
    // Delayed-Rejection variant (first entry in mc_variants) — run for each cluster_type and tolerance value
    const std::string dr_variant = mc_variants[0];
    const std::string dr_variant_expand = mc_variants[1];
    
    //for(int i = 0; i < 10; i++) {
      //  std::vector<int> properties = build_cluster_tree_("BJ=-0.44_NewMC", 0.00f, 0.1f, 6, 1000, 1, true);
    //}






    for (int cluster_type : {1}) {
        for (float tol = 0.50f; tol <= 1.0f; tol += 1.0f) {
            for (float scp = 0.50f; scp <= 1.00f; scp += 1.0f) {
            //    std::cout << "\nRunning " << dr_variant << " | cluster_type=" << cluster_type
            //             << " | tol=" << tol << " | scp=" << scp << std::endl;

               all_filenames.push_back(run_and_save_delayed_rejection(dr_variant, iterations, cluster_type, max_num_swaps, true, scp, 0.0f, tol));   
           //    all_filenames.push_back(run_and_save_delayed_rejection(dr_variant_expand, iterations, cluster_type, max_num_swaps, true, scp, 0.5f, tol));       
           //    all_filenames.push_back(run_and_save(base_variant, iterations, 2, true, scp, 0.0f, 0.2f));
           //    all_filenames.push_back(run_and_save(base_variant_expand, iterations, 2, true, scp, 0.9f, 0.2f));    
                                  
            }
        }
    }
    // Output all filenames in JSON-style array format
    std::cout << "\nAll results files (no .csv):" << std::endl;
    std::cout << "[";
    for (size_t i = 0; i < all_filenames.size(); ++i) {
        std::cout << "\"" << all_filenames[i] << "\"";
        if (i + 1 < all_filenames.size()) std::cout << ",";
    }
    std::cout << "]" << std::endl;



    return 0;
}
#endif

// Lightweight test driver for subset utilities
int main() {
    std::cout << "\n============ COMBINATORICS TESTING SUITE ============" << std::endl;
    std::cout << "Testing subset and permutation utilities from globals.hpp/cpp\n" << std::endl;
    
    // Test 1: Bit manipulation method with unique elements
    std::cout << "=== TEST 1: Bit Manipulation Method ===" << std::endl;
    std::vector<int> unique_nums = {1, 2, 3};
    auto subsets1 = getSubsetsBitManipulation(unique_nums);
    printSubsets(subsets1, "Subsets from Bit Manipulation (Input: {1, 2, 3})");
    
    // Test 2: Backtracking method with duplicate elements
    std::cout << "=== TEST 2: Backtracking Method (handles duplicates) ===" << std::endl;
    std::vector<int> duplicate_nums = {1, 2, 2};
    auto subsets2 = getSubsetsBacktracking(duplicate_nums);
    printSubsets(subsets2, "Subsets from Backtracking (Input: {1, 2, 2})");
    
    // Test 3: Permutations of subset collections
    std::cout << "=== TEST 3: Permutations of Subset Collections ===" << std::endl;
    std::vector<std::vector<int>> subset_set = {{1}, {2, 3}, {}};
    auto permutations = getAllPermutations(subset_set);
    printPermutations(permutations, "Permutations of Subsets (Input: {{1}, {2, 3}, {}})");
    
    // Test 4: Larger example with bit manipulation
    std::cout << "=== TEST 4: Larger Example (4 elements) ===" << std::endl;
    std::vector<int> four_elements = {0, 1, 2, 3};
    auto subsets4 = getSubsetsBitManipulation(four_elements);
    printSubsets(subsets4, "All Subsets of {0, 1, 2, 3} (16 total)");
    
    // Test 5: Small permutation example
    std::cout << "=== TEST 5: Small Permutation Example ===" << std::endl;
    std::vector<std::vector<int>> small_subsets = {{1}, {2}};
    auto small_perms = getAllPermutations(small_subsets);
    printPermutations(small_perms, "Permutations of {{1}, {2}}");
    
    // Test 6: Edge case - empty subset in collection
    std::cout << "=== TEST 6: Edge Case with Empty Subset ===" << std::endl;
    std::vector<int> single_element = {42};
    auto single_subsets = getSubsetsBitManipulation(single_element);
    printSubsets(single_subsets, "Subsets of {42} (includes empty set)");
    
    std::cout << "\n============ COMBINATORICS TESTING COMPLETE ============" << std::endl;
    std::cout << "All subset generation and permutation functions working correctly!" << std::endl;
    
    return 0;
}

