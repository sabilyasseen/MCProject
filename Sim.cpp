#include "Sim.hpp"
#include "globals.hpp"
#include <random>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <sys/stat.h>
#ifdef _WIN32
    #include <direct.h>
    #define MKDIR(dir) _mkdir(dir)
#else
    #include <sys/stat.h>
    #define MKDIR(dir) mkdir(dir, 0755)
#endif

void Sim::initialize(const SimConditions& sim_conditions) {
    //std::cout << "  Initializing simulation with conditions:" << std::endl;
    //std::cout << "    J: " << sim_conditions.J << std::endl;
    //std::cout << "    Beta: " << sim_conditions.Beta << std::endl;
    //std::cout << "    Tolerance: " << sim_conditions.tolerance << std::endl;
    
    conditions = sim_conditions;
    
    //std::cout << "  Initializing random number generators" << std::endl;
    gen = std::mt19937(rd());

    //std::cout << "  Calculating initial selection probabilities" << std::endl;
    calculate_initial_selection_probs_weighted(false);
    variables.backup();
    
    //std::cout << "  Computing initial system energy" << std::endl;
    system_energy = compute_total_energy();
    
    // Initialize grid configuration history with initial state
    // Clear any existing history and set the initial configuration
    grid_config_starting_history.clear();
    grid_config_potential_history.clear();
    possibly_changed_history.clear();
    proposed_atoms_history.clear();
    accepted_atoms_history.clear();
    accepted_possibly_changed_history.clear();
    
    //std::cout << "  Initialization complete" << std::endl;
}

void Sim::cluster_conditions(int cluster_condition) {
    switch(cluster_condition) {
        case 0:
            cluster_conditions_0();
            break;
        case 1:
            cluster_conditions_1();
            break;
        case 2:
            cluster_conditions_2();
            break;
        case 3:
            cluster_conditions_3();
            break;
        default:
            cluster_conditions_2();
            break;
    }
}

void Sim::cluster_conditions_0() {
    // Sample configuration - Basic binary patterns
    std::vector<std::vector<int>> property_matrix = {
        {-1}   // Cluster 1: patterns 0 and 1
    };
    cluster_conditions_from_indices_matrix(property_matrix, 1);
}
void Sim::cluster_conditions_1() {
    // Sample configuration - Alternating patterns
    std::vector<std::vector<int>> property_matrix = {
        {0,1},      // Cluster 1: patterns 0 and 1
        {4,5}          // Cluster 3: pattern 4 only
    };
    cluster_conditions_from_indices_matrix(property_matrix, 1);
}

void Sim::cluster_conditions_2() {
    // Sample configuration - Complex patterns  
    std::vector<std::vector<int>> property_matrix = {
        {1},   // Cluster 1: patterns 0, 4, and 8
        {0},      // Cluster 2: patterns 1 and 5
        {2,4},   // Cluster 3: patterns 2, 6, and 9
        {3,5}       // Cluster 4: patterns 3 and 7
    };
    cluster_conditions_from_indices_matrix(property_matrix, 1);
}

void Sim::cluster_conditions_3() {
    // Sample configuration - Dense pattern set
    std::vector<std::vector<int>> property_matrix = {
        {-1},
        {0, 2, 4, 6}, // Cluster 1: patterns 0, 2, 4, 6
        {1, 3, 5},    // Cluster 2: patterns 1, 3, 5
        {7, 8, 9},  // Cluster 3: patterns 7, 8, 9
        {10,11}
    };
    cluster_conditions_from_indices_matrix_uninverted(property_matrix, 2);
}

void Sim::reset_cluster_conditions(){
    if(grid.clusters.size()>0){
        for(int i=0;i<grid.Size*grid.Size;i++){
            grid.remove_atom(i);
        }
        grid.clusters.clear();
    }
    // grid.clusters.resize(0); // removed redundant resize to avoid default ctor requirement
    grid.num_clusters=0;
}

SwapType Sim::determine_swap_type(int atom1, int atom2) {
    // Get the cells corresponding to the two atoms
    Cell& cell1 = grid.array[atom1];
    Cell& cell2 = grid.array[atom2];

    // Check if the atoms are boundary or interior
    bool is_boundary1 = cell1.cell_type==CellType::BOUNDARY;
    bool is_boundary2 =cell2.cell_type==CellType::BOUNDARY;

    // Determine the swap type based on their statuses
    if (is_boundary1 && is_boundary2) {
        return SwapType::BB; // Both atoms are boundary atoms
    } else if (!is_boundary1 && !is_boundary2) {
        return SwapType::II; // Both atoms are interior atoms
    } else {
        return SwapType::BI; // One atom is boundary, the other is interior
    }
}

// Function to calculate the probability of selecting a specific swap type and pair of atoms
float Sim::swap_type_probability(Swap input_swap) {
    // Get the probability of selecting the swap type
   float swap_type_prob = 0.0f;
   float num_same_cluster;
   float num_different_cluster;
   float prob_same_cluster_base=conditions.prob_same_cluster_base;
   float prob_swap_type_base=0.0f;
    //initial swap type probabilities will be calculated based on the total number of possible swaps of said type over the total number of possible swaps
    //these terms will be biased via the prob_bb_base float values, but since they will always be normalized, we can simply use the following form: only the base values are needed here
    float prob_same_cluster_initial=0.0f;
    float prob_diff_cluster_initial=0.0f;
    switch (input_swap.this_swap_type) {
        case SwapType::BB:
            swap_type_prob = variables.prob_bb_initial;
            num_same_cluster = variables.num_bb_same_cluster;
            num_different_cluster = variables.num_bb_diff_cluster;
            prob_same_cluster_initial = variables.prob_bb_same_cluster_initial;
            prob_diff_cluster_initial = variables.prob_bb_diff_cluster_initial;
            prob_swap_type_base = conditions.prob_bb_base;
            break;
        case SwapType::BI:
            num_same_cluster = variables.num_bi_same_cluster;
            num_different_cluster = variables.num_bi_diff_cluster;
            prob_same_cluster_initial = variables.prob_bi_same_cluster_initial;
            prob_diff_cluster_initial = variables.prob_bi_diff_cluster_initial;
            prob_swap_type_base = conditions.prob_bi_base;
            break;
        case SwapType::II:
            num_same_cluster = variables.num_ii_same_cluster;
            num_different_cluster = variables.num_ii_diff_cluster;
            prob_same_cluster_initial = variables.prob_ii_same_cluster_initial;
            prob_diff_cluster_initial = variables.prob_ii_diff_cluster_initial;
            prob_swap_type_base = conditions.prob_ii_base;
            break;
        default:
            std::cout<<"Error: Invalid SwapType provided." << std::endl;
            std::cerr << "Error: Invalid SwapType provided." << std::endl;
            return 0.0f;
    }
    if(input_swap.same_cluster){
        swap_type_prob=prob_same_cluster_base*prob_swap_type_base;
    } else {
        swap_type_prob= (1.0f-prob_same_cluster_base)*prob_swap_type_base;
    }
    if(input_swap.primary_atom_clusterID>=0 && input_swap.secondary_atom_clusterID>=0){ 
    swap_type_prob*=grid.clusters[input_swap.primary_atom_clusterID].cluster_weight*grid.clusters[input_swap.secondary_atom_clusterID].cluster_weight;
    }
    else{return 0.0f;}
    // Return the combined probability
    return (swap_type_prob/variables.Norm);
}

float Sim::compute_total_energy() {
    float total_energy = 0.0f;
    float config_energy_sum = 0.0f;

    // Iterate only over the currently allowed lattice sites
    for (int i_idx = 0; i_idx < static_cast<int>(allowed_sites.size()); ++i_idx) {
        int i = allowed_sites[i_idx];
        std::vector<int> config = grid.array[i].local_config;
        if (config.empty() || config.size() < 5) {
            std::cout << "[DEBUG] Invalid configuration found for cell " << i << std::endl;
            continue; // Skip cells with invalid configurations
        }
        
        float energy_contribution = conditions.J * Energy.at(config);
        config_energy_sum += energy_contribution;
    }
    
    // We divide by 2 to avoid double-counting interactions
    total_energy = config_energy_sum ;
    
    return 2*total_energy;
}
void Sim::save_data(const Move& move) {
    try {
        // Save energy and rolling acceptance total
        energy_history.push_back(system_energy);
        total_acceptances_history.push_back(total_acceptances);
        // Save move data
        forward_prob_history.push_back(move.final_forward_probability);
        reverse_prob_history.push_back(move.final_reverse_probability);
        num_atoms_swapped_history.push_back(move.num_changed);
        delta_e_history.push_back(move.DeltaE);
        // Save Norm after each move
        norm_history.push_back(variables.Norm);
        
        // Extract possibly_changed atoms from the move
        std::vector<int> possibly_changed_atoms = move.possibly_changed_vec;
        possibly_changed_history.push_back(possibly_changed_atoms);
        
        // Extract primary_atoms from the move (flatten all tiers)
        std::vector<int> primary_atoms_flattened;
        for (const auto& tier : move.primary_atoms) {
            for (int atom : tier) {
                primary_atoms_flattened.push_back(atom);
            }
        }
        proposed_atoms_history.push_back(primary_atoms_flattened);
        
        // Save cluster sizes
        std::vector<int> current_boundary_sizes;
        std::vector<int> current_interior_sizes;
        if (grid.clusters.empty()) {
            throw std::runtime_error("No clusters found in grid");
        }
        for (const auto& cluster : grid.clusters) {
            if (cluster.num_species <= 0) {
                throw std::runtime_error("Invalid number of species in cluster");
            }
            int total_boundary = 0;
            int total_interior = 0;
            // Sum sizes across all species (now 0 and 1)
            for (int i = 0; i < cluster.num_species; i++) {
                if (i < cluster.boundary_cells.size() && i < cluster.interior_cells.size()) {
                    total_boundary += cluster.boundary_cells[i].size();
                    total_interior += cluster.interior_cells[i].size();
                }
            }
            current_boundary_sizes.push_back(total_boundary);
            current_interior_sizes.push_back(total_interior);
        }
        boundary_sizes_history.push_back(current_boundary_sizes);
        interior_sizes_history.push_back(current_interior_sizes);
    } catch (const std::exception& e) {
        std::cerr << "Error in save_data: " << e.what() << std::endl;
        // Continue execution but log the error
    }
}

void Sim::accept_move(Move& move) {
    // Update system energy
    system_energy += 2 * move.DeltaE;

    // Finalise the cluster modifications
    grid.apply_cluster_change(move.possibly_changed_vec);

    // Commit variable backup and statistics
    variables.backup();
    total_acceptances++;
    
    // For accepted moves, the accepted atoms are the primary_atoms from the move
    std::vector<int> accepted_atoms_flattened;
    for (const auto& tier : move.primary_atoms) {
        for (int atom : tier) {
            accepted_atoms_flattened.push_back(atom);
        }
    }
    accepted_atoms_history.push_back(accepted_atoms_flattened);
    
    // For accepted moves, also save the possibly_changed atoms
    accepted_possibly_changed_history.push_back(move.possibly_changed_vec);
    
    save_data(move);
}

void Sim::undo_move(Move& move) {
    //std::cout << "[DEBUG] undo_move called: DeltaE=" << move.DeltaE << ", num_changed=" << move.num_changed << std::endl;
    variables.restore();
    grid.undo_cluster_changes(move.possibly_changed_vec);
    
    // For rejected moves, no atoms were accepted
    std::vector<int> accepted_atoms; // Empty vector
    accepted_atoms_history.push_back(accepted_atoms);
    
    // For rejected moves, no possibly_changed atoms were accepted either
    std::vector<int> accepted_possibly_changed; // Empty vector
    accepted_possibly_changed_history.push_back(accepted_possibly_changed);
    
    // Save data
    save_data(move);
}

float Sim::expansion_probability(SwapType swap_type) {
    switch (swap_type) {
        case SwapType::BB:
            return conditions.prob_bb_expansion;
        case SwapType::BI:
            return conditions.prob_bi_expansion;
        case SwapType::II:
            return conditions.prob_ii_expansion;
        default:
            std::cerr << "Error: Invalid SwapType provided." << std::endl;
            return 0.0f;
    }
}

Swap Sim::calculate_reverse_swap(const Move& move) {
    // Use build_swap to construct the reverse swap from the current grid state
    return build_swap(move.forward.primary_atom, move.forward.secondary_atom,false);
}

void Sim::save_results_to_csv(const std::string& filename) {
    // Create outputs directory if it doesn't exist
    size_t pos = filename.find('/');
    if (pos != std::string::npos) {
        std::string dir = filename.substr(0, pos);
        struct stat st;
        if (stat(dir.c_str(), &st) != 0) {
            // Directory doesn't exist, create it
            if (MKDIR(dir.c_str()) != 0) {
                std::cerr << "Warning: Could not create directory " << dir << std::endl;
            }
        }
    }
    
    // Validate input data before attempting to save
    std::vector<std::string> empty_vectors;
    if (energy_history.empty()) empty_vectors.push_back("energy_history");
    if (total_acceptances_history.empty()) empty_vectors.push_back("total_acceptances_history");
    if (forward_prob_history.empty()) empty_vectors.push_back("forward_prob_history");
    if (reverse_prob_history.empty()) empty_vectors.push_back("reverse_prob_history");
    if (num_atoms_swapped_history.empty()) empty_vectors.push_back("num_atoms_swapped_history");
    if (delta_e_history.empty()) empty_vectors.push_back("delta_e_history");
    if (norm_history.empty()) empty_vectors.push_back("norm_history");
    if (possibly_changed_history.empty()) empty_vectors.push_back("possibly_changed_history");
    if (proposed_atoms_history.empty()) empty_vectors.push_back("proposed_atoms_history");
    if (accepted_atoms_history.empty()) empty_vectors.push_back("accepted_atoms_history");
    if (accepted_possibly_changed_history.empty()) empty_vectors.push_back("accepted_possibly_changed_history");
    if (grid_config_starting_history.empty()) empty_vectors.push_back("grid_config_starting_history");
    if (grid_config_potential_history.empty()) empty_vectors.push_back("grid_config_potential_history");
    if (boundary_sizes_history.empty()) empty_vectors.push_back("boundary_sizes_history");
    if (interior_sizes_history.empty()) empty_vectors.push_back("interior_sizes_history");
    if (!empty_vectors.empty()) {
        std::cerr << "Error: The following data vectors are empty: ";
        for (size_t i = 0; i < empty_vectors.size(); ++i) {
            std::cerr << empty_vectors[i];
            if (i + 1 < empty_vectors.size()) std::cerr << ", ";
        }
        std::cerr << std::endl;
        throw std::runtime_error("Cannot save results - missing data");
    }

    // Verify all vectors have same size
    size_t expected_size = energy_history.size();
    if (total_acceptances_history.size() != expected_size ||
        forward_prob_history.size() != expected_size ||
        reverse_prob_history.size() != expected_size ||
        num_atoms_swapped_history.size() != expected_size ||
        delta_e_history.size() != expected_size ||
        norm_history.size() != expected_size ||
        possibly_changed_history.size() != expected_size ||
        proposed_atoms_history.size() != expected_size ||
        accepted_atoms_history.size() != expected_size ||
        accepted_possibly_changed_history.size() != expected_size ||
        grid_config_starting_history.size() != expected_size ||
        grid_config_potential_history.size() != expected_size ||
        boundary_sizes_history.size() != expected_size ||
        interior_sizes_history.size() != expected_size) {
        std::cerr << "Error: Data vector size mismatch" << std::endl;
        throw std::runtime_error("Cannot save results - inconsistent data sizes");
    }

    // Open file for writing
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    // Write header
    file << "Energy,Total_Acceptances,ForwardProb,ReverseProb,NumAtomsSwapped,DeltaE,Norm,PossiblyChanged,ProposedAtoms,AcceptedAtoms,AcceptedPossiblyChanged,GridConfigStarting,GridConfigPotential";
    for (size_t i = 0; i < grid.clusters.size(); ++i) {
        file << ",BoundarySize_Cluster" << i << ",InteriorSize_Cluster" << i;
    }
    file << "\n";

    // Write data
    for (size_t i = 0; i < energy_history.size(); ++i) {
        file << energy_history[i] << ","
             << total_acceptances_history[i] << ","
             << forward_prob_history[i] << ","
             << reverse_prob_history[i] << ","
             << num_atoms_swapped_history[i] << ","
             << delta_e_history[i] << ","
             << norm_history[i] << ",";

        // Write possibly_changed atoms as comma-separated list enclosed in quotes
        file << "\"";
        if (i < possibly_changed_history.size()) {
            for (size_t j = 0; j < possibly_changed_history[i].size(); ++j) {
                file << static_cast<int>(possibly_changed_history[i][j]);
                if (j < possibly_changed_history[i].size() - 1) {
                    file << ",";
                }
            }
        }
        file << "\",";

        // Write proposed atoms (primary_atoms) as comma-separated list enclosed in quotes
        file << "\"";
        if (i < proposed_atoms_history.size()) {
            for (size_t j = 0; j < proposed_atoms_history[i].size(); ++j) {
                file << static_cast<int>(proposed_atoms_history[i][j]);
                if (j < proposed_atoms_history[i].size() - 1) {
                    file << ",";
                }
            }
        }
        file << "\",";

        // Write accepted atoms as comma-separated list enclosed in quotes
        file << "\"";
        if (i < accepted_atoms_history.size()) {
            for (size_t j = 0; j < accepted_atoms_history[i].size(); ++j) {
                file << static_cast<int>(accepted_atoms_history[i][j]);
                if (j < accepted_atoms_history[i].size() - 1) {
                    file << ",";
                }
            }
        }
        file << "\",";

        // Write accepted possibly changed atoms as comma-separated list enclosed in quotes
        file << "\"";
        if (i < accepted_possibly_changed_history.size()) {
            for (size_t j = 0; j < accepted_possibly_changed_history[i].size(); ++j) {
                file << static_cast<int>(accepted_possibly_changed_history[i][j]);
                if (j < accepted_possibly_changed_history[i].size() - 1) {
                    file << ",";
                }
            }
        }
        file << "\",";

        // Write starting grid configuration enclosed in quotes
        file << "\"";
        if (i < grid_config_starting_history.size()) {
            file << grid_config_starting_history[i];
        }
        file << "\",";

        // Write potential grid configuration enclosed in quotes
        file << "\"";
        if (i < grid_config_potential_history.size()) {
            file << grid_config_potential_history[i];
        }
        file << "\"";

        // Write cluster sizes
        for (size_t j = 0; j < grid.clusters.size(); ++j) {
            if (j < boundary_sizes_history[i].size() && j < interior_sizes_history[i].size()) {
                file << "," << static_cast<int>(boundary_sizes_history[i][j]) << ","
                     << static_cast<int>(interior_sizes_history[i][j]);
            } else {
                file << ",NA,NA";
            }
        }
        file << "\n";
    }

    file.close();
    if (!file) {
        throw std::runtime_error("Error occurred while closing file");
    }

    std::cout << "Successfully saved results to " << filename << std::endl;
}

void Sim::save_grid_configuration_to_csv(const std::string& filename) {
    // Create outputs directory if it doesn't exist
    size_t pos = filename.find('/');
    if (pos != std::string::npos) {
        std::string dir = filename.substr(0, pos);
        struct stat st;
        if (stat(dir.c_str(), &st) != 0) {
            // Directory doesn't exist, create it
            if (MKDIR(dir.c_str()) != 0) {
                std::cerr << "Warning: Could not create directory " << dir << std::endl;
            }
        }
    }
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        //std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    // Write header
    file << "x,y,species,cluster_id,cell_type\n";

    // Write data for each cell
    for (int y = 0; y < grid.Size; y++) {
        for (int x = 0; x < grid.Size; x++) {
            int index = y * grid.Size + x;
            const Cell& cell = grid.array[index];
            
            // Convert cell type to string
            std::string cell_type_str;
            switch (cell.cell_type) {
                case CellType::BOUNDARY:
                    cell_type_str = "BOUNDARY";
                    break;
                case CellType::INTERIOR:
                    cell_type_str = "INTERIOR";
                    break;
                default:
                    cell_type_str = "ANY";
            }

            file << x << "," << y << "," 
                 << cell.cell_species << "," 
                 << cell.cell_cluster_ID << "," 
                 << cell_type_str << "\n";
        }
    }

    file.close();
}

std::string Sim::grid_config_to_string() const {
    std::ostringstream config_stream;
    for (int i = 0; i < grid.Size * grid.Size; i++) {
        const Cell& cell = grid.array[i];
        // Format: species|cluster_id|cell_type
        char cell_type_char;
        switch (cell.cell_type) {
            case CellType::BOUNDARY: cell_type_char = 'B'; break;
            case CellType::INTERIOR: cell_type_char = 'I'; break;
            default: cell_type_char = 'A'; break;
        }
        config_stream << cell.cell_species << "|" << cell.cell_cluster_ID << "|" << cell_type_char;
        if (i < grid.Size * grid.Size - 1) {
            config_stream << ";";
        }
    }
    return config_stream.str();
}

Sim::Sim() {
    conditions = SimConditions(1.0f, 1.0f, 1.0f, 0.5f, 0.5f, 0.5f, 1.0f, 1.0f, 0.1f);
    gen = std::mt19937(rd());
    total_acceptances = 0;
    grid.minsize5_boundary=conditions.minsize5_boundary_setting;

    // Initialise the allowed_sites vector to cover the full lattice
    allowed_sites.reserve(grid.Size * grid.Size);
    for (int i = 0; i < grid.Size * grid.Size; ++i) {
        allowed_sites.push_back(i);
    }
}

Swap Sim::select_initial_atoms_weighted(bool debug) {
    if (debug) {
        std::cout << "\n=== SELECT INITIAL ATOMS WEIGHTED (DEBUG) ===" << std::endl;
    }

    // ------------------------------------------------------------
    // STEP-0 : choose the primary cluster according to cluster_prob
    // ------------------------------------------------------------
    // Use UNNORMALISED cluster probabilities (raw pair counts) for the draw
    double total_raw_weight = static_cast<double>(variables.Norm);
    if (total_raw_weight <= 0.0) {
        std::cerr << "[ERROR] select_initial_atoms_weighted: Norm is zero" << std::endl;
        return build_swap(-1, -1, debug);
    }

    std::uniform_real_distribution<double> dist01_raw(0.0, 1.0);
    double target_primary = dist01_raw(gen);
    double cum_raw        = 0.0;

    int primary_cluster = -1;
    for (int c = 0; c < grid.num_clusters; ++c) {
        double raw_w = static_cast<double>(grid.clusters[c].cluster_prob);
        cum_raw += raw_w;
        if (target_primary < cum_raw) {
            primary_cluster = c;
            break;
        }
    }
    if (primary_cluster == -1) {
        std::cerr << "[ERROR] select_initial_atoms_weighted: failed to pick primary cluster" << std::endl;
        return build_swap(-1, -1, debug);
    }
    Cluster *A = &grid.clusters[primary_cluster];

    if (debug) {
        std::cout << "  Picked primary cluster " << primary_cluster << " (r=" << target_primary << ")" << std::endl;
    }

    // ----------------------------------------------------------------------
    // STEP-1 : compute the total weight of *all* valid (A,B) ordered pairs
    //          Weight = (#pairs) × base probabilities × cluster weights
    // ----------------------------------------------------------------------
    
    double total_pair_weight = grid.clusters[primary_cluster].cluster_prob*variables.Norm;

  
    if (total_pair_weight <= 0.0) {
        std::cerr << "[ERROR] select_initial_atoms_weighted: total_pair_weight == 0" << std::endl;
        return build_swap(-1, -1, debug);
    }

    // ----------------------------------------------------------------------
    // STEP-2 : draw a second *independent* uniform random number to pick one
    //          concrete pair according to the weights computed above
    // ----------------------------------------------------------------------
    float total_pairs=800.0f*800.0f;
    std::uniform_real_distribution<double> dist_total(0.0, total_pair_weight);
    double target = dist_total(gen);
    double acc    = 0.0;

    for (int b = 0; b <= primary_cluster; ++b) {
        Cluster *B = &grid.clusters[b];

        double cluster_factor = ((b == primary_cluster)
                                 ? conditions.prob_same_cluster_base
                                 : (1.0 - conditions.prob_same_cluster_base));
        cluster_factor *= A->cluster_weight * B->cluster_weight;

        for (int locA = 0; locA < 2; ++locA) {
            for (int locB = 0; locB < 2; ++locB) {
                for (int sp = 0; sp < 2; ++sp) {
                    if (b == primary_cluster && locA == locB && sp == 1) continue; // avoid BB/II double count
                    if(b==primary_cluster && locA!=locB && sp==1)continue;

                    int sizeA = (locA == 0) ? A->boundary_cells[sp].size()
                                            : A->interior_cells[sp].size();
                    int sizeB = (locB == 0) ? B->boundary_cells[1 - sp].size()
                                            : B->interior_cells[1 - sp].size();
                    if (sizeA == 0 || sizeB == 0) continue;

                    double swap_factor = cluster_factor;
                    if      (locA == 0 && locB == 0) swap_factor *= conditions.prob_bb_base;
                    else if (locA == 0 || locB == 0) swap_factor *= conditions.prob_bi_base;
                    else                             swap_factor *= conditions.prob_ii_base;

                    double contrib = static_cast<double>(sizeA) * sizeB * swap_factor/total_pairs;
                    if (contrib >= target) {
                        // We found the interval – sample one concrete pair uniformly
                        if (A->requires_refresh) A->refresh_vectors();
                        if (B->requires_refresh) B->refresh_vectors();

                        std::uniform_int_distribution<int> pick(0, sizeA * sizeB - 1);
                        int idx  = pick(gen);
                        int posA = idx % sizeA;
                        int posB = idx / sizeA;

                        const std::vector<int>& vecA = (locA == 0) ? A->boundary_atoms_vector[sp]
                                                                     : A->interior_atoms_vector[sp];
                        const std::vector<int>& vecB = (locB == 0) ? B->boundary_atoms_vector[1 - sp]
                                                                     : B->interior_atoms_vector[1 - sp];

                        return build_swap(vecA[posA], vecB[posB], debug);
                    }
                    target-=contrib;;
                }
            }
        }
    }

    // Should never get here
    std::cerr << "[ERROR] select_initial_atoms_weighted: exhaustive search failed" << std::endl;
    return build_swap(-1, -1, debug);
}

void Sim::copy_grid(Grid& grid_copy){
  
  if(grid.Size != grid_copy.Size){
    std::cerr << "[ERROR] copy_grid: grid sizes do not match" << std::endl;
    return;
  }
  for(int i=0; i< grid.Size * grid.Size; i++){
    grid.remove_atom(i);
    grid.array[i].cell_species = grid_copy.array[i].cell_species;
    grid.array[i].local_config=grid_copy.array[i].local_config;
  }

  // rebuild local_config for all cells to match new species layout
  for(int i=0;i<grid.Size*grid.Size;i++){
    grid.refresh_local_config(i);
  }

  system_energy = compute_total_energy();
  }
void Sim::reapply_cluster_weights(float cluster_bias, int tier_size) {
    if(tier_size <= 0) tier_size = 1; // safeguard
    for(int i = 0; i < grid.num_clusters; ++i) {
        if(i < tier_size) {
            grid.clusters[i].cluster_weight = 0.0f;
        }
        else if(cluster_bias != 0.0f) {
            float tier_index = static_cast<float>(i) / static_cast<float>(tier_size);
            grid.clusters[i].cluster_weight = std::exp(conditions.J * conditions.Beta * (tier_index) * (1.0f / cluster_bias));
        }
        else {
            grid.clusters[i].cluster_weight = 1.0f;
        }
    }
}

// -----------------------------------------------------------------------------
// Build cluster set from a vector of indices referencing the
// GLOBAL_PATTERN_VECTORS table (includes inversions)
// Each entry spawns a new min-size-5 cluster with that property.
// -----------------------------------------------------------------------------
void Sim::cluster_conditions_from_indices(const std::vector<int>& property_indices, bool catch_all) {
    // 1. Reset any existing clusters / assignments
    reset_cluster_conditions();

    // 2. Always add the original minsize-1 cluster with {0} and {1} properties first (0-weight catch-all)
    grid.add_cluster(1);
    grid.clusters.back().add_property({0});
    grid.clusters.back().add_property({1});

    // 3. Conditionally add an additional weighted catch-all cluster
    if (catch_all) {
        grid.add_cluster(1);
        grid.clusters.back().add_property({0});
        grid.clusters.back().add_property({1});
    }

    if (property_indices.empty()) {
        reapply_cluster_weights(conditions.cluster_bias);
        grid.assign_initial_clusters(allowed_sites);
        return;
    }

    // 4. For each index add a cluster with the specified property pattern
    for (int idx : property_indices) {
        if (idx >= 0 && idx < static_cast<int>(GLOBAL_PATTERN_VECTORS.size())) {
            // Regular minsize-5 cluster with the full pattern (rotations+inversions)
            grid.add_cluster(5);
            const auto &pattern = GLOBAL_PATTERN_VECTORS[idx];
            grid.clusters.back().add_property(pattern);
            grid.clusters.back().property.generate_inversions();
        } else {
            std::cerr << "[cluster_conditions_from_indices] Invalid property index " << idx << std::endl;
            continue;
        }
    }

    // 5. Apply cluster weights (tier_size defaults to 1)
    reapply_cluster_weights(conditions.cluster_bias, tier_size);

    // 6. Keep the original cluster at index 0 with 0 weight, but set the weighted catch-all appropriately
    if (catch_all && grid.num_clusters > 2) {
        grid.clusters[0].cluster_weight = 0.0f; // Original 0-weight catch-all
        grid.clusters[1].cluster_weight = grid.clusters.back().cluster_weight * 0.1f; // Weighted catch-all
    } else if (!catch_all && grid.num_clusters > 1) {
        grid.clusters[0].cluster_weight = 0.0f; // Original 0-weight catch-all
    }

    // 7. Finalise assignment
    grid.assign_initial_clusters(allowed_sites);
}

bool Sim::iterate_improved(bool debug) {
    // Recalculate selection probabilities before each attempt
    calculate_initial_selection_probs_weighted(false);

    // 1. Pick an initial swap according to the current probabilities
    Swap initial_swap = select_initial_atoms_weighted(debug);
    if (initial_swap.primary_atom == -1 || initial_swap.secondary_atom == -1) {
        // No valid swap found – treat as rejection, do nothing
        return false;
    }

    // 1.5. Capture the starting grid state before any changes
    grid_config_starting_history.push_back(grid_config_to_string());

    // 2. Build a Move deterministically from this swap
    Move move(conditions.J, conditions.Beta, conditions.tolerance);
    bool build_ok = move.deterministic_build_from_swap(initial_swap, grid, debug);
    if(!build_ok) {
        // Invalid move – revert any tentative changes and treat as rejection
        // Add empty potential grid (since move failed)
        grid_config_potential_history.push_back("");
        undo_move(move);
        return false;
    }

    // 2.5. Capture the potential grid state after move is applied but before accept/reject
    grid_config_potential_history.push_back(grid_config_to_string());

    // 3. Construct the corresponding reverse swap and obtain acceptance prob
    Swap reverse_swap = calculate_reverse_swap(move);
    float acceptance_prob = move.get_final_acceptance_deterministic(grid, reverse_swap, debug);

    // 4. Metropolis criterion
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    float rand_num = dist(gen);

    if (rand_num < acceptance_prob) {
        accept_move(move);
        return true;
    } else {
        undo_move(move);
        return false;
    }
}

// -----------------------------------------------------------------------------
// Build cluster set from indices referencing the rotation-only library
// (GLOBAL_PATTERN_VECTORS_ROT).  No inversions are added automatically.
// -----------------------------------------------------------------------------
void Sim::cluster_conditions_from_indices_uninverted(const std::vector<int>& property_indices, int tier_size, bool catch_all) {
    // 1. Reset clusters
    reset_cluster_conditions();

    // 2. Always add the original universal catch-all cluster ({0}+{1}) at index 0 (0-weight)
    grid.add_cluster(1);
    grid.clusters.back().add_property({0});
    grid.clusters.back().add_property({1});

    // 3. Conditionally add an additional weighted catch-all cluster
    if (catch_all) {
        grid.add_cluster(1);
        grid.clusters.back().add_property({0});
        grid.clusters.back().add_property({1});
    }

    // 4. Early exit if no additional indices supplied
    if(property_indices.empty()) {
        reapply_cluster_weights(conditions.cluster_bias, tier_size);
        if (catch_all && grid.num_clusters > 1) {
            grid.clusters[0].cluster_weight = 0.0f; // Original 0-weight catch-all
            grid.clusters[1].cluster_weight = 1.0f; // Weighted catch-all
        } else {
            grid.clusters[0].cluster_weight = 0.0f; // Original 0-weight catch-all
        }
        grid.assign_initial_clusters(allowed_sites);
        return;
    }

    // 5. Process indices
    for(int idx : property_indices) {
        if(idx >= 0 && idx < static_cast<int>(GLOBAL_PATTERN_VECTORS_ROT.size())) {
            grid.add_cluster(5);
            const auto &pattern = GLOBAL_PATTERN_VECTORS_ROT[idx];
            grid.clusters.back().add_property(pattern); // rotations handled by ClusterProperty
            // << no inversions >>
        } else {
            std::cerr << "[cluster_conditions_from_indices_uninverted] Invalid property index " << idx << std::endl;
        }
    }

    // 6. Apply weighted tiers
    reapply_cluster_weights(conditions.cluster_bias, tier_size);
    if(catch_all && grid.num_clusters > 2) {
        grid.clusters[0].cluster_weight = 0.0f; // Original 0-weight catch-all
        grid.clusters[1].cluster_weight = grid.clusters.back().cluster_weight * 0.1f; // Weighted catch-all
    } else if (!catch_all && grid.num_clusters > 1) {
        grid.clusters[0].cluster_weight = 0.0f; // Original 0-weight catch-all
    }

    // 7. Finalise
    grid.assign_initial_clusters(allowed_sites);
}

void Sim::calculate_initial_selection_probs_weighted(bool debug){
    //debug=true;
    variables.Norm = 0.0;
    // Clear any existing data
    variables.num_bb = 0;
    variables.num_bi = 0;
    variables.num_ii = 0;

    variables.num_bb_same_cluster = 0;
    variables.num_bi_same_cluster = 0;
    variables.num_ii_same_cluster = 0;

    variables.num_bb_diff_cluster = 0;
    variables.num_bi_diff_cluster = 0;
    variables.num_ii_diff_cluster = 0;

    variables.prob_bb_same_cluster_initial=0.0;
    variables.prob_bi_same_cluster_initial=0.0;
    variables.prob_ii_same_cluster_initial=0.0;
    variables.prob_bb_diff_cluster_initial=0.0;
    variables.prob_bi_diff_cluster_initial=0.0;
    variables.prob_ii_diff_cluster_initial=0.0;
    
    
    if (debug) {
        std::cout << "\n=== CALCULATE INITIAL SELECTION PROBABILITIES (DEBUG) ===" << std::endl;
        std::cout << "Step 1: Clearing existing data and initializing counters" << std::endl;
    }
   variables.num_bb = grid.total_boundary_atoms[0] * grid.total_boundary_atoms[1];
    
    // Count boundary-interior pairs between different species
    variables.num_bi = grid.total_boundary_atoms[0] * grid.total_interior_atoms[1]
                    + grid.total_interior_atoms[0] * grid.total_boundary_atoms[1];
    
    // Count interior-interior pairs between species 0 and 1
    variables.num_ii = grid.total_interior_atoms[0] * grid.total_interior_atoms[1];


    float total_pairs=800.0f*800.0f;
    float total_pairs_grid=total_pairs;
    float total_swaps_checked=0.0f;
    double norm_accum = 0.0;
    for (int i = 0; i < grid.num_clusters; ++i) {
        Cluster &A = grid.clusters[i];
        A.cluster_prob = 0.0f;

        for (int j = 0; j <= i; ++j) {
            Cluster &B = grid.clusters[j];

            double cluster_factor = (i == j ? conditions.prob_same_cluster_base
                                            : (1.0 - conditions.prob_same_cluster_base));
            cluster_factor *= A.cluster_weight * B.cluster_weight;

            for (int locA = 0; locA < 2; ++locA) {
                for (int locB = 0; locB < 2; ++locB) {
                    for (int sp = 0; sp < 2; ++sp) {
                        // Avoid double-counting when A and B are the same cluster and the two
                        // locations are the same (BB or II).  In that situation, the pair with
                        // sp==1 is the mirror image of the pair with sp==0.
                        if (i == j && locA == locB && sp == 1) continue;
                        if(i==j && locA!=locB && sp==1)continue;
                        int nA = (locA == 0) ? A.boundary_cells[sp].size() : A.interior_cells[sp].size();
                        int nB = (locB == 0) ? B.boundary_cells[1 - sp].size() : B.interior_cells[1 - sp].size();
                        if (nA == 0 || nB == 0) continue;

                        double swap_factor = cluster_factor;
                        if      (locA == 0 && locB == 0) swap_factor *= conditions.prob_bb_base;
                        else if (locA == 0 || locB == 0) swap_factor *= conditions.prob_bi_base;
                        else                             swap_factor *= conditions.prob_ii_base;

                        total_swaps_checked+=nA*nB;

                        double contrib = static_cast<double>(nA) *nB * swap_factor/total_pairs;
                        if(debug){
                            std::cout<< "i=" << i << ", j=" << j << std::endl;
                            std::cout << "nA=" << nA << ", nB=" << nB << ", swap_factor=" << swap_factor << std::endl;}
                        if (i == j) {
                            if      ((locA == 0) && (locB == 0)) { variables.prob_bb_same_cluster_initial += contrib; }
                            else if ((locA == 0) || (locB == 0)) { variables.prob_bi_same_cluster_initial += contrib; }
                            else                                  { variables.prob_ii_same_cluster_initial += contrib; }
                        } else {
                            if      ((locA == 0) && (locB == 0)) { variables.prob_bb_diff_cluster_initial += contrib; }
                            else if ((locA == 0) || (locB == 0)) { variables.prob_bi_diff_cluster_initial += contrib; }
                            else                                  { variables.prob_ii_diff_cluster_initial += contrib; }
                        }

                        A.cluster_prob += contrib;
                    }
                }
            }
        }
        norm_accum += A.cluster_prob;
    }
    if(debug){
        std::cout << "total_swaps_checked=" << total_swaps_checked << std::endl;
    }
    // normalise
    variables.Norm = norm_accum;
    if (variables.Norm > 0.0f) {
        for (auto &cl : grid.clusters) cl.cluster_prob /= variables.Norm;

        variables.prob_bb_same_cluster_initial /= variables.Norm;
        variables.prob_bi_same_cluster_initial /= variables.Norm;
        variables.prob_ii_same_cluster_initial /= variables.Norm;
        variables.prob_bb_diff_cluster_initial /= variables.Norm;
        variables.prob_bi_diff_cluster_initial /= variables.Norm;
        variables.prob_ii_diff_cluster_initial /= variables.Norm;
    }

    variables.prob_bb_initial = variables.prob_bb_same_cluster_initial + variables.prob_bb_diff_cluster_initial;
    variables.prob_bi_initial = variables.prob_bi_same_cluster_initial + variables.prob_bi_diff_cluster_initial;
    variables.prob_ii_initial = variables.prob_ii_same_cluster_initial + variables.prob_ii_diff_cluster_initial;


    
    if (debug) {
        std::cout << "\nCluster probabilities:" << std::endl;
        for (int i = 0; i < grid.num_clusters; i++) {
            std::cout << "  Cluster " << i << " probability: " << grid.clusters[i].cluster_prob << std::endl;
        }
        std::cout << "\nStep 4: Final probability calculations" << std::endl;
        std::cout << "Individual probabilities:" << std::endl;
        std::cout << "  prob_bb_same_cluster_initial: " << variables.prob_bb_same_cluster_initial << std::endl;
        std::cout << "  prob_bb_diff_cluster_initial: " << variables.prob_bb_diff_cluster_initial << std::endl;
        std::cout << "  prob_bi_same_cluster_initial: " << variables.prob_bi_same_cluster_initial << std::endl;
        std::cout << "  prob_bi_diff_cluster_initial: " << variables.prob_bi_diff_cluster_initial << std::endl;
        std::cout << "  prob_ii_same_cluster_initial: " << variables.prob_ii_same_cluster_initial << std::endl;
        std::cout << "  prob_ii_diff_cluster_initial: " << variables.prob_ii_diff_cluster_initial << std::endl;
        
        std::cout << "\nCombined swap type probabilities:" << std::endl;
        std::cout << "  prob_bb_initial: " << variables.prob_bb_initial << std::endl;
        std::cout << "  prob_bi_initial: " << variables.prob_bi_initial << std::endl;
        std::cout << "  prob_ii_initial: " << variables.prob_ii_initial << std::endl;
        
        float total_prob = variables.prob_bb_initial + variables.prob_bi_initial + variables.prob_ii_initial;
        std::cout << "  Total probability (should be 1.0): " << total_prob << std::endl;
        std::cout << "=== END CALCULATE INITIAL SELECTION PROBABILITIES ===" << std::endl;
    }
}

// Implementation for build_swap
Swap Sim::build_swap(int primary_atom, int secondary_atom, bool debug) {
    Swap swap;
    
    if (primary_atom < 0 || primary_atom >= grid.Size * grid.Size ||
        secondary_atom < 0 || secondary_atom >= grid.Size * grid.Size) {
        // Invalid atoms
        swap.primary_atom = -1;
        swap.secondary_atom = -1;
        return swap;
    }
    
    // Fill in swap details
    swap.primary_atom = primary_atom;
    swap.secondary_atom = secondary_atom;
    swap.primary_atom_species = grid.array[primary_atom].cell_species;
    swap.secondary_atom_species = grid.array[secondary_atom].cell_species;
    swap.primary_atom_clusterID = grid.array[primary_atom].cell_cluster_ID;
    swap.secondary_atom_clusterID = grid.array[secondary_atom].cell_cluster_ID;
    swap.primary_atom_type = grid.array[primary_atom].cell_type;
    swap.secondary_atom_type = grid.array[secondary_atom].cell_type;
    
    // Determine swap type
    swap.this_swap_type = determine_swap_type(primary_atom, secondary_atom);
    
    // Check if same cluster
    swap.same_cluster = (swap.primary_atom_clusterID == swap.secondary_atom_clusterID);
    
    // Calculate selection probability
    swap.initial_selection_prob = swap_type_probability(swap);
    
    // Set expansion probability based on swap type
    swap.expansion_prob = expansion_probability(swap.this_swap_type);
    
    if (debug) {
        std::cout << "[DEBUG] Built swap: " << primary_atom << " <-> " << secondary_atom << std::endl;
        std::cout << "  Species: " << swap.primary_atom_species << " <-> " << swap.secondary_atom_species << std::endl;
        std::cout << "  Clusters: " << swap.primary_atom_clusterID << " <-> " << swap.secondary_atom_clusterID << std::endl;
        std::cout << "  Same cluster: " << (swap.same_cluster ? "YES" : "NO") << std::endl;
        std::cout << "  Selection prob: " << swap.initial_selection_prob << std::endl;
    }
    
    return swap;
}

bool Sim::delayed_rejection_iteration(int max_num_swaps, bool debug){
    calculate_initial_selection_probs_weighted(debug);
     //debug=true;
    if (debug) std::cout << "[DEBUG][Sim] Starting delayed_rejection_iteration with max_num_swaps=" << max_num_swaps << std::endl;
    std::unordered_set<int> visited;
    Move move(conditions.J, conditions.Beta, conditions.tolerance);
    move.DeltaE=0.0f;
    if (debug) std::cout << "[DEBUG][Sim] Selecting random swaps..." << std::endl;
    std::vector<Swap> random_swaps = select_random_swaps(max_num_swaps, false);
    if (random_swaps.empty()) {
        std::cout << "[DEBUG][Sim] No valid swaps found, saving dummy move data." << std::endl;
        // Capture starting grid state even for failed moves
        grid_config_starting_history.push_back(grid_config_to_string());
        grid_config_potential_history.push_back(""); // Empty potential grid for failed move
        
        // Fill a dummy move with default values
        Move dummy_move(conditions.J, conditions.Beta, conditions.tolerance);
        dummy_move.DeltaE = 0.0f;
        dummy_move.num_changed = 0;
        dummy_move.final_forward_probability = 0.0f;
        dummy_move.final_reverse_probability = 0.0f;
        dummy_move.forward.this_swap_type = SwapType::BB;
        
        // For delayed rejection with no valid swaps, we still need to populate proposed atoms from random_swaps
        // but since random_swaps is empty, the proposed atoms will be empty, which is correct
        undo_move(dummy_move);
        return false;
    }
    
    // Capture starting grid state before any changes
    grid_config_starting_history.push_back(grid_config_to_string());
    
    move.max_num_atoms=random_swaps.size()*2;
    if (random_swaps.back().primary_atom_species == random_swaps.back().secondary_atom_species) {

        std::cout << "[DEBUG][Sim] Move rejected, same species" << std::endl;
        
        // Add empty potential grid for failed move
        grid_config_potential_history.push_back("");
        undo_move(move);
        return false;
    }

    if (debug) std::cout << "[DEBUG][Sim] Building move from swaps..." << std::endl;
    move.deterministic_build_from_swap_delayed_rejection(random_swaps, grid, debug);
    if (debug) std::cout << "[DEBUG][Sim] Evaluating candidates..." << std::endl;
    move.evaluate_candidates_dr(grid, debug);
    move.calculate_cluster_changes(grid, debug);
    
    // Capture potential grid state after move is applied but before accept/reject
    grid_config_potential_history.push_back(grid_config_to_string());
    
    if (debug) std::cout << "[DEBUG][Sim] Calculating reverse swaps..." << std::endl;
    std::vector<Swap> reverse_swaps = calculate_reverse_swaps(random_swaps, debug);
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    if (debug) std::cout << "[DEBUG][Sim] Calculating acceptance probability..." << std::endl;
    float acceptance_prob = move.calculate_reverse_delayed_rejection(reverse_swaps, grid, debug);
    float rand_num = dist(gen);
    int total_expansions=0;
    if (debug) {std::cout << "[DEBUG][Sim] Random number for acceptance: " << rand_num << " (acceptance_prob=" << acceptance_prob << ")" << std::endl;}

    if(rand_num < acceptance_prob){
        if (debug) std::cout << "[DEBUG][Sim] Move accepted on first try." << std::endl;
        accept_move(move);
        return true;
    }
    else{



        if (debug) std::cout << "[DEBUG][Sim] Move rejected, starting expansion loop..." << std::endl;
        // Begin expansion loop after rejection
        while (total_expansions<1 && random_swaps.begin()->expansion_prob>0.0f){
            std::vector<Swap> new_swaps=expand_swap_dr_greater_structure(grid, *random_swaps.begin(), total_expansions, debug);
            total_expansions++;
            float deltaE=move.DeltaE;
            move.swap_history.push_back(new_swaps);
            move.evaluate_candidates_dr(grid, debug);

            float acceptance_y_star = std::min(1.0f,std::exp(-conditions.Beta*(-deltaE)));
            float ratio1= (1.0f-acceptance_y_star)/(1.0f-acceptance_prob);
            float ratio2=std::exp(-conditions.Beta*move.DeltaE);
            float new_acceptance_prob=ratio1*ratio2;
            if (debug) {
                std::cout << "[DEBUG][Expansion] new_acceptance_prob=" << new_acceptance_prob << ", ratio1=" << ratio1 << ", ratio2=" << ratio2 << std::endl;
            }
            float rand_num = dist(gen);
            if(rand_num < new_acceptance_prob){
                move.calculate_cluster_changes(grid, debug);
                accept_move(move);
                return true;
            }
            else{
                acceptance_prob=new_acceptance_prob;
            }
        }
        if (debug) std::cout << "[DEBUG][Sim] No further expansions possible, move rejected. Final DeltaE=" << move.DeltaE << std::endl;
        // If we reach here, the move was rejected after all expansions. Save the data for the last attempted move.
        undo_move(move);

        if (debug) std::cout << "[DEBUG][Sim] Exiting delayed_rejection_iteration with REJECT after all expansions." << std::endl;
        return false;
    }
}

// Implementation for select_random_swaps
std::vector<Swap> Sim::select_random_swaps(int max_num_swaps, bool debug) {
    std::vector<Swap> swaps;
    std::unordered_set<int> used_atoms;
    
    for (int i = 0; i < max_num_swaps; i++) {
        // Try to find a valid swap
        int attempts = 0;
        const int max_attempts = 100;
        
        while (attempts < max_attempts) {
            Swap candidate_swap = select_initial_atoms_weighted(false);
            
            if (candidate_swap.primary_atom == -1 || candidate_swap.secondary_atom == -1) {
                attempts++;
                continue;
            }
            
            // Check if atoms are already used
            if (used_atoms.find(candidate_swap.primary_atom) != used_atoms.end() ||
                used_atoms.find(candidate_swap.secondary_atom) != used_atoms.end()) {
                attempts++;
                continue;
            }
            
            // Valid swap found
            used_atoms.insert(candidate_swap.primary_atom);
            used_atoms.insert(candidate_swap.secondary_atom);
            swaps.push_back(candidate_swap);
            break;
        }
        
        if (attempts >= max_attempts) {
            if (debug) std::cout << "Failed to find swap " << i << " after " << max_attempts << " attempts" << std::endl;
            break;
        }
    }
    
    if (debug) {
        std::cout << "Selected " << swaps.size() << " swaps out of " << max_num_swaps << " requested" << std::endl;
    }
    
    return swaps;
}

// Implementation for calculate_reverse_swaps
std::vector<Swap> Sim::calculate_reverse_swaps(std::vector<Swap> input_swaps, bool debug) {
    std::vector<Swap> reverse_swaps;
    
    for (const Swap& forward_swap : input_swaps) {
        Swap reverse_swap = build_swap(forward_swap.secondary_atom, forward_swap.primary_atom, debug);
        reverse_swaps.push_back(reverse_swap);
    }
    
    return reverse_swaps;
}

std::vector<Swap> Sim::expand_swap_dr_greater_structure(Grid& grid, Swap swap, int iteration, bool debug){
    std::vector<Swap> expanded_swaps;
    std::unordered_set<int> visited;
    visited.insert(swap.primary_atom);
    visited.insert(swap.secondary_atom);
    std::vector<int> primary_neighbors=grid.array[swap.primary_atom].neighbor_indexes;
    std::vector<int> secondary_neighbors={grid.array[swap.secondary_atom].neighbor_indexes[2],grid.array[swap.secondary_atom].neighbor_indexes[3],grid.array[swap.secondary_atom].neighbor_indexes[0],grid.array[swap.secondary_atom].neighbor_indexes[1]};


    for(int i=0;i<primary_neighbors.size();i++){
        if(visited.find(primary_neighbors[i])==visited.end() && visited.find(secondary_neighbors[i])==visited.end()){
            visited.insert(primary_neighbors[i]);
            visited.insert(secondary_neighbors[i]);
            expanded_swaps.push_back(build_swap(primary_neighbors[i], secondary_neighbors[i],debug));

        }
        std::vector<int> new_primary_neighbors=grid.array[primary_neighbors[i]].neighbor_indexes;
        std::vector<int> new_secondary_neighbors={grid.array[secondary_neighbors[i]].neighbor_indexes[2],grid.array[secondary_neighbors[i]].neighbor_indexes[3],grid.array[secondary_neighbors[i]].neighbor_indexes[0],grid.array[secondary_neighbors[i]].neighbor_indexes[1]};
        if(iteration>0){
            for(int j=0;j<new_primary_neighbors.size();j++){
                if(visited.find(new_primary_neighbors[j])==visited.end() && visited.find(new_secondary_neighbors[j])==visited.end()){
                    visited.insert(new_primary_neighbors[j]);
                    visited.insert(new_secondary_neighbors[j]); 
                    expanded_swaps.push_back(build_swap(new_primary_neighbors[j], new_secondary_neighbors[j],debug));
                }
            }
        }
    }
    return expanded_swaps;
}

// Build clusters from a matrix where each row represents a cluster
// Each element in a row is an index into GLOBAL_PATTERN_VECTORS (includes inversions)
// Special case: if first element of first row is -1, add weighted catch-all cluster
void Sim::cluster_conditions_from_indices_matrix(const std::vector<std::vector<int>>& property_matrix, int tier_size) {
    // 1. Reset any existing clusters / assignments
    reset_cluster_conditions();

    // 2. Always add the original minsize-1 cluster with {0} and {1} properties first (0-weight catch-all)
    grid.add_cluster(1);
    grid.clusters.back().add_property({0});
    grid.clusters.back().add_property({1});

    // 3. Check if we should add weighted catch-all cluster (first element of first row is -1)
    bool use_catch_all = false;
    if (!property_matrix.empty() && !property_matrix[0].empty() && property_matrix[0][0] == -1) {
        use_catch_all = true;
        grid.add_cluster(1);
        grid.clusters.back().add_property({0});
        grid.clusters.back().add_property({1});
    }

    if (property_matrix.empty()) {
        reapply_cluster_weights(conditions.cluster_bias, tier_size);
        grid.clusters[0].cluster_weight = 0.0f; // Original 0-weight catch-all
        grid.assign_initial_clusters(allowed_sites);
        return;
    }

    // 4. For each row in the matrix, create a minsize-5 cluster
    for (const auto& cluster_indices : property_matrix) {
        if (cluster_indices.empty()) {
            continue; // Skip empty rows
        }
        
        // Skip the special -1 catch-all indicator row
        if (cluster_indices.size() == 1 && cluster_indices[0] == -1) {
            continue;
        }
        
        // Add a new minsize-5 cluster
        grid.add_cluster(5);
        
        // Add each pattern index as a property to this cluster
        for (int idx : cluster_indices) {
            if (idx >= 0 && idx < static_cast<int>(GLOBAL_PATTERN_VECTORS.size())) {
                const auto &pattern = GLOBAL_PATTERN_VECTORS[idx];
                grid.clusters.back().add_property(pattern);
            } else if (idx != -1) { // Don't warn about -1 as it's a special case
                std::cerr << "[cluster_conditions_from_indices_matrix] Invalid property index " << idx << std::endl;
            }
        }
        
        // Generate inversions for this cluster
        grid.clusters.back().property.generate_inversions();
    }

    // 5. Apply cluster weights
    reapply_cluster_weights(conditions.cluster_bias, tier_size);

    // 6. Set appropriate weights for catch-all clusters
    if (use_catch_all && grid.num_clusters > 2) {
        grid.clusters[0].cluster_weight = 0.0f; // Original 0-weight catch-all
        grid.clusters[1].cluster_weight = grid.clusters.back().cluster_weight * 0.1f; // Weighted catch-all
    } else {
        grid.clusters[0].cluster_weight = 0.0f; // Original 0-weight catch-all
    }

    // 7. Finalise assignment
    grid.assign_initial_clusters(allowed_sites);
}

// Build clusters from a matrix using rotation-only patterns (no inversions)
// Each row represents a cluster, each element is an index into GLOBAL_PATTERN_VECTORS_ROT
// Special case: if first element of first row is -1, add weighted catch-all cluster
void Sim::cluster_conditions_from_indices_matrix_uninverted(const std::vector<std::vector<int>>& property_matrix, int tier_size) {
    // 1. Reset clusters
    reset_cluster_conditions();

    // 2. Always add the original universal catch-all cluster ({0}+{1}) at index 0 (0-weight)
    grid.add_cluster(1);
    grid.clusters.back().add_property({0});
    grid.clusters.back().add_property({1});

    // 3. Check if we should add weighted catch-all cluster (first element of first row is -1)
    bool use_catch_all = false;
    if (!property_matrix.empty() && !property_matrix[0].empty() && property_matrix[0][0] == -1) {
        use_catch_all = true;
        grid.add_cluster(1);
        grid.clusters.back().add_property({0});
        grid.clusters.back().add_property({1});
    }

    // 4. Early exit if no matrix provided
    if(property_matrix.empty()) {
        reapply_cluster_weights(conditions.cluster_bias, tier_size);
        grid.clusters[0].cluster_weight = 0.0f; // Original 0-weight catch-all
        grid.assign_initial_clusters(allowed_sites);
        return;
    }

    // 5. Process each row as a cluster
    for(const auto& cluster_indices : property_matrix) {
        if (cluster_indices.empty()) {
            continue; // Skip empty rows
        }
        
        // Skip the special -1 catch-all indicator row
        if (cluster_indices.size() == 1 && cluster_indices[0] == -1) {
            continue;
        }
        
        // Add a new minsize-5 cluster
        grid.add_cluster(5);
        
        // Add each pattern index as a property to this cluster
        for(int idx : cluster_indices) {
            if(idx >= 0 && idx < static_cast<int>(GLOBAL_PATTERN_VECTORS_ROT.size())) {
                const auto &pattern = GLOBAL_PATTERN_VECTORS_ROT[idx];
                grid.clusters.back().add_property(pattern); // rotations handled by ClusterProperty
                // << no inversions >>
            } else if (idx != -1) { // Don't warn about -1 as it's a special case
                std::cerr << "[cluster_conditions_from_indices_matrix_uninverted] Invalid property index " << idx << std::endl;
            }
        }
    }

    // 6. Apply weighted tiers
    reapply_cluster_weights(conditions.cluster_bias, tier_size);
    
    // 7. Set appropriate weights for catch-all clusters
    if (use_catch_all && grid.num_clusters > 2) {
        grid.clusters[0].cluster_weight = 0.0f; // Original 0-weight catch-all
        grid.clusters[1].cluster_weight = grid.clusters.back().cluster_weight * 0.1f; // Weighted catch-all
    } else {
        grid.clusters[0].cluster_weight = 0.0f; // Original 0-weight catch-all
    }

    grid.assign_initial_clusters(allowed_sites);
}

// ===============================================
// NEW: OPTIMIZATION SYSTEM IMPLEMENTATION
// ===============================================

// CHUNK 1: MAIN OPTIMIZATION ORCHESTRATOR
void Sim::run_cluster_seed_optimization(int iterations, bool debug) {
    if (debug) {
        std::cout << "\n========== CLUSTER SEED OPTIMIZATION STARTED ==========" << std::endl;
        std::cout << "[DEBUG] Target iterations per sim: " << iterations << std::endl;
        std::cout << "[DEBUG] Original system energy: " << system_energy << std::endl;
    }
    
    // Step 1: Generate all possible cluster combinations
    if (debug) std::cout << "\n--- CHUNK 1: GENERATE COMBINATIONS ---" << std::endl;
    std::vector<std::vector<std::vector<int>>> combinations = generate_cluster_combinations(debug);
    
    if (combinations.empty()) {
        if (debug) std::cout << "[DEBUG] No combinations generated, optimization aborted" << std::endl;
        return;
    }
    
    // Step 2: Initialize optimizer sims with combinations
    if (debug) std::cout << "\n--- CHUNK 2: INITIALIZE SIMS ---" << std::endl;
    initialize_optimizer_sims(combinations, debug);
    
    // Step 3: Run tests on all optimizer sims
    if (debug) std::cout << "\n--- CHUNK 3: RUN TESTS ---" << std::endl;
    run_optimizer_tests(iterations, debug);
    
    // Step 4: Select best performing sim
    if (debug) std::cout << "\n--- CHUNK 4: SELECT BEST ---" << std::endl;
    int best_index = select_best_optimizer_sim(debug);
    
    // Step 5: Apply best cluster seed to main sim
    if (debug) std::cout << "\n--- CHUNK 5: APPLY BEST ---" << std::endl;
    apply_best_cluster_seed(debug);
    
    if (debug) {
        std::cout << "\n========== OPTIMIZATION COMPLETED ==========" << std::endl;
        debug_print_optimization_summary(best_index, debug);
    }
}

// CHUNK 1: COMBINATION GENERATION
std::vector<std::vector<std::vector<int>>> Sim::generate_cluster_combinations(bool debug) {
    if (debug) std::cout << "[DEBUG] Starting combination generation..." << std::endl;
    
    // Step 1: Vectorize current cluster seed
    std::vector<std::vector<int>> vectorized_seed = vectorize_cluster_seed(debug);
    
    if (vectorized_seed.empty()) {
        if (debug) std::cout << "[DEBUG] Warning: Empty cluster seed vectorized" << std::endl;
        return std::vector<std::vector<std::vector<int>>>();
    }
    
    // Step 2: Determine if inversions are enabled
    bool use_inversions = true;  // Default assumption
    
    // Check if any cluster uses inversions by looking at cluster properties
    if (grid.num_clusters > 0) {
        // Check the first non-catch-all cluster to see if inversions are being used
        for (int i = 1; i < grid.num_clusters && i < static_cast<int>(grid.clusters.size()); i++) {
            ClusterProperty& prop = grid.clusters[i].property;
            // If cluster has multiple patterns that look like inversions, inversions are on
            if (prop.local_configs.size() > 12) { // Heuristic: more than base rotations
                use_inversions = true;
                break;
            }
        }
    }
    
    if (debug) {
        std::cout << "[DEBUG] Use inversions: " << (use_inversions ? "true" : "false") << std::endl;
        std::cout << "[DEBUG] Target cluster count: " << (use_inversions ? "4" : "2") << std::endl;
    }
    
    // Step 3: Create combinations from vectorized entries
    std::vector<std::vector<std::vector<int>>> combinations = 
        create_combinations_from_entries(vectorized_seed, use_inversions, debug);
    
    // Step 4: Distribute properties evenly
    distribute_properties_evenly(combinations, debug);
    
    if (debug) {
        std::cout << "[DEBUG] Generated " << combinations.size() << " total combinations" << std::endl;
    }
    
    return combinations;
}

// Helper: Vectorize cluster seed
std::vector<std::vector<int>> Sim::vectorize_cluster_seed(bool debug) {
    if (debug) std::cout << "[DEBUG] Vectorizing cluster seed..." << std::endl;
    
    std::vector<std::vector<int>> vectorized_seed;
    
    // Extract property indices from each cluster (skip catch-all cluster 0)
    for (int cluster_id = 1; cluster_id < grid.num_clusters && cluster_id < static_cast<int>(grid.clusters.size()); cluster_id++) {
        std::vector<int> cluster_properties;
        ClusterProperty& prop = grid.clusters[cluster_id].property;
        
        // Find which global pattern indices this cluster represents
        // This is a simplified approach - we'll collect the first few unique patterns
        for (size_t i = 0; i < prop.local_configs.size() && cluster_properties.size() < 3; i++) {
            bool found_match = false;
            std::vector<int> config = prop.local_configs[i];
            
            // Search in GLOBAL_PATTERN_VECTORS to find matching index
            for (size_t global_idx = 0; global_idx < GLOBAL_PATTERN_VECTORS.size(); global_idx++) {
                if (GLOBAL_PATTERN_VECTORS[global_idx] == config) {
                    // Check if we already have this index
                    if (std::find(cluster_properties.begin(), cluster_properties.end(), static_cast<int>(global_idx)) == cluster_properties.end()) {
                        cluster_properties.push_back(static_cast<int>(global_idx));
                        found_match = true;
                        break;
                    }
                }
            }
            
            if (found_match && cluster_properties.size() >= 3) break;
        }
        
        // If we couldn't find matches, create some default entries
        if (cluster_properties.empty()) {
            cluster_properties.push_back(cluster_id - 1); // Simple fallback
        }
        
        vectorized_seed.push_back(cluster_properties);
        
        if (debug) {
            std::cout << "[DEBUG] Cluster " << cluster_id << " vectorized to: [";
            for (size_t i = 0; i < cluster_properties.size(); i++) {
                std::cout << cluster_properties[i];
                if (i < cluster_properties.size() - 1) std::cout << ",";
            }
            std::cout << "]" << std::endl;
        }
    }
    
    if (debug) {
        std::cout << "[DEBUG] Vectorized " << vectorized_seed.size() << " clusters" << std::endl;
    }
    
    return vectorized_seed;
}

// Helper: Create combinations from entries
std::vector<std::vector<std::vector<int>>> Sim::create_combinations_from_entries(
    const std::vector<std::vector<int>>& entries, bool use_inversions, bool debug) {
    
    if (debug) {
        std::cout << "[DEBUG] Creating combinations from " << entries.size() << " entry groups" << std::endl;
    }
    
    std::vector<std::vector<std::vector<int>>> combinations;
    
    // Flatten all entries into a single pool
    std::vector<int> all_entries;
    for (const auto& entry_group : entries) {
        for (int entry : entry_group) {
            all_entries.push_back(entry);
        }
    }
    
    if (debug) {
        std::cout << "[DEBUG] Total entries to distribute: " << all_entries.size() << std::endl;
    }
    
    // Determine target number of clusters
    int target_clusters = use_inversions ? 4 : 2;
    
    if (all_entries.empty()) {
        if (debug) std::cout << "[DEBUG] No entries to create combinations from" << std::endl;
        return combinations;
    }
    
    // Generate all possible ways to distribute entries across target_clusters
    // This is a complex combinatorial problem, so we'll use a simplified approach
    
    // Method: Generate some representative combinations
    int max_combinations = 50; // Limit to prevent explosion
    int combinations_generated = 0;
    
    // Basic combination: distribute entries round-robin
    std::vector<std::vector<int>> basic_combination(target_clusters);
    for (size_t i = 0; i < all_entries.size(); i++) {
        basic_combination[i % target_clusters].push_back(all_entries[i]);
    }
    combinations.push_back(basic_combination);
    combinations_generated++;
    
    if (debug) {
        std::cout << "[DEBUG] Generated basic round-robin combination" << std::endl;
        debug_print_combination(basic_combination, 0, debug);
    }
    
    // Generate some permutations by rotating the distribution
    for (int shift = 1; shift < target_clusters && combinations_generated < max_combinations; shift++) {
        std::vector<std::vector<int>> shifted_combination(target_clusters);
        for (size_t i = 0; i < all_entries.size(); i++) {
            shifted_combination[(i + shift) % target_clusters].push_back(all_entries[i]);
        }
        combinations.push_back(shifted_combination);
        combinations_generated++;
        
        if (debug) {
            std::cout << "[DEBUG] Generated shifted combination (shift=" << shift << ")" << std::endl;
        }
    }
    
    // Generate combinations with different clustering strategies
    if (all_entries.size() >= 4 && combinations_generated < max_combinations) {
        // Strategy: Group consecutive entries
        std::vector<std::vector<int>> grouped_combination(target_clusters);
        int entries_per_cluster = all_entries.size() / target_clusters;
        int remainder = all_entries.size() % target_clusters;
        
        int entry_index = 0;
        for (int cluster = 0; cluster < target_clusters; cluster++) {
            int cluster_size = entries_per_cluster + (cluster < remainder ? 1 : 0);
            for (int i = 0; i < cluster_size && entry_index < static_cast<int>(all_entries.size()); i++) {
                grouped_combination[cluster].push_back(all_entries[entry_index++]);
            }
        }
        combinations.push_back(grouped_combination);
        combinations_generated++;
        
        if (debug) {
            std::cout << "[DEBUG] Generated grouped combination" << std::endl;
        }
    }
    
    if (debug) {
        std::cout << "[DEBUG] Generated " << combinations_generated << " combinations total" << std::endl;
    }
    
    return combinations;
}

// Helper: Distribute properties evenly
void Sim::distribute_properties_evenly(std::vector<std::vector<std::vector<int>>>& combinations, bool debug) {
    if (debug) std::cout << "[DEBUG] Distributing properties evenly across combinations..." << std::endl;
    
    // This function ensures properties are distributed evenly within each combination
    for (size_t combo_idx = 0; combo_idx < combinations.size(); combo_idx++) {
        auto& combination = combinations[combo_idx];
        
        if (debug) {
            std::cout << "[DEBUG] Processing combination " << combo_idx << std::endl;
        }
        
        // Count total properties in this combination
        int total_properties = 0;
        for (const auto& cluster : combination) {
            total_properties += cluster.size();
        }
        
        if (total_properties == 0) continue;
        
        int target_clusters = combination.size();
        int properties_per_cluster = total_properties / target_clusters;
        int remainder = total_properties % target_clusters;
        
        if (debug) {
            std::cout << "[DEBUG] Total properties: " << total_properties 
                      << ", Target per cluster: " << properties_per_cluster 
                      << ", Remainder: " << remainder << std::endl;
        }
        
        // Collect all properties and redistribute
        std::vector<int> all_properties;
        for (const auto& cluster : combination) {
            for (int prop : cluster) {
                all_properties.push_back(prop);
            }
        }
        
        // Clear and redistribute
        for (auto& cluster : combination) {
            cluster.clear();
        }
        
        int prop_index = 0;
        for (int cluster_idx = 0; cluster_idx < target_clusters; cluster_idx++) {
            int cluster_size = properties_per_cluster + (cluster_idx < remainder ? 1 : 0);
            for (int i = 0; i < cluster_size && prop_index < static_cast<int>(all_properties.size()); i++) {
                combination[cluster_idx].push_back(all_properties[prop_index++]);
            }
        }
        
        if (debug) {
            std::cout << "[DEBUG] Redistributed combination " << combo_idx << ":" << std::endl;
            debug_print_combination(combination, combo_idx, debug);
        }
    }
    
    if (debug) {
        std::cout << "[DEBUG] Property distribution completed for " << combinations.size() << " combinations" << std::endl;
    }
}

// CHUNK 2: SIMULATOR INITIALIZATION
void Sim::initialize_optimizer_sims(const std::vector<std::vector<std::vector<int>>>& combinations, bool debug) {
    if (debug) {
        std::cout << "[DEBUG] Initializing " << combinations.size() << " optimizer sims..." << std::endl;
    }
    
    // Clear any existing optimizer sims
    optimizer_sims.clear();
    optimizer_sims.reserve(combinations.size());
    
    for (size_t i = 0; i < combinations.size(); i++) {
        if (debug) {
            std::cout << "[DEBUG] Initializing optimizer sim " << i << std::endl;
        }
        
        // Create new sim object
        Sim new_sim;
        
        // Copy parent properties (conditions, allowed sites, etc.)
        copy_parent_properties(new_sim, debug);
        
        // Apply the specific cluster combination
        try {
            // Reset cluster conditions first
            new_sim.reset_cluster_conditions();
            
            // Apply the combination using the matrix method
            if (debug) {
                std::cout << "[DEBUG] Applying cluster combination " << i << ":" << std::endl;
                debug_print_combination(combinations[i], i, debug);
            }
            
            // Determine which method to use based on whether inversions are used
            bool use_inversions = true;  // This should match the logic from generate_cluster_combinations
            
            if (use_inversions) {
                new_sim.cluster_conditions_from_indices_matrix(combinations[i], 1);
            } else {
                new_sim.cluster_conditions_from_indices_matrix_uninverted(combinations[i], 2);
            }
            
            // Initialize the new sim
            new_sim.initialize(new_sim.conditions);
            
            if (debug) {
                std::cout << "[DEBUG] Sim " << i << " initialized successfully. Energy: " 
                          << new_sim.system_energy << std::endl;
            }
            
        } catch (const std::exception& e) {
            if (debug) {
                std::cout << "[DEBUG] Error initializing sim " << i << ": " << e.what() << std::endl;
            }
            // Create a dummy sim with high energy to mark as failed
            new_sim.system_energy = 1e9f;
        }
        
        optimizer_sims.push_back(std::move(new_sim));
    }
    
    if (debug) {
        std::cout << "[DEBUG] Successfully initialized " << optimizer_sims.size() << " optimizer sims" << std::endl;
    }
}

// CHUNK 3: TEST EXECUTION
void Sim::run_optimizer_tests(int iterations, bool debug) {
    if (debug) {
        std::cout << "[DEBUG] Running " << iterations << " iterations on " 
                  << optimizer_sims.size() << " optimizer sims..." << std::endl;
    }
    
    for (size_t sim_idx = 0; sim_idx < optimizer_sims.size(); sim_idx++) {
        if (debug) {
            std::cout << "[DEBUG] Running test on sim " << sim_idx 
                      << " (initial energy: " << optimizer_sims[sim_idx].system_energy << ")" << std::endl;
        }
        
        Sim& test_sim = optimizer_sims[sim_idx];
        
        // Skip if sim failed to initialize properly
        if (test_sim.system_energy > 1e8f) {
            if (debug) {
                std::cout << "[DEBUG] Skipping sim " << sim_idx << " (failed initialization)" << std::endl;
            }
            continue;
        }
        
        float initial_energy = test_sim.system_energy;
        int accepted_moves = 0;
        
        // Run the iterations
        for (int iter = 0; iter < iterations; iter++) {
            try {
                bool accepted = test_sim.iterate_improved(debug && (iter % 20 == 0)); // Debug every 20th iteration
                if (accepted) {
                    accepted_moves++;
                }
            } catch (const std::exception& e) {
                if (debug) {
                    std::cout << "[DEBUG] Error during iteration " << iter << " on sim " << sim_idx 
                              << ": " << e.what() << std::endl;
                }
                break; // Stop iterations for this sim on error
            }
        }
        
        float final_energy = test_sim.system_energy;
        float energy_change = final_energy - initial_energy;
        
        if (debug) {
            std::cout << "[DEBUG] Sim " << sim_idx << " completed: " 
                      << accepted_moves << "/" << iterations << " accepted, "
                      << "Energy: " << initial_energy << " -> " << final_energy 
                      << " (Δ=" << energy_change << ")" << std::endl;
        }
    }
    
    if (debug) {
        std::cout << "[DEBUG] All optimizer tests completed" << std::endl;
    }
}

// CHUNK 4: BEST SELECTION
int Sim::select_best_optimizer_sim(bool debug) {
    if (debug) {
        std::cout << "[DEBUG] Selecting best optimizer sim from " << optimizer_sims.size() << " candidates..." << std::endl;
    }
    
    if (optimizer_sims.empty()) {
        if (debug) std::cout << "[DEBUG] No optimizer sims available for selection" << std::endl;
        return -1;
    }
    
    int best_index = 0;
    float best_energy = optimizer_sims[0].system_energy;
    
    if (debug) {
        std::cout << "[DEBUG] Energy comparison:" << std::endl;
    }
    
    for (size_t i = 0; i < optimizer_sims.size(); i++) {
        float energy = optimizer_sims[i].system_energy;
        
        if (debug) {
            std::cout << "[DEBUG]   Sim " << i << ": " << energy;
            if (energy < best_energy) std::cout << " <- NEW BEST";
            std::cout << std::endl;
        }
        
        if (energy < best_energy) {
            best_energy = energy;
            best_index = i;
        }
    }
    
    if (debug) {
        std::cout << "[DEBUG] Selected sim " << best_index << " with energy " << best_energy << std::endl;
    }
    
    return best_index;
}

// CHUNK 5: APPLY BEST RESULT
void Sim::apply_best_cluster_seed(bool debug) {
    if (debug) {
        std::cout << "[DEBUG] Applying best cluster seed to main simulation..." << std::endl;
    }
    
    int best_index = select_best_optimizer_sim(debug);
    
    if (best_index < 0 || best_index >= static_cast<int>(optimizer_sims.size())) {
        if (debug) {
            std::cout << "[DEBUG] Invalid best index " << best_index << ", cannot apply" << std::endl;
        }
        return;
    }
    
    Sim& best_sim = optimizer_sims[best_index];
    
    if (debug) {
        std::cout << "[DEBUG] Copying configuration from sim " << best_index << std::endl;
        std::cout << "[DEBUG] Original energy: " << system_energy << std::endl;
        std::cout << "[DEBUG] Best energy: " << best_sim.system_energy << std::endl;
    }
    
    // Store original energy for comparison
    float original_energy = system_energy;
    
    // Reset current cluster conditions
    reset_cluster_conditions();
    
    // Copy the grid state from the best sim
    try {
        // Copy grid configuration
        copy_grid(best_sim.grid);
        
        // Copy cluster configuration
        grid.clusters = best_sim.grid.clusters;
        grid.num_clusters = best_sim.grid.num_clusters;
        
        // Update system energy
        system_energy = best_sim.system_energy;
        
        // Reinitialize with current conditions to ensure consistency
        variables.backup();
        calculate_initial_selection_probs_weighted(debug);
        
        if (debug) {
            std::cout << "[DEBUG] Cluster seed application completed" << std::endl;
            std::cout << "[DEBUG] Energy change: " << original_energy << " -> " << system_energy 
                      << " (improvement: " << (original_energy - system_energy) << ")" << std::endl;
        }
        
    } catch (const std::exception& e) {
        if (debug) {
            std::cout << "[DEBUG] Error applying best cluster seed: " << e.what() << std::endl;
        }
        // Restore original energy if copy failed
        system_energy = original_energy;
    }
}

// UTILITY FUNCTIONS

void Sim::copy_parent_properties(Sim& child_sim, bool debug) {
    if (debug) {
        std::cout << "[DEBUG] Copying parent properties to child sim..." << std::endl;
    }
    
    // Copy simulation conditions
    child_sim.conditions = conditions;
    
    // Copy allowed sites
    child_sim.allowed_sites = allowed_sites;
    
    // Copy grid size and basic structure
    child_sim.grid.Size = grid.Size;
    
    // Initialize random number generator
    child_sim.gen = std::mt19937(child_sim.rd());
    
    if (debug) {
        std::cout << "[DEBUG] Copied conditions, allowed_sites (" << allowed_sites.size() 
                  << " sites), and grid size (" << grid.Size << "x" << grid.Size << ")" << std::endl;
    }
}

void Sim::debug_print_combination(const std::vector<std::vector<int>>& combination, int index, bool debug) {
    if (!debug) return;
    
    std::cout << "[DEBUG]   Combination " << index << ":" << std::endl;
    for (size_t cluster_idx = 0; cluster_idx < combination.size(); cluster_idx++) {
        std::cout << "[DEBUG]     Cluster " << cluster_idx << ": [";
        for (size_t prop_idx = 0; prop_idx < combination[cluster_idx].size(); prop_idx++) {
            std::cout << combination[cluster_idx][prop_idx];
            if (prop_idx < combination[cluster_idx].size() - 1) std::cout << ",";
        }
        std::cout << "]" << std::endl;
    }
}

void Sim::debug_print_optimization_summary(int best_index, bool debug) {
    if (!debug) return;
    
    std::cout << "\n[DEBUG] ========== OPTIMIZATION SUMMARY ==========" << std::endl;
    std::cout << "[DEBUG] Total combinations tested: " << optimizer_sims.size() << std::endl;
    std::cout << "[DEBUG] Best performing sim: " << best_index << std::endl;
    
    if (best_index >= 0 && best_index < static_cast<int>(optimizer_sims.size())) {
        std::cout << "[DEBUG] Best energy: " << optimizer_sims[best_index].system_energy << std::endl;
        std::cout << "[DEBUG] Best sim cluster count: " << optimizer_sims[best_index].grid.num_clusters << std::endl;
    }
    
    std::cout << "[DEBUG] Final main sim energy: " << system_energy << std::endl;
    std::cout << "[DEBUG] =============================================" << std::endl;
}