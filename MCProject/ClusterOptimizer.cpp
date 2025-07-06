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
#include <unordered_set>
#include <unordered_map>
#include <mutex>
#include <fstream>
#include <chrono>
#include <limits>
#include <omp.h>

// =============================================================================
// CLUSTER OPTIMIZER: Successor to Simtest.cpp
// Goals:
// 1. Check every possible cluster property combination for fastest convergence
// 2. Implement restricted partitioning schemes
// =============================================================================

// Global containers for optimization results
std::unordered_map<std::vector<std::vector<int>>, float, VectorVectorHash> tested_configurations;
std::mutex config_mutex;

// Custom hash function for cluster configurations
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

// Standardize cluster configuration for consistent comparison
std::vector<std::vector<int>> standardize_configuration(const std::vector<std::vector<int>>& config) {
    if (config.empty()) return config;
    
    std::vector<std::vector<int>> result = config;
    
    // Sort and remove duplicates within each cluster
    for (auto& cluster : result) {
        if (!cluster.empty()) {
            std::sort(cluster.begin(), cluster.end());
            cluster.erase(std::unique(cluster.begin(), cluster.end()), cluster.end());
        }
    }
    
    // Sort clusters by first element for consistent ordering
    std::sort(result.begin(), result.end(), [](const std::vector<int>& a, const std::vector<int>& b) {
        if (a.empty() && b.empty()) return false;
        if (a.empty()) return true;
        if (b.empty()) return false;
        return a[0] < b[0];
    });
    
    return result;
}

// Test a single cluster configuration and return final energy
float test_cluster_configuration(const std::vector<std::vector<int>>& config, 
                                const std::string& conditions, int iterations) {
    
    // Check if already tested
    std::vector<std::vector<int>> std_config = standardize_configuration(config);
    {
        std::lock_guard<std::mutex> lock(config_mutex);
        auto it = tested_configurations.find(std_config);
        if (it != tested_configurations.end()) {
            return it->second;
        }
    }
    
    // Setup and run simulation
    Sim sim;
    sim.conditions = SimConditions(SIM_CONDITIONS.at(conditions));
    sim.cluster_conditions(2);
    sim.cluster_conditions_from_indices_matrix(config, 1);
    sim.initialize(sim.conditions);
    
    // Check if simulation is viable
    if (sim.variables.Norm == 0.0) {
        std::lock_guard<std::mutex> lock(config_mutex);
        tested_configurations[std_config] = 1e6f;
        return 1e6f;
    }
    
    // Run simulation
    for (int iter = 0; iter < iterations; ++iter) {
        sim.iterate_improved();
    }
    
    float final_energy = sim.system_energy;
    
    // Store result
    {
        std::lock_guard<std::mutex> lock(config_mutex);
        tested_configurations[std_config] = final_energy;
    }
    
    return final_energy;
}
    
// Find optimal cluster configuration (Goal 1)
std::vector<std::vector<int>> find_optimal_cluster_configuration(const std::string& conditions, 
                                                               int iterations, 
                                                               int max_clusters, 
                                                               int max_properties_per_cluster) {
    std::cout << "[GOAL 1] Finding optimal cluster configuration..." << std::endl;
    
    int num_properties = static_cast<int>(GLOBAL_PATTERN_VECTORS.size());
    std::vector<std::vector<int>> best_config;
    float best_energy = std::numeric_limits<float>::max();
    int configs_tested = 0;
    
    // Generate and test all possible configurations
    std::function<void(std::vector<std::vector<int>>&, std::vector<int>&)> test_configs = 
        [&](std::vector<std::vector<int>>& current_config, std::vector<int>& unused_props) {
            
            // Test current configuration if it has at least 1 cluster
            if (!current_config.empty()) {
                float energy = test_cluster_configuration(current_config, conditions, iterations);
                configs_tested++;
                
                if (energy < best_energy) {
                    best_energy = energy;
                    best_config = current_config;
                    std::cout << "[INFO] New best config found with energy: " << energy << std::endl;
                }
            }
            
            // Stop if at maximum clusters
            if (current_config.size() >= max_clusters) return;
            
            // Try adding new clusters
            for (size_t i = 0; i < unused_props.size(); ++i) {
                std::vector<std::vector<int>> new_config = current_config;
                new_config.push_back({unused_props[i]});
                
                std::vector<int> new_unused = unused_props;
                new_unused.erase(new_unused.begin() + i);
                
                test_configs(new_config, new_unused);
            }
            
            // Try adding properties to existing clusters
            for (size_t cluster_idx = 0; cluster_idx < current_config.size(); ++cluster_idx) {
                if (current_config[cluster_idx].size() < max_properties_per_cluster) {
                    for (size_t prop_idx = 0; prop_idx < unused_props.size(); ++prop_idx) {
                        std::vector<std::vector<int>> new_config = current_config;
                        new_config[cluster_idx].push_back(unused_props[prop_idx]);
                        
                        std::vector<int> new_unused = unused_props;
                        new_unused.erase(new_unused.begin() + prop_idx);
                        
                        test_configs(new_config, new_unused);
                    }
                }
            }
        };
    
    // Start with all properties available
    std::vector<int> all_properties;
    for (int i = 0; i < std::min(num_properties, 8); ++i) { // Limit to first 8 properties for efficiency
        all_properties.push_back(i);
    }
    
    std::vector<std::vector<int>> empty_config;
    test_configs(empty_config, all_properties);
    
    std::cout << "[INFO] Tested " << configs_tested << " configurations" << std::endl;
    std::cout << "[INFO] Best energy: " << best_energy << std::endl;
    
    return best_config;
}

// Run simulation with restricted partitioning (Goal 2)
std::string run_partitioned_simulation(const std::string& conditions, 
                                     const std::vector<std::vector<int>>& cluster_config,
                                     int iterations, int num_partitions) {
    std::cout << "[GOAL 2] Running " << num_partitions << "-partition simulation..." << std::endl;
    
    // Generate filename
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << "ClusterOptimizer_PART" << num_partitions << "_" << conditions 
       << "_iter" << iterations << "_" << time_t << ".csv";
    std::string filename = ss.str();
    
    // Create spatial partitions
    int grid_size = 40;
    int partition_size = grid_size / static_cast<int>(std::sqrt(num_partitions));
    std::vector<std::vector<int>> partitions(num_partitions);
    
    for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            int site_index = i * grid_size + j;
            int partition_row = i / partition_size;
            int partition_col = j / partition_size;
            int partition_id = partition_row * (grid_size / partition_size) + partition_col;
            
            if (partition_id < num_partitions) {
                partitions[partition_id].push_back(site_index);
            }
        }
    }
    
    // Run simulations on each partition
    std::vector<float> partition_energies(num_partitions);
    
    #pragma omp parallel for
    for (int p = 0; p < num_partitions; ++p) {
        Sim sim;
        sim.conditions = SimConditions(SIM_CONDITIONS.at(conditions));
        sim.cluster_conditions(2);
        sim.allowed_sites = partitions[p];
        
        if (!cluster_config.empty()) {
            sim.cluster_conditions_from_indices_matrix(cluster_config, 1);
        }
        
        sim.initialize(sim.conditions);
        
        // Run simulation
        for (int iter = 0; iter < iterations / num_partitions; ++iter) {
            sim.iterate_improved();
        }
        
        partition_energies[p] = sim.system_energy;
    }
    
    // Save results
    std::ofstream file(filename);
    if (file.is_open()) {
        file << "partition_id,final_energy\n";
        for (int p = 0; p < num_partitions; ++p) {
            file << p << "," << partition_energies[p] << "\n";
        }
        float avg_energy = 0.0f;
        for (float e : partition_energies) avg_energy += e;
        avg_energy /= num_partitions;
        file << "average," << avg_energy << "\n";
        file.close();
    }
    
    std::cout << "[INFO] Partitioned simulation saved to: " << filename << std::endl;
    return filename;
}

// Standard simulation runner (borrowed from original)
std::string run_standard_simulation(const std::string& conditions, int iterations, 
                                  const std::vector<std::vector<int>>& cluster_config) {
    
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << "ClusterOptimizer_" << conditions << "_iter" << iterations << "_" << time_t << ".csv";
    std::string filename = ss.str();
    
    Sim sim;
    sim.conditions = SimConditions(SIM_CONDITIONS.at(conditions));
    sim.cluster_conditions(2);
    
    if (!cluster_config.empty()) {
        sim.cluster_conditions_from_indices_matrix(cluster_config, 1);
    }
    
    sim.initialize(sim.conditions);
    
    for (int iter = 0; iter < iterations; ++iter) {
        sim.iterate_improved();
    }
    
    sim.save_results_to_csv(filename);
    std::cout << "[INFO] Standard simulation saved to: " << filename << std::endl;
    return filename;
}

// Delayed rejection simulation runner (borrowed from original)
std::string run_delayed_rejection_simulation(const std::string& conditions, int iterations, 
                                           int max_swaps, const std::vector<std::vector<int>>& cluster_config) {
    
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << "ClusterOptimizer_DR_" << conditions << "_swaps" << max_swaps 
       << "_iter" << iterations << "_" << time_t << ".csv";
    std::string filename = ss.str();
    
    Sim sim;
    sim.conditions = SimConditions(SIM_CONDITIONS.at(conditions));
    sim.cluster_conditions(2);
    
    if (!cluster_config.empty()) {
        sim.cluster_conditions_from_indices_matrix(cluster_config, 1);
    }
    
    sim.initialize(sim.conditions);
    
    for (int iter = 0; iter < iterations; ++iter) {
        sim.delayed_rejection_iteration(max_swaps, false);
    }
    
    sim.save_results_to_csv(filename);
    std::cout << "[INFO] Delayed rejection simulation saved to: " << filename << std::endl;
    return filename;
}

// Main program
int main() {
    // Initialize global libraries
    initialize_energy();
    initialize_pattern_library();
    
    std::cout << "=============================================================================\n";
    std::cout << "                            CLUSTER OPTIMIZER\n";
    std::cout << "=============================================================================\n";
    
    // Parameters
    const std::string base_conditions = "BJ=-0.44_DR";
    const int iterations = 500;
    std::vector<std::string> output_files;
    
    // =============================================================================
    // GOAL 1: FIND OPTIMAL CLUSTER CONFIGURATION
    // =============================================================================
    
    std::vector<std::vector<int>> optimal_config = find_optimal_cluster_configuration(
        base_conditions, iterations, 3, 2);
    
    // Run simulation with optimal configuration
    std::string optimal_file = run_standard_simulation(base_conditions, iterations, optimal_config);
    output_files.push_back(optimal_file);
    
    // =============================================================================
    // GOAL 2: RESTRICTED PARTITIONING SCHEMES
    // =============================================================================
    
    // Test different partition counts
    for (int partitions : {4, 9, 16}) {
        std::string part_file = run_partitioned_simulation(base_conditions, optimal_config, 
                                                         iterations, partitions);
        output_files.push_back(part_file);
    }
    
    // =============================================================================
    // COMPARATIVE BENCHMARKS
    // =============================================================================
    
    // Kawasaki benchmark
    std::string kawasaki_file = run_standard_simulation("BJ=-0.44_Kawasaki", iterations, {});
    output_files.push_back(kawasaki_file);
    
    // Delayed rejection variants
    for (int swaps : {1, 2, 3}) {
        std::string dr_file = run_delayed_rejection_simulation(base_conditions, iterations, 
                                                             swaps, optimal_config);
        output_files.push_back(dr_file);
    }
    
    // =============================================================================
    // OUTPUT RESULTS
    // =============================================================================
    
    std::cout << "\n=============================================================================\n";
    std::cout << "                            RESULTS SUMMARY\n";
    std::cout << "=============================================================================\n";
    
    // Output filenames in JSON format
    std::cout << "\nOutput files:\n[";
    for (size_t i = 0; i < output_files.size(); ++i) {
        std::cout << "\"" << output_files[i] << "\"";
        if (i + 1 < output_files.size()) std::cout << ",";
    }
    std::cout << "]\n";
    
    // Summary
    std::cout << "\nSummary:\n";
    std::cout << "- Total configurations tested: " << tested_configurations.size() << "\n";
    std::cout << "- Output files generated: " << output_files.size() << "\n";
    
    if (!optimal_config.empty()) {
        std::cout << "- Optimal configuration: ";
        for (size_t i = 0; i < optimal_config.size(); ++i) {
            std::cout << "{";
            for (size_t j = 0; j < optimal_config[i].size(); ++j) {
                std::cout << optimal_config[i][j];
                if (j + 1 < optimal_config[i].size()) std::cout << ",";
            }
            std::cout << "}";
            if (i + 1 < optimal_config.size()) std::cout << ",";
        }
        std::cout << "\n";
    }
    
    return 0;
}