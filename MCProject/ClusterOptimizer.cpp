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
#include <thread>
#include <future>
#include <chrono>
#include <omp.h>

// =============================================================================
// CLUSTER OPTIMIZER: Successor to Simtest.cpp
// Goals:
// 1. System to check every possible set/subset of cluster properties and 
//    combinations at a given temperature to find fastest convergence
// 2. Implement restricted partitioning schemes breaking down main simulation
//    into sub-children
// =============================================================================

// Global hash function for cluster configurations
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

// Global containers for optimization
std::unordered_map<std::vector<std::vector<int>>, float, VectorVectorHash> tested_configurations;
std::mutex config_mutex;

// =============================================================================
// UTILITY FUNCTIONS
// =============================================================================

// Standardize cluster configuration matrix
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
    
    // Remove duplicate clusters
    for (size_t i = 0; i < result.size(); ++i) {
        for (size_t j = i + 1; j < result.size(); ) {
            if (result[i] == result[j]) {
                result.erase(result.begin() + j);
            } else {
                ++j;
            }
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

// Calculate convergence rate based on energy stabilization
float calculate_convergence_rate(const std::vector<float>& energy_history) {
    if (energy_history.size() < 10) return 0.0f;
    
    size_t n = energy_history.size();
    size_t window = std::min(n / 4, size_t(100)); // Use last 25% or 100 steps
    
    float final_energy = energy_history.back();
    float variance = 0.0f;
    
    for (size_t i = n - window; i < n; ++i) {
        float diff = energy_history[i] - final_energy;
        variance += diff * diff;
    }
    
    variance /= window;
    return 1.0f / (1.0f + variance); // Higher value = faster convergence
}

// =============================================================================
// CLUSTER PROPERTY OPTIMIZATION
// =============================================================================

struct OptimizationResult {
    std::vector<std::vector<int>> configuration;
    float final_energy;
    float convergence_rate;
    int iterations_to_convergence;
    std::string filename;
    
    OptimizationResult() : final_energy(0.0f), convergence_rate(0.0f), iterations_to_convergence(0) {}
};

class ClusterPropertyOptimizer {
private:
    std::string base_conditions;
    float temperature;
    int max_iterations;
    int min_clusters;
    int max_clusters;
    int min_properties_per_cluster;
    int max_properties_per_cluster;
    bool use_rotation_only;
    
public:
    ClusterPropertyOptimizer(const std::string& conditions, float temp, int max_iter = 1000,
                           int min_clust = 1, int max_clust = 4, 
                           int min_prop = 1, int max_prop = 3,
                           bool rotation_only = false) 
        : base_conditions(conditions), temperature(temp), max_iterations(max_iter),
          min_clusters(min_clust), max_clusters(max_clust),
          min_properties_per_cluster(min_prop), max_properties_per_cluster(max_prop),
          use_rotation_only(rotation_only) {}
    
    // Test a specific cluster configuration
    OptimizationResult test_configuration(const std::vector<std::vector<int>>& config) {
        OptimizationResult result;
        result.configuration = config;
        
        // Check if already tested
        std::vector<std::vector<int>> std_config = standardize_configuration(config);
        {
            std::lock_guard<std::mutex> lock(config_mutex);
            auto it = tested_configurations.find(std_config);
            if (it != tested_configurations.end()) {
                result.final_energy = it->second;
                return result;
            }
        }
        
        // Setup and run simulation
        Sim sim;
        sim.conditions = SimConditions(SIM_CONDITIONS.at(base_conditions));
        sim.cluster_conditions(2);
        
        // Apply cluster configuration
        if (use_rotation_only) {
            sim.cluster_conditions_from_indices_matrix_uninverted(config, 2);
        } else {
            sim.cluster_conditions_from_indices_matrix(config, 1);
        }
        
        sim.initialize(sim.conditions);
        
        // Check if simulation is viable
        if (sim.variables.Norm == 0.0) {
            result.final_energy = 1e6f;
            result.convergence_rate = 0.0f;
            result.iterations_to_convergence = max_iterations;
            return result;
        }
        
        // Run simulation
        auto start_time = std::chrono::high_resolution_clock::now();
        
        for (int iter = 0; iter < max_iterations; ++iter) {
            sim.iterate_improved();
            
            // Check for early convergence
            if (iter > 100 && iter % 50 == 0) {
                float convergence = calculate_convergence_rate(sim.energy_history);
                if (convergence > 0.95f) {
                    result.iterations_to_convergence = iter;
                    break;
                }
            }
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        result.final_energy = sim.system_energy;
        result.convergence_rate = calculate_convergence_rate(sim.energy_history);
        if (result.iterations_to_convergence == 0) {
            result.iterations_to_convergence = max_iterations;
        }
        
        // Store result
        {
            std::lock_guard<std::mutex> lock(config_mutex);
            tested_configurations[std_config] = result.final_energy;
        }
        
        return result;
    }
    
    // Generate all possible cluster configurations
    std::vector<std::vector<std::vector<int>>> generate_all_configurations() {
        std::vector<std::vector<std::vector<int>>> all_configs;
        
        int num_properties = use_rotation_only ? 
            static_cast<int>(GLOBAL_PATTERN_VECTORS_ROT.size()) : 
            static_cast<int>(GLOBAL_PATTERN_VECTORS.size());
        
        // Generate configurations recursively
        std::function<void(std::vector<std::vector<int>>&, std::vector<int>&)> generate_configs = 
            [&](std::vector<std::vector<int>>& current_config, std::vector<int>& unused_properties) {
                // Add current configuration if valid
                if (current_config.size() >= min_clusters && current_config.size() <= max_clusters) {
                    bool valid = true;
                    for (const auto& cluster : current_config) {
                        if (cluster.size() < min_properties_per_cluster || 
                            cluster.size() > max_properties_per_cluster) {
                            valid = false;
                            break;
                        }
                    }
                    if (valid) {
                        all_configs.push_back(current_config);
                    }
                }
                
                // Stop if at maximum clusters
                if (current_config.size() >= max_clusters) return;
                
                // Try adding new clusters
                for (size_t i = 0; i < unused_properties.size(); ++i) {
                    std::vector<std::vector<int>> new_config = current_config;
                    new_config.push_back({unused_properties[i]});
                    
                    std::vector<int> new_unused = unused_properties;
                    new_unused.erase(new_unused.begin() + i);
                    
                    generate_configs(new_config, new_unused);
                }
                
                // Try adding properties to existing clusters
                for (size_t cluster_idx = 0; cluster_idx < current_config.size(); ++cluster_idx) {
                    if (current_config[cluster_idx].size() < max_properties_per_cluster) {
                        for (size_t prop_idx = 0; prop_idx < unused_properties.size(); ++prop_idx) {
                            std::vector<std::vector<int>> new_config = current_config;
                            new_config[cluster_idx].push_back(unused_properties[prop_idx]);
                            
                            std::vector<int> new_unused = unused_properties;
                            new_unused.erase(new_unused.begin() + prop_idx);
                            
                            generate_configs(new_config, new_unused);
                        }
                    }
                }
            };
        
        std::vector<int> all_properties;
        for (int i = 0; i < num_properties; ++i) {
            all_properties.push_back(i);
        }
        
        std::vector<std::vector<int>> empty_config;
        generate_configs(empty_config, all_properties);
        
        return all_configs;
    }
    
    // Find optimal cluster configuration
    OptimizationResult find_optimal_configuration() {
        std::cout << "[INFO] Generating all possible cluster configurations..." << std::endl;
        auto all_configs = generate_all_configurations();
        std::cout << "[INFO] Generated " << all_configs.size() << " configurations to test" << std::endl;
        
        OptimizationResult best_result;
        best_result.final_energy = std::numeric_limits<float>::max();
        
        // Test configurations in parallel
        #pragma omp parallel for
        for (size_t i = 0; i < all_configs.size(); ++i) {
            OptimizationResult result = test_configuration(all_configs[i]);
            
            // Update best result (thread-safe)
            #pragma omp critical
            {
                if (result.final_energy < best_result.final_energy ||
                    (result.final_energy == best_result.final_energy && 
                     result.convergence_rate > best_result.convergence_rate)) {
                    best_result = result;
                }
            }
        }
        
        std::cout << "[INFO] Optimal configuration found with energy: " << best_result.final_energy << std::endl;
        return best_result;
    }
};

// =============================================================================
// RESTRICTED PARTITIONING SCHEMES
// =============================================================================

class RestrictedPartitioningScheme {
private:
    std::string base_conditions;
    int grid_size;
    int num_partitions;
    
public:
    RestrictedPartitioningScheme(const std::string& conditions, int grid_sz = 40, int num_part = 4)
        : base_conditions(conditions), grid_size(grid_sz), num_partitions(num_part) {}
    
    // Generate spatial partitions of the grid
    std::vector<std::vector<int>> generate_spatial_partitions() {
        std::vector<std::vector<int>> partitions(num_partitions);
        
        int partition_size = grid_size / static_cast<int>(std::sqrt(num_partitions));
        int partitions_per_row = grid_size / partition_size;
        
        for (int i = 0; i < grid_size; ++i) {
            for (int j = 0; j < grid_size; ++j) {
                int site_index = i * grid_size + j;
                int partition_row = i / partition_size;
                int partition_col = j / partition_size;
                int partition_id = partition_row * partitions_per_row + partition_col;
                
                if (partition_id < num_partitions) {
                    partitions[partition_id].push_back(site_index);
                }
            }
        }
        
        return partitions;
    }
    
    // Run simulation on a specific partition
    OptimizationResult run_partition_simulation(const std::vector<int>& partition_sites, 
                                               const std::vector<std::vector<int>>& cluster_config,
                                               int iterations) {
        OptimizationResult result;
        
        Sim sim;
        sim.conditions = SimConditions(SIM_CONDITIONS.at(base_conditions));
        sim.cluster_conditions(2);
        
        // Set allowed sites to the partition
        sim.allowed_sites = partition_sites;
        
        // Apply cluster configuration
        sim.cluster_conditions_from_indices_matrix(cluster_config, 1);
        sim.initialize(sim.conditions);
        
        // Run simulation
        for (int iter = 0; iter < iterations; ++iter) {
            sim.iterate_improved();
        }
        
        result.final_energy = sim.system_energy;
        result.convergence_rate = calculate_convergence_rate(sim.energy_history);
        
        return result;
    }
    
    // Test cluster configuration using restricted partitioning
    std::vector<OptimizationResult> test_with_partitioning(
        const std::vector<std::vector<int>>& cluster_config, int iterations_per_partition) {
        
        auto partitions = generate_spatial_partitions();
        std::vector<OptimizationResult> results(num_partitions);
        
        // Run simulations on each partition in parallel
        #pragma omp parallel for
        for (int i = 0; i < num_partitions; ++i) {
            results[i] = run_partition_simulation(partitions[i], cluster_config, iterations_per_partition);
        }
        
        return results;
    }
    
    // Combine results from multiple partitions
    OptimizationResult combine_partition_results(const std::vector<OptimizationResult>& partition_results) {
        OptimizationResult combined;
        
        combined.final_energy = 0.0f;
        combined.convergence_rate = 0.0f;
        combined.iterations_to_convergence = 0;
        
        for (const auto& result : partition_results) {
            combined.final_energy += result.final_energy;
            combined.convergence_rate += result.convergence_rate;
            combined.iterations_to_convergence += result.iterations_to_convergence;
        }
        
        combined.final_energy /= partition_results.size();
        combined.convergence_rate /= partition_results.size();
        combined.iterations_to_convergence /= partition_results.size();
        
        return combined;
    }
};

// =============================================================================
// SIMULATION RUNNER AND FILE MANAGEMENT
// =============================================================================

class SimulationRunner {
private:
    std::string base_name;
    std::vector<std::string> output_files;
    
public:
    SimulationRunner(const std::string& name) : base_name(name) {}
    
    // Run standard simulation and save results
    std::string run_standard_simulation(const std::string& conditions, int iterations, 
                                      int cluster_type, const std::vector<std::vector<int>>& cluster_config) {
        
        // Generate filename with timestamp
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        std::stringstream ss;
        ss << base_name << "_" << conditions << "_cluster" << cluster_type 
           << "_iter" << iterations << "_" << time_t << ".csv";
        std::string filename = ss.str();
        
        // Setup and run simulation
        Sim sim;
        sim.conditions = SimConditions(SIM_CONDITIONS.at(conditions));
        sim.cluster_conditions(cluster_type);
        
        if (!cluster_config.empty()) {
            sim.cluster_conditions_from_indices_matrix(cluster_config, 1);
        }
        
        sim.initialize(sim.conditions);
        
        // Run simulation
        for (int iter = 0; iter < iterations; ++iter) {
            sim.iterate_improved();
        }
        
        // Save results
        sim.save_results_to_csv(filename);
        output_files.push_back(filename);
        
        return filename;
    }
    
    // Run delayed rejection simulation
    std::string run_delayed_rejection_simulation(const std::string& conditions, int iterations,
                                                int cluster_type, int max_num_swaps,
                                                const std::vector<std::vector<int>>& cluster_config) {
        
        // Generate filename
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        std::stringstream ss;
        ss << base_name << "_DR_" << conditions << "_cluster" << cluster_type 
           << "_swaps" << max_num_swaps << "_iter" << iterations << "_" << time_t << ".csv";
        std::string filename = ss.str();
        
        // Setup and run simulation
        Sim sim;
        sim.conditions = SimConditions(SIM_CONDITIONS.at(conditions));
        sim.cluster_conditions(cluster_type);
        
        if (!cluster_config.empty()) {
            sim.cluster_conditions_from_indices_matrix(cluster_config, 1);
        }
        
        sim.initialize(sim.conditions);
        
        // Run delayed rejection simulation
        for (int iter = 0; iter < iterations; ++iter) {
            sim.delayed_rejection_iteration(max_num_swaps, false);
        }
        
        // Save results
        sim.save_results_to_csv(filename);
        output_files.push_back(filename);
        
        return filename;
    }
    
    // Run partitioned simulation
    std::string run_partitioned_simulation(const std::string& conditions, int iterations,
                                         const std::vector<std::vector<int>>& cluster_config,
                                         int num_partitions) {
        
        // Generate filename
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        std::stringstream ss;
        ss << base_name << "_PART" << num_partitions << "_" << conditions 
           << "_iter" << iterations << "_" << time_t << ".csv";
        std::string filename = ss.str();
        
        // Setup partitioning scheme
        RestrictedPartitioningScheme partitioner(conditions, 40, num_partitions);
        
        // Run partitioned simulation
        auto partition_results = partitioner.test_with_partitioning(cluster_config, iterations / num_partitions);
        auto combined_result = partitioner.combine_partition_results(partition_results);
        
        // Save combined results (create a simple CSV with the combined metrics)
        std::ofstream file(filename);
        if (file.is_open()) {
            file << "partition_id,final_energy,convergence_rate,iterations_to_convergence\n";
            for (size_t i = 0; i < partition_results.size(); ++i) {
                file << i << "," << partition_results[i].final_energy << ","
                     << partition_results[i].convergence_rate << ","
                     << partition_results[i].iterations_to_convergence << "\n";
            }
            file << "combined," << combined_result.final_energy << ","
                 << combined_result.convergence_rate << ","
                 << combined_result.iterations_to_convergence << "\n";
            file.close();
        }
        
        output_files.push_back(filename);
        return filename;
    }
    
    // Get all output filenames
    std::vector<std::string> get_output_files() const {
        return output_files;
    }
};

// =============================================================================
// MAIN PROGRAM
// =============================================================================

int main() {
    // Initialize global libraries
    initialize_energy();
    initialize_pattern_library();
    initialize_pattern_library_rot();
    
    std::cout << "=============================================================================\n";
    std::cout << "                            CLUSTER OPTIMIZER\n";
    std::cout << "=============================================================================\n";
    
    // Simulation parameters
    const std::string base_conditions = "BJ=-0.44_DR";
    const float temperature = -0.44f;
    const int base_iterations = 1000;
    
    // Create simulation runner
    SimulationRunner runner("ClusterOptimizer");
    
    // =============================================================================
    // GOAL 1: FIND OPTIMAL CLUSTER PROPERTIES FOR FASTEST CONVERGENCE
    // =============================================================================
    
    std::cout << "\n[GOAL 1] Finding optimal cluster properties for fastest convergence...\n";
    
    // Test different cluster optimization configurations
    std::vector<std::tuple<int, int, int, int, bool>> optimization_configs = {
        {2, 3, 1, 2, false},  // 2-3 clusters, 1-2 properties each, with inversions
        {2, 4, 1, 3, false},  // 2-4 clusters, 1-3 properties each, with inversions
        {3, 4, 2, 2, true},   // 3-4 clusters, 2 properties each, rotation only
        {2, 3, 2, 3, true}    // 2-3 clusters, 2-3 properties each, rotation only
    };
    
    std::vector<OptimizationResult> optimization_results;
    
    for (const auto& config : optimization_configs) {
        int min_clusters, max_clusters, min_props, max_props;
        bool rotation_only;
        std::tie(min_clusters, max_clusters, min_props, max_props, rotation_only) = config;
        
        std::cout << "\n[INFO] Testing configuration: " << min_clusters << "-" << max_clusters 
                  << " clusters, " << min_props << "-" << max_props << " properties"
                  << (rotation_only ? " (rotation only)" : " (with inversions)") << "\n";
        
        ClusterPropertyOptimizer optimizer(base_conditions, temperature, base_iterations,
                                         min_clusters, max_clusters, min_props, max_props, rotation_only);
        
        OptimizationResult result = optimizer.find_optimal_configuration();
        optimization_results.push_back(result);
        
        // Run full simulation with optimal configuration
        std::string filename = runner.run_standard_simulation(base_conditions, base_iterations, 2, result.configuration);
        result.filename = filename;
        
        std::cout << "[INFO] Optimal configuration saved to: " << filename << "\n";
    }
    
    // =============================================================================
    // GOAL 2: RESTRICTED PARTITIONING SCHEMES
    // =============================================================================
    
    std::cout << "\n[GOAL 2] Testing restricted partitioning schemes...\n";
    
    // Use the best configuration from optimization
    std::vector<std::vector<int>> best_config;
    if (!optimization_results.empty()) {
        auto best_it = std::min_element(optimization_results.begin(), optimization_results.end(),
            [](const OptimizationResult& a, const OptimizationResult& b) {
                return a.final_energy < b.final_energy;
            });
        best_config = best_it->configuration;
        std::cout << "[INFO] Using best configuration from optimization with energy: " << best_it->final_energy << "\n";
    } else {
        // Default configuration if optimization failed
        best_config = {{0, 1}, {2, 3}};
        std::cout << "[INFO] Using default configuration: {{0,1}, {2,3}}\n";
    }
    
    // Test different partitioning schemes
    std::vector<int> partition_counts = {4, 9, 16};
    
    for (int num_partitions : partition_counts) {
        std::cout << "\n[INFO] Testing " << num_partitions << "-partition scheme...\n";
        
        std::string filename = runner.run_partitioned_simulation(base_conditions, base_iterations, 
                                                               best_config, num_partitions);
        
        std::cout << "[INFO] Partitioned simulation (" << num_partitions << " partitions) saved to: " << filename << "\n";
    }
    
    // =============================================================================
    // COMPARATIVE BENCHMARKS
    // =============================================================================
    
    std::cout << "\n[BENCHMARK] Running comparative benchmarks...\n";
    
    // Standard Kawasaki dynamics
    std::string kawasaki_file = runner.run_standard_simulation("BJ=-0.44_Kawasaki", base_iterations, 0, {});
    std::cout << "[INFO] Kawasaki benchmark saved to: " << kawasaki_file << "\n";
    
    // Delayed rejection with different swap counts
    for (int max_swaps : {1, 2, 3}) {
        std::string dr_file = runner.run_delayed_rejection_simulation(base_conditions, base_iterations, 
                                                                    1, max_swaps, best_config);
        std::cout << "[INFO] Delayed rejection (" << max_swaps << " swaps) saved to: " << dr_file << "\n";
    }
    
    // =============================================================================
    // RESULTS SUMMARY
    // =============================================================================
    
    std::cout << "\n=============================================================================\n";
    std::cout << "                            RESULTS SUMMARY\n";
    std::cout << "=============================================================================\n";
    
    // Output all filenames in JSON format
    auto all_files = runner.get_output_files();
    
    std::cout << "\nAll output files in JSON format:\n";
    std::cout << "[";
    for (size_t i = 0; i < all_files.size(); ++i) {
        std::cout << "\"" << all_files[i] << "\"";
        if (i + 1 < all_files.size()) std::cout << ",";
    }
    std::cout << "]\n";
    
    // Summary statistics
    std::cout << "\nOptimization Summary:\n";
    std::cout << "- Total configurations tested: " << tested_configurations.size() << "\n";
    std::cout << "- Optimization runs completed: " << optimization_results.size() << "\n";
    std::cout << "- Partitioning schemes tested: " << partition_counts.size() << "\n";
    std::cout << "- Total output files: " << all_files.size() << "\n";
    
    if (!optimization_results.empty()) {
        auto best_result = *std::min_element(optimization_results.begin(), optimization_results.end(),
            [](const OptimizationResult& a, const OptimizationResult& b) {
                return a.final_energy < b.final_energy;
            });
        
        std::cout << "\nBest Configuration Found:\n";
        std::cout << "- Final energy: " << best_result.final_energy << "\n";
        std::cout << "- Convergence rate: " << best_result.convergence_rate << "\n";
        std::cout << "- Iterations to convergence: " << best_result.iterations_to_convergence << "\n";
        std::cout << "- Configuration: ";
        for (size_t i = 0; i < best_result.configuration.size(); ++i) {
            std::cout << "{";
            for (size_t j = 0; j < best_result.configuration[i].size(); ++j) {
                std::cout << best_result.configuration[i][j];
                if (j + 1 < best_result.configuration[i].size()) std::cout << ",";
            }
            std::cout << "}";
            if (i + 1 < best_result.configuration.size()) std::cout << ",";
        }
        std::cout << "\n";
    }
    
    std::cout << "\n=============================================================================\n";
    std::cout << "                        CLUSTER OPTIMIZER COMPLETE\n";
    std::cout << "=============================================================================\n";
    
    return 0;
}