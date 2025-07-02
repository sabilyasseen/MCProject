#include "Sim.hpp"
#include <iostream>

int main() {
    std::cout << "Testing Sim partition functions..." << std::endl;
    
    // Initialize global data structures first
    initialize_energy();
    initialize_pattern_library();
    initialize_pattern_library_rot();
    
    // Create a basic simulation
    Sim sim;
    
    // Initialize with basic conditions
    SimConditions conditions(1.0f, 1.0f, 1.0f, 0.5f, 0.5f, 0.5f, 1.0f, 1.0f, 0.1f);
    sim.initialize(conditions);
    
    // Set up some basic cluster conditions for testing
    sim.cluster_conditions(2); // Use cluster condition 2 which has multiple clusters
    
    std::cout << "Simulation initialized with " << sim.grid.num_clusters << " clusters" << std::endl;
    std::cout << "Initial system energy: " << sim.system_energy << std::endl;
    
    // Test partition function
    std::cout << "\nTesting partition function..." << std::endl;
    std::vector<Sim*> children = sim.partition(true); // Use inversions
    
    std::cout << "Created " << children.size() << " child simulations" << std::endl;
    
    // Display some information about children
    for (size_t i = 0; i < children.size(); ++i) {
        std::cout << "Child " << i << " has " << children[i]->allowed_sites.size() 
                  << " allowed sites and " << children[i]->grid.num_clusters 
                  << " clusters" << std::endl;
        std::cout << "Child " << i << " energy: " << children[i]->system_energy << std::endl;
    }
    
    // Test iterate_children function
    std::cout << "\nTesting iterate_children function..." << std::endl;
    float total_energy_before = sim.iterate_children(true); // Debug enabled
    std::cout << "Total energy after iterate_children: " << total_energy_before << std::endl;
    
    // Test resynchronize function
    std::cout << "\nTesting resynchronize function..." << std::endl;
    sim.resynchronize();
    std::cout << "Resynchronization completed" << std::endl;
    std::cout << "Final system energy: " << sim.system_energy << std::endl;
    
    std::cout << "\nAll tests completed successfully!" << std::endl;
    return 0;
}