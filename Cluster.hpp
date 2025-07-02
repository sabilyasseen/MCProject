#ifndef _CLUSTER_HPP
#define _CLUSTER_HPP
#include <omp.h>
#include <mutex>
#include "Cell.hpp"
#include "ClusterProperty.hpp"
#include <unordered_set>
#include <vector>
#include <stdexcept>

class Cluster {
public:
    Cluster(int ID, int minsize);
    std::mutex *mutex_boundary;  // Array of mutexes for boundary atoms of each species
    std::mutex *mutex_interior;  // Array of mutexes for interior atoms of each species
    // Core cluster properties
    int cluster_id;
    int cluster_size;
    bool requires_refresh;
    int num_species=2;
    float cluster_weight=1.0f;
    float cluster_prob=0.0f;
    // Cluster data structures
    ClusterProperty property;
    std::vector<std::unordered_map<int,int>> boundary_cells;  // Vector of maps for each species
    std::vector<std::unordered_map<int,int>> interior_cells;  // Vector of maps for each species

    // Optimized vectors for quick access
    std::vector<std::vector<int>> interior_atoms_vector;  // Vector of vectors for each species
    std::vector<std::vector<int>> boundary_atoms_vector;  // Vector of vectors for each species


    // Cell management functions

    void add_property(std::vector<int> config);

    // Property checking
    bool satisfies_properties(std::vector<int> config);
    void add_atom_interior(int atom, int species);
    void add_atom_boundary(int atom, int species);
    void remove_atom_interior(int atom, int species);
    void remove_atom_boundary(int atom, int species);
    void refresh_vectors();

};

#endif
