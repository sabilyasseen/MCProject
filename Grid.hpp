#ifndef _LATTICE_HPP
#define _LATTICE_HPP

#include "Cluster.hpp"
#include "Cell.hpp"
#include "ClusterProperty.hpp"
#include <vector>
#include <iostream>
#include <random>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <ctime>
#include <unordered_set>


// This 'Grid' class represents the lattice or the entire simulation.
// It sorts the atoms into different clusters and performs the Monte Carlo simulation on those clusters.

class Grid {
public:
    Grid(int N, bool minsize5_boundary_bool=false);
    ~Grid();
    int *total_boundary_atoms;
    int *total_interior_atoms;
    // Core grid properties
    int Size;
    int num_clusters;
    Cell* array;
    std::vector<Cluster> clusters;
    bool minsize5_boundary;
    // Lattice site filtering sets
    std::vector<std::unordered_set<int>> interior_lattice_sites;
    std::vector<std::unordered_set<int>> boundary_lattice_sites;
    // Random number generation
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_int_distribution<> species_dist;
    std::vector<Cluster> restricted_clusters; // List of forbidden clusters

    // Grid operations
    void add_cluster(int minsize);
    void assign_initial_clusters();
    void assign_initial_clusters(const std::vector<int>& allowed_indices);
    
    void apply_cluster_change(std::vector<int> atom_list, int filter_index = -1);
    void undo_cluster_changes(std::vector<int> atom_list, int filter_index = -1);
    void change_specie(int atom);
    void refresh_local_config(int atom);
    
    // Unified atom management functions
    void add_atom(int atom, int cluster_index, CellType type);
    void remove_atom(int atom);
    
    bool reclassify_cluster_status(std::vector<int> atom_list, int filter_index = -1);
    void reclassify_boundary_status_minsize5(std::vector<int> atom_list, int filter_index = -1);
    void reclassify_boundary_status_minsize2(std::vector<int> atom_list, int filter_index = -1);
    void reclassify_boundary_status(std::vector<int> atom_list, int filter_index = -1);
    void adjust_total_boundary(int specie, int increment);
    void adjust_total_interior(int specie, int increment);
    bool any_neighbors_property(int atom,int cluster_id);
    bool all_neighbors_property(int atom,int cluster_id);
    bool any_neighbors_id(int atom,int cluster_id);
    bool all_neighbors_id(int atom,int cluster_id);
    // Add a restricted cluster property (to be prohibited during reclassification)
    void add_restricted_cluster(int minsize);
};

#endif
