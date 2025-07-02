#include "Cluster.hpp"
#include <algorithm>
#include <iostream>

Cluster::Cluster(int ID, int minsize) {
    if (ID < 0) {
        throw std::invalid_argument("Invalid cluster ID");
    }
    if (minsize <= 0) {
        throw std::invalid_argument("Invalid minimum size");
    }
    requires_refresh = true;
    property.min_size = minsize;
    cluster_id = ID;
    num_species=2;

    mutex_interior = new std::mutex[num_species];
    mutex_boundary = new std::mutex[num_species];
    
    // Initialize vectors of sets for each species
    interior_cells.resize(num_species);
    boundary_cells.resize(num_species);
    
    // Initialize vectors of vectors for each species
    interior_atoms_vector.resize(num_species);
    boundary_atoms_vector.resize(num_species);
    
    // Pre-allocate vectors
    for (int i = 0; i < num_species; i++) {
        interior_atoms_vector[i].reserve(800);
        boundary_atoms_vector[i].reserve(800);
    }
}

void Cluster::add_property(std::vector<int> config) {
    property.add_config(config);
}

bool Cluster::satisfies_properties(std::vector<int> config) {
    if (config.empty()) {
        throw std::invalid_argument("Empty configuration");
    }
    return property.adhere(config);
}

void Cluster::refresh_vectors() {
    // Convert sets to vectors for each species
   // for (int i = 0; i < num_species; i++) {
     //   interior_atoms_vector[i] = std::vector<int>(interior_cells[i].begin(), interior_cells[i].end());
     //   boundary_atoms_vector[i] = std::vector<int>(boundary_cells[i].begin(), boundary_cells[i].end());
    //}
    requires_refresh = false;
}

void Cluster::add_atom_interior(int atom, int species) {
    if (species < 0 || species >= num_species) {
        std::cout << "[DEBUG][Cluster] Invalid species ID: " << species << " num_species: " << num_species << std::endl;
        throw std::invalid_argument("Invalid species ID");
    }
    mutex_interior[species].lock();
    if(interior_cells[species].find(atom)==interior_cells[species].end()){
    interior_cells[species][atom]=interior_cells[species].size();
    interior_atoms_vector[species].push_back(atom);
    }
    mutex_interior[species].unlock();
}

void Cluster::add_atom_boundary(int atom, int species) {
    if (species < 0 || species >= num_species) {
        std::cout << "[DEBUG][Cluster] Invalid species ID: " << species << " num_species: " << num_species << std::endl;
        throw std::invalid_argument("Invalid species ID: " + std::to_string(species));
    }
    mutex_boundary[species].lock();
    if(boundary_cells[species].find(atom)==boundary_cells[species].end()){
    boundary_cells[species][atom]=boundary_cells[species].size();
    boundary_atoms_vector[species].push_back(atom);
    }
    mutex_boundary[species].unlock();
}

void Cluster::remove_atom_interior(int atom, int species) {
    if (species < 0 || species >= num_species) {
        std::cout << "[DEBUG][Cluster] Invalid species ID: " << species << " num_species: " << num_species << std::endl;
        throw std::invalid_argument("Invalid species ID: " + std::to_string(species));
    }
    mutex_interior[species].lock();
    auto it = interior_cells[species].find(atom);
    if(it != interior_cells[species].end()){
        int index = it->second;
        int last_atom = interior_atoms_vector[species].back();
        interior_atoms_vector[species][index] = last_atom;
        interior_cells[species][last_atom] = index;
        interior_atoms_vector[species].pop_back();    
        interior_cells[species].erase(it);
    }
    mutex_interior[species].unlock();
}

void Cluster::remove_atom_boundary(int atom, int species) {
    if (species < 0 || species >= num_species) {
        std::cout << "[DEBUG][Cluster] Invalid species ID: " << species << " num_species: " << num_species << std::endl;
        throw std::invalid_argument("Invalid species ID: " + std::to_string(species));
    }
    mutex_boundary[species].lock();
    auto it = boundary_cells[species].find(atom);
    if(it != boundary_cells[species].end()){
        int index = it->second;
        int last_atom = boundary_atoms_vector[species].back();
        boundary_atoms_vector[species][index] = last_atom;
        boundary_cells[species][last_atom] = index;
        boundary_atoms_vector[species].pop_back();    
        boundary_cells[species].erase(it);
    }
    mutex_boundary[species].unlock();
}