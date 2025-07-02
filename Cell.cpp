#include "Cell.hpp"

#include <vector>
#include <iostream>
#include <unordered_map>
#include <functional>

// Constructors & Destructor
Cell::Cell() : index(-1), cell_species(0), cell_type(CellType::CELL_ANY), cell_cluster_ID(-1) {}

Cell::Cell(int xcord, int ycord, int N, int specie) {
    if (xcord < 0 || ycord < 0 || N <= 0) {
        throw std::invalid_argument("Invalid grid coordinates or size");
    }
    cell_cord.cordx = xcord;
    cell_cord.cordy = ycord;
    index = ycord * N + xcord;
    cell_species = specie;
    cell_cluster_ID = -1;
    cell_type = CellType::CELL_ANY;
}

Cell::~Cell() {}

// Setters and Getters
void Cell::set_atom(int xcord, int ycord, int N, int specie) {
    cell_cord.cordx = xcord;
    cell_cord.cordy = ycord;
    index = ycord * N + xcord;
    cell_species = specie;
    prior_species = specie;
    cell_type = CellType::CELL_ANY;
    prior_type = CellType::CELL_ANY;
    cell_cluster_ID = -1;
    prior_cluster_ID = -1;
    min_size = 0;
    prior_min_size = 0;
}

void Cell::revert_cell() {
    // Revert all cell properties back to their previous state
    cell_species = prior_species;
    cell_type = prior_type;
    cell_cluster_ID = prior_cluster_ID;     
    min_size = prior_min_size;
    local_config = prior_local_config;
}

void Cell::accept_cell() {
    // Update previous state to match current state
    prior_species = cell_species;
    prior_type = cell_type;
    prior_cluster_ID = cell_cluster_ID;
    prior_min_size = min_size;
    prior_local_config = local_config;
}
