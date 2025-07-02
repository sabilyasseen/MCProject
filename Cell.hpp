#ifndef _CELL_HPP
#define _CELL_HPP
#include "globals.hpp"
#include <vector>
struct Cell_Info{
    int index;
    int cell_species;
    CellType cell_type;
    int cell_cluster_ID;
    int min_size;
};
class Cell {
public:
    // Constructors & Destructor
    Cell();
    Cell(int xcord, int ycord, int N, int specie);
    ~Cell();

    // --- Data Members ---
    int index;
    cord cell_cord;
    int cell_species;
    int prior_species;
    CellType cell_type;
    CellType prior_type;

    int cell_cluster_ID;
    int prior_cluster_ID;

    int min_size;
    int prior_min_size;

    std::vector<int> local_config;  // [center, up, right, down, left]
    std::vector<int> prior_local_config;
    std::vector<int> neighbor_indexes;
    std::vector<int> expansion_indexes;

    // --- Setters and Getters ---
    void set_atom(int xcord, int ycord, int N, int specie);
    void revert_cell();
    void accept_cell();
    
    // --- New Member Functions ---
};

#endif