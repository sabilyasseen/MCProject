#ifndef _MOVE_HPP
#define _MOVE_HPP
#include "Grid.hpp"
#include "globals.hpp"
#include <unordered_set>
#include <random>
#include <omp.h> 
#include <mutex>
#include <algorithm>

class Move {
public: 
    Move(float j, float beta, float tolerance) : J(j), Beta(beta), tolerance(tolerance), gen(rd()) {}
    int num_changed=0;
    int tolerant_forward_atoms_accepted=0;
    std::vector<int> tolerant_forward_atoms_accepted_tiers;
    int tolerant_forward_atoms_rejected=0;
    std::vector<int> tolerant_forward_atoms_rejected_tiers;
    std::vector<std::vector<Swap>> swap_history;
    std::vector<std::vector<Swap>> reverse_swap_history;
    int tier=0;
    int max_num_atoms=0;
    bool cluster_distinction=true;
    bool celltype_distinction=true;
    // Core simulation parameters
    float J;
    float Beta;  // Inverse temperature parameter
    float tolerance;
    // Forward and reverse swap information
    Swap forward;
    Swap reverse;
    float DeltaE;
    std::vector<float> DeltaE_teirs;
    bool any_new=false;
    // Expansion history and probabilities
    std::vector<std::vector<expansion_candidate>> primary_expansion_history;
    std::vector<std::vector<expansion_candidate>> secondary_expansion_history;
    std::vector<float> forward_expansion_teir_probabilities;
    std::vector<float> reverse_expansion_teir_probabilities;
    float final_forward_probability = 1.0f;
    float final_reverse_probability = 1.0f;
    float reverse_expansion_prob = 1.0f;
    // Selected atoms for the move
    std::vector<std::vector<int>> primary_atoms;
    std::vector<std::vector<int>> secondary_atoms;

    // Tracking visited atoms
    std::unordered_set<int> possibly_changed;
    std::vector<int> possibly_changed_vec;
    std::unordered_set<int> visited;

    // Random number generation
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<float> dist{0.0f, 1.0f};

    // Core functionality
    bool evaluate_deterministic_candidates(Grid& grid, bool debug);
    bool evaluate_candidates_dr(Grid& grid, bool debug);
    bool calculate_cluster_changes(Grid& grid, bool debug);
    void possibly_changed_atoms(int atom, Grid& grid);
    int get_primary_expansion_index_candidates(Grid& grid, bool debug);
    int get_secondary_expansion_index_candidates(Grid& grid, bool debug);
    std::vector<expansion_candidate> deterministic_expand(Grid& grid, bool debug);
    bool deterministic_build_from_swap(Swap input_swap, Grid& grid, bool debug);
    void deterministic_build_from_swap_delayed_rejection(std::vector<Swap> input_swaps, Grid& grid, bool debug);
    float calculate_reverse_expansion_prob(std::vector<Swap> input_swaps, Grid& grid, bool debug);
    void evaluate_reverse_expansion_deterministic(Grid& grid, bool debug);
    float get_final_acceptance_deterministic(Grid& grid, Swap reverse_swap, bool debug);
    bool is_valid_swap(Grid& grid, int primary_candidate, int secondary_candidate, Swap reference_swap, bool cluster_match, bool celltype_match, bool debug);
    bool reclassify_singles(Grid& grid, std::vector<int> main_atoms);
    float calculate_reverse_delayed_rejection(std::vector<Swap> reverse_swaps, Grid& grid, bool debug);
    bool expand_swap_dr(Grid& grid, bool debug);
    bool expand_swap_dr_greater_structure(Grid& grid, bool debug);
};

#endif

extern std::unordered_map<std::vector<int>, float, VectorHash> Energy;
extern std::unordered_map<std::vector<int>, float, VectorHash> delta_energy;
