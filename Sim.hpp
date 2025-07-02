    #ifndef _SIM_HPP
    #define _SIM_HPP
    #include "Move.hpp"
#include "globals.hpp"
#include <unordered_map>
#include <vector>
#include <functional>
#include <random>
#include <memory>

    // Custom hash function for std::vector<int>


    class Sim {

    public:
        Sim();  // Declaration only
        Grid grid{40};  // Fixed: Use uniform initialization
        SimConditions conditions;
        SimVariables variables;

        float system_energy;

        // Random number generation
        std::random_device rd;
        std::mt19937 gen;

        std::unordered_map<int,int> expansion_candidates_primary;
        std::unordered_map<int,int> expansion_candidates_secondary;
        std::vector<expansion_candidate> primary_candidate;
        std::unordered_set<int> visited;
        Swap swap;
        std::vector<std::vector<expansion_candidate>> expansion_history;

        // Energy and acceptance tracking
        std::vector<float> energy_history;
        std::vector<int> total_acceptances_history;  // Changed to track rolling total
        int total_acceptances;
        // New vectors to store additional data
        std::vector<SwapType> swap_type_history;
        std::vector<float> forward_prob_history;
        std::vector<float> reverse_prob_history;
        std::vector<int> num_atoms_swapped_history;
        std::vector<float> delta_e_history;

        // Vectors to store cluster sizes
        std::vector<std::vector<int>> boundary_sizes_history;  // Changed to store sizes only
        std::vector<std::vector<int>> interior_sizes_history;  // Changed to store sizes only

        // New vector to store Norm after each move
        std::vector<float> norm_history;

        // New vectors to store atom data for each iteration
        std::vector<std::vector<int>> possibly_changed_history;  // possibly_changed_vec from move
        std::vector<std::vector<int>> proposed_atoms_history;     // primary_atoms from move
        std::vector<std::vector<int>> accepted_atoms_history;     // primary_atoms from move but only if accepted
        std::vector<std::vector<int>> accepted_possibly_changed_history;  // possibly_changed_vec from move but only if accepted
        
        // Vectors to store grid configuration history as strings
        std::vector<std::string> grid_config_starting_history;    // Grid state before move is attempted
        std::vector<std::string> grid_config_potential_history;   // Grid state after move is applied (before accept/reject)

            // ------------------------------
    // New: allowed lattice sites
    // ------------------------------
    // Stores the indices of lattice sites that are currently active in the
    // simulation.  All energy calculations and initial cluster assignment
    // routines will be restricted to this subset.  By default this is the
    // full lattice (0..N*N-1).
    std::vector<int> allowed_sites;

    // ------------------------------
    // Child simulations for partitioning
    // ------------------------------
    std::vector<std::unique_ptr<Sim>> children_simulations;

        // Function declarations
        void cluster_conditions(int cluster_condition );
        void cluster_conditions_0();
        void cluster_conditions_1();
        void cluster_conditions_2();
        void cluster_conditions_3();
        void initialize(const SimConditions& sim_conditions);
        SwapType determine_swap_type(int atom1, int atom2);
        bool iterate_improved(bool debug = false);
        float swap_type_probability(Swap input_swap);
        float expansion_probability(SwapType);
        std::vector<expansion_candidate> evaluate_expansion(std::vector<expansion_candidate> primary_set, std::vector<expansion_candidate> secondary_set);
        float compute_total_energy();
        void accept_move(Move& move);
        void undo_move(Move& move);
        Swap calculate_reverse_swap(const Move& move);
        void save_data(const Move& move);  // New function to save move data
        void save_results_to_csv(const std::string& filename);
        std::string grid_config_to_string() const;  // Helper function to convert grid config to string
        void save_grid_configuration_to_csv(const std::string& filename);
        Swap build_swap(int primary_atom, int secondary_atom,bool debug);

        std::vector<Swap> expand_swap_dr_greater_structure(Grid& grid, Swap swap, int iteration, bool debug);

        ///Functionality for multiple initial swaps
        std::vector<Swap> select_random_swaps(int max_num_swaps, bool debug);
        std::vector<Swap> calculate_reverse_swaps(std::vector<Swap> input_swaps, bool debug);
        bool delayed_rejection_iteration(int max_num_swaps, bool debug);

         Swap select_initial_atoms_weighted(bool debug);
         void calculate_initial_selection_probs_weighted(bool debug);
         void reset_cluster_conditions();
         void copy_grid(Grid& grid_copy);
         void reapply_cluster_weights(float cluster_bias, int tier_size = 1);
         // Build clusters using a list of indices into GLOBAL_CLUSTER_PROPERTIES (includes inversions)
         void cluster_conditions_from_indices(const std::vector<int>& property_indices, bool catch_all = true);
         // Build clusters using rotation-only pattern set (no inversions)
         void cluster_conditions_from_indices_uninverted(const std::vector<int>& property_indices, int tier_size = 1, bool catch_all = true);
         // Build clusters from a matrix where each row represents a cluster (includes inversions)
         void cluster_conditions_from_indices_matrix(const std::vector<std::vector<int>>& property_matrix, int tier_size = 1);
         // Build clusters from a matrix using rotation-only patterns (no inversions)
         void cluster_conditions_from_indices_matrix_uninverted(const std::vector<std::vector<int>>& property_matrix, int tier_size = 2);
         
         // Partitioning and child simulation functions
         std::vector<Sim*> partition(bool use_inversions = true);
         float iterate_children(bool debug = false);
         void resynchronize();
    };



    #endif