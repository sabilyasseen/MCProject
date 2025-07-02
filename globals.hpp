#ifndef Globals_HPP
#define Globals_HPP
#include <unordered_map>
#include <vector>
#include <string>
#include <random>
#include <stdexcept>
#include <atomic>

enum CellType { CELL_ANY, BOUNDARY, INTERIOR };
enum Species { SPECIES_ANY, ONE, TWO };
enum direction {UP,DOWN};

struct cord {
    int cordx;
    int cordy;
};

struct VectorEqual {
    bool operator()(const std::vector<int>& a, const std::vector<int>& b) const {
        if (a.size() != b.size()) return false;
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] != b[i]) return false;
        }
        return true;
    }
};

struct VectorHash {
    std::size_t operator()(const std::vector<int>& v) const {
        std::size_t seed = 0;
        for (int x : v) {
            seed ^= std::hash<int>{}(x) + 0x9e3779b9 + (seed << 6) + (seed >> 2); // boost::hash_combine
        }
        return seed;
    }
};
enum SwapType { BB, BI, II };

struct Swap {
    int primary_atom;
    int primary_atom_species;
    int primary_atom_clusterID;
    CellType primary_atom_type;

    int secondary_atom;
    int secondary_atom_species;
    int secondary_atom_clusterID;
    CellType secondary_atom_type;

    SwapType this_swap_type;
    bool same_cluster;
    float initial_selection_prob;
    float expansion_prob;
};

struct expansion_candidate{
    
    int position=-1;
    int num_neighbors=0;
    float prob_denial=1.0f;
    float base_prob_denial=1.0f;
    bool acceptance=false;

};

struct SimConditions {
    float prob_bb_base;      // Base probability for boundary-boundary swaps
    float prob_bi_base;      // Base probability for boundary-interior swaps
    float prob_ii_base;      // Base probability for interior-interior swaps
    float prob_same_cluster_base; // Base probability for same cluster swaps

    float prob_bb_expansion; // Expansion probability for boundary-boundary swaps
    float prob_bi_expansion; // Expansion probability for boundary-interior swaps
    float prob_ii_expansion; // Expansion probability for interior-interior swaps
    float J;                 // Coupling constant
    float Beta;              // Inverse temperature parameter
    bool minsize5_boundary_setting;
    float tolerance;         // Tolerance parameter for expansion probabilities
    float cluster_bias;

    // Default constructor
    SimConditions() 
        : prob_bb_base(1.0f/3.0f), 
          prob_bi_base(1.0f/3.0f), 
          prob_ii_base(1.0f/3.0f),
          prob_same_cluster_base(0.5f),
          prob_bb_expansion(0.0f), 
          prob_bi_expansion(0.0f), 
          prob_ii_expansion(0.0f),
          J(10.0f), 
          Beta(0.01f),
          tolerance(0.1f),
          cluster_bias(0.0f),
          minsize5_boundary_setting(false) {} 

    // Parameterized constructor
    SimConditions(float bb_base, float bi_base, float ii_base, float same_cluster_base,
                 float bb_exp, float bi_exp, float ii_exp,
                 float j, float beta = 1.0f, float tol = 0.1f, float cluster_bias = 0.0f, bool minsize5_boundary_setting_bool = false)
        : prob_bb_base(bb_base), prob_bi_base(bi_base), prob_ii_base(ii_base), prob_same_cluster_base(same_cluster_base),
          prob_bb_expansion(bb_exp), prob_bi_expansion(bi_exp), prob_ii_expansion(ii_exp),
          J(j), Beta(beta), tolerance(tol), cluster_bias(cluster_bias), minsize5_boundary_setting(false) {}
};

struct SimVariables {
    // Current state variables
    int total_boundary_species_0;
    int total_boundary_species_1;
    int total_interior_species_0;
    int total_interior_species_1;

    float num_bb;           // Number of boundary-boundary pairs
    float num_bb_same_cluster;
    float num_bb_diff_cluster;
    float num_bi;
    float num_bi_same_cluster;
    float num_bi_diff_cluster;
    float num_ii;           // Number of interior-interior pairs
    float num_ii_same_cluster;
    float num_ii_diff_cluster;

    float prob_bb_initial;  // Initial probability for BB swaps
    float prob_bi_initial;  // Initial probability for BI swaps
    float prob_ii_initial;  // Initial probability for II swaps
    float prob_bb_same_cluster_initial;
    float prob_bi_same_cluster_initial;
    float prob_ii_same_cluster_initial;
    float prob_bb_diff_cluster_initial;
    float prob_bi_diff_cluster_initial;
    float prob_ii_diff_cluster_initial;
    double Norm_BB;          // Normalization factor for BB swaps
    double Norm_BI;          // Normalization factor for BI swaps
    double Norm_II;          // Normalization factor for II swaps
    double Norm;             // Overall normalization factor
    double system_energy;    // Current system energy
    
    // Backup variables
    float backup_num_bb;
    float backup_num_bi;
    float backup_num_ii;
    float backup_prob_bb_initial;
    float backup_prob_bi_initial;
    float backup_prob_ii_initial;
    float backup_prob_bb_same_cluster_initial;
    float backup_prob_bi_same_cluster_initial;
    float backup_prob_ii_same_cluster_initial;
    float backup_prob_bb_diff_cluster_initial;
    float backup_prob_bi_diff_cluster_initial;
    float backup_prob_ii_diff_cluster_initial;
    double backup_Norm_BB;
    double backup_Norm_BI;
    double backup_Norm_II;
    double backup_Norm;
    float backup_num_bb_same_cluster;
    float backup_num_bi_same_cluster;
    float backup_num_ii_same_cluster;
    float backup_num_bb_diff_cluster;
    float backup_num_bi_diff_cluster;
    float backup_num_ii_diff_cluster;

    // Default constructor
    SimVariables()
        : num_bb(0.0f), num_bb_same_cluster(0.0f), num_bb_diff_cluster(0.0f),
          num_bi(0.0f), num_bi_same_cluster(0.0f), num_bi_diff_cluster(0.0f), 
          num_ii(0.0f), num_ii_same_cluster(0.0f), num_ii_diff_cluster(0.0f),
          prob_bb_initial(0.0f), prob_bi_initial(0.0f), prob_ii_initial(0.0f),
          prob_bb_same_cluster_initial(0.0f), prob_bi_same_cluster_initial(0.0f), prob_ii_same_cluster_initial(0.0f),
          prob_bb_diff_cluster_initial(0.0f), prob_bi_diff_cluster_initial(0.0f), prob_ii_diff_cluster_initial(0.0f),
          Norm_BB(0.0f), Norm_BI(0.0f), Norm_II(0.0f), Norm(0.0f),
          system_energy(0.0f),
          backup_num_bb(0.0f), backup_num_bi(0.0f), backup_num_ii(0.0f),
          backup_prob_bb_initial(0.0f), backup_prob_bi_initial(0.0f), backup_prob_ii_initial(0.0f),
          backup_prob_bb_same_cluster_initial(0.0f), backup_prob_bi_same_cluster_initial(0.0f), backup_prob_ii_same_cluster_initial(0.0f),
          backup_prob_bb_diff_cluster_initial(0.0f), backup_prob_bi_diff_cluster_initial(0.0f), backup_prob_ii_diff_cluster_initial(0.0f),
          backup_Norm_BB(0.0f), backup_Norm_BI(0.0f), backup_Norm_II(0.0f), backup_Norm(0.0f),
          backup_num_bb_same_cluster(0.0f), backup_num_bi_same_cluster(0.0f), backup_num_ii_same_cluster(0.0f),
          backup_num_bb_diff_cluster(0.0f), backup_num_bi_diff_cluster(0.0f), backup_num_ii_diff_cluster(0.0f) {}

    // Method to backup current state
    void backup() {
        backup_num_bb = num_bb;
        backup_num_bi = num_bi;
        backup_num_ii = num_ii;
        backup_prob_bb_initial = prob_bb_initial;
        backup_prob_bi_initial = prob_bi_initial;
        backup_prob_ii_initial = prob_ii_initial;
        backup_prob_bb_same_cluster_initial = prob_bb_same_cluster_initial;

        backup_prob_bi_same_cluster_initial = prob_bi_same_cluster_initial;
        backup_prob_ii_same_cluster_initial = prob_ii_same_cluster_initial;

        backup_prob_bb_diff_cluster_initial = prob_bb_diff_cluster_initial;
        backup_prob_bi_diff_cluster_initial = prob_bi_diff_cluster_initial;
        backup_prob_ii_diff_cluster_initial = prob_ii_diff_cluster_initial;
        backup_Norm_BB = Norm_BB;
        backup_Norm_BI = Norm_BI;
        backup_Norm_II = Norm_II;
        backup_Norm = Norm;
        backup_num_bb_same_cluster = num_bb_same_cluster;
        backup_num_bi_same_cluster = num_bi_same_cluster;
        backup_num_ii_same_cluster = num_ii_same_cluster;
        backup_num_bb_diff_cluster = num_bb_diff_cluster;
        backup_num_bi_diff_cluster = num_bi_diff_cluster;
        backup_num_ii_diff_cluster = num_ii_diff_cluster;
    }

    // Method to restore from backup
    void restore() {
        num_bb = backup_num_bb;
        num_bi = backup_num_bi;
        num_ii = backup_num_ii;
        prob_bb_initial = backup_prob_bb_initial;
        prob_bi_initial = backup_prob_bi_initial;
        prob_ii_initial = backup_prob_ii_initial;
        prob_bb_same_cluster_initial = backup_prob_bb_same_cluster_initial;
        prob_bi_same_cluster_initial = backup_prob_bi_same_cluster_initial;
        prob_ii_same_cluster_initial = backup_prob_ii_same_cluster_initial;
        prob_bb_diff_cluster_initial = backup_prob_bb_diff_cluster_initial;
        prob_bi_diff_cluster_initial = backup_prob_bi_diff_cluster_initial;
        prob_ii_diff_cluster_initial = backup_prob_ii_diff_cluster_initial;
        Norm_BB = backup_Norm_BB;
        Norm_BI = backup_Norm_BI;
        Norm_II = backup_Norm_II;
        Norm = backup_Norm;
        num_bb_same_cluster = backup_num_bb_same_cluster;
        num_bi_same_cluster = backup_num_bi_same_cluster;
        num_ii_same_cluster = backup_num_ii_same_cluster;
        num_bb_diff_cluster = backup_num_bb_diff_cluster;
        num_bi_diff_cluster = backup_num_bi_diff_cluster;
        num_ii_diff_cluster = backup_num_ii_diff_cluster;
    }
};



extern std::unordered_map<std::vector<int>, float, VectorHash> Energy;
extern std::unordered_map<std::vector<int>, float, VectorHash> delta_energy;
extern const std::unordered_map<std::string, SimConditions> SIM_CONDITIONS;

// Global repository of basic 5-length binary patterns
extern std::vector<std::vector<int>> GLOBAL_PATTERN_VECTORS;
// Lookup from pattern -> index within GLOBAL_PATTERN_VECTORS
extern std::unordered_map<std::vector<int>, int, VectorHash, VectorEqual> PATTERN_LOOKUP;

extern std::vector<std::vector<int>> GLOBAL_PATTERN_VECTORS_ROT; // Patterns without inversions (rotations only)
extern std::unordered_map<std::vector<int>, int, VectorHash, VectorEqual> PATTERN_LOOKUP_ROT; // lookup for rotation-only patterns

void initialize_energy();
void initialize_delta_energy();
void initialize_pattern_library();
void initialize_pattern_library_rot();

struct SwapSelector {
    // Current probabilities and counts
    float prob_bb;
    float prob_bi;
    float prob_ii;
    int count_bb;
    int count_bi;
    int count_ii;

    // Random number generation
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<float> dist;

    SwapSelector() : gen(rd()), dist(0.0f, 1.0f) {}

    // Initialize with current state
    void initialize(const SimVariables& vars) {
        prob_bb = vars.prob_bb_initial;
        prob_bi = vars.prob_bi_initial;
        prob_ii = vars.prob_ii_initial;
        count_bb = static_cast<int>(vars.num_bb);
        count_bi = static_cast<int>(vars.num_bi);
        count_ii = static_cast<int>(vars.num_ii);
    }

    // Select a random swap type based on current probabilities
    SwapType select_swap_type() {
        float r = dist(gen);
        if (r < prob_bb) return SwapType::BB;
        if (r < prob_bb + prob_bi) return SwapType::BI;
        return SwapType::II;
    }

    // Select random atoms for a given swap type
    std::pair<int, int> select_atoms(SwapType type, const std::vector<int>& boundary_atoms, 
                                    const std::vector<int>& interior_atoms) {
        switch (type) {
            case SwapType::BB: {
                std::uniform_int_distribution<> dist(0, boundary_atoms.size() - 1);
                int idx1 = dist(gen);
                int idx2;
                do { idx2 = dist(gen); } while (idx2 == idx1 && boundary_atoms.size() > 1);
                return {boundary_atoms[idx1], boundary_atoms[idx2]};
            }
            case SwapType::II: {
                std::uniform_int_distribution<> dist(0, interior_atoms.size() - 1);
                int idx1 = dist(gen);
                int idx2;
                do { idx2 = dist(gen); } while (idx2 == idx1 && interior_atoms.size() > 1);
                return {interior_atoms[idx1], interior_atoms[idx2]};
            }
            case SwapType::BI: {
                std::uniform_int_distribution<> dist_b(0, boundary_atoms.size() - 1);
                std::uniform_int_distribution<> dist_i(0, interior_atoms.size() - 1);
                return {boundary_atoms[dist_b(gen)], interior_atoms[dist_i(gen)]};
            }
            default:
                throw std::runtime_error("Invalid swap type");
        }
    }

    // Calculate the probability of selecting a specific swap
    float calculate_selection_probability(SwapType type) {
        switch (type) {
            case SwapType::BB:
                return count_bb > 1 ? prob_bb / (count_bb * (count_bb - 1) / 2) : 0.0f;
            case SwapType::II:
                return count_ii > 1 ? prob_ii / (count_ii * (count_ii - 1) / 2) : 0.0f;
            case SwapType::BI:
                return count_bi > 0 ? prob_bi / count_bi : 0.0f;
            default:
                return 0.0f;
        }
    }
};

#endif