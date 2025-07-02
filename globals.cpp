#include "globals.hpp"
#include <array>
#include <algorithm>

std::unordered_map<std::vector<int>, float, VectorHash> Energy;
std::unordered_map<std::vector<int>, float, VectorHash> delta_energy;
const std::unordered_map<std::string, SimConditions> SIM_CONDITIONS = {
    {"BJ=0.1_Kawasaki", SimConditions(1.0/3.0f, 1.0/3.0f, 1.0/3.0f, 0.5f, 0.0f, 0.0f, 0.0f, 10.0f, 0.01f, 0.0f)}, //BJ = 10 * 0.01 = 0.1
    {"BJ=0.1_NewMC", SimConditions(0.55f, 0.31f, 0.14f,0.2f, 0.0f, 0.0f, 0.0f, 10.0f, 0.01f, 0.0f)}, //BJ = 10 * 0.01 = 0.1
    {"BJ=0.1Expand_NewMC", SimConditions(0.55f, 0.31f, 0.14f,0.2f, 0.5f, 0.5f, 0.5f, 10.0f, 0.01f, 0.0f)}, //BJ = 10 * 0.01 = 0.1

    {"BJ=0.44_Kawasaki", SimConditions(1.0/3.0f, 1.0/3.0f, 1.0/3.0f, 0.5f, 0.0f, 0.0f, 0.0f, 10.0f, 0.044f, 0.0f)}, //BJ = 10 * 0.044 = 0.44
    {"BJ=0.44_NewMC", SimConditions(0.86f, 0.13f, 0.01f, 0.2f, 0.0f, 0.0f, 0.0f, 10.0f, 0.044f, 0.0f)}, //BJ = 10 * 0.044 = 0.44
    {"BJ=0.44Expand_NewMC", SimConditions(0.86f, 0.13f, 0.01f, 0.2f, 0.5f, 0.5f, 0.5f, 10.0f, 0.044f, 0.0f)}, //BJ = 10 * 0.044 = 0.44
   
    {"BJ=1.5_Kawasaki", SimConditions(1.0/3.0f, 1.0/3.0f, 1.0/3.0f, 0.5f, 0.0f, 0.0f, 0.0f, 10.0f, 0.15f, 0.0f)}, //BJ = 10 * 0.15 = 1.5
    {"BJ=1.5_NewMC", SimConditions(0.984f, 0.013f, 0.003f, 0.2f, 0.0f, 0.0f, 0.0f, 10.0f, 0.15f, 0.0f)}, //BJ = 10 * 0.15 = 1.5
    {"BJ=1.5Expand_NewMC", SimConditions(0.984f, 0.013f, 0.003f, 0.2f, 0.5f, 0.5f, 0.5f, 10.0f, 0.15f, 0.0f)}, //BJ = 10 * 0.15 = 1.5

    {"BJ=-0.1_Kawasaki", SimConditions(1.0/3.0f, 1.0/3.0f, 1.0/3.0f, 0.5f, 0.0f, 0.0f, 0.0f, -10.0f, 0.01f, 0.0f, 0.0f, false)}, //BJ = -10 * 0.01 = -0.1
    {"BJ=-0.1_DR", SimConditions(1.0/3.0f, 1.0/3.0f, 1.0/3.0f, 0.2f, 0.0f, 0.0f, 0.0f, -10.0f, 0.01f, 0.0f, 0.5f, false)}, //BJ = -10 * 0.01 = -0.1
    {"BJ=-0.1_DR_Expand", SimConditions(1.0/3.0f, 1.0/3.0f, 1.0/3.0f, 0.2f, 0.5f, 0.5f, 0.5f, -10.0f, 0.01f, 0.0f, 0.5f, false)}, //BJ = -10 * 0.01 = -0.1
    {"BJ=-0.1_NewMC", SimConditions(0.55f, 0.31f, 0.14f, 0.2f, 0.0f, 0.0f, 0.0f, -10.0f, 0.01f, 0.0f, 0.5f, true)}, //BJ = -10 * 0.01 = -0.1
    {"BJ=-0.1_NewMC_Expand", SimConditions(0.55f, 0.31f, 0.14f, 0.5f, 0.5f, 0.5f, 0.5f, -10.0f, 0.01f, 0.2f, 0.5f, true)}, //BJ = -10 * 0.01 = -0.1

    {"BJ=-0.44_Kawasaki", SimConditions(1.0/3.0f, 1.0/3.0f, 1.0/3.0f, 0.5f, 0.0f, 0.0f, 0.0f, -10.0f, 0.044f, 0.0f, 0.0f, false)}, //BJ = -10 * 0.044 = -0.44
    {"BJ=-0.44_DR", SimConditions(1.0/3.0f, 1.0/3.0f, 1.0/3.0f, 0.5f, 0.0f, 0.0f, 0.0f, -10.0f, 0.044f, 0.0f, 0.5f, false)}, //BJ = -10 * 0.044 = -0.44
    {"BJ=-0.44_DR_Expand", SimConditions(1.0/3.0f, 1.0/3.0f, 1.0/3.0f, 0.5f, 0.5f, 0.5f, 0.5f, -10.0f, 0.044f, 0.0f, 0.5f, false)}, //BJ = -10 * 0.044 = -0.44
    {"BJ=-0.44_NewMC", SimConditions(0.86f, 0.13f, 0.01f, 0.5f, 0.0f, 0.0f, 0.0f, -10.0f, 0.044f, 0.0f, 0.5f, true)}, //BJ = -10 * 0.044 = -0.44
    {"BJ=-0.44_NewMC_Expand", SimConditions(0.86f, 0.13f, 0.01f, 0.5f, 0.5f, 0.5f, 0.5f, -10.0f, 0.044f, 0.2f, 0.5f, true)}, //BJ = -10 * 0.044 = -0.44

    {"BJ=-1.5_Kawasaki", SimConditions(1.0/3.0f, 1.0/3.0f, 1.0/3.0f, 0.5f, 0.0f, 0.0f, 0.0f, -10.0f, 0.15f, 0.0f, 0.0f, false)}, //BJ = -10 * 0.15 = -1.5
    {"BJ=-1.5_DR", SimConditions(1.0/3.0f, 1.0/3.0f, 1.0/3.0f, 0.2f, 0.0f, 0.0f, 0.0f, -10.0f, 0.15f, 0.0f, 0.5f, false)}, //BJ = -10 * 0.15 = -1.5
    {"BJ=-1.5_DR_Expand", SimConditions(1.0/3.0f, 1.0/3.0f, 1.0/3.0f, 0.2f, 0.5f, 0.5f, 0.5f, -10.0f, 0.15f, 0.0f, 0.5f, false)}, //BJ = -10 * 0.15 = -1.5
    {"BJ=-1.5_NewMC", SimConditions(0.984f, 0.013f, 0.003f, 0.2f, 0.0f, 0.0f, 0.0f, -10.0f, 0.15f, 0.0f, 0.5f, true)}, //BJ = -10 * 0.15 = -1.5
    {"BJ=-1.5_NewMC_Expand", SimConditions(0.984f, 0.013f, 0.003f, 0.5f, 0.5f, 0.5f, 0.5f, -10.0f, 0.15f, 0.2f, 0.5f, true)}, //BJ = -10 * 0.15 = -1.5



};

// ---------------- ClusterProperty global tables ----------------
std::vector<std::vector<int>> GLOBAL_PATTERN_VECTORS;
std::unordered_map<std::vector<int>, int, VectorHash, VectorEqual> PATTERN_LOOKUP;

void initialize_pattern_library() {
    GLOBAL_PATTERN_VECTORS.clear();
    PATTERN_LOOKUP.clear();

    auto invert = [](const std::vector<int> &vec){
        std::vector<int> inv(vec.size());
        for(size_t i=0;i<vec.size();++i) inv[i] = vec[i] == 0 ? 1 : 0;
        return inv;
    };

    auto insert_with_variants = [&](const std::vector<int>& base, int current_index){
        std::vector<int> temp = base;
        for(int r=0;r<5;++r){
            // rotation
            PATTERN_LOOKUP[temp] = current_index;
            PATTERN_LOOKUP[invert(temp)] = current_index;
            // rotate keeping first entry fixed (rotate positions 1..4)
            std::rotate(temp.begin()+1, temp.begin()+2, temp.end());
        }
    };
    
    for (int a = 0; a <= 1; ++a)
        for (int b = 0; b <= 1; ++b)
            for (int c = 0; c <= 1; ++c)
                for (int d = 0; d <= 1; ++d)
                    for (int e = 0; e <= 1; ++e) {
                        std::vector<int> pattern = {a,b,c,d,e};
                        if(PATTERN_LOOKUP.find(pattern) != PATTERN_LOOKUP.end()) {
                            continue; // pattern (or equivalent) already stored
                        }
                        int idx = static_cast<int>(GLOBAL_PATTERN_VECTORS.size());
                        GLOBAL_PATTERN_VECTORS.push_back(pattern);
                        insert_with_variants(pattern, idx);
                    }
}

void initialize_energy(){
    // Initialize energy map
    for (int center=0; center<=1; center++){
        for (int up=0; up<=1; up++){
            for (int right=0; right<=1; right++){
                for (int down=0; down<=1; down++){
                    for (int left=0; left<=1; left++){
                        std::vector<int>config={center,up,right,down,left};
                        float sum=0.0f;
                        float centerval=1.0f;
                        std::vector<int>shell={up,right,down,left};
                        if(center==0){centerval*=-1.0f;};
                        for (int bond : shell){
                            if(bond==0){sum+=centerval*-1.0f;}
                            else{sum+=centerval*1.0f;}
                        }
                        Energy[config]=-0.5f*sum;
                    }
                }
            }
        }
    }
    
    // After initializing Energy, initialize delta_energy
    initialize_delta_energy();
}

void initialize_delta_energy() {
    // Clear any existing entries
    delta_energy.clear();
    
    // Initialize delta_energy map
    for (const auto& entry : Energy) {
        const std::vector<int>& post_config = entry.first;
        float post_energy = entry.second;
        
        // Create pre_config by switching the first value (center) to the other species
        std::vector<int> pre_config = post_config;
        pre_config[0] = (pre_config[0] == 0) ? 1 : 0;
        
        // Get the pre-energy from the Energy map
        float pre_energy = Energy[pre_config];
        
        // Calculate and store delta_energy = post_energy - pre_energy
        delta_energy[post_config] = 2*(post_energy - pre_energy);
    }
}

std::vector<std::vector<int>> GLOBAL_PATTERN_VECTORS_ROT;
std::unordered_map<std::vector<int>, int, VectorHash, VectorEqual> PATTERN_LOOKUP_ROT;

void initialize_pattern_library_rot() {
    // Clear previous content
    GLOBAL_PATTERN_VECTORS_ROT.clear();
    PATTERN_LOOKUP_ROT.clear();

    // Utility lambda: rotate positions 1..end keeping first entry fixed
    auto rotate_keep_first = [](std::vector<int>& vec) {
        if (vec.size() > 1) {
            std::rotate(vec.begin() + 1, vec.begin() + 2, vec.end());
        }
    };

    // Helper to insert a pattern and its 4 rotations (total 5 including base)
    auto insert_with_rotations = [&](const std::vector<int>& base, int index) {
        std::vector<int> temp = base;
        for (int r = 0; r < static_cast<int>(base.size()); ++r) {
            PATTERN_LOOKUP_ROT[temp] = index;
            rotate_keep_first(temp);
        }
    };

    // Generate all 5-length binary patterns (no inversions – each pattern treated independently)
    int current_index = 0;

    for (int a = 0; a <= 1; ++a)
        for (int b = 0; b <= 1; ++b)
            for (int c = 0; c <= 1; ++c)
                for (int d = 0; d <= 1; ++d)
                    for (int e = 0; e <= 1; ++e) {
                        std::vector<int> pattern = {a, b, c, d, e};
                        // Skip patterns already represented by a rotation of an earlier pattern
                        if (PATTERN_LOOKUP_ROT.find(pattern) != PATTERN_LOOKUP_ROT.end()) {
                            continue;
                        }
                        // Store the canonical pattern
                        GLOBAL_PATTERN_VECTORS_ROT.push_back(pattern);
                        insert_with_rotations(pattern, current_index);
                        ++current_index;
                    }
} 