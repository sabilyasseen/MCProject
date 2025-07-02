#include "Move.hpp"
#include <algorithm>

void Move::deterministic_build_from_swap_delayed_rejection(std::vector<Swap> input_swaps, Grid& grid, bool debug){
    if (debug) std::cout << "[DEBUG][Move] Entering deterministic_build_from_swap_delayed_rejection with " << input_swaps.size() << " swaps." << std::endl;
    swap_history.push_back(input_swaps);
    any_new=false;
    tier=0;
    DeltaE=0.0f;
    DeltaE_teirs.push_back(0.0f);
    final_forward_probability=1.0f;
    final_reverse_probability=1.0f;

    if (debug) std::cout << "[DEBUG][Move] Finished deterministic_build_from_swap_delayed_rejection." << std::endl;
}

bool Move::evaluate_candidates_dr(Grid& grid, bool debug){
    if (debug) std::cout << "[DEBUG][Move] Entering evaluate_candidates_dr for swap_history size " << swap_history.size() << std::endl;
    primary_atoms.push_back(std::vector<int>(swap_history.back().size()));
    secondary_atoms.push_back(std::vector<int>(swap_history.back().size()));
    forward_expansion_teir_probabilities.push_back(1.0f);
    reverse_expansion_teir_probabilities.push_back(1.0f);
    DeltaE_teirs.push_back(0.0f);
    for(int i=0; i<swap_history.back().size(); i++){
        if(visited.find(swap_history.back()[i].primary_atom) != visited.end() || visited.find(swap_history.back()[i].secondary_atom) != visited.end()) {
            swap_history.back().erase(swap_history.back().begin() + i);
            primary_atoms.back().resize(primary_atoms.back().size() - 1);
            secondary_atoms.back().resize(secondary_atoms.back().size() - 1);
            continue;
        }
        Swap swap = swap_history.back()[i];
        final_forward_probability*=swap.initial_selection_prob;
        forward_expansion_teir_probabilities.back()*=swap.initial_selection_prob;
        visited.insert(swap.primary_atom);
        primary_atoms.back()[i]=swap.primary_atom;
        visited.insert(swap.secondary_atom);
        secondary_atoms.back()[i]=swap.secondary_atom;
        possibly_changed_atoms(swap.primary_atom, grid);
        possibly_changed_atoms(swap.secondary_atom, grid);
        num_changed+=2;
        // Only process if species are different
        if(grid.array[swap.primary_atom].cell_species != grid.array[swap.secondary_atom].cell_species) {
            if(debug) std::cout << "[DEBUG][Move] Swapping species for atoms " << swap.primary_atom << " and " << swap.secondary_atom << std::endl;
            // Remove from cluster before switching species
            grid.remove_atom(swap.primary_atom);
            grid.remove_atom(swap.secondary_atom);
            grid.change_specie(swap.primary_atom);
            DeltaE_teirs.back()+=J*delta_energy[grid.array[swap.primary_atom].local_config];
            grid.change_specie(swap.secondary_atom);
            DeltaE_teirs.back()+=J*delta_energy[grid.array[swap.secondary_atom].local_config]; 
        } else {
            if(debug) std::cout << "[DEBUG][Move] Skipping swap for atoms " << swap.primary_atom << " and " << swap.secondary_atom << " (same species)" << std::endl;
        }
    }
    DeltaE+=DeltaE_teirs.back();
    if (debug) std::cout << "[DEBUG][Move] Finished evaluate_candidates_dr." << std::endl;
    return true;
}

int Move::get_primary_expansion_index_candidates(Grid& grid, bool debug){   
    return 0;
}

int Move::get_secondary_expansion_index_candidates(Grid& grid, bool debug){
    return 0;
}   

float Move::calculate_reverse_expansion_prob(std::vector<Swap> input_swaps, Grid& grid, bool debug){
    if (debug) std::cout << "[DEBUG][Move] Entering calculate_reverse_expansion_prob with " << input_swaps.size() << " swaps." << std::endl;
    Swap placeholder_swap;
    placeholder_swap.primary_atom=0;
    placeholder_swap.secondary_atom=0;
    placeholder_swap.primary_atom_species=0;
    placeholder_swap.secondary_atom_species=1;
    reverse_expansion_prob=1.0f;
    std::unordered_set<int> swapping_into;
    for(Swap swap : input_swaps){
        if(swapping_into.find(swap.secondary_atom) == swapping_into.end()){
            swapping_into.insert(swap.secondary_atom);
        }
        else if(swapping_into.find(swap.primary_atom) == swapping_into.end()){
            swapping_into.insert(swap.primary_atom);
        }
        else{
            std::cout << "[ERROR] Repeated swap in expansion teir" << std::endl;
            std::cout << "Swap: " << swap.primary_atom << " " << swap.secondary_atom << std::endl;
            std::cout << "Swapping into: ";
            for(int atom : swapping_into){
                std::cout << atom << " ";
            }
            std::cout << std::endl;
            return 0.0f;
        }
    }
    std::vector<int> new_primary;
    std::vector<int> new_secondary;
    std::unordered_set<int> seen;
    for(Swap swap : swap_history.back()){
        for(auto expansion_index : grid.array[swap.primary_atom].expansion_indexes){
            if(swapping_into.find(expansion_index) == swapping_into.end() && seen.find(expansion_index) == seen.end()){
                new_primary.push_back(expansion_index);
                seen.insert(expansion_index);
            }
        }
        for(auto expansion_index : grid.array[swap.secondary_atom].expansion_indexes){
            if(swapping_into.find(expansion_index) == swapping_into.end() && seen.find(expansion_index) == seen.end()){
                new_secondary.push_back(expansion_index);
                seen.insert(expansion_index);
            }
        }
        // Filter by species - completely skip processing if species don't match
        bool species_match = (grid.array[swap.secondary_atom].cell_species == placeholder_swap.primary_atom_species) &&
                            (grid.array[swap.primary_atom].cell_species == placeholder_swap.secondary_atom_species);
        
        if(species_match) {
            bool valid_swap=is_valid_swap(grid, swap.secondary_atom, swap.primary_atom, placeholder_swap, false, false, debug);
            if(valid_swap){
                reverse_expansion_prob*=1.0f-tolerance;
            }
            else{
                reverse_expansion_prob*=tolerance; 
            }
        }
        // If species don't match, completely omit from consideration - no probability modification
    }
    int min_size = std::min(new_primary.size(), new_secondary.size());
    std::random_device rd;
    std::mt19937 rng(rd());
    std::shuffle(new_primary.begin(), new_primary.end(), rng);
    std::shuffle(new_secondary.begin(), new_secondary.end(), rng);


    int max_value=0;
    int N=min_size+ input_swaps.size()-1;
    for(int i=0;i<input_swaps.size();i++){
        int random_value=std::uniform_int_distribution<int>(0, N)(rng);
        if(random_value>max_value){
            max_value=random_value;
        }
    }
    int new_size = max_value-input_swaps.size();
    if(new_size<0){
        new_size=0;
    }
    new_primary.resize(new_size); 
    new_secondary.resize(new_size);

    for(int i=0; i<new_size; i++){
        // Filter by species - completely skip processing if species don't match
        bool species_match = (grid.array[new_primary[i]].cell_species == placeholder_swap.primary_atom_species) &&
                            (grid.array[new_secondary[i]].cell_species == placeholder_swap.secondary_atom_species);
        
        if(species_match) {
            bool valid_swap=false;
            valid_swap=is_valid_swap(grid, new_primary[i], new_secondary[i], placeholder_swap, false, false, debug);
            if(valid_swap){
                reverse_expansion_prob*=tolerance;
            }
            else{   
                reverse_expansion_prob*=1.0f-tolerance;
            }
        }
        // If species don't match, completely omit from consideration - no probability modification
    }
    if (debug) std::cout << "[DEBUG][Move] Finished calculate_reverse_expansion_prob. Result: " << reverse_expansion_prob << std::endl;
    return reverse_expansion_prob;
}
float Move::calculate_reverse_delayed_rejection(std::vector<Swap> reverse_swaps, Grid& grid, bool debug){
    if (debug) std::cout << "[DEBUG][Move] Entering calculate_reverse_delayed_rejection with " << reverse_swaps.size() << " swaps." << std::endl;
    
    final_reverse_probability=1.0f;
    if (debug) std::cout << "[DEBUG][Move] Initial final_reverse_probability: " << final_reverse_probability << std::endl;
    
    for(Swap swap : reverse_swaps){
        if (debug) std::cout << "[DEBUG][Move] Processing swap with initial_selection_prob: " << swap.initial_selection_prob << std::endl;
        final_reverse_probability*=swap.initial_selection_prob;
        if (debug) std::cout << "[DEBUG][Move] Updated final_reverse_probability: " << final_reverse_probability << std::endl;
    }
    
    if (debug) std::cout << "[DEBUG][Move] final_forward_probability: " << final_forward_probability << std::endl;
    float selection_ratio=final_reverse_probability/final_forward_probability;
    if (debug) std::cout << "[DEBUG][Move] Selection ratio (reverse/forward): " << selection_ratio << std::endl;
    
    if (debug) std::cout << "[DEBUG][Move] Beta: " << Beta << ", DeltaE: " << DeltaE << std::endl;
    float boltzmann_factor = std::exp(-Beta*DeltaE);
    if (debug) std::cout << "[DEBUG][Move] Boltzmann factor: " << boltzmann_factor << std::endl;
    
    float acceptance_prob = selection_ratio*boltzmann_factor;

    if (debug) {
        std::cout << "[DEBUG][Move] Final probabilities breakdown:" << std::endl;
        std::cout << "  - Final reverse probability: " << final_reverse_probability << std::endl;
        std::cout << "  - Final forward probability: " << final_forward_probability << std::endl;
        std::cout << "  - Selection ratio: " << selection_ratio << std::endl;
        std::cout << "  - Boltzmann factor: " << boltzmann_factor << std::endl;
        std::cout << "  - Final acceptance probability: " << acceptance_prob << std::endl;
    }
    acceptance_prob=std::min(1.0f,acceptance_prob);
    return acceptance_prob;
}
bool Move::expand_swap_dr(Grid& grid, bool debug){
    if(debug){
        std::cout << "[DEBUG][Move] Expanding swap" << std::endl;
    }
    Swap placeholder_swap;
    placeholder_swap.primary_atom=0;
    placeholder_swap.secondary_atom=0;
    placeholder_swap.primary_atom_species=0;
    placeholder_swap.secondary_atom_species=1;
    swap_history.push_back(std::vector<Swap>(0));
    bool pairs_found=false;
    std::vector<int> new_primary;
    std::vector<int> new_secondary;
    std::unordered_set<int> new_primary_set;
    std::unordered_set<int> new_secondary_set;
    int num_atoms = 0;
    for(int i=0; i<primary_atoms.back().size(); i++){
            for(auto expansion_index : grid.array[primary_atoms.back()[i]].expansion_indexes){
                if(visited.find(expansion_index) == visited.end() && new_primary_set.find(expansion_index) == new_primary_set.end()){
                    new_primary.push_back(expansion_index);
                    new_primary_set.insert(expansion_index);
                }
            }
            for(auto expansion_index : grid.array[secondary_atoms.back()[i]].expansion_indexes){
                if(visited.find(expansion_index) == visited.end() && new_secondary_set.find(expansion_index) == new_secondary_set.end()){
                    new_secondary.push_back(expansion_index);
                    new_secondary_set.insert(expansion_index);
                }
            }
    }
    int min_size = std::min(new_primary.size(), new_secondary.size());
    min_size=std::min(min_size,max_num_atoms/2);
    std::shuffle(new_primary.begin(), new_primary.end(), gen);
    std::shuffle(new_secondary.begin(), new_secondary.end(), gen);
    new_primary.resize(min_size);
    new_secondary.resize(min_size);
    reverse_expansion_prob=1.0f;
    int num_swapped=0;
    for(int i=0; i<min_size; i++){
        if(num_swapped>=max_num_atoms){
            break;
        }
        // Filter by species - completely skip this iteration if species don't match
        bool species_match = (grid.array[new_primary[i]].cell_species == placeholder_swap.primary_atom_species) &&
                            (grid.array[new_secondary[i]].cell_species == placeholder_swap.secondary_atom_species);
        
        if(!species_match) {
            // Species don't match - completely omit this pair from consideration
            continue;
        }
        
        bool valid_swap=is_valid_swap(grid, new_primary[i], new_secondary[i], placeholder_swap, false, false, debug);
        float r = dist(gen);
        if(valid_swap){
            if(r<tolerance){
                valid_swap=false;
                reverse_expansion_prob*=tolerance;
                reverse_expansion_teir_probabilities.back()*=tolerance;
            }
            else{
                pairs_found=true;
                reverse_expansion_prob*=1.0f-tolerance;
                reverse_expansion_teir_probabilities.back()*=1.0f-tolerance;
            }
        }
        else{
            if(r<tolerance){    
                pairs_found=true;
                valid_swap=true;
                reverse_expansion_prob*=tolerance;
                reverse_expansion_teir_probabilities.back()*=tolerance;
            }
            else{   
                reverse_expansion_prob*=1.0f-tolerance; 
                reverse_expansion_teir_probabilities.back()*=1.0f-tolerance;
            }
        }
        if(pairs_found){
            Swap swap;
            swap.primary_atom=new_primary[i];
            swap.secondary_atom=new_secondary[i];
            swap.initial_selection_prob=1.0f;
            swap_history.back().push_back(swap);
            num_swapped+=2;
        }
    }
    if(min_size>0){
        pairs_found=true;
    }


    final_forward_probability*=reverse_expansion_prob;
    if (debug) std::cout << "[DEBUG][Move] Finished expand_swap_dr. pairs_found=" << pairs_found << std::endl;
    return pairs_found;
}

bool Move::reclassify_singles(Grid& grid, std::vector<int> main_atoms){
    std::unordered_set<int> main_atoms_set;
    std::unordered_set<int> neighbor_atoms_set;
    std::unordered_set<int> neighbor_neighbors_atoms_set;
    std::vector<int> main_atoms_vec;
    std::vector<int> neighbor_atoms_vec;
    std::vector<int> neighbor_neighbors_atoms_vec;

    // First tier: Add main atoms and collect their neighbors
    for (int atom : main_atoms) {
        if (main_atoms_set.insert(atom).second) {
            main_atoms_vec.push_back(atom);
        }
    }

    for (int atom : main_atoms) {
        if(neighbor_atoms_set.insert(atom).second){
            main_atoms_vec.push_back(atom);
        }
        for (int neighbor : grid.array[atom].neighbor_indexes) {
            if (neighbor_atoms_set.insert(neighbor).second) {
                neighbor_atoms_vec.push_back(neighbor);
            }
        }
    }

    for (int neighbor : neighbor_atoms_vec) {
        if(neighbor_neighbors_atoms_set.insert(neighbor).second){
            neighbor_neighbors_atoms_vec.push_back(neighbor);
        }
        for (int neighbor_neighbor : grid.array[neighbor].neighbor_indexes) {
            if (neighbor_neighbors_atoms_set.insert(neighbor_neighbor).second) {
                neighbor_neighbors_atoms_vec.push_back(neighbor_neighbor);
            }
        }
    }


    // neighbor_neighbors_atoms_vec now includes all atoms (tier 2)
    // neighbor_atoms_vec includes main atoms and their direct neighbors (tier 1)
    // main_atoms_vec contains only the original atoms (tier 0)

    // Reclassify all affected atoms, starting from the outermost tier
    grid.reclassify_cluster_status(neighbor_neighbors_atoms_vec);

    grid.reclassify_cluster_status(neighbor_atoms_vec);

    grid.reclassify_cluster_status(main_atoms_vec);

    return true;
}



bool Move::calculate_cluster_changes(Grid& grid, bool debug) {
    if(debug){  
        std::cout << "[DEBUG][Move] Calculating cluster changes" << std::endl;
    }
    if(debug){
        std::cout << "[DEBUG][Move] Reclassifying cluster status..." << std::endl;
    }
    if(!grid.reclassify_cluster_status(possibly_changed_vec)) {
        if(debug) {
            std::cout << "[DEBUG][Move] Restricted pattern encountered during reclassification – aborting move." << std::endl;
        }
        return false;
    }
    if(debug){
        std::cout << "[DEBUG][Move] Reclassifying Boundary status minsize5: " << possibly_changed_vec.size() << std::endl;
    }
    grid.reclassify_boundary_status_minsize5(possibly_changed_vec);
    if(debug){
        std::cout << "[DEBUG][Move] Reclassifying Boundary status minsize2: " << possibly_changed_vec.size() << std::endl;
    }
    grid.reclassify_boundary_status_minsize2(possibly_changed_vec);
    if(debug){
        std::cout << "[DEBUG][Move] Reclassifying Boundary status: " << possibly_changed_vec.size() << std::endl;
    }
    grid.reclassify_boundary_status(possibly_changed_vec);

    return true;
}




void Move::possibly_changed_atoms(int atom, Grid& grid) {
    // Add the atom itself
    if (possibly_changed.insert(atom).second) {
        possibly_changed_vec.push_back(atom);
    }

    // Add all neighbors
    for (int neighbor : grid.array[atom].neighbor_indexes) {
        if (possibly_changed.insert(neighbor).second) {
            possibly_changed_vec.push_back(neighbor);
        }
        
        // Add all neighbors' neighbors
        for (int neighbor_neighbor : grid.array[neighbor].neighbor_indexes) {
            if (possibly_changed.insert(neighbor_neighbor).second) {
                possibly_changed_vec.push_back(neighbor_neighbor);
            }
            
            // Add neighbors' neighbors' neighbors
            for (int neighbor_neighbor_neighbor : grid.array[neighbor_neighbor].neighbor_indexes) {
                if (possibly_changed.insert(neighbor_neighbor_neighbor).second) {
                    possibly_changed_vec.push_back(neighbor_neighbor_neighbor);
                }
            }

        }
    }
}

std::vector<expansion_candidate> Move::deterministic_expand(Grid& grid, bool debug) {
    any_new=false;
    if(debug){
        std::cout << "[DEBUG] Deterministic expand" << std::endl;
    }
    std::vector<int> primary_seeds = primary_atoms.back();
    std::vector<int> secondary_seeds = secondary_atoms.back();
    // Ensure these vectors are large enough for the number of expansion indexes
    // (not just the number of tiers)
    // This prevents out-of-bounds access below
    tolerant_forward_atoms_accepted_tiers.push_back(0);
    tolerant_forward_atoms_rejected_tiers.push_back(0);
    
    std::vector<int> primary_expansion_indexes;
    std::vector<int> secondary_expansion_indexes;
    
    // Determine max number of seeds between primary and secondary
    int max_seeds = std::max(primary_seeds.size(), secondary_seeds.size());

    // 50% chance to determine which set claims first

    // Roll once to determine order
    bool check_primary_first = (dist(gen) < 0.5f);
    if(debug){
        std::cout << "[DEBUG] Check primary first: " << check_primary_first << std::endl;
    }
    // For each potential seed index
    for(int i = 0; i < max_seeds; i++) {
        if(check_primary_first) {
            // Check primary seed first if it exists
            if(i < primary_seeds.size()) {
                for(int expansion_index : grid.array[primary_seeds[i]].expansion_indexes) {
                    if(visited.find(expansion_index) == visited.end()) {
                        primary_expansion_indexes.push_back(expansion_index);
                        if(debug){
                            std::cout << "[DEBUG] Primary expansion index: " << expansion_index << std::endl;
                        }
                        visited.insert(expansion_index);
                    }
                }
            }

            // Then check secondary seed if it exists
            if(i < secondary_seeds.size()) {
                for(int expansion_index : grid.array[secondary_seeds[i]].expansion_indexes) {
                    if(visited.find(expansion_index) == visited.end()) {
                        secondary_expansion_indexes.push_back(expansion_index);
                        if(debug){
                            std::cout << "[DEBUG] Secondary expansion index: " << expansion_index << std::endl;
                        }
                        visited.insert(expansion_index);
                    }
                }
            }
        } else {
            // Check secondary seed first if it exists
            if(i < secondary_seeds.size()) {
                for(int expansion_index : grid.array[secondary_seeds[i]].expansion_indexes) {
                    if(visited.find(expansion_index) == visited.end()) {
                        secondary_expansion_indexes.push_back(expansion_index);
                        visited.insert(expansion_index);
                        if(debug){
                            std::cout << "[DEBUG] Secondary expansion index: " << expansion_index << std::endl;
                        }
                    }
                }
            }

            // Then check primary seed if it exists
            if(i < primary_seeds.size()) {
                for(int expansion_index : grid.array[primary_seeds[i]].expansion_indexes) {
                    if(visited.find(expansion_index) == visited.end()) {
                        primary_expansion_indexes.push_back(expansion_index);
                        visited.insert(expansion_index);
                        if(debug){
                            std::cout << "[DEBUG] Primary expansion index: " << expansion_index << std::endl;
                        }
                    }
                }
            }
        }
    }
    std::shuffle(primary_expansion_indexes.begin(), primary_expansion_indexes.end(), gen);
    std::shuffle(secondary_expansion_indexes.begin(), secondary_expansion_indexes.end(), gen);
        
    //They should be the same size, but just in case they arent then they are forced to be after shuffling
    if(primary_expansion_indexes.size() != secondary_expansion_indexes.size()) {
        int min_size = std::min(primary_expansion_indexes.size(), secondary_expansion_indexes.size());
        primary_expansion_indexes.resize(min_size);
        secondary_expansion_indexes.resize(min_size);
    }
    std::vector<expansion_candidate> primary_expansion_candidates(primary_expansion_indexes.size());
    std::vector<expansion_candidate> secondary_expansion_candidates(secondary_expansion_indexes.size());

    for(int i=0; i<primary_expansion_indexes.size(); i++) {
        // Filter by species - completely skip processing if species don't match
        bool species_match = (grid.array[primary_expansion_indexes[i]].cell_species == forward.primary_atom_species) &&
                            (grid.array[secondary_expansion_indexes[i]].cell_species == forward.secondary_atom_species);
        
        primary_expansion_candidates[i].position = primary_expansion_indexes[i];
        primary_expansion_candidates[i].num_neighbors = 1;
        secondary_expansion_candidates[i].position = secondary_expansion_indexes[i];
        secondary_expansion_candidates[i].num_neighbors = 1;

        if(!species_match) {
            // Species don't match - mark as invalid and skip tolerance processing
            primary_expansion_candidates[i].prob_denial = 1.0f;
            primary_expansion_candidates[i].acceptance = false;
            primary_expansion_candidates[i].base_prob_denial = 1.0f;
            secondary_expansion_candidates[i].prob_denial = 1.0f;
            secondary_expansion_candidates[i].acceptance = false;
            secondary_expansion_candidates[i].base_prob_denial = 1.0f;
            continue;
        }
        
        bool valid_swap = is_valid_swap(grid, primary_expansion_indexes[i], secondary_expansion_indexes[i], forward, cluster_distinction, celltype_distinction, debug);

        if(!valid_swap){
            float r = dist(gen);

            if (r < tolerance) {
                valid_swap = true;
                tolerant_forward_atoms_accepted++;
                tolerant_forward_atoms_accepted_tiers.back()++;
                final_forward_probability*=tolerance;
                forward_expansion_teir_probabilities.back()*=tolerance;
            }
            else{
                tolerant_forward_atoms_rejected++;
                tolerant_forward_atoms_rejected_tiers.back()++;
                final_forward_probability*=1.0f-tolerance;
                forward_expansion_teir_probabilities.back()*=1.0f-tolerance;
            }
        }
        else if(valid_swap){
            float r = dist(gen);
            if (r < tolerance) {
                valid_swap = false;
                tolerant_forward_atoms_accepted++;
                tolerant_forward_atoms_accepted_tiers.back()++;
                final_forward_probability*=tolerance;
                forward_expansion_teir_probabilities.back()*=tolerance;
            }
            else{   
                tolerant_forward_atoms_rejected++;
                tolerant_forward_atoms_rejected_tiers.back()++;
                final_forward_probability*=1.0f-tolerance;
                forward_expansion_teir_probabilities.back()*=1.0f-tolerance;
            }
        }
        if(valid_swap){
            any_new=true;
            primary_expansion_candidates[i].prob_denial = 0.0f;
            primary_expansion_candidates[i].acceptance = true;
            primary_expansion_candidates[i].base_prob_denial = 0.0f;

            secondary_expansion_candidates[i].prob_denial = 0.0f;
            secondary_expansion_candidates[i].acceptance = true;
            secondary_expansion_candidates[i].base_prob_denial = 0.0f;

        }
        else{
            primary_expansion_candidates[i].prob_denial = 1.0f;
            primary_expansion_candidates[i].acceptance = false;
            primary_expansion_candidates[i].base_prob_denial = 1.0f;
            secondary_expansion_candidates[i].prob_denial = 1.0f;
            secondary_expansion_candidates[i].acceptance = false;
            secondary_expansion_candidates[i].base_prob_denial = 1.0f;
            
        }
    }
    primary_expansion_history.push_back(primary_expansion_candidates);
    secondary_expansion_history.push_back(secondary_expansion_candidates);
    return primary_expansion_candidates;
}


bool Move::evaluate_deterministic_candidates(Grid& grid, bool debug){
    // Make sure primary_atoms and secondary_atoms vectors have enough slots for this tier
    // and the next tier (expansion_teir+1)
    DeltaE_teirs.push_back(0.0f);
    primary_atoms.push_back(std::vector<int>());
    secondary_atoms.push_back(std::vector<int>());



    primary_atoms.back().resize(primary_expansion_history.back().size());
    secondary_atoms.back().resize(secondary_expansion_history.back().size());

    if (debug) {
        std::cout << "\n[DEBUG] --- evaluate_candidates (tier " << tier << ") ---\n";
    }


    
    size_t min_size = std::min(primary_expansion_history.back().size(), secondary_expansion_history.back().size());
    
  //  std::shuffle(primary_expansion_history[expansion_teir].begin(), primary_expansion_history[expansion_teir].end(), gen);
  //  std::shuffle(secondary_expansion_history[expansion_teir].begin(), secondary_expansion_history[expansion_teir].end(), gen);
    


    for (size_t idx = 0; idx < min_size; ++idx) {
        auto& primary = primary_expansion_history.back()[idx];
        auto& secondary = secondary_expansion_history.back()[idx];



        if (primary.acceptance && secondary.acceptance) {
        visited.insert(primary.position);
        visited.insert(secondary.position);
        any_new=true;
            primary_atoms.back().push_back(primary.position);
            secondary_atoms.back().push_back(secondary.position);
            // Only now update probability for accepted pairs
            // Check if atoms are from different species
                if (debug) {
                    std::cout << "[DEBUG] Swapping atoms: " << primary.position << " and " << secondary.position << std::endl;
                }

            if(grid.array[primary.position].cell_species != grid.array[secondary.position].cell_species){
                grid.remove_atom(primary.position);
                grid.remove_atom(secondary.position);
                grid.change_specie(primary.position);
                float primary_deltaE=J * delta_energy[grid.array[primary.position].local_config];
                grid.change_specie(secondary.position);
                float secondary_deltaE=J * delta_energy[grid.array[secondary.position].local_config];
                DeltaE += primary_deltaE + secondary_deltaE;
                DeltaE_teirs.back() += primary_deltaE;
                DeltaE_teirs.back() += secondary_deltaE;
                possibly_changed_atoms(primary.position, grid);
                possibly_changed_atoms(secondary.position, grid);
                num_changed += 2;

            }
        }  

    }

 
    if (debug) {
        std::cout << "[DEBUG] Number of atoms swapped: " << num_changed << std::endl;
        std::cout << "[DEBUG] Primary atoms: ";
        for (auto a : primary_atoms.back()) std::cout << a << " ";
        std::cout << "\n[DEBUG] Secondary atoms: ";
        for (auto a : secondary_atoms.back()) std::cout << a << " ";
        std::cout << std::endl;
    }
    return true;
}


bool Move::deterministic_build_from_swap(Swap input_swap, Grid& grid, bool debug){
    if (debug) std::cout << "[DEBUG] Entering deterministic_build_from_swap" << std::endl;
    // Clear any existing data
    primary_atoms.clear();
    secondary_atoms.clear();
    primary_expansion_history.clear();
    secondary_expansion_history.clear();
    forward_expansion_teir_probabilities.clear();
    reverse_expansion_teir_probabilities.clear();
    
    final_forward_probability = input_swap.initial_selection_prob;
    expansion_candidate primary_candidate;
    primary_candidate.position = input_swap.primary_atom;
    primary_candidate.num_neighbors = 1;
    primary_candidate.prob_denial = 0.0f;
    primary_candidate.acceptance = true;
    primary_candidate.base_prob_denial=0.0f;
    primary_expansion_history.push_back({primary_candidate});
    
    expansion_candidate secondary_candidate;
    secondary_candidate.position = input_swap.secondary_atom;
    secondary_candidate.num_neighbors = 1;
    secondary_candidate.prob_denial = 0.0f;
    secondary_candidate.acceptance = true;
    secondary_candidate.base_prob_denial=0.0f;
    secondary_expansion_history.push_back({secondary_candidate});
    
    visited.insert(input_swap.primary_atom);
    visited.insert(input_swap.secondary_atom);
    possibly_changed_atoms(input_swap.primary_atom, grid);
    possibly_changed_atoms(input_swap.secondary_atom, grid);
    forward = input_swap;
    forward.expansion_prob=input_swap.expansion_prob;

    DeltaE = 0.0f;

    // Ensure vectors are initialized
    primary_atoms.resize(1);
    secondary_atoms.resize(1);
    final_forward_probability=forward.initial_selection_prob;
    forward_expansion_teir_probabilities.resize(1,1.0f);
    forward_expansion_teir_probabilities[0]=forward.initial_selection_prob;
    evaluate_deterministic_candidates(grid, debug); //Uses the regular evaluate function to properly track everything
    std::uniform_real_distribution<float> dist_type(0.0f, 1.0f);
    ///Any_new is true if the most recent expansion had any valid swaps in it. 
    ///Expand is true if r<forward.expansion_prob
    //The first 'expansion tier' is the swap itself
    tier=0;
    float r = dist_type(gen);   
    bool expand = (r<forward.expansion_prob);
    while(expand){
        if(debug){
            std::cout << "[DEBUG] Expanding" << std::endl;
        }
        forward_expansion_teir_probabilities.push_back(forward.expansion_prob);
        final_forward_probability*=forward.expansion_prob;
        deterministic_expand(grid, debug);
        evaluate_deterministic_candidates(grid, debug);

        if(any_new && num_changed<=6){
            tier++;
            r=dist_type(gen);
            expand=(r<forward.expansion_prob);
        }
        else{
            break;
        }

    }



    if(any_new){
        forward_expansion_teir_probabilities.back()*=1.0f-forward.expansion_prob;
        final_forward_probability*=1.0f-forward.expansion_prob;
        ///This means that the swap could have expanded, but it failed due to epansion prob
    }

    bool success = calculate_cluster_changes(grid, debug);
    return success;
}
void Move::evaluate_reverse_expansion_deterministic(Grid& grid, bool debug){
    reverse_expansion_teir_probabilities.clear();
    reverse_expansion_teir_probabilities.resize(forward_expansion_teir_probabilities.size(), reverse.expansion_prob);
    reverse_expansion_teir_probabilities[0]=reverse.initial_selection_prob;
    final_reverse_probability=reverse.initial_selection_prob;

    if(any_new){
        reverse_expansion_teir_probabilities.back()*=1.0f-reverse.expansion_prob;
        final_reverse_probability*=1.0f-reverse.expansion_prob;
    }
    final_reverse_probability*=std::pow(reverse.expansion_prob, forward_expansion_teir_probabilities.size()-1);


    for(int i=1; i<primary_expansion_history.size(); i++){
        
        for(int j=0; j<primary_expansion_history[i].size(); j++){

            // Filter by species - completely skip processing if species don't match
            bool species_match = (grid.array[primary_expansion_history[i][j].position].cell_species == reverse.primary_atom_species) &&
                                (grid.array[secondary_expansion_history[i][j].position].cell_species == reverse.secondary_atom_species);
            
            if(!species_match) {
                // Species don't match - completely omit from consideration
                continue;
            }
            
            bool valid_swap = is_valid_swap(grid, primary_expansion_history[i][j].position, secondary_expansion_history[i][j].position, reverse, cluster_distinction, celltype_distinction, debug);


            if(valid_swap && primary_expansion_history[i][j].acceptance && secondary_expansion_history[i][j].acceptance){
                reverse_expansion_teir_probabilities[i]*=tolerance;
                final_reverse_probability *= 1.0f-tolerance;
            }
            else if (!valid_swap && !primary_expansion_history[i][j].acceptance && !secondary_expansion_history[i][j].acceptance){
                reverse_expansion_teir_probabilities[i]*=1.0f-tolerance;
                final_reverse_probability *= 1.0f-tolerance;
            }
            else{
                reverse_expansion_teir_probabilities[i]*=tolerance;
                final_reverse_probability *= tolerance;
            }

        }
    }

}

float Move::get_final_acceptance_deterministic(Grid& grid, Swap reverse_swap, bool debug){
    
    reverse = reverse_swap;
    final_reverse_probability = reverse.initial_selection_prob;


    evaluate_reverse_expansion_deterministic(grid, debug);

    if(final_forward_probability <= 1e-20){
        if (debug) {
            std::cout << "[ERROR] final_forward_probability is zero or too small, setting acceptance_prob to 0." << std::endl;
        }
        return 0.0f;
    }
    
    float ratio = final_reverse_probability / final_forward_probability;
    float boltzmann_factor = std::exp(-Beta * DeltaE);
    float acceptance_prob = std::min(1.0f, ratio * boltzmann_factor);
    if (std::isnan(acceptance_prob)) {
        acceptance_prob = 0.0f;
    }
    if (debug) {
        std::cout << "[DEBUG] Reverse probability: " << final_reverse_probability << std::endl;
        std::cout << "[DEBUG] Forward probability: " << final_forward_probability << std::endl;
        std::cout << "[DEBUG] Ratio: " << ratio << std::endl;
        std::cout << "[DEBUG] Boltzmann factor: " << boltzmann_factor << std::endl;
        std::cout << "[DEBUG] Acceptance probability: " << acceptance_prob << std::endl;
    }
    return acceptance_prob;
}

bool Move::is_valid_swap(Grid& grid, int primary_candidate, int secondary_candidate, Swap reference_swap, bool cluster_match, bool celltype_match, bool debug){
    if(debug){
        std::cout << "[DEBUG] Checking if swap is valid" << std::endl;
    }
    int primary_atom=reference_swap.primary_atom;
    int secondary_atom=reference_swap.secondary_atom;
    
    if(cluster_match){
        if(grid.array[primary_candidate].cell_cluster_ID != reference_swap.primary_atom_clusterID){
            return false;
        }
        if(grid.array[secondary_candidate].cell_cluster_ID != reference_swap.secondary_atom_clusterID){
            return false;
        }
        if(debug){
            std::cout << "[DEBUG] Cluster match" << std::endl;
        }
    }
    if(celltype_match){
        if(grid.array[primary_candidate].cell_type != reference_swap.primary_atom_type){
            return false;
        }
        if(grid.array[secondary_candidate].cell_type != reference_swap.secondary_atom_type){
            return false;
        }
        if(debug){
            std::cout << "[DEBUG] Cell type match" << std::endl;
        }
    }
    return true;
}
    bool deterministic_build_from_swap(Swap input_swap, Grid& grid, bool debug);
