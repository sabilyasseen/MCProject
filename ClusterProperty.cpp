#include "ClusterProperty.hpp"
#include <algorithm>
#include <unordered_set>

ClusterProperty::ClusterProperty() {}

ClusterProperty::ClusterProperty(int minsize) {
    min_size = minsize;
}

int ClusterProperty::get_minsize() {
    return min_size;
}

void ClusterProperty::add_config(std::vector<int> property) {
    if (property.size() != min_size) {
        return;
    }

    // If the base configuration (or any of its derived forms) is already present,
    // it will already exist in the hash_set. Just check the original vector.
    if (hash_set.find(property) != hash_set.end()) {
        return; // nothing new to add
    }

    // Not present yet – generate rotations/inversions and store all unique ones
    generate_configurations(property);
}

void ClusterProperty::generate_configurations(std::vector<int> property) {
    auto invert = [](const std::vector<int>& v){
        std::vector<int> inv(v.size());
        for (size_t i = 0; i < v.size(); ++i) inv[i] = (v[i] == 0 ? 1 : 0);
        return inv;
    };

    auto insert_unique = [this](const std::vector<int>& cfg){
        if (hash_set.insert(cfg).second) {
            properties.push_back(cfg);
        }
    };

    if (min_size == 5) {
        // Generate 5 rotations of the original pattern
        std::vector<int> current = property;
        for (int r = 0; r < 5; ++r) {
            insert_unique(current);
            std::rotate(current.begin() + 1, current.begin() + 2, current.end());
        }
    } else if (min_size == 2) {
        insert_unique(property);
    } else if (min_size == 1) {
        insert_unique(property);
    }

    // Update homogeneity flag
    homogenous = true;
    for (int v : property) {
        if (v != property[0]) { homogenous = false; break; }
    }
}

// Generate inversions for every stored property (runs once manually)
void ClusterProperty::generate_inversions() {
    auto invert = [](const std::vector<int>& v){
        std::vector<int> inv(v.size());
        for (size_t i = 0; i < v.size(); ++i) inv[i] = (v[i] == 0 ? 1 : 0);
        return inv;
    };

    // Work on a snapshot to avoid iteration issues when pushing new elements
    std::vector<std::vector<int>> original_properties = properties;
    for (const auto& cfg : original_properties) {
        std::vector<int> inv_cfg = invert(cfg);
        if (hash_set.insert(inv_cfg).second) {
            properties.push_back(inv_cfg);
        }
    }
}

// Placeholder implementation – to be fleshed out later
std::vector<std::vector<int>> ClusterProperty::get_children_properties() {
    std::vector<std::vector<int>> results;
    // Helper to ensure uniqueness across results
    std::unordered_set<std::vector<int>, VectorHash, VectorEqual> seen;

    auto mark_rotations_into_seen = [&seen](const std::vector<int>& vec){
        std::vector<int> rot = vec;
        size_t len = rot.size();
        for (size_t r=0; r<len; ++r) {
            seen.insert(rot);
            if(len > 1) {
                std::rotate(rot.begin() + 1, rot.begin() + 2, rot.end());
            }
        }
    };

    if (min_size == 5) {
        for (const auto &prop : properties) {
            if (prop.size() != 5) continue;
            for (int idx = 1; idx <= 4; ++idx) {
                std::vector<int> child = {prop[0], prop[idx]};
                if (seen.find(child) == seen.end()) {
                    results.push_back(child);
                    mark_rotations_into_seen(child);
                }
            }
        }
    } else if (min_size == 2) {
        for (const auto &prop : properties) {
            if (prop.size() != 2) continue;
            std::vector<int> child = {prop[0]};
            if (seen.find(child) == seen.end()) {
                results.push_back(child);
                mark_rotations_into_seen(child);
            }
        }
    }
    return results;
}

// Check if the local configuration adheres to the property
bool ClusterProperty::adhere(std::vector<int> local_config) {
    if (min_size == 5) {
        // For min_size 5, check if the full configuration exists in the hash set
        return hash_set.find(local_config) != hash_set.end();
    } else if (min_size == 1) {
        // For min_size 1, check if the first value matches
        if (local_config.empty()) return false;
        int first_value = local_config[0];
        std::vector<int> single_value = {first_value};
        return (hash_set.find(single_value) != hash_set.end()) ;
    }
    else if (min_size == 2) {
        // For min_size 2, check if the first value matches
        if (local_config.empty()) return false;
        int first_value = local_config[0];
        int second_value = local_config[1];
        int third_value = local_config[2];
        int fourth_value = local_config[3];
        int fifth_value = local_config[4];
        std::vector<int> first_config = {first_value, second_value};
        std::vector<int> second_config = {first_value, third_value};
        std::vector<int> third_config = {first_value, fourth_value};
        std::vector<int> fourth_config = {first_value, fifth_value};
        return (hash_set.find(first_config) != hash_set.end()) || (hash_set.find(second_config) != hash_set.end()) || (hash_set.find(third_config) != hash_set.end()) || (hash_set.find(fourth_config) != hash_set.end());
    }
    return false;
}
