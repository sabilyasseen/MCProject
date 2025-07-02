#ifndef _CLUSTERPROPERTY_HPP
#define _CLUSTERPROPERTY_HPP
#include "globals.hpp"
#include "Cell.hpp"
#include <vector>
#include <unordered_set>
#include <functional> // For std::hash

class ClusterProperty {
public:
    ClusterProperty();
    ClusterProperty(int minsize);
    ~ClusterProperty() = default;
    bool homogenous = true;

    std::vector<std::vector<int>> properties;
    int min_size;

    int get_minsize();
    void add_config(std::vector<int> property);
    bool adhere(std::vector<int> local_config);
    void generate_configurations(std::vector<int> property);
    void generate_inversions();
    std::vector<std::vector<int>> get_children_properties();

    std::unordered_set<std::vector<int>, VectorHash, VectorEqual> hash_set;
};

#endif