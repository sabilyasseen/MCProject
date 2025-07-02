cmake_minimum_required(VERSION 3.10)
project(ResearchCpp)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add source files
set(SOURCES
    main.cpp
    Sim.cpp
    Grid.cpp
    Cell.cpp
    Cluster.cpp
    ClusterProperty.cpp
    Move.cpp
    SimTest.cpp
)

# Add header files
set(HEADERS
    Sim.hpp
    Grid.hpp
    Cell.hpp
    Cluster.hpp
    ClusterProperty.hpp
    Move.hpp
    globals.hpp
    SimTest.hpp
)

# Create executable
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

# Include directories
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}) 