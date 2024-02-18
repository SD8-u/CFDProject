cmake_minimum_required(VERSION 3.14)
project(bloodFlowCFD)

set(CMAKE_CXX_COMPILER "/usr/bin/mpic++")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set the directory containing the source files
set(SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Set the directory for the build output
set(OUTDIR ${CMAKE_CURRENT_SOURCE_DIR}/build)

# Set the directory for the module output
set(MODDIR ${CMAKE_CURRENT_SOURCE_DIR}/scripts)

# Collect all source files
file(GLOB SOURCES ${SRCDIR}/*.cpp)

# Create the executable target
add_executable(main ${SOURCES})

set(CMAKE_PREFIX_PATH "venv/lib/python3.10/site-packages/pybind11/share/cmake/pybind11")
find_package(pybind11 REQUIRED)

# Set the include directories for the c++ executable
target_include_directories(main PRIVATE
    /usr/include/petsc
    /usr/include
    /usr/include/python3.10
    ${CMAKE_CURRENT_SOURCE_DIR}/venv/lib/python3.10/site-packages/pybind11/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Link the required libraries for the c++ executable
target_link_libraries(main PRIVATE
    petsc
    stdc++
    gmsh
    python3.10
)

# Set compiler flags for the c++ executable
target_compile_options(main PRIVATE
    -fdiagnostics-color=always
    -g
)

# Set output directory for the c++ executable
set_target_properties(main PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${OUTDIR}
)

# Optional: Create directory for build output
file(MAKE_DIRECTORY ${OUTDIR})

# Create the python bloodflow module target
pybind11_add_module(bloodflow ${SOURCES})

# Set include directories for the bloodflow module
target_include_directories(bloodflow PRIVATE
    /usr/include/petsc
    /usr/include
    /usr/include/python3.10
    ${CMAKE_CURRENT_SOURCE_DIR}/venv/lib/python3.10/site-packages/pybind11/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Link the required libraries for the bloodflow module
target_link_libraries(bloodflow PRIVATE
    petsc
    stdc++
    gmsh
    python3.10
)

# Set output directory for the bloodflow module
set_target_properties(bloodflow PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${MODDIR}
)

# Locate GTest
find_package(GTest REQUIRED)
include_directories(
    ${GTEST_INCLUDE_DIRS}
    /usr/include/petsc
    /usr/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

add_executable(runTests 
    tests/globalbuildertest.cpp
    tests/solvertest.cpp 
    src/globalbuilder.cpp 
    src/localbuilder.cpp 
    src/mesh.cpp 
    src/solver.cpp
)

target_link_libraries(runTests PRIVATE 
    petsc
    stdc++
    gmsh
    ${GTEST_LIBRARIES} 
)