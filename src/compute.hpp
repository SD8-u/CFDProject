#include <gmsh.h>
#include <petsc.h>
#include <iostream>
#include <vector>
#include <mesh.hpp>

using namespace std;

Mat computeMassMatrix(size_t elementTag);