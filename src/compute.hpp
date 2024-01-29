#include <gmsh.h>
#include <petsc.h>
#include <iostream>
#include <vector>
#include <mesh.hpp>

using namespace std;

Mat computeMassMatrix(size_t elementTag);
Mat computeViscosityMatrix(size_t elementTag);
Mat computeConvectionMatrix(vector<double> elemVec, size_t elementTag);
Mat computeGradientMatrix(size_t elementTag);
Vec computeFirstStep(Mesh *msh, size_t elementTag);