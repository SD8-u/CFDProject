#include <gmsh.h>
#include <petsc.h>
#include <iostream>
#include <vector>

#ifndef MESH_HPP
#define MESH_HPP
#include "mesh.hpp"
#endif

using namespace std;

Mat computeMassMatrix(size_t elementTag, bool inverse, bool full);
Mat computeViscosityMatrix(size_t elementTag);
Mat computeConvectionMatrix(size_t elementTag);
Mat computeGradientMatrix(size_t elementTag);
Mat computeFinalMatrix(size_t elementTag, double dt);