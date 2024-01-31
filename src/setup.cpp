#include <iostream>
#include <vector>
#include <string>
#include <petsc.h>
#include <gmsh.h>
#include <compute.hpp>
#include <solver.hpp>
#include <pybind11/pybind11.h>

using namespace std;

// Load Gmsh script and generate square mesh
Mesh* generateMesh(int refinement){
   gmsh::merge("geometry/example.geo");
   gmsh::model::geo::synchronize();
   gmsh::model::mesh::generate(2);

   //Refine mesh uniformly
   for(int x = 0; x < refinement; x++)
      gmsh::model::mesh::refine();
   
   //Construct quadratic triangle elements
   gmsh::model::mesh::setOrder(2);

   gmsh::write("geometry/example.msh");

   return new Mesh("geometry/example.msh");
}

void computeFlow(int refinement){
   gmsh::initialize();

   Mesh *msh = generateMesh(refinement);
   
   //Vec vint = computeFirstStep(msh, msh->elementTags[0][100]);

   Solver* solver = new Solver(msh);
   solver->assembleMatrices();

   gmsh::fltk::run();
   gmsh::finalize();
}

PYBIND11_MODULE(MeshExtension, m) {
    m.def("computeFlow", &computeFlow, "A function to compute flow");
}

//Basic Petsc Code
int main(int argc, char **argv)
{
   PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
   computeFlow(1);

   PetscInt n = 5;
   Vec x;

   VecCreate(PETSC_COMM_WORLD, &x);
   VecSetSizes(x, PETSC_DECIDE, n);
   VecSetFromOptions(x);

   PetscReal value = 1.0;
   VecSet(x, value);
   PetscPrintf(PETSC_COMM_WORLD, "Vector x:\n");
   VecView(x, PETSC_VIEWER_STDOUT_WORLD);
   VecDestroy(&x);
   PetscFinalize();
   return 0;
}