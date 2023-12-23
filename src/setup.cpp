#include <iostream>
#include <vector>
#include <string>
#include <petsc.h>
#include <gmsh.h>
#include <pybind11/pybind11.h>

using namespace std;

// Load Gmsh script and generate mesh
void generateMesh(){
   gmsh::initialize();
   gmsh::merge("geometry/example.geo");
   gmsh::model::geo::synchronize();
   gmsh::model::mesh::generate(2);
   gmsh::write("geometry/example.msh");
   gmsh::finalize();
}

PYBIND11_MODULE(MeshExtension, m) {
    m.def("generateMesh", &generateMesh, "A function to generate a basic mesh");
}

//Basic Petsc Code
int main(int argc, char **argv)
{
   PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

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