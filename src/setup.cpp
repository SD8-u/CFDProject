#include "setup.hpp"

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

void computeFlow(int refinement, int steps){
   gmsh::initialize();
   PetscInitializeNoArguments();

   Mesh *msh = generateMesh(refinement);
   Solver* solver = new Solver(msh);
   solver->assembleMatrices();
   vector<vector<double>> fluid = solver->computeTimeStep(steps);

   gmsh::fltk::run();
   PetscFinalize();
   gmsh::finalize();
}

PYBIND11_MODULE(bloodflow, m) {
    m.def("computeFlow", &computeFlow, "A function to compute flow");
}

//Basic Petsc Code
int main(int argc, char **argv)
{
   PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
   computeFlow(3, 200);

   PetscFinalize();
   return 0;
}