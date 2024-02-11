#include "setup.hpp"

using namespace std;

// Load Gmsh script and generate square mesh
Mesh* generateMesh(int refinement, double vel){
   gmsh::merge("geometry/Aneurysm15.vtk");
   gmsh::model::geo::synchronize();
   gmsh::model::mesh::generate(2);

   //Refine mesh uniformly
   for(int x = 0; x < refinement; x++)
      gmsh::model::mesh::refine();

   //Construct quadratic triangle elements
   gmsh::model::mesh::setOrder(2);
   gmsh::write("geometry/example.msh");

   return new Mesh("geometry/example.msh", vel);
}

pybind11::tuple computeFlow(int refinement, int steps, double vel, double dt, double visc){
   gmsh::initialize();
   PetscInitializeNoArguments();

   Mesh *msh = generateMesh(refinement, vel);
   Solver* solver = new Solver(msh, dt, visc);
   solver->assembleMatrices();
   vector<vector<double>> fluid = solver->computeTimeStep(steps);
   vector<vector<double>> coord = msh->getNodeCoords();

   pybind11::array_t<double> np_X(coord[0].size(), coord[0].data());
   pybind11::array_t<double> np_Y(coord[1].size(), coord[1].data());   

   pybind11::array_t<double> np_U(fluid[0].size(), fluid[0].data());
   pybind11::array_t<double> np_V(fluid[1].size(), fluid[1].data());
   pybind11::array_t<double> np_P(fluid[2].size(), fluid[2].data());

   pybind11::tuple result(5);
   result[0] = np_X;
   result[1] = np_Y;
   result[2] = np_U;
   result[3] = np_V;
   result[4] = np_P;

   //gmsh::fltk::run();
   PetscFinalize();
   gmsh::finalize();
   return result;
}

PYBIND11_MODULE(bloodflow, m) {
   m.doc() = "computes flow given (refinement, steps, velocity, dt, and visc)";
   m.def("computeFlow", &computeFlow, "A function to compute flow");
}

int main(int argc, char **argv)
{
   gmsh::initialize();
   generateMesh(1, 1);
   gmsh::fltk::run();
   gmsh::finalize();

   return 0;
}