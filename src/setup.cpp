#include "setup.hpp"
#include <omp.h>

using namespace std;

pybind11::tuple computeFlow(int refinement, int steps, double vel, double dt, double visc){
   gmsh::initialize();
   PetscInitializeNoArguments();

   Mesh *msh = new Mesh("geometry/example.geo", refinement, vel);
   Solver* solver = new Solver(msh, dt, visc);

   solver->computeTimeStep(steps);
   vector<vector<double>> solData = solver->interpolateSolution(0.03);

   pybind11::array_t<double> np_X(solData[3].size(), solData[3].data());
   pybind11::array_t<double> np_Y(solData[4].size(), solData[4].data());   

   pybind11::array_t<double> np_U(solData[1].size(), solData[1].data());
   pybind11::array_t<double> np_V(solData[2].size(), solData[2].data());
   pybind11::array_t<double> np_P(solData[0].size(), solData[0].data());

   pybind11::tuple result(5);
   result[0] = np_X;
   result[1] = np_Y;
   result[2] = np_U;
   result[3] = np_V;
   result[4] = np_P;

   delete(solver);

   gmsh::clear();
   return result;
}

PYBIND11_MODULE(bloodflow, m) {
   m.doc() = "computes flow given (refinement, steps, velocity, dt, and visc)";
   m.def("computeFlow", &computeFlow, "A function to compute flow");
}

void computeFlowC(int refinement, int steps, double vel, double dt, double visc){
   PetscInitializeNoArguments();
   gmsh::initialize();

   Mesh *msh = new Mesh("geometry/example.geo", refinement, vel);
   Solver* solver = new Solver(msh, dt, visc);

   solver->computeTimeStep(steps);
   vector<vector<double>> solData = solver->interpolateSolution(0.03);

   delete(solver);
   gmsh::finalize();
   PetscFinalize();
}

int main(int argc, char **argv)
{
   computeFlowC(5, 1, 100, 0.001, 10);
   return 0;
}