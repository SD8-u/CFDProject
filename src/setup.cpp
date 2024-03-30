#include "setup.hpp"

pybind11::tuple computeFlow(int refinement, int steps, double vel, double dt,
                            double visc, string file) {
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("Processor Rank %d out of %d\n", rank, size);

  pybind11::tuple result(5);
  result[0] = NULL;
  result[1] = NULL;
  result[2] = NULL;
  result[3] = NULL;
  result[4] = NULL;

  if (rank == 0) Mesh::generateMesh(file, refinement);

  MPI_Barrier(MPI_COMM_WORLD);

  double start = MPI_Wtime();
  Mesh *msh = new Mesh("geometry/temp.msh", vel);
  Solver *solver = new Solver(msh, dt, visc);
  solver->computeTimeSteps(steps);
  double stop = MPI_Wtime() - start;
  vector<vector<double>> solData = solver->interpolateSolution(0.005, rank);

  if (rank == 0) {
    pybind11::array_t<double> np_X(solData[3].size(), solData[3].data());
    pybind11::array_t<double> np_Y(solData[4].size(), solData[4].data());
    pybind11::array_t<double> np_U(solData[1].size(), solData[1].data());
    pybind11::array_t<double> np_V(solData[2].size(), solData[2].data());
    pybind11::array_t<double> np_P(solData[0].size(), solData[0].data());

    result[0] = np_X;
    result[1] = np_Y;
    result[2] = np_U;
    result[3] = np_V;
    result[4] = np_P;
  }

  delete (solver);
  delete (msh);
  gmsh::clear();
  cout << "Time: " << stop << "\n";
  return result;
}

void startUp() {
  gmsh::initialize();
  PetscInitializeNoArguments();
}

void cleanUp() {
  gmsh::finalize();
  PetscFinalize();
}

PYBIND11_MODULE(cfd, m) {
  m.doc() = "computes flow given (refinement, steps, velocity, dt, and visc)";
  m.def("computeFlow", &computeFlow, "A function to compute flow");
  m.def("startUp", &startUp, "A function to allocate resources");
  m.def("cleanUp", &cleanUp, "A function to free resources");
}

void computeFlowC(int refinement, int steps, double vel, double dt,
                  double visc) {
  PetscInitializeNoArguments();
  gmsh::initialize();

  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("Processor Rank %d out of %d\n", rank, size);

  if (rank == 0) Mesh::generateMesh("geometry/lidcavity.geo", refinement);

  MPI_Barrier(MPI_COMM_WORLD);

  double start = MPI_Wtime();
  Mesh *msh = new Mesh("geometry/temp.msh", vel);
  Solver *solver = new Solver(msh, dt, visc);

  solver->computeTimeSteps(steps);
  double stop = MPI_Wtime() - start;
  vector<vector<double>> solData = solver->interpolateSolution(0.02, rank);

  delete (solver);
  delete (msh);

  cout << "Time: " << stop << "\n";
  gmsh::finalize();
  PetscFinalize();
}

int main(int argc, char **argv) {
  computeFlowC(5, 1, 100, 0.001, 10);
  return 0;
}
