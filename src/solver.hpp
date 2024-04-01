#ifndef SRC_SOLVER_HPP_
#define SRC_SOLVER_HPP_

#include <petsc.h>

#include <vector>

#include "globalbuilder.hpp"
#include "mesh.hpp"

class Solver {
 private:
  double dt, viscosity;
  KSP stp1Solver, stp2Solver;
  IS vecMapping;
  VecScatter vecScatter1, vecScatter2;
  Mesh *msh;
  GlobalBuilder *globalBuild;

  void initialiseSolver();
  void applyDirichletConditions(Mat *m, Vec *v, bool full);
  void interpolateValues(vector<double> *coord, vector<vector<double>> *solData,
                         vector<size_t> *nodeTags, Vec *solVec);
  void updateVectors(Vec *vec1, Vec *vec2, bool vel);
  void bdf2(Mat *tempMat, Vec *tempVec, bool euler);
  void computeFirstStep(bool euler);
  void computeSecondStep();

 public:
  Solver(Mesh *msh, double dt, double visc);
  ~Solver();

  void computeTimeSteps(int steps);
  vector<vector<double>> interpolateSolution(double resolution, int rank);
};

#endif  // SRC_SOLVER_HPP_
