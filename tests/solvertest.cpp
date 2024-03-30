#include <gtest/gtest.h>
#include <omp.h>

#include "solver.hpp"

class SolverTest : public testing::Test {
 protected:
  Mesh *msh;
  Solver *solver;

  static void SetUpTestSuite() { ::testing::internal::CaptureStdout(); }

  void SetUp() override {
    Mesh::generateMesh("geometry/lidcavity.geo", 4);
    msh = new Mesh("geometry/temp.msh", 1);
    solver = new Solver(msh, 0.0001, 0.1);
  }

  void TearDown() {
    delete (msh);
    delete (solver);
    gmsh::clear();
  }

  static void TearDownTestSuite() {
    gmsh::finalize();
    PetscFinalize();

    ::testing::internal::GetCapturedStdout();
  }
};

TEST_F(SolverTest, SolverDataSize) {
  solver->computeTimeSteps(2);
  vector<vector<double>> sol = solver->interpolateSolution(0.1, 0);
  ASSERT_EQ(sol.size(), 5);
}

TEST_F(SolverTest, SolverBoundaryVelX) {
  solver->computeTimeSteps(2);
  vector<vector<double>> sol = solver->interpolateSolution(0.01, 0);
  int avg = 0;
  for (double vel : sol[1]) avg += vel;
  avg /= sol[1].size();
  ASSERT_LE(avg, 1);
}

TEST_F(SolverTest, SolverBoundaryVelY) {
  solver->computeTimeSteps(2);
  vector<vector<double>> sol = solver->interpolateSolution(0.01, 0);
  int avg = 0;
  for (double vel : sol[2]) avg += vel;
  avg /= sol[1].size();
  ASSERT_LE(avg, 1);
}

int main(int argc, char **argv) {
  printf("Running all tests...\n");
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
