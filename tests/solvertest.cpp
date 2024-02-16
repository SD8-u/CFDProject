#include <gtest/gtest.h>
#include "solver.hpp"

class SolverTest : public testing::Test {
    protected:
        Mesh *msh;
        Solver *solver;

        static void SetUpTestSuite() {
            //::testing::internal::CaptureStdout();

            gmsh::initialize();
            PetscInitializeNoArguments();

        }

        void SetUp() override {
            msh = new Mesh("geometry/example.geo", 4, 1);
            solver = new Solver(msh, 0.0001, 0.1);
        }

        void TearDown() {
            delete(msh);
            delete(solver);
            gmsh::clear();
        }

        static void TearDownTestSuite(){
            gmsh::finalize();
            PetscFinalize();

            //::testing::internal::GetCapturedStdout();
        }
};

TEST_F(SolverTest, SolverDataSize) {
    solver->computeTimeStep(2);
    vector<vector<double>> sol = solver->interpolateSolution(0.1);
    ASSERT_EQ(sol.size(), 5);
}

TEST_F(SolverTest, SolverBoundaryVelX) {
    solver->computeTimeStep(2);
    vector<vector<double>> sol = solver->interpolateSolution(0.01);
    int avg = 0;
    for(double vel : sol[1])
        avg+=vel;
    avg/=sol[1].size();
    ASSERT_LE(avg, 1);
}

TEST_F(SolverTest, SolverBoundaryVelY) {
    solver->computeTimeStep(2);
    vector<vector<double>> sol = solver->interpolateSolution(0.01);
    int avg = 0;
    for(double vel : sol[2])
        avg+=vel;
    avg/=sol[1].size();
    ASSERT_LE(avg, 1);
}

int main(int argc, char **argv) {
  printf("Running solver tests from %s\n", __FILE__);
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
