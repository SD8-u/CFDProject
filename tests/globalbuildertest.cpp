#include <gtest/gtest.h>

#include "solver.hpp"

class GlobalBuilderTest : public testing::Test {
 protected:
  Mesh *msh;
  GlobalBuilder *build;

  static void SetUpTestSuite() { ::testing::internal::CaptureStdout(); }

  // Initialise mesh and Global Builder
  void SetUp() override {
    Mesh::generateMesh("geometry/lidcavity.geo", 4);
    msh = new Mesh("geometry/temp.msh", 1);
    build = new GlobalBuilder(2, 0.001, 0.1, msh);
  }

  // Clean up test objects
  void TearDown() {
    delete (msh);
    delete (build);
    gmsh::clear();
  }

  static void TearDownTestSuite() { ::testing::internal::GetCapturedStdout(); }
};

// Global Builder Unit Tests

TEST_F(GlobalBuilderTest, GlobalBuilderMatSize) {
  build->assembleMatrices();
  PetscInt row, col;
  MatGetSize(build->globalMassMat, &row, &col);
  ASSERT_EQ(row, msh->p2Size() * 2);
  ASSERT_EQ(col, msh->p2Size() * 2);
  MatGetSize(build->globalViscMat, &row, &col);
  ASSERT_EQ(row, msh->p2Size() * 2);
  ASSERT_EQ(col, msh->p2Size() * 2);
  MatGetSize(build->globalFullMat, &row, &col);
  ASSERT_EQ(row, msh->p2Size() * 2 + msh->p1Size());
  ASSERT_EQ(col, msh->p2Size() * 2 + msh->p1Size());
}

TEST_F(GlobalBuilderTest, GlobalBuilderConvSize) {
  build->assembleVectors();
  build->assembleConvectionMatrix();
  PetscInt row, col;
  MatGetSize(build->globalConvMat, &row, &col);
  ASSERT_EQ(row, msh->p2Size() * 2);
  ASSERT_EQ(col, msh->p2Size() * 2);
}

TEST_F(GlobalBuilderTest, GlobalBuilderVecSize) {
  build->assembleVectors();
  PetscInt size;
  VecGetSize(build->currVelVec, &size);
  ASSERT_EQ(size, msh->p2Size() * 2);
  VecGetSize(build->fullVec, &size);
  ASSERT_EQ(size, msh->p2Size() * 2 + msh->p1Size());
}

TEST_F(GlobalBuilderTest, GlobalBuilderVecLid) {
  build->assembleVectors();
  PetscReal max;
  VecMax(build->currVelVec, NULL, &max);
  ASSERT_EQ(max, 1);
  ASSERT_EQ(max, 1);
}

TEST_F(GlobalBuilderTest, GlobalBuilderVecWall) {
  build->assembleVectors();
  PetscReal min;
  VecMin(build->currVelVec, NULL, &min);
  ASSERT_EQ(min, 0);
  ASSERT_EQ(min, 0);
}

TEST_F(GlobalBuilderTest, GlobalBuilderNonZeroMat) {
  build->assembleVectors();
  build->assembleMatrices();
  build->assembleConvectionMatrix();
  double massVal, viscVal, convVal, fullVal;
  for (int i = 0; i < msh->p2Size() * 2 + msh->p1Size(); i++) {
    for (int j = 0; j < msh->p2Size() * 2 + msh->p1Size(); j++) {
      PetscScalar val;
      if (i < msh->p2Size() * 2 && j < msh->p2Size() * 2) {
        MatGetValue(build->globalMassMat, i, j, &val);
        massVal += val;
        MatGetValue(build->globalViscMat, i, j, &val);
        viscVal += val;
        MatGetValue(build->globalConvMat, i, j, &val);
        convVal += val;
      }
      MatGetValue(build->globalFullMat, i, j, &val);
      fullVal += val;
    }
  }
  ASSERT_NE(massVal, 0);
  ASSERT_NE(viscVal, 0);
  ASSERT_NE(convVal, 0);
  ASSERT_NE(fullVal, 0);
}

TEST_F(GlobalBuilderTest, GlobalBuilderNonZeroVec) {
  build->assembleVectors();
  double velVal, nodalVal;
  for (int i = 0; i < msh->p2Size() * 2 + msh->p1Size(); i++) {
    PetscScalar val;
    if (i < msh->p2Size() * 2) {
      VecGetValues(build->currVelVec, 1, &i, &val);
      velVal += val;
    }
    VecGetValues(build->fullVec, 1, &i, &val);
    nodalVal += val;
  }
  ASSERT_NE(velVal, 0);
  ASSERT_NE(nodalVal, 0);
}
