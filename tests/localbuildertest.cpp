#include <gtest/gtest.h>
#include "solver.hpp"

class LocalBuilderTest : public testing::Test {
    protected:
        Mesh *msh;
        LocalBuilder *build;

        static void SetUpTestSuite() {
            ::testing::internal::CaptureStdout();
        }

        void SetUp() override {
            msh = new Mesh("geometry/example.geo", 4, 1);
        }

        void TearDown() {
            delete(msh);
            delete(build);
            gmsh::clear();
        }

        static void TearDownTestSuite(){
            ::testing::internal::GetCapturedStdout();
        }
};

TEST_F(LocalBuilderTest, LocalBuilderMatSize) {
    build = new LocalBuilder(0.01);
    build->assembleMatrices(msh->elementTags[0][0], 0);
    PetscInt row, col;
    MatGetSize(build->localMassMat, &row, &col);
    ASSERT_EQ(row, 12);
    ASSERT_EQ(col, 12);
    MatGetSize(build->localViscMat, &row, &col);
    ASSERT_EQ(row, 12);
    ASSERT_EQ(col, 12);
    MatGetSize(build->localFullMat, &row, &col);
    ASSERT_EQ(row, 15);
    ASSERT_EQ(col, 15);
}

TEST_F(LocalBuilderTest, LocalBuilderConvSize) {
    build = new LocalBuilder();
    Vec velVec;
    VecCreate(PETSC_COMM_WORLD, &velVec);
    VecSetSizes(velVec, PETSC_DECIDE, 12);
    VecSetFromOptions(velVec);
    VecZeroEntries(velVec);
    VecAssemblyBegin(velVec);
    VecAssemblyEnd(velVec);

    build->computeConvectionMatrix(msh->elementTags[0][0], &velVec);
    PetscInt row, col;
    MatGetSize(build->localConvMat, &row, &col);
    ASSERT_EQ(row, 12);
    ASSERT_EQ(col, 12);
    VecDestroy(&velVec);
}

TEST_F(LocalBuilderTest, LocalBuilderNonZeroMat) {
    build = new LocalBuilder(0.01);
    build->assembleMatrices(msh->elementTags[0][0], 0);
    double massVal = 0, viscVal = 0, fullVal = 0;
    for(int i = 0; i < 15; i++){
        for(int j = 0; j < 15; j++){
            PetscScalar val;
            if(i < 12 && j < 12){
                MatGetValue(build->localMassMat, i, j, &val);
                massVal += val;
                MatGetValue(build->localViscMat, i, j, &val);
                viscVal += val;
            }
            MatGetValue(build->localFullMat, i, j, &val);
            fullVal += val;
        }
    }
    MatView(build->localViscMat, PETSC_VIEWER_STDOUT_WORLD);
    ASSERT_NE(massVal, 0);
    ASSERT_NE(viscVal, 0);
    ASSERT_NE(fullVal, 0);
}


TEST_F(LocalBuilderTest, LocalBuilderNonZeroConvMat) {
    build = new LocalBuilder();
    Vec velVec;
    VecCreate(PETSC_COMM_WORLD, &velVec);
    VecSetSizes(velVec, PETSC_DECIDE, 12);
    VecSetFromOptions(velVec);
    for(int i = 0; i < 12; i++){
        VecSetValue(velVec, i, 1, INSERT_VALUES);
    }
    VecAssemblyBegin(velVec);
    VecAssemblyEnd(velVec);
    build->computeConvectionMatrix(msh->elementTags[0][0], &velVec);
    double convVal = 0;
    for(int i = 0; i < 12; i++){
        for(int j = 0; j < 12; j++){
            PetscScalar val;
            MatGetValue(build->localConvMat, i, j, &val);
            convVal += val;
        }
    }
    ASSERT_NE(convVal, 0);
    VecDestroy(&velVec);
}

