#include <gtest/gtest.h>
#include "solver.hpp"

class MeshTest : public testing::Test {
    protected:
        Mesh *msh;

        static void SetUpTestSuite() {
            gmsh::initialize();
            PetscInitializeNoArguments();
            ::testing::internal::CaptureStdout();
        }

        void TearDown() {
            delete(msh);
            gmsh::clear();
        }

        static void TearDownTestSuite(){
            ::testing::internal::GetCapturedStdout();
        }
};

TEST_F(MeshTest, MeshElementSize1) {
    msh = new Mesh("geometry/example.geo", 1, 1);
    ASSERT_EQ(msh->elementSize, 16);
}

TEST_F(MeshTest, MeshElementSize2) {
    msh = new Mesh("geometry/example.geo", 2, 1);
    ASSERT_EQ(msh->elementSize, 64);
}

TEST_F(MeshTest, MeshElementSize3) {
    msh = new Mesh("geometry/example.geo", 3, 1);
    ASSERT_EQ(msh->elementSize, 256);
}

TEST_F(MeshTest, MeshElementSize4) {
    msh = new Mesh("geometry/example.geo", 4, 1);
    ASSERT_EQ(msh->elementSize, 1024);
}

TEST_F(MeshTest, MeshElementSize5) {
    msh = new Mesh("geometry/example.geo", 5, 1);
    ASSERT_EQ(msh->elementSize, 4096);
}

TEST_F(MeshTest, MeshElementSize6) {
    msh = new Mesh("geometry/example.geo", 6, 1);
    ASSERT_EQ(msh->elementSize, 16384);
}

TEST_F(MeshTest, NodeTest) {
    msh = new Mesh("geometry/example.geo", 4, 1);
    vector<double> coord;
    vector<double> paramCoord;
    for(int i = 0; i < msh->nNodes; i++){
        gmsh::model::mesh::getNode(msh->nodeIds[i], coord, paramCoord);
        ASSERT_EQ(coord[0], msh->nodes[msh->nodeIds[i]].x);
        ASSERT_EQ(coord[1], msh->nodes[msh->nodeIds[i]].y);
    }
}

TEST_F(MeshTest, ElementTest) {
    msh = new Mesh("geometry/example.geo", 4, 1);
    for(int i = 0; i < msh->elements.size(); i++){
        vector<size_t> nodeTags;
        int elementType;
        gmsh::model::mesh::getElement(msh->elementTags[0][i], elementType, nodeTags);
        for(int j = 0; j < nodeTags.size(); j++) {
            ASSERT_EQ(msh->elements[msh->elementTags[0][i]][j], nodeTags[j]);
        } 
    }
}
