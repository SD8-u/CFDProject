#include <gtest/gtest.h>

#include "solver.hpp"

class MeshTest : public testing::Test {
 protected:
  Mesh *msh;

  // Initialise petsc and gmsh for Mesh
  static void SetUpTestSuite() {
    gmsh::initialize();
    PetscInitializeNoArguments();
    ::testing::internal::CaptureStdout();
  }

  // Clean up test objects
  void TearDown() {
    delete (msh);
    gmsh::clear();
  }

  static void TearDownTestSuite() { ::testing::internal::GetCapturedStdout(); }
};

// Mesh Unit Tests

TEST_F(MeshTest, MeshElementSize1) {
  Mesh::generateMesh("geometry/lidcavity.geo", 1);
  msh = new Mesh("geometry/temp.msh", 1);
  ASSERT_EQ(msh->elementSize(), 16);
}

TEST_F(MeshTest, MeshElementSize2) {
  Mesh::generateMesh("geometry/lidcavity.geo", 2);
  msh = new Mesh("geometry/temp.msh", 1);
  ASSERT_EQ(msh->elementSize(), 64);
}

TEST_F(MeshTest, MeshElementSize3) {
  Mesh::generateMesh("geometry/lidcavity.geo", 3);
  msh = new Mesh("geometry/temp.msh", 1);
  ASSERT_EQ(msh->elementSize(), 256);
}

TEST_F(MeshTest, MeshElementSize4) {
  Mesh::generateMesh("geometry/lidcavity.geo", 4);
  msh = new Mesh("geometry/temp.msh", 1);
  ASSERT_EQ(msh->elementSize(), 1024);
}

TEST_F(MeshTest, MeshElementSize5) {
  Mesh::generateMesh("geometry/lidcavity.geo", 5);
  msh = new Mesh("geometry/temp.msh", 1);
  ASSERT_EQ(msh->elementSize(), 4096);
}

TEST_F(MeshTest, MeshElementSize6) {
  Mesh::generateMesh("geometry/lidcavity.geo", 6);
  msh = new Mesh("geometry/temp.msh", 1);
  ASSERT_EQ(msh->elementSize(), 16384);
}

TEST_F(MeshTest, NodeTest) {
  Mesh::generateMesh("geometry/lidcavity.geo", 4);
  msh = new Mesh("geometry/temp.msh", 1);
  vector<double> coord;
  vector<double> paramCoord;
  vector<size_t> nodeTags;
  gmsh::model::mesh::getNodes(nodeTags, coord, paramCoord);
  for (int i = 0; i < nodeTags.size(); i++) {
    Node node = msh->getNode(nodeTags[i]);
    ASSERT_EQ(coord[i * 3], node.x);
    ASSERT_EQ(coord[i * 3 + 1], node.y);
  }
}

TEST_F(MeshTest, ElementTest) {
  Mesh::generateMesh("geometry/lidcavity.geo", 4);
  msh = new Mesh("geometry/temp.msh", 1);

  for (int i = 0; i < msh->elementSize(); i++) {
    vector<size_t> nodeTags;
    int elementType;
    gmsh::model::mesh::getElement(msh->getElement(i), elementType, nodeTags);

    for (int j = 0; j < nodeTags.size(); j++) {
      Node n1 = msh->getNode(msh->getElement(i), j);
      Node n2 = msh->getNode(nodeTags[j]);
      ASSERT_EQ(n1.id, n2.id);
    }
  }
}
