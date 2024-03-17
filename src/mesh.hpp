#ifndef SRC_MESH_HPP_
#define SRC_MESH_HPP_

#include <gmsh.h>
#include <petsc.h>

#include <iostream>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>

using std::cout;
using std::map;
using std::string;
using std::unordered_set;
using std::vector;

// Node values for fluid domain
struct Node {
  int id;
  int pid;
  bool boundary;
  bool inlet;
  double force[2];
  double velocity[2];
  double pressure;
  double x, y;  // 2 Dimensional
};

class Mesh {
 private:
  int nP1;
  int nP2;
  int nElements;
  int nDirichlet;
  double boundaryVel = 1.0;

  map<size_t, Node> nodes;
  map<int, size_t> nodeIds;
  vector<PetscInt> dirichletIds;
  map<size_t, vector<size_t>> elements;

  void getNodes(int* nId, vector<size_t> nodeTags, vector<double> nodeCoords,
                double boundaryVel, bool boundary, bool inlet);
  void getElementConnectivity();

 public:
  vector<vector<size_t>> elementTags;

  static void generateMesh(string filePath, int refinement);
  Mesh(string filePath, double boundaryVel);

  int p1Size();
  int p2Size();
  int elementSize();
  int dirichletSize();

  Node getNode(int i);
  Node getNode(size_t nodeTag);
  Node getNode(size_t elementTag, int i);

  size_t getElement(int e);
  vector<size_t> getElements();

  PetscInt* getDirichlet();
};

#endif  // SRC_MESH_HPP_
