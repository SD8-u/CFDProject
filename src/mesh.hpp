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
  double boundaryVel = 1.0;
  void getNodes(int *nId, vector<size_t> nodeTags, vector<double> nodeCoords,
                double boundaryVel, bool boundary, bool inlet);
  void getElementConnectivity();

 public:
  int elementSize;
  int nNodes;
  int nLinear;

  map<size_t, Node> nodes;
  map<int, size_t> nodeIds;
  vector<vector<size_t>> elementTags;
  map<size_t, vector<size_t>> elements;
  vector<PetscInt> dirichletIds;

  Mesh(string filePath, double boundaryVel);
  static void generateMesh(string filePath, int refinement);
};

#endif  // SRC_MESH_HPP_
