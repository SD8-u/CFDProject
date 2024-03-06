#include <gmsh.h>
#include <petsc.h>
#include <vector>
#include <string>
#include <map>
#include <iostream>

using namespace std;

//Node values for fluid domain
struct Node {
    int id;
    int pid;
    bool boundary;
    bool inlet;
    double force[2];
    double velocity[2];
    double pressure;
    double x, y; //2 Dimensional
};

class Mesh {
    public:
        int elementSize;
        int nNodes;
        int nLinear;
        double boundaryVel = 1.0;
        map<size_t, Node> nodes;
        map<int, size_t> nodeIds;
        vector<vector<size_t>> elementTags;
        map<size_t, vector<size_t>> elements;
        vector<PetscInt> dirichletIds;

        Mesh(string filePath, double boundaryVel);
        static void generateMesh(string filePath, int refinement);
        vector<vector<double>> getNodeCoords();
};