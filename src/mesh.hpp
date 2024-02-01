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
        int nLinear;
        map<size_t, Node> nodes;
        vector<vector<size_t>> elementTags;
        map<size_t, vector<size_t>> elements;

        Mesh(string filePath);
        void getElementVector(size_t element, PetscScalar* elemVec, bool pressure);
        void getForceVector(size_t element, PetscScalar* force);
};