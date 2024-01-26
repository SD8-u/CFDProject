#include <gmsh.h>
#include <vector>
#include <string>
#include <map>
#include <iostream>

using namespace std;

//Node values for fluid domain
struct Node {
    size_t id;
    bool boundary;
    bool inlet;
    double force[2];
    double velocity[2];
    double pressure;
    double x, y; //2 Dimensional
};

class Mesh {
    private:
        map<size_t, Node> nodes;
        vector<vector<size_t>> elements;
    public:
        int elementSize;
        vector<vector<size_t>> elementTags;
        Mesh(string filePath);
};