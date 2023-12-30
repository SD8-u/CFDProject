#include <gmsh.h>
#include <vector>
#include <string>
#include <iostream>

using namespace std;

//Node values for fluid domain
struct Node {
    size_t id;
    bool velocity;
    double value;
    double coords[2]; //2 dimensional
};

//Local system triangle elements
class TriangleElement{
    public:
        Node nodes[9];
};

class Mesh{
    private:
        vector<vector<Node>> nodes;
    public:
        int elementSize;
        Mesh(string filePath){
            gmsh::initialize();
            gmsh::open(filePath);
            
            vector<size_t> nodeIds;
            vector<double> nodeCoords;
            vector<double> paramCoords;
            gmsh::model::mesh::getNodes(nodeIds, nodeCoords, paramCoords, 2, -1 , true, true);

            for(int node = 0; node < nodeCoords.size()/3; node++){
                cout << "Node ID: " << nodeIds[node] << endl;
                cout << "Coord: " << nodeCoords[node * 3] << " : " << nodeCoords[node * 3 + 1] << endl;
                cout << "Size: " << nodeCoords.size() << endl;
            }

            vector<int> elementTypes;
            vector<vector<size_t>> elementTags, elementNodeTags;
            vector<double> nodeCoords1, paramCoords1;
            gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags, 2, -1);

            gmsh::finalize();
        }
};
