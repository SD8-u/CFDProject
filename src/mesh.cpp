#include "mesh.hpp"
#include <unordered_set>

// Load Gmsh script and generate msh file
void generateMesh(string filePath, int refinement, double vel){
    gmsh::merge(filePath);
    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate(2);
    
    //Refine mesh uniformly
    for(int x = 0; x < refinement; x++)
        gmsh::model::mesh::refine();
    
    //Construct quadratic triangle elements
    gmsh::model::mesh::setOrder(2);
    
    gmsh::write("geometry/example.msh");
    gmsh::open("geometry/example.msh");
}

Mesh::Mesh(string filePath, int refinement, double boundaryVel){

    generateMesh(filePath, refinement, boundaryVel);
    this->boundaryVel = boundaryVel;

    vector<size_t> nodeTags;
    vector<double> nodeCoords;
    vector<double> paramCoords;
    int nId = 0;
    int pId = 0;

    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, paramCoords);
    nNodes = nodeTags.size();
    gmsh::vectorpair dimTags;
    gmsh::model::getPhysicalGroups(dimTags);

    for(auto dimTag : dimTags){
        string name;
        gmsh::model::getPhysicalName(dimTag.first, dimTag.second, name);
        cout << "Name: " << name << " Tag: " << dimTag.second << " Dim: " << dimTag.first << "\n";
    }

    //Retrieve nodes at inlet
    int inletTag = 2;
    gmsh::model::mesh::getNodesForPhysicalGroup(1, inletTag, nodeTags, nodeCoords);
    for(int node = 0; node < nodeCoords.size()/3; node++){
        dirichletIds.push_back(nId);
        dirichletIds.push_back(nId + nNodes);
        nodeIds[nId] = nodeTags[node];
        nodes[nodeTags[node]] = Node{nId++, -1, false, true, 
        {0}, {node == 0 || node == ((nodeCoords.size() - 3)/3) ? 0 : boundaryVel, 0}, 
        0, nodeCoords[node * 3], nodeCoords[node * 3 + 1]};
    }

    //Retreive nodes in boundary
    int boundaryTag = 1;
    gmsh::model::mesh::getNodesForPhysicalGroup(1, boundaryTag, nodeTags, nodeCoords);
    for(int node = 0; node < nodeCoords.size()/3; node++){
        if(nodes.find(nodeTags[node]) == nodes.end()){
            nodeIds[nId] = nodeTags[node];
            dirichletIds.push_back(nId);
            dirichletIds.push_back(nId + nNodes);
            nodes[nodeTags[node]] = Node{nId++, -1, true, false, 
            {0}, {0, 0}, 0, nodeCoords[node * 3], nodeCoords[node * 3 + 1]};
        }
    }
    
    //Retrieve all other nodes
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, paramCoords);
    for(int node = 0; node < nodeCoords.size()/3; node++){
        if(nodes.find(nodeTags[node]) == nodes.end()){
            nodeIds[nId] = nodeTags[node];
            nodes[nodeTags[node]] = Node{nId++, -1, false, false, 
            {0}, {0, 0}, 0, nodeCoords[node * 3], nodeCoords[node * 3 + 1]};
        }
    }
    
    //Encode element connectivity
    vector<int> elementTypes;
    vector<vector<size_t>> elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags, 2, -1);
    unordered_set<size_t> pNodes;

    this->elementSize = elementTags[0].size();
    for(int element = 0; element < elementTags[0].size(); element++){
        int type;
        vector<size_t> nodeTags;
        gmsh::model::mesh::getElement(elementTags[0][element], type, nodeTags);

        for(int node = 0; node < nodeTags.size(); node++){
            if(node < 3){
                if(pNodes.find(nodeTags[node]) == pNodes.end()){
                    nodes[nodeTags[node]].pid = pId++;
                    pNodes.emplace(nodeTags[node]);
                }
            }
            elements[elementTags[0][element]].push_back(nodeTags[node]);
        }
    }
    nLinear = pId;
    cout << "Linear Nodes: " << nLinear << "\n";
}

//Returns vector of coordinates for each node
vector<vector<double>> Mesh::getNodeCoords(){
    vector<vector<double>> coord = vector<vector<double>>(2);
    for(int i = 0; i < nNodes; i++){
        coord[0].push_back(nodes[nodeIds[i]].x);
        coord[1].push_back(nodes[nodeIds[i]].y);
    }
    return coord;
}

