#include <mesh.hpp>

Mesh::Mesh(string filePath){
    gmsh::open(filePath);
    
    vector<size_t> nodeIds;
    vector<double> nodeCoords;
    vector<double> paramCoords;
    //Retreive nodes in boundary
    for(int boundaryTag = 2; boundaryTag < 5; boundaryTag++){
        gmsh::model::mesh::getNodesForPhysicalGroup(2, boundaryTag, nodeIds, nodeCoords);
        for(int node = 0; node < nodeCoords.size()/3; node++){
            nodes[nodeIds[node]] = Node{nodeIds[node], true, false, 
            {0}, {0, 0}, 0, nodeCoords[node * 3], nodeCoords[node * 3 + 1]};
        }
    }
    
    //Retrieve nodes at inlet
    int inletTag = 1;
    gmsh::model::mesh::getNodesForPhysicalGroup(2, inletTag, nodeIds, nodeCoords);
    for(int node = 0; node < nodeCoords.size()/3; node++){
        nodes[nodeIds[node]] = Node{nodeIds[node], false, true, 
        {0}, {1, 0}, 0, nodeCoords[node * 3], nodeCoords[node * 3 + 1]};
    }
    
    //Retrieve fluid nodes
    gmsh::model::mesh::getNodes(nodeIds, nodeCoords, paramCoords, 2);
    for(int node = 0; node < nodeCoords.size()/3; node++){
        nodes[nodeIds[node]] = Node{nodeIds[node], false, false, 
        {0}, {1, 0}, 0, nodeCoords[node * 3], nodeCoords[node * 3 + 1]};
    
    }
    
    //Encode element connectivity
    vector<int> elementTypes;
    vector<vector<size_t>> elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags, 2, -1);
    //elements = vector<vector<size_t>>(elementTags[0].size());
    for(int element = 0; element < elementTags[0].size(); element++){
        int type;
        vector<size_t> nodeTags;
        gmsh::model::mesh::getElement(elementTags[0][element], type, nodeTags);

        for(int node = 0; node < nodeTags.size(); node++){
            elements[elementTags[0][element]].push_back(nodeTags[node]);
        }
    }
}

//Returns vector of velocity and pressure for given element
void Mesh::getElementVector(size_t element, PetscScalar* elemVec, bool pressure){
    vector<double> xVel;
    vector<double> yVel;
    vector<double> p;

    int i = 0;

    for(size_t node : elements[element]){
        xVel.push_back(nodes[node].velocity[0]);
        yVel.push_back(nodes[node].velocity[1]);
        p.push_back(nodes[node].pressure);
    }
    for(double velocity : xVel){
        elemVec[i++] = velocity;
    }
    for(double velocity: yVel){
        elemVec[i++] = velocity;
    }
    if(pressure){
        for(double pressure : p){
            elemVec[i++] = pressure;
        }
    }
}

void Mesh::getForceVector(size_t element, PetscScalar* force){
    vector<double> xForce;
    vector<double> yForce;
    int i = 0;

    for(size_t node: elements[element]){
        xForce.push_back(nodes[node].force[0]);
        yForce.push_back(nodes[node].force[1]);
    }
    for(double f : xForce){
        force[i++] = f;
    }
    for(double f : yForce){
        force[i++] = f;
    }
}

