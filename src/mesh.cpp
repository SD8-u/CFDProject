#include <mesh.hpp>
#include <unordered_set>

Mesh::Mesh(string filePath){
    gmsh::open(filePath);
    
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
        {0}, {300, 0}, 0, nodeCoords[node * 3], nodeCoords[node * 3 + 1]};
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
    
    //Retrieve fluid nodes
    int fluidTag = 3;
    gmsh::model::mesh::getNodesForPhysicalGroup(2, fluidTag, nodeTags, nodeCoords);
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

