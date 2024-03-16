#include "mesh.hpp"

// Load Gmsh script and generate msh file
void Mesh::generateMesh(string filePath, int refinement) {
  gmsh::merge(filePath);
  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate(2);

  // Refine mesh uniformly
  for (int x = 0; x < refinement; x++) gmsh::model::mesh::refine();

  // Construct quadratic triangle elements
  gmsh::model::mesh::setOrder(2);

  gmsh::write("geometry/temp.msh");
}

Mesh::Mesh(string filePath, double boundaryVel) {
  gmsh::open(filePath);
  this->boundaryVel = boundaryVel;

  vector<size_t> nodeTags;
  vector<double> nodeCoords;
  vector<double> paramCoords;
  int nId = 0;

  gmsh::model::mesh::getNodes(nodeTags, nodeCoords, paramCoords);
  nNodes = nodeTags.size();

  int inletTag = 2;
  gmsh::model::mesh::getNodesForPhysicalGroup(1, inletTag, nodeTags,
                                              nodeCoords);
  getNodes(&nId, nodeTags, nodeCoords, boundaryVel, false, true);

  int boundaryTag = 1;
  gmsh::model::mesh::getNodesForPhysicalGroup(1, boundaryTag, nodeTags,
                                              nodeCoords);
  getNodes(&nId, nodeTags, nodeCoords, 0, true, false);

  // Retrieve all other nodes
  gmsh::model::mesh::getNodes(nodeTags, nodeCoords, paramCoords);
  getNodes(&nId, nodeTags, nodeCoords, 0, false, false);

  // Retrieve element connectivity
  getElementConnectivity();
}

void Mesh::getNodes(int *nId, vector<size_t> nodeTags,
                    vector<double> nodeCoords, double boundaryVel,
                    bool boundary, bool inlet) {
  // Retrieve nodes
  for (int node = 0; node < nodeCoords.size() / 3; node++) {
    if (nodes.find(nodeTags[node]) == nodes.end() || inlet) {
      if (boundary || inlet) {
        dirichletIds.push_back(*nId);
        dirichletIds.push_back(*nId + nNodes);
      }
      nodeIds[*nId] = nodeTags[node];

      nodes[nodeTags[node]] = Node{(*nId)++,
                                   -1,
                                   boundary,
                                   inlet,
                                   {0},
                                   {boundaryVel, 0},
                                   0,
                                   nodeCoords[node * 3],
                                   nodeCoords[node * 3 + 1]};
    }
  }
}

// Encode Element Connectivity
void Mesh::getElementConnectivity() {
  vector<int> elementTypes;
  vector<vector<size_t>> elementNodeTags;
  gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags, 2,
                                 -1);
  unordered_set<size_t> pNodes;
  int pId = 0;

  this->elementSize = elementTags[0].size();
  for (int element = 0; element < elementTags[0].size(); element++) {
    int type;
    vector<size_t> nodeTags;
    gmsh::model::mesh::getElement(elementTags[0][element], type, nodeTags);

    for (int node = 0; node < nodeTags.size(); node++) {
      if (node < 3) {
        if (pNodes.find(nodeTags[node]) == pNodes.end()) {
          nodes[nodeTags[node]].pid = pId++;
          pNodes.emplace(nodeTags[node]);
        }
      }
      elements[elementTags[0][element]].push_back(nodeTags[node]);
    }
  }
  nLinear = pId;
}
