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

// Retrieve nodes from Gmsh
void Mesh::getNodes(int* nId, vector<size_t> nodeTags,
                    vector<double> nodeCoords, double boundaryVel,
                    bool boundary, bool inlet) {
  // Enumerate over given nodes
  for (int node = 0; node < nodeCoords.size() / 3; node++) {
    if (nodes.find(nodeTags[node]) == nodes.end()) {
      // Add to Dirichlet indices if part of the boundary
      if (boundary || inlet) {
        dirichletIds.push_back(*nId);
        dirichletIds.push_back(*nId + nP2);
      }
      nodeIds[*nId] = nodeTags[node];
      // Initialise node struct and add to set of nodes
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

  this->nElements = elementTags[0].size();
  // Enumerate all elements
  for (int element = 0; element < elementTags[0].size(); element++) {
    int type;
    vector<size_t> nodeTags;
    gmsh::model::mesh::getElement(elementTags[0][element], type, nodeTags);

    // Add nodes owned by current element
    for (int node = 0; node < nodeTags.size(); node++) {
      if (node < 3) {
        // Determine linear nodes (p1) for pressure
        if (pNodes.find(nodeTags[node]) == pNodes.end()) {
          nodes[nodeTags[node]].pid = pId++;
          pNodes.emplace(nodeTags[node]);
        }
      }
      elements[elementTags[0][element]].push_back(nodeTags[node]);
    }
  }
  nP1 = pId;
}

Mesh::Mesh(string filePath, double boundaryVel) {
  gmsh::open(filePath);
  this->boundaryVel = boundaryVel;

  vector<size_t> nodeTags;
  vector<double> nodeCoords;
  vector<double> paramCoords;
  int nId = 0;

  gmsh::model::mesh::getNodes(nodeTags, nodeCoords, paramCoords);
  nP2 = nodeTags.size();

  // Retrieve nodes at boundary
  int boundaryTag = 1;
  gmsh::model::mesh::getNodesForPhysicalGroup(1, boundaryTag, nodeTags,
                                              nodeCoords);
  getNodes(&nId, nodeTags, nodeCoords, 0, true, false);

  // Retrieve nodes at inlet
  int inletTag = 2;
  gmsh::model::mesh::getNodesForPhysicalGroup(1, inletTag, nodeTags,
                                              nodeCoords);
  getNodes(&nId, nodeTags, nodeCoords, boundaryVel, false, true);

  // Retrieve all other nodes
  gmsh::model::mesh::getNodes(nodeTags, nodeCoords, paramCoords);
  getNodes(&nId, nodeTags, nodeCoords, 0, false, false);

  this->nDirichlet = dirichletIds.size();

  // Retrieve element connectivity
  getElementConnectivity();
}

// Number of linear (p1) and quadratic (p2) node points
int Mesh::p1Size() { return nP1; }
int Mesh::p2Size() { return nP2; }

// Number of elements and dirichlet nodes
int Mesh::elementSize() { return nElements; }
int Mesh::dirichletSize() { return nDirichlet; }

// Node retrieving methods
Node Mesh::getNode(int i) { return nodes[nodeIds[i]]; }
Node Mesh::getNode(size_t nodeTag) { return nodes[nodeTag]; }
Node Mesh::getNode(size_t elementTag, int i) {
  return nodes[elements[elementTag][i]];
}

// Element retrieving methods
size_t Mesh::getElement(int e) { return elementTags[0][e]; }
vector<size_t> Mesh::getElements() { return elementTags[0]; }

// Retrieve Dirichlet index node
PetscInt* Mesh::getDirichlet() { return dirichletIds.data(); }
