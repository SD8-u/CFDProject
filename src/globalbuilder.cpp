#include "globalbuilder.hpp"

// Get subset of fluid domain for current processor
void GlobalBuilder::getDomain() {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  domainSize = (msh->elementSize() / size);
  domainStart = rank * domainSize;
}

// Initialise memory for matrices
void GlobalBuilder::initialiseMats() {
  MatCreate(PETSC_COMM_WORLD, &globalMassMat);
  MatSetSizes(globalMassMat, PETSC_DECIDE, PETSC_DECIDE, msh->p2Size() * 2,
              msh->p2Size() * 2);
  MatSetFromOptions(globalMassMat);
  MatSetUp(globalMassMat);

  MatCreate(PETSC_COMM_WORLD, &globalFullMat);
  MatCreate(PETSC_COMM_WORLD, &globalFullMat);
  MatSetSizes(globalFullMat, PETSC_DECIDE, PETSC_DECIDE,
              msh->p2Size() * 2 + msh->p1Size(),
              msh->p2Size() * 2 + msh->p1Size());
  MatSetFromOptions(globalFullMat);
  MatSetUp(globalFullMat);

  MatDuplicate(globalMassMat, MAT_DO_NOT_COPY_VALUES, &globalViscMat);
  MatDuplicate(globalMassMat, MAT_DO_NOT_COPY_VALUES, &globalConvMat);

  MatZeroEntries(globalMassMat);
  MatZeroEntries(globalViscMat);
  MatZeroEntries(globalConvMat);
  MatZeroEntries(globalFullMat);
}

// Initialise memory for vectors
void GlobalBuilder::initialiseVecs() {
  VecCreate(PETSC_COMM_WORLD, &currVelVec);
  VecSetSizes(currVelVec, PETSC_DECIDE, msh->p2Size() * 2);
  VecSetFromOptions(currVelVec);
  VecDuplicate(currVelVec, &prevVelVec);

  VecCreate(PETSC_COMM_WORLD, &fullVec);
  VecSetSizes(fullVec, PETSC_DECIDE, msh->p2Size() * 2 + msh->p1Size());
  VecSetFromOptions(fullVec);
}

GlobalBuilder::GlobalBuilder(int dim, double dt, double visc, Mesh *msh) {
  this->dt = 1 / dt;
  this->viscosity = visc;
  this->msh = msh;

  getDomain();
  initialiseMats();
  initialiseVecs();
}

// Destroy GlobalBuilder
GlobalBuilder::~GlobalBuilder() {
  MatDestroy(&globalMassMat);
  MatDestroy(&globalViscMat);
  MatDestroy(&globalConvMat);
  MatDestroy(&globalFullMat);
  VecDestroy(&currVelVec);
  VecDestroy(&fullVec);
}

// Add entries of local matrix to their global matrix entries
void GlobalBuilder::localToGlobalMat(size_t elementTag, Mat *localMat,
                                     Mat *globalMat, bool final = false) {
  int add = final ? 3 : 0;
  // Enumerate entries of local matrix
  for (int i = 0; i < 12 + add; i++) {
    // Get global row index of node
    Node iNode = msh->getNode(elementTag, i % 6);
    int x, y;
    x = iNode.id;
    if (i > 5) {
      x += msh->p2Size();
    }
    if (i > 11) {
      x = iNode.pid;
      x += msh->p2Size() * 2;
    }
    for (int j = 0; j < 12 + add; j++) {
      // Get global column index of node
      Node jNode = msh->getNode(elementTag, j % 6);
      y = jNode.id;
      if (j > 5) {
        y += msh->p2Size();
      }
      if (j > 11) {
        y = jNode.pid;
        y += msh->p2Size() * 2;
      }
      PetscScalar matVal;
      // Add value from local matrix to corresponding global entry
      MatGetValue(*localMat, i, j, &matVal);
      MatSetValue(*globalMat, x, y, matVal, ADD_VALUES);
    }
  }
}

// Constuct the global velocity/pressure vectors from initial conditions
void GlobalBuilder::localToGlobalVec(bool full) {
  Vec *vec = &currVelVec;
  if (full) {
    vec = &fullVec;
  }

  // Determine processors ownership domain of global vector
  int high, low;
  VecGetOwnershipRange(*vec, &low, &high);

  // Enumerate nodes of all elements
  vector<size_t> elements = msh->getElements();
  for (size_t element : elements) {
    for (int i = 0; i < 6; i++) {
      Node n = msh->getNode(element, i);
      int x = n.id;
      int y = n.pid;

      // Set inlet/initial conditions if processor owns vector entry
      if (n.inlet) {
        if (x >= low && x < high) {
          VecSetValue(*vec, x, n.velocity[0], INSERT_VALUES);
        }
        if (x + msh->p2Size() >= low && x + msh->p2Size() < high) {
          VecSetValue(*vec, x + msh->p2Size(), n.velocity[1], INSERT_VALUES);
        }
      } else {
        if (x >= low && x < high) {
          VecSetValue(*vec, x, 0.0, INSERT_VALUES);
        }
        if (x + msh->p2Size() >= low && x + msh->p2Size() < high) {
          VecSetValue(*vec, x + msh->p2Size(), 0.0, INSERT_VALUES);
        }
      }
      if (i < 3 && full) {
        if (y + 2 * msh->p2Size() >= low && y + 2 * msh->p2Size() < high) {
          VecSetValue(*vec, 2 * msh->p2Size() + y, n.pressure, INSERT_VALUES);
        }
      }
    }
  }
  VecAssemblyBegin(*vec);
  VecAssemblyEnd(*vec);
}

// Gather parallel distributed vector entries into a local vector
void GlobalBuilder::globalToLocalVec(size_t elementTag, Vec *localVec) {
  PetscInt gblIndices[12];
  for (int node = 0; node < 6; node++) {
    PetscInt ix = msh->getNode(elementTag, node).id;
    PetscInt iy = ix + msh->p2Size();
    PetscScalar vx, vy;
    gblIndices[node] = ix;
    gblIndices[node + 6] = iy;
  }
  IS is;
  VecScatter vs;

  // Create vector scatter to obtain entries from multiple processors
  ISCreateGeneral(PETSC_COMM_WORLD, 12, gblIndices, PETSC_COPY_VALUES, &is);
  VecScatterCreate(currVelVec, is, *localVec, NULL, &vs);
  VecScatterBegin(vs, currVelVec, *localVec, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(vs, currVelVec, *localVec, INSERT_VALUES, SCATTER_FORWARD);

  VecScatterDestroy(&vs);
  ISDestroy(&is);
  VecAssemblyBegin(*localVec);
  VecAssemblyEnd(*localVec);
}

// Assemble all unchanged global matrices
void GlobalBuilder::assembleMatrices() {
  // Use LocalBuilder for each element
  localBuild = new LocalBuilder(dt, viscosity);

  int domainEnd = domainStart + domainSize;
  if (msh->elementSize() - domainEnd < domainSize) {
    domainEnd += msh->elementSize() - domainEnd;
  }

  // Enumerate assigned domain elements and construct local matrices
  for (int e = domainStart; e < domainEnd; e++) {
    size_t element = msh->getElement(e);
    localBuild->assembleMatrices(element);
    localToGlobalMat(element, &localBuild->localMassMat, &globalMassMat);
    localToGlobalMat(element, &localBuild->localViscMat, &globalViscMat);
    localToGlobalMat(element, &localBuild->localFullMat, &globalFullMat, true);
  }
  delete (localBuild);

  // Synchronise all matrix updates between processors
  MatAssemblyBegin(globalMassMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(globalMassMat, MAT_FINAL_ASSEMBLY);

  MatAssemblyBegin(globalViscMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(globalViscMat, MAT_FINAL_ASSEMBLY);

  MatAssemblyBegin(globalFullMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(globalFullMat, MAT_FINAL_ASSEMBLY);

  // Scale mass and viscosity matrices
  MatScale(globalViscMat, viscosity);
  MatScale(globalMassMat, dt);
}

// Assemble global convection matrix (changes on each time step)
void GlobalBuilder::assembleConvectionMatrix() {
  MatZeroEntries(globalConvMat);
  Vec localVelVec;
  VecCreate(PETSC_COMM_SELF, &localVelVec);
  VecSetSizes(localVelVec, PETSC_DECIDE, 12);
  VecSetFromOptions(localVelVec);

  // Use LocalBuilder for each element
  localBuild = new LocalBuilder();

  int domainEnd = domainStart + domainSize;

  // Enumerate assigned domain elements and construct local matrices
  for (int e = domainStart; e < domainEnd; e++) {
    size_t elementTag = msh->getElement(e);
    globalToLocalVec(elementTag, &localVelVec);
    localBuild->assembleConvectionMatrix(elementTag, &localVelVec);
    localToGlobalMat(elementTag, &localBuild->localConvMat, &globalConvMat);
  }

  // Synchronise all matrix updates between processors
  MatAssemblyBegin(globalConvMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(globalConvMat, MAT_FINAL_ASSEMBLY);

  delete (localBuild);
  VecDestroy(&localVelVec);
}

// Assemble velocity/pressure vectors for solution
void GlobalBuilder::assembleVectors() {
  localToGlobalVec(false);
  localToGlobalVec(true);
}
