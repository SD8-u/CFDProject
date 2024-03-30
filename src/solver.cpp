#include "solver.hpp"

// Update distibuted vectors using vector scatters
void Solver::updateVectors(Vec *vec1, Vec *vec2, bool vel) {
  VecScatter *vecScatter = vel ? &vecScatter1 : &vecScatter2;
  VecScatterBegin(*vecScatter, *vec1, *vec2, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(*vecScatter, *vec1, *vec2, INSERT_VALUES, SCATTER_FORWARD);
  VecAssemblyBegin(*vec2);
  VecAssemblyEnd(*vec2);
}

// Initialise the KSP GMRES solver instance
void Solver::initialiseSolver() {
  PC preConditioner;
  KSPCreate(PETSC_COMM_WORLD, &stp2Solver);
  KSPSetType(stp2Solver, KSPGMRES);
  KSPSetOperators(stp2Solver, globalBuild->globalFullMat,
                  globalBuild->globalFullMat);

  KSPGetPC(stp2Solver, &preConditioner);
  PCSetType(preConditioner, PCNONE);
  KSPSetFromOptions(stp2Solver);
}

Solver::Solver(Mesh *msh, double dt, double viscosity) {
  const int kVelSize = msh->p2Size() * 2;
  this->msh = msh;

  // Initialise system using GlobalBuilder
  Vec tempVec;
  globalBuild = new GlobalBuilder(2, dt, viscosity, msh);
  globalBuild->assembleMatrices();
  globalBuild->assembleVectors();

  VecCreate(PETSC_COMM_WORLD, &tempVec);
  VecSetSizes(tempVec, PETSC_DECIDE, kVelSize + msh->p1Size());
  VecSetFromOptions(tempVec);

  applyDirichletConditions(&globalBuild->globalFullMat, &tempVec, false);
  VecDestroy(&tempVec);

  initialiseSolver();

  // Construct scatter mappings for distributed vectors
  PetscInt vecIndices[kVelSize];
  for (int i = 0; i < kVelSize; i++) vecIndices[i] = i;

  ISCreateGeneral(PETSC_COMM_WORLD, kVelSize, vecIndices, PETSC_COPY_VALUES,
                  &vecMapping);
  VecScatterCreate(globalBuild->fullVec, vecMapping, globalBuild->currVelVec,
                   vecMapping, &vecScatter1);
  VecScatterCreate(globalBuild->currVelVec, vecMapping, globalBuild->fullVec,
                   vecMapping, &vecScatter2);
}

// Destroy solver
Solver::~Solver() {
  KSPDestroy(&stp2Solver);
  VecScatterDestroy(&vecScatter1);
  VecScatterDestroy(&vecScatter2);
  ISDestroy(&vecMapping);
  delete (globalBuild);
}

// Apply Dirichlet conditions to system
void Solver::applyDirichletConditions(Mat *m, Vec *v, bool expl) {
  PetscInt *rows = msh->getDirichlet();
  int high, low;
  VecGetOwnershipRange(*v, &low, &high);
  for (int i = 0; i < msh->p2Size(); i++) {
    Node n = msh->getNode(i);
    if (n.boundary || n.inlet) {
      if (i >= low && i < high) {
        VecSetValue(*v, i, n.velocity[0], INSERT_VALUES);
      }
      if (i + msh->p2Size() >= low && i + msh->p2Size() < high) {
        VecSetValue(*v, i + msh->p2Size(), n.velocity[1], INSERT_VALUES);
      }
    }
  }

  VecAssemblyBegin(*v);
  VecAssemblyEnd(*v);

  if (!expl) {
    MatZeroRows(*m, msh->dirichletSize(), rows, 1.0, *v, *v);
  }
  MatAssemblyBegin(*m, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*m, MAT_FINAL_ASSEMBLY);
}

// Assemble matrix and RHS for Crank Nicolson time stepping
void Solver::crankNicolson(Mat *tempMat, Vec *tempVec) {
  MatConvert(globalBuild->globalMassMat, MATSAME, MAT_INITIAL_MATRIX, tempMat);
  // Subtract half viscous + convection terms from right hand side
  MatAXPY(*tempMat, -1 / 2, globalBuild->globalViscMat,
          DIFFERENT_NONZERO_PATTERN);
  MatAXPY(*tempMat, -1 / 2, globalBuild->globalConvMat,
          DIFFERENT_NONZERO_PATTERN);

  MatMult(*tempMat, globalBuild->currVelVec, *tempVec);
  // Reapply viscous and convection terms to system of equations
  MatAXPY(*tempMat, 1.0, globalBuild->globalViscMat, DIFFERENT_NONZERO_PATTERN);
  MatAXPY(*tempMat, 1.0, globalBuild->globalConvMat, DIFFERENT_NONZERO_PATTERN);
}

// Compute intermediate velocity in Chorin-Temam
void Solver::computeFirstStep() {
  Mat tempMat;
  Vec tempVec;
  PC preConditioner;

  VecCreate(PETSC_COMM_WORLD, &tempVec);
  VecSetSizes(tempVec, PETSC_DECIDE, msh->p2Size() * 2);
  VecSetFromOptions(tempVec);

  // Assemble system
  globalBuild->assembleConvectionMatrix();
  crankNicolson(&tempMat, &tempVec);

  // Impose Dirichlet Conditions
  applyDirichletConditions(&tempMat, &tempVec, false);

  // Initialise KSP and Solve system
  KSPCreate(PETSC_COMM_WORLD, &stp1Solver);
  KSPSetType(stp1Solver, KSPGMRES);
  KSPSetOperators(stp1Solver, tempMat, tempMat);

  KSPGetPC(stp1Solver, &preConditioner);
  PCSetType(preConditioner, PCNONE);
  KSPSetFromOptions(stp1Solver);

  KSPSolve(stp1Solver, tempVec, globalBuild->currVelVec);
  applyDirichletConditions(&tempMat, &globalBuild->currVelVec, true);

  // Cleanup
  VecDestroy(&tempVec);
  MatDestroy(&tempMat);
  KSPDestroy(&stp1Solver);
}

// Compute final velocity and pressure in Chorin-Temam
void Solver::computeSecondStep() {
  Vec tempVec;
  Vec solVec;
  Mat tempMat;

  VecCreate(PETSC_COMM_WORLD, &tempVec);
  VecCreate(PETSC_COMM_WORLD, &solVec);
  VecSetSizes(tempVec, PETSC_DECIDE, msh->p2Size() * 2);
  VecSetSizes(solVec, PETSC_DECIDE, msh->p2Size() * 2 + msh->p1Size());

  VecSetFromOptions(tempVec);
  VecSetFromOptions(solVec);

  // Compute RHS of second linear system
  MatMult(globalBuild->globalMassMat, globalBuild->currVelVec, tempVec);
  updateVectors(&tempVec, &solVec, false);

  // Apply Dirichlet conditions and solve the system
  applyDirichletConditions(&globalBuild->globalFullMat, &solVec, true);
  KSPSolve(stp2Solver, solVec, globalBuild->fullVec);
  applyDirichletConditions(&globalBuild->globalFullMat, &globalBuild->fullVec,
                           true);

  // Update the velocity vector
  updateVectors(&globalBuild->fullVec, &globalBuild->currVelVec, true);

  // Cleanup
  VecDestroy(&solVec);
  VecDestroy(&tempVec);

  PetscReal max;
  VecMax(globalBuild->currVelVec, NULL, &max);
  cout << "Max (Instability Metric): " << max << "\n";
}

// Perform time marching for transient simulation
void Solver::computeTimeSteps(int steps) {
  for (int x = 0; x < steps; x++) {
    this->computeFirstStep();
    this->computeSecondStep();
    cout << "Step: " << x << "\n";
  }
}

// Interpolate specific velocity/pressure values given coords
void Solver::interpolateValues(vector<double> *coord,
                               vector<vector<double>> *solData,
                               vector<size_t> *nodeTags, Vec *solVec) {
  int nComp, nOrien;
  double pressure = 0, xv = 0, yv = 0;
  double nodePre, nodeVx, nodeVy;
  vector<double> basisFuncsVel, basisFuncsPre;

  gmsh::model::mesh::getBasisFunctions(2, *coord, "Lagrange", nComp,
                                       basisFuncsPre, nOrien);
  gmsh::model::mesh::getBasisFunctions(9, *coord, "Lagrange", nComp,
                                       basisFuncsVel, nOrien);

  for (int n = 0; n < (*nodeTags).size(); n++) {
    Node node = msh->getNode((*nodeTags)[n]);
    PetscInt pi = node.pid + msh->p2Size() * 2;
    PetscInt vxi = node.id;
    PetscInt vyi = vxi + msh->p2Size();
    VecGetValues(*solVec, 1, &vxi, &nodeVx);
    VecGetValues(*solVec, 1, &vyi, &nodeVy);
    VecGetValues(*solVec, 1, &pi, &nodePre);
    if (n < 3) {
      pressure += basisFuncsPre[n] * nodePre;
    }

    if (basisFuncsVel[n] <= 1) {
      xv += basisFuncsVel[n] * nodeVx;
      yv += basisFuncsVel[n] * nodeVy;
    }
  }

  (*solData)[0].push_back(pressure);
  (*solData)[1].push_back(xv);
  (*solData)[2].push_back(yv);
  (*coord).clear();
}

// Interpolate the velocity/pressure fields via basis functions
vector<vector<double>> Solver::interpolateSolution(double resolution,
                                                   int rank) {
  Vec solVec;
  VecScatter vs;
  vector<vector<double>> solData = vector<vector<double>>(5);

  VecScatterCreateToZero(globalBuild->fullVec, &vs, &solVec);
  VecScatterBegin(vs, globalBuild->fullVec, solVec, INSERT_VALUES,
                  SCATTER_FORWARD);
  VecScatterEnd(vs, globalBuild->fullVec, solVec, INSERT_VALUES,
                SCATTER_FORWARD);

  if (rank != 0) {
    VecDestroy(&solVec);
    return solData;
  }

  size_t elementTag;
  int elementType;
  vector<size_t> nodeTags;
  double u, v, w;
  double xmax, ymax, zmax;
  double xmin, ymin, zmin;
  vector<double> coord;

  gmsh::model::getBoundingBox(-1, -1, xmin, ymin, zmin, xmax, ymax, zmax);
  for (double x = xmin; x < xmax; x += resolution) {
    for (double y = ymin; y < ymax; y += resolution) {
      solData[3].push_back(x);
      solData[4].push_back(y);

      try {
        gmsh::model::mesh::getElementByCoordinates(
            x, y, 0, elementTag, elementType, nodeTags, u, v, w, -1, true);
      } catch (...) {
        solData[0].push_back(0);
        solData[1].push_back(0);
        solData[2].push_back(0);
        continue;
      }

      coord.push_back(u);
      coord.push_back(v);
      coord.push_back(w);

      interpolateValues(&coord, &solData, &nodeTags, &solVec);
    }
  }
  VecDestroy(&solVec);
  return solData;
}
