#include "localbuilder.hpp"

// Destroy sets of matrices
void cleanUp(vector<Mat> matrices) {
  for (Mat m : matrices) {
    MatDestroy(&m);
  }
}

void LocalBuilder::setUp() {
  // Initialise vector sets
  inverseJacobian = vector<Mat>(6);
  basisGradMats = vector<Mat>(6);
  basisMats = vector<Mat>(6);
  j = vector<double>();
  jdets = vector<double>();
  basisFuncs = vector<double>();
  basisFuncsGrad = vector<double>();

  int nComponents;
  int nOrientations;

  // Initialise memory for matrices
  for (int m = 0; m < 6; m++) {
    MatCreate(PETSC_COMM_SELF, &inverseJacobian[m]);
    MatSetSizes(inverseJacobian[m], PETSC_DECIDE, PETSC_DECIDE, 2, 2);
    MatSetFromOptions(inverseJacobian[m]);
    MatSetUp(inverseJacobian[m]);

    MatCreate(PETSC_COMM_SELF, &basisGradMats[m]);
    MatSetSizes(basisGradMats[m], PETSC_DECIDE, PETSC_DECIDE, 2, 12);
    MatSetFromOptions(basisGradMats[m]);
    MatSetUp(basisGradMats[m]);

    MatCreate(PETSC_COMM_SELF, &basisMats[m]);
    MatSetSizes(basisMats[m], PETSC_DECIDE, PETSC_DECIDE, 12, 2);
    MatSetFromOptions(basisMats[m]);
    MatSetUp(basisMats[m]);
  }

  // Obtain basis functions from Gmsh for provided gauss points
  gmsh::model::mesh::getBasisFunctions(9, gaussPoints, "Lagrange", nComponents,
                                       basisFuncs, nOrientations);
  gmsh::model::mesh::getBasisFunctions(2, gaussPoints, "Lagrange", nComponents,
                                       basisFuncsPres, nOrientations);
  gmsh::model::mesh::getBasisFunctions(9, gaussPoints, "GradLagrange",
                                       nComponents, basisFuncsGrad,
                                       nOrientations);
  buildBasisMatrix();
}

// Constructor for local convection matrices
LocalBuilder::LocalBuilder() {
  this->conv = true;
  setUp();
  MatCreate(PETSC_COMM_SELF, &localConvMat);
  MatSetSizes(localConvMat, PETSC_DECIDE, PETSC_DECIDE, 12, 12);
  MatSetFromOptions(localConvMat);
  MatSetUp(localConvMat);
}

// Constructor for local unchanged matrices
LocalBuilder::LocalBuilder(double dt, double viscosity) {
  this->dt = dt;
  this->viscosity = viscosity;
  this->conv = false;

  setUp();

  // Initialise memory for local matrices
  MatCreate(PETSC_COMM_SELF, &localMassMat);
  MatSetSizes(localMassMat, PETSC_DECIDE, PETSC_DECIDE, 12, 12);
  MatSetFromOptions(localMassMat);
  MatSetUp(localMassMat);
  MatDuplicate(localMassMat, MAT_DO_NOT_COPY_VALUES, &localViscMat);

  MatCreate(PETSC_COMM_SELF, &localGradMat);
  MatSetSizes(localGradMat, PETSC_DECIDE, PETSC_DECIDE, 12, 3);
  MatSetFromOptions(localGradMat);
  MatSetUp(localGradMat);

  MatCreate(PETSC_COMM_SELF, &localFullMat);
  MatCreate(PETSC_COMM_SELF, &localFullMat);
  MatSetSizes(localFullMat, PETSC_DECIDE, PETSC_DECIDE, 15, 15);
  MatSetFromOptions(localFullMat);
  MatSetUp(localFullMat);
}

// Destroy LocalBuilder
LocalBuilder::~LocalBuilder() {
  if (!conv) {
    MatDestroy(&localGradMat);
    MatDestroy(&localMassMat);
    MatDestroy(&localViscMat);
    MatDestroy(&localFullMat);
  } else {
    MatDestroy(&localConvMat);
  }
  cleanUp(basisMats);
  cleanUp(inverseJacobian);
  cleanUp(basisGradMats);
}

// Build matrices containing basis functions
void LocalBuilder::buildBasisMatrix() {
  int m = 0;
  // Add a matrix for each integration point
  for (int point = 0; point < basisFuncs.size(); point += 6) {
    // Construct each matrix with evaluation of basis at Gauss points
    for (int i = 0; i < 12; i++) {
      for (int j = 0; j < 2; j++) {
        if (j == 0 && i < 6) {
          MatSetValue(basisMats[m], i, j, basisFuncs[m * 6 + i], INSERT_VALUES);
        } else if (j == 1 && i > 5) {
          MatSetValue(basisMats[m], i, j, basisFuncs[m * 6 + (i - 6)],
                      INSERT_VALUES);
        } else {
          MatSetValue(basisMats[m], i, j, 0, INSERT_VALUES);
        }
      }
    }
    MatAssemblyBegin(basisMats[m], MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(basisMats[m], MAT_FINAL_ASSEMBLY);
    m++;
  }
}

void LocalBuilder::buildInverseJacobian() {
  // Inverse jacobian matrix for each integration point
  int m = 0;
  for (int point = 0; point < j.size(); point += 9) {
    MatZeroEntries(inverseJacobian[m]);
    MatSetValue(inverseJacobian[m], 0, 0, j[point + 4] * 1 / jdets[m],
                INSERT_VALUES);
    MatSetValue(inverseJacobian[m], 0, 1, -j[point + 1] * 1 / jdets[m],
                INSERT_VALUES);
    MatSetValue(inverseJacobian[m], 1, 0, -j[point + 3] * 1 / jdets[m],
                INSERT_VALUES);
    MatSetValue(inverseJacobian[m], 1, 1, j[point] * 1 / jdets[m],
                INSERT_VALUES);

    MatAssemblyBegin(inverseJacobian[m], MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(inverseJacobian[m], MAT_FINAL_ASSEMBLY);
    m++;
  }
}

// Compute basis function gradient matrix w.r.t spatial coordinates
void LocalBuilder::buildBasisGradMatrix() {
  int m = 0;
  // Compute matrix for each integration point
  for (int point = 0; point < basisFuncsGrad.size(); point += 18) {
    Mat temp;
    MatCreate(PETSC_COMM_SELF, &temp);
    MatSetSizes(temp, PETSC_DECIDE, PETSC_DECIDE, 2, 12);
    MatSetFromOptions(temp);
    MatSetUp(temp);

    // Construct local basis gradient matrix
    for (int i = 0; i < 2; i++) {
      int j = 0;
      for (int rep = 0; rep < 2; rep++) {
        for (int basis = 0; basis < 18; basis += 3) {
          MatSetValue(temp, i, j, basisFuncsGrad[point + basis + i],
                      INSERT_VALUES);
          j++;
        }
      }
    }

    MatAssemblyBegin(temp, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(temp, MAT_FINAL_ASSEMBLY);

    // Perform conversion from local coordinates to physical via inverse
    // Jacobian
    MatDestroy(&basisGradMats[m]);
    MatMatMult(inverseJacobian[m], temp, MAT_INITIAL_MATRIX, PETSC_DEFAULT,
               &basisGradMats[m]);
    MatDestroy(&temp);

    MatAssemblyBegin(basisGradMats[m], MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(basisGradMats[m], MAT_FINAL_ASSEMBLY);
    m++;
  }
}

void LocalBuilder::computeMassMatrix() {
  MatZeroEntries(localMassMat);

  // Compute each entry for mass matrix (product of basis functions)
  for (int i = 0; i < 12; i++) {
    PetscScalar massVal = 0.0;
    for (int j = 0; j < 12; j++) {
      int w = 0;
      //   Gauss quadrature (integration)
      for (int gp = 0; gp < 36; gp += 6) {
        massVal += basisFuncs[(i % 6) + gp] * basisFuncs[(j % 6) + gp] *
                   gaussWeights[w] * jdets[w++];
      }
    }
    MatSetValue(localMassMat, i, i, massVal, INSERT_VALUES);
  }

  MatAssemblyBegin(localMassMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(localMassMat, MAT_FINAL_ASSEMBLY);
}

void LocalBuilder::computeViscosityMatrix() {
  MatZeroEntries(localViscMat);
  vector<Mat> viscMats = vector<Mat>(basisGradMats.size());

  // Multiply basis function gradient matrices to obtain second order
  for (int m = 0; m < basisGradMats.size(); m++) {
    MatTransposeMatMult(basisGradMats[m], basisGradMats[m], MAT_INITIAL_MATRIX,
                        PETSC_DEFAULT, &viscMats[m]);
  }

  // Compute each entry for viscosity matrix
  for (int i = 0; i < 12; i++) {
    for (int j = 0; j < 12; j++) {
      PetscScalar viscVal = 0.0;
      int w = 0;
      // Gauss quadrature (integration)
      for (int gp = 0; gp < 36; gp += 6) {
        PetscScalar matVal;
        MatGetValue(viscMats[w], i, j, &matVal);
        viscVal += matVal * gaussWeights[w] * jdets[w++];
      }
      MatSetValue(localViscMat, i, j, viscVal, INSERT_VALUES);
    }
  }

  MatAssemblyBegin(localViscMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(localViscMat, MAT_FINAL_ASSEMBLY);

  // Clean up matrices
  cleanUp(viscMats);
}

void LocalBuilder::computeConvectionMatrix(PetscScalar *v) {
  for (int i = 0; i < 12; i++) {
    for (int j = 0; j < 12; j++) {
      PetscScalar convVal = 0.0;
      int w = 0;
      // Approximate integration using gauss points
      for (int gp = 0; gp < 36; gp += 6) {
        PetscScalar basisVali, basisValj, gradXi, gradYi, gradXj, gradYj;
        MatGetValue(basisMats[w], i % 6, 0, &basisVali);
        MatGetValue(basisMats[w], j % 6, 0, &basisValj);
        MatGetValue(basisGradMats[w], 0, i, &gradXi);
        MatGetValue(basisGradMats[w], 1, i, &gradYi);
        MatGetValue(basisGradMats[w], 0, j, &gradXj);
        MatGetValue(basisGradMats[w], 1, j, &gradYj);

        // Apply basis functions to current velocity
        double velU = 0, velV = 0;
        for (int x = 0; x < 12; x++) {
          if (x < 6) {
            velU += v[x] * basisFuncs[gp + x % 6];
          } else {
            velV += v[x] * basisFuncs[gp + x % 6];
          }
        }

        PetscScalar matVal;
        // Add convection term at current gauss point
        matVal = basisVali * (velU * gradXj + velV * gradYj);
        convVal += (matVal * gaussWeights[w] * jdets[w++]);
      }

      MatSetValue(localConvMat, i, j, (convVal), INSERT_VALUES);
    }
  }
}

void LocalBuilder::computeGradientMatrix() {
  MatZeroEntries(localGradMat);

  // Compute each entry for gradient operator
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 12; j++) {
      int w = 0;
      PetscScalar gradVal = 0.0;
      int x = (j > 5);
      // Gauss quadrature (integration)
      for (int gp = 0; gp < 18; gp += 3) {
        PetscScalar matVal = 0.0;
        MatGetValue(basisGradMats[w], x, j, &matVal);
        gradVal +=
            basisFuncsPres[gp + i] * matVal * gaussWeights[w] * jdets[w++];
      }
      MatSetValue(localGradMat, j, i, gradVal * -1, INSERT_VALUES);
    }
  }

  MatAssemblyBegin(localGradMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(localGradMat, MAT_FINAL_ASSEMBLY);
}

// Compute matrix for pressure and final velocity
void LocalBuilder::computeFinalMatrix() {
  MatZeroEntries(localFullMat);

  for (int i = 0; i < 12; i++) {
    for (int j = 0; j < 12; j++) {
      PetscScalar massVal;
      MatGetValue(localMassMat, i, j, &massVal);
      MatSetValue(localFullMat, i, j, massVal * dt, INSERT_VALUES);
    }
  }

  for (int i = 0; i < 12; i++) {
    for (int j = 12; j < 15; j++) {
      PetscScalar matVal;
      MatGetValue(localGradMat, i, j - 12, &matVal);
      MatSetValue(localFullMat, i, j, matVal, INSERT_VALUES);
    }
  }

  for (int i = 12; i < 15; i++) {
    for (int j = 0; j < 12; j++) {
      PetscScalar matVal;
      MatGetValue(localGradMat, j, i - 12, &matVal);
      MatSetValue(localFullMat, i, j, matVal, INSERT_VALUES);
    }
  }

  for (int i = 12; i < 15; i++) {
    for (int j = 12; j < 15; j++) {
      MatSetValue(localFullMat, i, j, 0, INSERT_VALUES);
    }
  }

  MatAssemblyBegin(localFullMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(localFullMat, MAT_FINAL_ASSEMBLY);
}

// Assemble local unchanged matrices
void LocalBuilder::assembleMatrices(size_t elementTag) {
  this->elementTag = elementTag;
  vector<double> coords;

  gmsh::model::mesh::getJacobian(elementTag, gaussPoints, j, jdets, coords);

  // Inverse Jacobian matrix and build basis function gradient matrix
  buildInverseJacobian();
  buildBasisGradMatrix();

  // Compute all local matrices
  computeMassMatrix();
  computeViscosityMatrix();
  computeGradientMatrix();
  computeFinalMatrix();
}

// Assemble local convection matrix (changes on each timestep)
void LocalBuilder::assembleConvectionMatrix(size_t elementTag,
                                            Vec *velocityVec) {
  // Obtain Jacobian from Gmsh using the elements tag
  this->elementTag = elementTag;
  vector<double> coords;
  gmsh::model::mesh::getJacobian(elementTag, gaussPoints, j, jdets, coords);

  // Inverse Jacobian matrix and build basis function gradient matrix for
  // convection term
  buildInverseJacobian();
  MatZeroEntries(localConvMat);
  buildBasisGradMatrix();

  PetscScalar v[12] = {0};
  for (int a = 0; a < 12; a++) {
    VecGetValues(*velocityVec, 1, &a, &v[a]);
  }

  computeConvectionMatrix(v);

  MatAssemblyBegin(localConvMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(localConvMat, MAT_FINAL_ASSEMBLY);
}
