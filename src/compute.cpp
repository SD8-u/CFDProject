#include <compute.hpp>

//Retrieve guass quadrature parameters
void getIntegrationPoints(int type, vector<double> &gaussPoints, vector<double> &gaussWeights){
    switch(type){
        //linear triangle
        case(1):
        //quadratic triangle
        case(2):
            gaussPoints = {0.44594849092, 0.44594849092, 0, 0.10810301817, 0.44594849092, 0, 
                           0.44594849092, 0.10810301817, 0, 0.09157621351, 0.09157621351, 0, 
                           0.81684757289, 0.09157621351, 0, 0.09157621351, 0.81684757289, 0};

            gaussWeights = {0.111690794839, 0.111690794839, 0.111690794839, 
                            0.054975871827, 0.054975871827, 0.054975871827};
            break;
    }
}

vector<Mat> getBasisMatrix(vector<double> gaussPoints, vector<double> basisFuncs){
    vector<Mat> basisMats = vector<Mat>(basisFuncs.size()/6);
    int m = 0;
    for(int point = 0; point < basisFuncs.size(); point+=6){
        
        MatCreate(PETSC_COMM_WORLD, &basisMats[m]);
        MatSetSizes(basisMats[m], PETSC_DECIDE, PETSC_DECIDE, 12, 2);
        MatSetFromOptions(basisMats[m]);
        MatSetUp(basisMats[m]);

        for(int i = 0; i < 12; i++){
            for(int j = 0; j < 2; j++){
                if(j == 0 && i < 6){
                    MatSetValue(basisMats[m], i, j, basisFuncs[m * 6 + i], INSERT_VALUES);
                }
                else if(j == 1 && i > 5){
                    MatSetValue(basisMats[m], i, j, basisFuncs[m * 6 + (i - 6)], INSERT_VALUES);
                }
                else{
                    MatSetValue(basisMats[m], i, j, 0, INSERT_VALUES);
                }
            }
        }

        MatAssemblyBegin(basisMats[m], MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(basisMats[m], MAT_FINAL_ASSEMBLY);
        m++;
    }

    return basisMats;
}


vector<Mat> computeInverseJacobian(vector<double> j, vector<double> dets){
    //Inverse jacobian matrix for each integration point
    vector<Mat> inverseJacobian = vector<Mat>(j.size()/9);
    int m = 0;
    for(int point = 0; point < j.size(); point+=9){

        MatCreate(PETSC_COMM_WORLD, &inverseJacobian[m]);
        MatSetSizes(inverseJacobian[m], PETSC_DECIDE, PETSC_DECIDE, 2, 2);
        MatSetFromOptions(inverseJacobian[m]);
        MatSetUp(inverseJacobian[m]);

        MatSetValue(inverseJacobian[m], 0, 0, j[point + 4] * 1/dets[m], INSERT_VALUES);
        MatSetValue(inverseJacobian[m], 0, 1, -j[point + 1] * 1/dets[m], INSERT_VALUES);
        MatSetValue(inverseJacobian[m], 1, 0, -j[point + 3] * 1/dets[m], INSERT_VALUES);
        MatSetValue(inverseJacobian[m], 1, 1, j[point] * 1/dets[m], INSERT_VALUES);

        MatAssemblyBegin(inverseJacobian[m], MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(inverseJacobian[m], MAT_FINAL_ASSEMBLY);
        m++;
    }
    return inverseJacobian;
}

//Compute basis function gradient matrix w.r.t spatial coordinates
vector<Mat> computeBasisGradMatrix(vector<double> basisFuncsGrad, vector<double> j, vector<double> dets){
    vector<Mat> basisGradMat = vector<Mat>(6);
    vector<Mat> inverseJacob = computeInverseJacobian(j, dets);
    int m = 0;
    //Compute matrix for each integration point
    for(int point = 0; point < basisFuncsGrad.size(); point+=18){
        Mat temp;
        MatCreate(PETSC_COMM_WORLD, &temp);
        MatCreate(PETSC_COMM_WORLD, &basisGradMat[m]);
        MatSetSizes(temp, PETSC_DECIDE, PETSC_DECIDE, 2, 12);
        MatSetSizes(basisGradMat[m], PETSC_DECIDE, PETSC_DECIDE, 2, 12);
        MatSetFromOptions(temp);
        MatSetFromOptions(basisGradMat[m]);
        MatSetUp(temp);
        MatSetUp(basisGradMat[m]);

        //Construct local basis gradient matrix
        for(int i = 0; i < 2; i++){
            int j = 0;
            for(int rep = 0; rep < 2; rep++){
                for(int basis = 0; basis < 18; basis+=3){
                    MatSetValue(temp, i, j, basisFuncsGrad[point+basis+i], INSERT_VALUES);
                    j++;
                }
            }
        }
            MatAssemblyBegin(temp, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(temp, MAT_FINAL_ASSEMBLY);

            //Perform conversion from local coordinates to physical via inverse jacobian
            MatMatMult(inverseJacob[m], temp, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &basisGradMat[m]);

            MatDestroy(&temp);
            MatAssemblyBegin(basisGradMat[m], MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(basisGradMat[m], MAT_FINAL_ASSEMBLY);

        m++;
    }

    return basisGradMat;
}

Mat computeMassMatrix(size_t elementTag){
    vector<double> coords;
    vector<double> gaussPoints;
    vector<double> gaussWeights;
    vector<double> basisFuncs;
    vector<double> j;
    vector<double> jdets;
    int nComponents;
    int nOrientations;

    //Initialise Mass Matrix (components are the integrals of basis functions)
    Mat massMatrix;
    MatCreate(PETSC_COMM_WORLD, &massMatrix);
    MatSetSizes(massMatrix, PETSC_DECIDE, PETSC_DECIDE, 12, 12);
    MatSetFromOptions(massMatrix);
    MatSetUp(massMatrix);

    getIntegrationPoints(2, gaussPoints, gaussWeights);

    gmsh::model::mesh::getBasisFunctions(9, gaussPoints, "Lagrange", nComponents, basisFuncs, nOrientations);
    gmsh::model::mesh::getJacobian(elementTag, gaussPoints, j, jdets, coords);

    PetscScalar massValues[144] = {0.0};
    PetscInt row[144];
    PetscInt col[144];

    //Enumerate the matrix
    for(int i = 0; i < 12; i++){
        for(int j = 0; j < 12; j++){
            PetscScalar massVal = 0.0;
            int w = 0;
            int offset = (j > 5 && i > 5) ? 6 : 0;
            
            if(j < 6 && i < 6 || j > 5 && i > 5){
                //Approximate integration using gauss points over the products of basis functions
                for(int gp = 0; gp < 36; gp+=6){
                    massVal += basisFuncs[(i-offset) + gp] * basisFuncs[(j-offset) + gp] 
                    * gaussWeights[w] * jdets[w++];
                }
            }
            MatSetValue(massMatrix, i, j, massVal, INSERT_VALUES);
        }
    }

    MatAssemblyBegin(massMatrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(massMatrix, MAT_FINAL_ASSEMBLY);

    // Print the mass matrix
    MatView(massMatrix, PETSC_VIEWER_STDOUT_WORLD);

    return massMatrix;
}

Mat computeViscosityMatrix(size_t elementTag){
    double viscosity = 1;
    vector<double> coords;
    vector<double> gaussPoints;
    vector<double> gaussWeights;
    vector<double> basisFuncsGrad;
    vector<double> j;
    vector<double> jdets;
    int nComponents;
    int nOrientations;

    //Initialise Viscosity Matrix (components are the integrals of basis function gradients * kinematic viscosity)
    Mat viscosityMatrix;
    MatCreate(PETSC_COMM_WORLD, &viscosityMatrix);
    MatSetSizes(viscosityMatrix, PETSC_DECIDE, PETSC_DECIDE, 12, 12);
    MatSetFromOptions(viscosityMatrix);
    MatSetUp(viscosityMatrix);

    getIntegrationPoints(2, gaussPoints, gaussWeights);

    gmsh::model::mesh::getBasisFunctions(9, gaussPoints, "GradLagrange", nComponents, basisFuncsGrad, nOrientations);
    gmsh::model::mesh::getJacobian(elementTag, gaussPoints, j, jdets, coords);

    vector<Mat> basisGradMats = computeBasisGradMatrix(basisFuncsGrad, j, jdets);
    vector<Mat> basisGradMatsT = vector<Mat>(basisGradMats.size());
    vector<Mat> viscMats = vector<Mat>(basisGradMats.size());

    for(int m = 0; m < basisGradMats.size(); m++){
        MatTranspose(basisGradMats[m], MAT_INITIAL_MATRIX, &basisGradMatsT[m]);
    }

    for(int m = 0; m < basisGradMats.size(); m++){
        MatMatMult(basisGradMatsT[m], basisGradMats[m], 
        MAT_INITIAL_MATRIX, PETSC_DEFAULT, &viscMats[m]);
    }

    //Enumerate the matrix
    for(int i = 0; i < 12; i++){
        for(int j = 0; j < 12; j++){
            PetscScalar viscVal = 0.0;
            int w = 0;
            //Approximate integration using gauss points
            for(int gp = 0; gp < 36; gp+=6){
                PetscScalar matVal;
                MatGetValue(viscMats[w], i, j, &matVal);
                viscVal +=  matVal * gaussWeights[w] * jdets[w++];
            }
            viscVal *= viscosity;
            MatSetValue(viscosityMatrix, i, j, viscVal, INSERT_VALUES);
        }
    }
    
    MatAssemblyBegin(viscosityMatrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(viscosityMatrix, MAT_FINAL_ASSEMBLY);

    // Print the viscosity matrix
    MatView(viscosityMatrix, PETSC_VIEWER_STDOUT_WORLD);

    return viscosityMatrix;
}

Mat computeConvectionMatrix(PetscScalar* elemVec, size_t elementTag){
    vector<double> coords;
    vector<double> gaussPoints;
    vector<double> gaussWeights;
    vector<double> basisFuncs;
    vector<double> basisFuncsGrad;
    vector<double> j;
    vector<double> jdets;
    int nComponents;
    int nOrientations;

    //Initialise Convection Matrix (components are the integrals of basis functions * velocity gradient)
    Mat convectionMatrix;
    MatCreate(PETSC_COMM_WORLD, &convectionMatrix);
    MatSetSizes(convectionMatrix, PETSC_DECIDE, PETSC_DECIDE, 12, 12);
    MatSetFromOptions(convectionMatrix);
    MatSetUp(convectionMatrix);

    getIntegrationPoints(2, gaussPoints, gaussWeights);

    gmsh::model::mesh::getBasisFunctions(9, gaussPoints, "Lagrange", nComponents, basisFuncs, nOrientations);
    gmsh::model::mesh::getBasisFunctions(9, gaussPoints, "GradLagrange", nComponents, basisFuncsGrad, nOrientations);
    gmsh::model::mesh::getJacobian(elementTag, gaussPoints, j, jdets, coords);

    vector<Mat> basisMats = getBasisMatrix(gaussPoints, basisFuncs);
    vector<Mat> basisGradMats = computeBasisGradMatrix(basisFuncsGrad, j, jdets);
    vector<Mat> convMats = vector<Mat>(basisMats.size());

    for(int m = 0; m < basisMats.size(); m++){
        MatMatMult(basisMats[m], basisGradMats[m], MAT_INITIAL_MATRIX, PETSC_DEFAULT, &convMats[m]);
    }

    for(int i = 0; i < 12; i++){
        for(int j = 0; j < 12; j++){
            PetscScalar convVal = 0.0;
            int w = 0;
            //Approximate integration using gauss points
            for(int gp = 0; gp < 36; gp+=6){
                PetscScalar matVal;
                MatGetValue(convMats[w], i, j, &matVal);
                convVal +=  matVal * elemVec[j] * gaussWeights[w] * jdets[w++];
            }
            MatSetValue(convectionMatrix, i, j, convVal, INSERT_VALUES);
        }
    }

    MatAssemblyBegin(convectionMatrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(convectionMatrix, MAT_FINAL_ASSEMBLY);

    // Print the convection matrix
    MatView(convectionMatrix, PETSC_VIEWER_STDOUT_WORLD);

    return convectionMatrix;
}

Mat computeGradientMatrix(size_t elementTag){
    vector<double> coords;
    vector<double> gaussPoints;
    vector<double> gaussWeights;
    vector<double> basisPres;
    vector<double> basisGradVelo;
    vector<double> j;
    vector<double> jdets;
    int nComponents;
    int nOrientations;

    //Initialise Gradient Matrix (components are the integrals of veloc basis grad * pressure basis)
    Mat gradientMatrix;
    MatCreate(PETSC_COMM_WORLD, &gradientMatrix);
    MatSetSizes(gradientMatrix, PETSC_DECIDE, PETSC_DECIDE, 12, 3);
    MatSetFromOptions(gradientMatrix);
    MatSetUp(gradientMatrix);

    getIntegrationPoints(2, gaussPoints, gaussWeights);

    gmsh::model::mesh::getBasisFunctions(2, gaussPoints, "Lagrange", nComponents, basisPres, nOrientations);
    gmsh::model::mesh::getBasisFunctions(9, gaussPoints, "GradLagrange", nComponents, basisGradVelo, nOrientations);
    gmsh::model::mesh::getJacobian(elementTag, gaussPoints, j, jdets, coords);

    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 12; j++){
            int w = 0;
            PetscScalar gradVal = 0.0;
            int x2 = (j > 5);
            for(int gp = 0; gp < 36; gp+=6){
                gradVal += basisPres[gp + i] * basisGradVelo[gp * 3 + j * 3 + x2] * gaussWeights[w] * jdets[w++];
            }
            MatSetValue(gradientMatrix, j, i, gradVal, INSERT_VALUES);
        }
    }

    MatAssemblyBegin(gradientMatrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(gradientMatrix, MAT_FINAL_ASSEMBLY);

    // Print the gradient matrix
    MatView(gradientMatrix, PETSC_VIEWER_STDOUT_WORLD);

    return gradientMatrix;
}

Vec computeFirstStep(Mesh *msh, size_t elementTag){
    double dt = 0.001;

    PetscScalar forceArr[12]; 
    msh->getForceVector(elementTag, forceArr);
    PetscScalar velVec[12];
    msh->getElementVector(elementTag, velVec, false);

    Mat massMat = computeMassMatrix(elementTag);
    Mat viscMat = computeViscosityMatrix(elementTag);
    Mat convMat = computeConvectionMatrix(velVec, elementTag);

    MatAXPY(convMat, 1.0, viscMat, SAME_NONZERO_PATTERN);

    Vec vint;
    Vec f;
    Vec v;
    Vec t;

    VecCreate(PETSC_COMM_WORLD, &vint);
    VecCreate(PETSC_COMM_WORLD, &f);
    VecCreate(PETSC_COMM_WORLD, &v);
    VecCreate(PETSC_COMM_WORLD, &t);
    VecSetSizes(vint, PETSC_DECIDE, 12);
    VecSetSizes(f, PETSC_DECIDE, 12);
    VecSetSizes(v, PETSC_DECIDE, 12);
    VecSetSizes(t, PETSC_DECIDE, 12);
    VecSetFromOptions(vint);
    VecSetFromOptions(f);
    VecSetFromOptions(v);
    VecSetFromOptions(t);

    PetscInt indices[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    VecSetValues(f, 12, indices, forceArr, INSERT_VALUES);
    VecSetValues(v, 12, indices, velVec, INSERT_VALUES);
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
    VecAssemblyBegin(v);
    VecAssemblyEnd(v);

    VecScale(v, 1/dt);
    MatMult(massMat, v, t);
    MatAYPX(massMat, 1/dt, convMat, SAME_NONZERO_PATTERN);
    VecAXPY(f, 1.0, t);

    KSP solver;
    KSPCreate(PETSC_COMM_WORLD, &solver);
    KSPSetOperators(solver, massMat, massMat);
    KSPSetType(solver, KSPGMRES);
    KSPSetFromOptions(solver);

    KSPSolve(solver, f, vint);

    cout << "First Step Solution: \n";
    VecView(vint, PETSC_VIEWER_STDOUT_WORLD);

    return vint;
}