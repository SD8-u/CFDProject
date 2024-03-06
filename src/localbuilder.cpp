#include "localbuilder.hpp"

void cleanUp(vector<Mat> matrices){
    for(Mat m : matrices){
        MatDestroy(&m);
    }
}

void LocalBuilder::setUp(){
    inverseJacobian = vector<Mat>(6);
    basisGradMats = vector<Mat>(6);
    basisMats = vector<Mat>(6);
    j = vector<double>();
    jdets = vector<double>();
    basisFuncs = vector<double>();
    basisFuncsGrad = vector<double>();

    int nComponents;
    int nOrientations;

    for(int m = 0; m < 6; m++){
        MatCreate(PETSC_COMM_WORLD, &inverseJacobian[m]);
        MatSetSizes(inverseJacobian[m], PETSC_DECIDE, PETSC_DECIDE, 2, 2);
        MatSetFromOptions(inverseJacobian[m]);
        MatSetUp(inverseJacobian[m]);

        MatCreate(PETSC_COMM_WORLD, &basisGradMats[m]);
        MatSetSizes(basisGradMats[m], PETSC_DECIDE, PETSC_DECIDE, 2, 12);
        MatSetFromOptions(basisGradMats[m]);
        MatSetUp(basisGradMats[m]);

        MatCreate(PETSC_COMM_WORLD, &basisMats[m]);
        MatSetSizes(basisMats[m], PETSC_DECIDE, PETSC_DECIDE, 12, 2);
        MatSetFromOptions(basisMats[m]);
        MatSetUp(basisMats[m]);
    }

    gmsh::model::mesh::getBasisFunctions(9, gaussPoints, "Lagrange", 
    nComponents, basisFuncs, nOrientations);
    gmsh::model::mesh::getBasisFunctions(2, gaussPoints, "Lagrange", 
    nComponents, basisFuncsPres, nOrientations);
    gmsh::model::mesh::getBasisFunctions(9, gaussPoints, "GradLagrange", 
    nComponents, basisFuncsGrad, nOrientations);
    buildBasisMatrix();
}

LocalBuilder::LocalBuilder(){
    this->conv = true;
    setUp();
    MatCreate(PETSC_COMM_WORLD, &localConvMat);
    MatSetSizes(localConvMat, PETSC_DECIDE, PETSC_DECIDE, 12, 12);
    MatSetFromOptions(localConvMat);
    MatSetUp(localConvMat);
}

LocalBuilder::LocalBuilder(double dt, double viscosity){
    this->dt = dt;
    this->viscosity = viscosity;
    this->conv = false;
    setUp();
    MatCreate(PETSC_COMM_WORLD, &localMassMat);
    MatSetSizes(localMassMat, PETSC_DECIDE, PETSC_DECIDE, 12, 12);
    MatSetFromOptions(localMassMat);
    MatSetUp(localMassMat);
    MatDuplicate(localMassMat, MAT_DO_NOT_COPY_VALUES, &localViscMat);

    MatCreate(PETSC_COMM_WORLD, &localGradMat);
    MatSetSizes(localGradMat, PETSC_DECIDE, PETSC_DECIDE, 12, 3);
    MatSetFromOptions(localGradMat);
    MatSetUp(localGradMat);

    MatCreate(PETSC_COMM_WORLD, &localFullMat);
    MatCreate(PETSC_COMM_WORLD, &localFullMat);
    MatSetSizes(localFullMat, PETSC_DECIDE, PETSC_DECIDE, 15, 15);
    MatSetFromOptions(localFullMat);
    MatSetUp(localFullMat);
}

LocalBuilder::~LocalBuilder(){
    if(!conv){
        MatDestroy(&localGradMat);
        MatDestroy(&localMassMat);
        MatDestroy(&localViscMat);
        MatDestroy(&localFullMat);
    }
    else{
        MatDestroy(&localConvMat);
    }
    cleanUp(basisMats);
    cleanUp(inverseJacobian);
    cleanUp(basisGradMats);
}

void LocalBuilder::buildBasisMatrix(){
    int m = 0;
    for(int point = 0; point < basisFuncs.size(); point+=6){
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
}

void LocalBuilder::buildInverseJacobian(){
    //Inverse jacobian matrix for each integration point
    int m = 0;
    for(int point = 0; point < j.size(); point+=9){
        MatZeroEntries(inverseJacobian[m]);
        MatSetValue(inverseJacobian[m], 0, 0, j[point + 4] * 1/jdets[m], INSERT_VALUES);
        MatSetValue(inverseJacobian[m], 0, 1, -j[point + 1] * 1/jdets[m], INSERT_VALUES);
        MatSetValue(inverseJacobian[m], 1, 0, -j[point + 3] * 1/jdets[m], INSERT_VALUES);
        MatSetValue(inverseJacobian[m], 1, 1, j[point] * 1/jdets[m], INSERT_VALUES);

        MatAssemblyBegin(inverseJacobian[m], MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(inverseJacobian[m], MAT_FINAL_ASSEMBLY);
        m++;
    }
}

//Compute basis function gradient matrix w.r.t spatial coordinates
void LocalBuilder::buildBasisGradMatrix(){

    int m = 0;
    //Compute matrix for each integration point
    for(int point = 0; point < basisFuncsGrad.size(); point+=18){
        Mat temp;
        MatCreate(PETSC_COMM_WORLD, &temp);
        MatSetSizes(temp, PETSC_DECIDE, PETSC_DECIDE, 2, 12);
        MatSetFromOptions(temp);
        MatSetUp(temp);

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
        MatDestroy(&basisGradMats[m]);
        MatMatMult(inverseJacobian[m], temp, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &basisGradMats[m]);
        MatDestroy(&temp);

        MatAssemblyBegin(basisGradMats[m], MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(basisGradMats[m], MAT_FINAL_ASSEMBLY);
        m++;
    }
}

void LocalBuilder::computeMassMatrix(){
    MatZeroEntries(localMassMat);

    //Enumerate the matrix
    for(int i = 0; i < 12; i++){
        for(int j = 0; j < 12; j++){
            PetscScalar massVal = 0.0;
            int w = 0;
            int offseti = (i > 5) ? 6 : 0;
            int offsetj = (j > 5) ? 6 : 0;
            
            if(j < 6 && i < 6 || j > 5 && i > 5){
                //Approximate integration using gauss points over the products of basis functions
                for(int gp = 0; gp < 36; gp+=6){
                    massVal += basisFuncs[(i-offseti) + gp] * basisFuncs[(j-offsetj) + gp] 
                    * gaussWeights[w] * jdets[w++];
                }
            }
            MatSetValue(localMassMat, i, j, massVal, INSERT_VALUES);
        }
    }

    MatAssemblyBegin(localMassMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(localMassMat, MAT_FINAL_ASSEMBLY);
}

void LocalBuilder::computeViscosityMatrix(){
    MatZeroEntries(localViscMat);

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
            MatSetValue(localViscMat, i, j, viscVal, INSERT_VALUES);
        }
    }
    
    MatAssemblyBegin(localViscMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(localViscMat, MAT_FINAL_ASSEMBLY);
    // Clean up matrices
    cleanUp(basisGradMatsT);
    cleanUp(viscMats);
}

void LocalBuilder::computeConvectionMatrix(size_t elementTag, Vec *velocityVec){
    this->elementTag = elementTag;
    vector<double> coords;
    gmsh::model::mesh::getJacobian(elementTag, gaussPoints, j, jdets, coords);
    buildInverseJacobian();
    MatZeroEntries(localConvMat);
    buildBasisGradMatrix();

    PetscScalar v[12] = {0};
    for(int a = 0; a < 12; a++){
        VecGetValues(*velocityVec, 1, &a, &v[a]);
    }

    for(int i = 0; i < 12; i++){
        for(int j = 0; j < 12; j++){
            PetscScalar convVal = 0.0;
            int w = 0;
            //Approximate integration using gauss points
            for(int gp = 0; gp < 36; gp+=6){
                PetscScalar basisVal, gradX, gradY;
                MatGetValue(basisMats[w], i % 6, 0, &basisVal);
                MatGetValue(basisGradMats[w], 0, j, &gradX);
                MatGetValue(basisGradMats[w], 1, j, &gradY);

                double velU = 0, velV = 0;
                velU = v[0] * basisFuncs[gp] + v[1] * basisFuncs[gp + 1] + 
                v[2] * basisFuncs[gp + 2] + v[3] * basisFuncs[gp + 3] + 
                v[4] * basisFuncs[gp + 4] + v[5] * basisFuncs[gp + 5];
    
                
                velV = v[6] * basisFuncs[gp] + v[7] * basisFuncs[gp + 1] + 
                v[8] * basisFuncs[gp + 2] + v[9] * basisFuncs[gp + 3] + 
                v[10] * basisFuncs[gp + 4] + v[11] * basisFuncs[gp + 5];
                
                PetscScalar matVal;
                matVal = basisVal * (velU * gradX + velV * gradY);
                convVal += (matVal * gaussWeights[w] * jdets[w++]);
            }
            
            MatSetValue(localConvMat, i, j, (convVal), INSERT_VALUES);
        }
    }

    MatAssemblyBegin(localConvMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(localConvMat, MAT_FINAL_ASSEMBLY);
}

void LocalBuilder::computeGradientMatrix(){
    MatZeroEntries(localGradMat);

    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 12; j++){
            int w = 0;
            PetscScalar gradVal = 0.0;
            int x = (j > 5);
            for(int gp = 0; gp < 18; gp+=3){
                PetscScalar matVal = 0.0;
                MatGetValue(basisGradMats[w], x, j, &matVal);
                gradVal += basisFuncsPres[gp + i] * matVal * gaussWeights[w] * jdets[w++];
            }
            MatSetValue(localGradMat, j, i, gradVal * -1, INSERT_VALUES);
        }
    }

    MatAssemblyBegin(localGradMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(localGradMat, MAT_FINAL_ASSEMBLY);
}

void LocalBuilder::computeFinalMatrix(){
    MatZeroEntries(localFullMat);
    MatScale(localMassMat, dt);
    //MatScale(localViscMat, viscosity);
    Mat localGradTMat;

    MatTranspose(localGradMat, MAT_INITIAL_MATRIX, &localGradTMat);

    for(int i = 0; i < 12; i++){
        for(int j = 0; j < 12; j++){
            PetscScalar massVal;
            PetscScalar viscVal;
            MatGetValue(localMassMat, i, j, &massVal);
            MatGetValue(localViscMat, i, j, &viscVal);
            MatSetValue(localFullMat, i, j, massVal, INSERT_VALUES);
        }
    }

    for(int i = 0; i < 12; i++){
        for(int j = 12; j < 15; j++){
            PetscScalar matVal;
            MatGetValue(localGradMat, i, j - 12, &matVal);
            MatSetValue(localFullMat, i, j, matVal, INSERT_VALUES);
        }
    }

    for(int i = 12; i < 15; i++){
        for(int j = 0; j < 12; j++){
            PetscScalar matVal;
            MatGetValue(localGradTMat, i - 12, j, &matVal);
            MatSetValue(localFullMat, i, j, matVal, INSERT_VALUES);
        }
    }

    for(int i = 12; i < 15; i++){
        for(int j = 12; j < 15; j++){
            MatSetValue(localFullMat, i, j, 0, INSERT_VALUES);
        }
    }

    MatAssemblyBegin(localFullMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(localFullMat, MAT_FINAL_ASSEMBLY);

    MatDestroy(&localGradTMat);
    MatScale(localMassMat, 1/dt);
    //MatScale(localViscMat, 1/viscosity);
}

void LocalBuilder::assembleMatrices(size_t elementTag){
    this->elementTag = elementTag;
    vector<double> coords;

    gmsh::model::mesh::getJacobian(elementTag, gaussPoints, j, jdets, coords);

    buildInverseJacobian();
    buildBasisGradMatrix();

    computeMassMatrix();
    computeViscosityMatrix();
    computeGradientMatrix();
    computeFinalMatrix();
}