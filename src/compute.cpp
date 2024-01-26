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
    MatSetSizes(massMatrix, PETSC_DECIDE, PETSC_DECIDE, 6, 6);
    MatSetFromOptions(massMatrix);
    MatSetUp(massMatrix);

    getIntegrationPoints(2, gaussPoints, gaussWeights);

    gmsh::model::mesh::getBasisFunctions(9, gaussPoints, "Lagrange", nComponents, basisFuncs, nOrientations);
    gmsh::model::mesh::getJacobian(elementTag, gaussPoints, j, jdets, coords);

    PetscScalar massValues[36] = {0.0};
    PetscInt row[36];
    PetscInt col[36];

    //Enumerate the matrix
    for(int i = 0; i < 6; i++){
        for(int j = 0; j < 6; j++){
            row[i * 6 + j] = i;
            col[i * 6 + j] = j;
            PetscScalar massVal = 0.0;
            int w = 0;
            //Approximate integration using gauss points over the products of basis functions
            for(int gp = 0; gp < 36; gp+=6){
                massVal += basisFuncs[i + gp] * basisFuncs[j + gp] * gaussWeights[w] * jdets[w++];
            }
            MatSetValue(massMatrix, i, j, massVal, INSERT_VALUES);
        }
    }

    MatAssemblyBegin(massMatrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(massMatrix, MAT_FINAL_ASSEMBLY);

    // Print the mass matrix
    MatView(massMatrix, PETSC_VIEWER_STDOUT_WORLD);
    MatDestroy(&massMatrix);

    return massMatrix;
}