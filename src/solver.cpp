#include <solver.hpp>
#include <compute.hpp>

Solver::Solver(Mesh* msh){
    this->msh = msh;
    nNodes = msh->nodes.size();

    MatCreate(PETSC_COMM_WORLD, &globalMassMat);
    MatCreate(PETSC_COMM_WORLD, &globalViscMat);
    MatCreate(PETSC_COMM_WORLD, &globalConvMat);
    MatCreate(PETSC_COMM_WORLD, &globalFullMat);

    MatSetSizes(globalMassMat, PETSC_DECIDE, PETSC_DECIDE, 
    nNodes * 2, nNodes * 2);
    MatSetSizes(globalConvMat, PETSC_DECIDE, PETSC_DECIDE, 
    nNodes * 2, nNodes * 2);
    MatSetSizes(globalViscMat, PETSC_DECIDE, PETSC_DECIDE, 
    nNodes * 2, nNodes * 2);
    MatSetSizes(globalFullMat, PETSC_DECIDE, PETSC_DECIDE, 
    nNodes * 2 + nNodes/2, nNodes * 2 + nNodes/2);

    MatSetFromOptions(globalMassMat);
    MatSetFromOptions(globalViscMat);
    MatSetFromOptions(globalConvMat);
    MatSetFromOptions(globalFullMat);

    MatSetUp(globalMassMat);
    MatSetUp(globalViscMat);
    MatSetUp(globalConvMat);
    MatSetUp(globalFullMat);

    MatZeroEntries(globalMassMat);
    MatZeroEntries(globalViscMat);
    MatZeroEntries(globalConvMat);
    MatZeroEntries(globalFullMat);
}

void Solver::localToGlobal(Mesh* msh, int type){
    for(size_t elementTag : msh->elementTags[0]){

        Mat localMat;
        Mat globalMat;
        switch(type){
            case 1:
                localMat = computeMassMatrix(elementTag);
                globalMat = globalMassMat;
            case 2:
                localMat = computeViscosityMatrix(elementTag);
                globalMat = globalViscMat;
            case 3:
                localMat = computeConvectionMatrix(elementTag);
                globalMat = globalConvMat;
        }
        for(int i = 0; i < 12; i++){
            int x, y;
            x = msh->nodes[msh->elements[elementTag][i % 6]].id;
            if(i > 5){
                x += nNodes;
            }
            for(int j = 0; j < 12; j++){
                y = msh->nodes[msh->elements[elementTag][j % 6]].id;
                if(j > 5){
                    y += nNodes;
                }
                PetscScalar matVal;

                MatGetValue(localMat, i, j, &matVal);
                MatSetValue(globalMat, x, y, matVal, ADD_VALUES);
            }
        }

        MatDestroy(&localMat);
    }
}

void Solver::assembleMatrices(){

    localToGlobal(msh, 1);
    MatScale(globalMassMat, 1/0.0001);
    
    MatAssemblyBegin(globalViscMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalViscMat, MAT_FINAL_ASSEMBLY);

    localToGlobal(msh, 2);
    MatAssemblyBegin(globalViscMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalViscMat, MAT_FINAL_ASSEMBLY);


    localToGlobal(msh, 3);
    MatAssemblyBegin(globalConvMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalConvMat, MAT_FINAL_ASSEMBLY);
}