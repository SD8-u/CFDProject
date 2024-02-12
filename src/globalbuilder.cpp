#include "globalbuilder.hpp"

GlobalBuilder::GlobalBuilder(int dim, double dt, double visc, Mesh* msh){
    this->dt = 1/dt;
    this->viscosity = visc;
    this->msh = msh;

    MatCreate(PETSC_COMM_WORLD, &globalMassMat);
    MatSetSizes(globalMassMat, PETSC_DECIDE, PETSC_DECIDE, 
    msh->nNodes * 2, msh->nNodes * 2);
    MatSetFromOptions(globalMassMat);
    MatSetUp(globalMassMat);

    MatCreate(PETSC_COMM_WORLD, &globalFullMat);
    MatSetSizes(globalFullMat, PETSC_DECIDE, PETSC_DECIDE, 
    msh->nNodes * 2 + msh->nLinear, msh->nNodes * 2 + msh->nLinear);
    MatSetFromOptions(globalFullMat);
    MatSetUp(globalFullMat);

    MatDuplicate(globalMassMat, MAT_DO_NOT_COPY_VALUES, &globalViscMat);
    MatDuplicate(globalMassMat, MAT_DO_NOT_COPY_VALUES, &globalConvMat);

    MatZeroEntries(globalMassMat);
    MatZeroEntries(globalViscMat);
    MatZeroEntries(globalConvMat);
    MatZeroEntries(globalFullMat);

    VecCreate(PETSC_COMM_WORLD, &velocityVec);
    VecSetSizes(velocityVec, PETSC_DECIDE, msh->nNodes * 2);
    VecSetFromOptions(velocityVec);

    VecCreate(PETSC_COMM_WORLD, &nodalVec);
    VecSetSizes(nodalVec, PETSC_DECIDE, msh->nNodes * 2 + msh->nLinear);
    VecSetFromOptions(nodalVec);
}

GlobalBuilder::~GlobalBuilder(){
    MatDestroy(&globalMassMat);
    MatDestroy(&globalViscMat);
    MatDestroy(&globalConvMat);
    MatDestroy(&globalFullMat);
    VecDestroy(&velocityVec);
    VecDestroy(&nodalVec);
}

void GlobalBuilder::localToGlobalMat(size_t elementTag, Mat *localMat, Mat *globalMat, bool final=false){
        int add = final ? 3 : 0;
        for(int i = 0; i < 12 + add; i++){
            int x, y;
            x = msh->nodes[msh->elements[elementTag][i % 6]].id;
            if(i > 5){
                x += msh->nNodes;
            }
            if(i > 11){
                x = msh->nodes[msh->elements[elementTag][i % 6]].pid;
                x += msh->nNodes * 2;
            }
            for(int j = 0; j < 12 + add; j++){
                y = msh->nodes[msh->elements[elementTag][j % 6]].id;
                if(j > 5){
                    y += msh->nNodes;
                }
                if(j > 11){
                    y = msh->nodes[msh->elements[elementTag][j % 6]].pid;
                    y += msh->nNodes * 2;
                }
                PetscScalar matVal;

                MatGetValue(*localMat, i, j, &matVal);
                MatSetValue(*globalMat, x, y, matVal, ADD_VALUES);
            }
        }
}

void GlobalBuilder::localToGlobalVec(bool full){
    Vec* vec = &velocityVec;
    if(full){
        vec = &nodalVec;
    }
    for(size_t elementTag : msh->elementTags[0]){
        for(int i = 0; i < 6; i++){
            Node n = msh->nodes[msh->elements[elementTag][i]];
            int x = n.id;
            int y = n.pid;
            if(n.inlet){
                VecSetValue(*vec, x, n.velocity[0], INSERT_VALUES);
                VecSetValue(*vec, x + msh->nNodes, n.velocity[1], INSERT_VALUES);
            }
            else{
                VecSetValue(*vec, x, 0.0, INSERT_VALUES);
                VecSetValue(*vec, x + msh->nNodes, 0.0, INSERT_VALUES);
            }
            if(i < 3 && full){
                VecSetValue(*vec, 2 * msh->nNodes + y, n.pressure, INSERT_VALUES);
            }
        }
    }
    VecAssemblyBegin(*vec);
    VecAssemblyEnd(*vec);
}

void GlobalBuilder::assembleMatrices(){
    for(size_t elementTag : msh->elementTags[0]){
        localBuild = new LocalBuilder(dt);
        localBuild->assembleMatrices(elementTag);
        localToGlobalMat(elementTag, &localBuild->localMassMat, &globalMassMat);
        localToGlobalMat(elementTag, &localBuild->localViscMat, &globalViscMat);
        localToGlobalMat(elementTag, &localBuild->localConvMat, &globalConvMat);
        localToGlobalMat(elementTag, &localBuild->localFullMat, &globalFullMat, true);
        delete(localBuild);
    }
    MatAssemblyBegin(globalMassMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalMassMat, MAT_FINAL_ASSEMBLY);

    MatAssemblyBegin(globalViscMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalViscMat, MAT_FINAL_ASSEMBLY);

    MatAssemblyBegin(globalConvMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalConvMat, MAT_FINAL_ASSEMBLY);

    MatAssemblyBegin(globalFullMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalFullMat, MAT_FINAL_ASSEMBLY);

    MatScale(globalViscMat, viscosity);
    MatScale(globalMassMat, dt);
}

void GlobalBuilder::assembleVectors(){
    localToGlobalVec(false);
    localToGlobalVec(true);
}

void GlobalBuilder::updateVelocity(){
    PetscScalar max = 0;
    for(int i = 0; i < msh->nNodes * 2; i++){
        PetscInt ind = i;
        PetscScalar val;
        VecGetValues(nodalVec, 1, &ind, &val);
        max = val > max ? val : max;
        VecSetValue(velocityVec, i, val, INSERT_VALUES);
    }
    cout << "MAX (instability metric): " << max << "\n";
    VecAssemblyBegin(velocityVec);
    VecAssemblyEnd(velocityVec);
}