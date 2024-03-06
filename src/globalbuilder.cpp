#include "globalbuilder.hpp"
#include <omp.h>

GlobalBuilder::GlobalBuilder(int dim, double dt, double visc, Mesh* msh){
    this->dt = 1/dt;
    this->viscosity = visc;
    this->msh = msh;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_split(MPI_COMM_WORLD, rank, 0, &this->comm);

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
        PetscInt m, n;
        MatGetOwnershipRange(*globalMat, &m, &n);
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
                if(y >= m && y <= n){
                    PetscScalar matVal;
                    MatGetValue(*localMat, i, j, &matVal);
                    MatSetValue(*globalMat, x, y, matVal, ADD_VALUES);
                }
            }
        }
}

void GlobalBuilder::localToGlobalVec(bool full){
    Vec* vec = &velocityVec;
    if(full){
        vec = &nodalVec;
    }

    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int domainSize = (msh->elementTags[0].size()/size);
    int domainStart = rank * domainSize;
    int domainEnd = domainStart + domainSize + 
    (msh->elementTags[0].size()%size) * (rank == size-1);

    for(int e = 0; e < domainEnd; e++){
        size_t elementTag = msh->elementTags[0][e];
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

void GlobalBuilder::globalToLocalVec(size_t elementTag, Vec *localVec){
    for(int node = 0; node < 6; node++){
        PetscInt ix = msh->nodes[msh->elements[elementTag][node]].id;
        PetscInt iy = ix + msh->nNodes;
        PetscScalar vx, vy;
        VecGetValues(velocityVec, 1, &ix, &vx);
        VecGetValues(velocityVec, 1, &iy, &vy);
        VecSetValue(*localVec, node, vx, INSERT_VALUES);
        VecSetValue(*localVec, node + 6, vy, INSERT_VALUES);
    }
    VecAssemblyBegin(*localVec);
    VecAssemblyEnd(*localVec);
}

void GlobalBuilder::assembleMatrices(){
    localBuild = new LocalBuilder(dt, viscosity, comm);
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int domainSize = (msh->elementTags[0].size()/size);
    int domainStart = rank * domainSize;
    int domainEnd = domainStart + domainSize + 
    (msh->elementTags[0].size()%size) * (rank == size-1);

    cout << "Size: "<< size << "\n";
    cout << "Rank: "<< rank << "\n";
    cout << "domain size: " << domainSize << "\n";
    cout << "domain start: " << domainStart << "\n";
    cout << "domain end: " << domainEnd << "\n";

    for(int e = 0; e < domainEnd; e++){

        size_t elementTag = msh->elementTags[0][e];
        localBuild->assembleMatrices(elementTag);

        localToGlobalMat(elementTag, &localBuild->localMassMat, &globalMassMat);
        localToGlobalMat(elementTag, &localBuild->localViscMat, &globalViscMat);
        localToGlobalMat(elementTag, &localBuild->localFullMat, &globalFullMat, true);
    }
    delete(localBuild);
    MatAssemblyBegin(globalMassMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalMassMat, MAT_FINAL_ASSEMBLY);

    MatAssemblyBegin(globalViscMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalViscMat, MAT_FINAL_ASSEMBLY);

    MatAssemblyBegin(globalFullMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalFullMat, MAT_FINAL_ASSEMBLY);

    MatScale(globalViscMat, viscosity);
    MatScale(globalMassMat, dt);
}

void GlobalBuilder::assembleConvectionMatrix(){
    MatZeroEntries(globalConvMat);
    Vec localVelVec;
    VecCreate(comm, &localVelVec);
    VecSetSizes(localVelVec, PETSC_DECIDE, 12);
    VecSetFromOptions(localVelVec);

    localBuild = new LocalBuilder(comm);

    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int domainSize = (msh->elementTags[0].size()/size);
    int domainStart = rank * domainSize;
    int domainEnd = domainStart + size + 
    (msh->elementTags[0].size()%size) * (rank == size-1);

    for(int e = 0; e < domainEnd; e++){
        size_t elementTag = msh->elementTags[0][e];

        globalToLocalVec(elementTag, &localVelVec);
        localBuild->computeConvectionMatrix(elementTag, &localVelVec);
        localToGlobalMat(elementTag, &localBuild->localConvMat, &globalConvMat);
    }

    MatAssemblyBegin(globalConvMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalConvMat, MAT_FINAL_ASSEMBLY);
    delete(localBuild);
    VecDestroy(&localVelVec);
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