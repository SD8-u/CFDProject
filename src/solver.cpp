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
    nNodes * 2 + msh->nLinear, nNodes * 2 + msh->nLinear);

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

    //VERY SLOW **FIX**
    for(int x = 0; x < nNodes * 2 + msh->nLinear; x++){
        for(int y = 0; y < nNodes * 2 + msh->nLinear; y++){
            MatSetValue(globalFullMat, x, y, 0, ADD_VALUES);
        }
    }

    VecCreate(PETSC_COMM_WORLD, &velocityVec);
    VecSetSizes(velocityVec, PETSC_DECIDE, nNodes * 2);
    VecSetFromOptions(velocityVec);

    VecCreate(PETSC_COMM_WORLD, &nodalVec);
    VecSetSizes(nodalVec, PETSC_DECIDE, nNodes * 2 + msh->nLinear);
    VecSetFromOptions(nodalVec);
}

void Solver::applyDirichletConditions(Mat *m, Vec *v, bool expl){
    for(int i = 0; i < nNodes; i++){
        Node n = msh->nodes[msh->nodeIds[i]];
        if(n.boundary || n.inlet) {
            if(!expl){
                for(int v = 0; v < nNodes + 1; v+=nNodes){
                    for(int j = 0; j < nNodes * 2; j++){
                        PetscScalar val = (j == i + v);
                        MatSetValue(*m, i + v, j, val, INSERT_VALUES);
                    }   
                }
            }
            VecSetValue(*v, i, 0.0, INSERT_VALUES);
            VecSetValue(*v, i + nNodes, 0.0, INSERT_VALUES);
        }
    }
    MatAssemblyBegin(*m, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*m, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(*v);
    VecAssemblyEnd(*v);
}

void Solver::localToGlobalVec(bool full){
    Vec* vec = &velocityVec;
    if(full){
        vec = &nodalVec;
    }
    for(size_t elementTag : msh->elementTags[0]){
        for(int i = 0; i < 6; i++){
            Node n = msh->nodes[msh->elements[elementTag][i]];
            int x = n.id;
            int y = n.pid;
            if(n.boundary){
                VecSetValue(*vec, x, 0, INSERT_VALUES);
                VecSetValue(*vec, x + nNodes, 0, INSERT_VALUES);
            }
            else{
                VecSetValue(*vec, x, n.velocity[0], INSERT_VALUES);
                VecSetValue(*vec, x + nNodes, n.velocity[1], INSERT_VALUES);
            }
            if(i < 3 && full){
                VecSetValue(*vec, 2 * nNodes + y, n.pressure, INSERT_VALUES);
            }
        }
    }
}

void Solver::localToGlobalMat(int type){
    for(size_t elementTag : msh->elementTags[0]){

        Mat localMat;
        Mat* globalMat;

        int add = 0;
        switch(type){
            case 1:
                localMat = computeMassMatrix(elementTag);
                globalMat = &globalMassMat;
                break;
            case 2:
                localMat = computeViscosityMatrix(elementTag);
                globalMat = &globalViscMat;
                break;
            case 3:
                localMat = computeConvectionMatrix(elementTag);
                globalMat = &globalConvMat;
                break;
            case 4:
                localMat = computeFinalMatrix(elementTag, dt);
                globalMat = &globalFullMat;
                add = 3;
                break;
        }
        for(int i = 0; i < 12 + add; i++){
            int x, y;
            x = msh->nodes[msh->elements[elementTag][i % 6]].id;
            if(i > 5){
                x += nNodes;
            }
            if(i > 11){
                x = msh->nodes[msh->elements[elementTag][i % 6]].pid;
                x += nNodes * 2;
            }
            for(int j = 0; j < 12 + add; j++){
                y = msh->nodes[msh->elements[elementTag][j % 6]].id;
                if(j > 5){
                    y += nNodes;
                }
                if(j > 11){
                    y = msh->nodes[msh->elements[elementTag][j % 6]].pid;
                    y += nNodes * 2;
                }
                PetscScalar matVal;

                MatGetValue(localMat, i, j, &matVal);
                MatSetValue(*globalMat, x, y, matVal, ADD_VALUES);
            }
        }

        MatDestroy(&localMat);
    }
}

void Solver::assembleMatrices(){
    cout << "ASSEMBLY\n";

    localToGlobalMat(1);
    MatAssemblyBegin(globalMassMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalMassMat, MAT_FINAL_ASSEMBLY);
    MatScale(globalMassMat, dt);

    localToGlobalMat(2);
    MatAssemblyBegin(globalViscMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalViscMat, MAT_FINAL_ASSEMBLY);

    localToGlobalMat(3);
    MatAssemblyBegin(globalConvMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalConvMat, MAT_FINAL_ASSEMBLY);

    localToGlobalMat(4);
    MatAssemblyBegin(globalFullMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalFullMat, MAT_FINAL_ASSEMBLY);

    localToGlobalVec(false);
    VecAssemblyBegin(velocityVec);
    VecAssemblyEnd(velocityVec);

    localToGlobalVec(true);
    VecAssemblyBegin(nodalVec);
    VecAssemblyEnd(nodalVec);
}

void Solver::computeFirstStep(){
    Mat tempMat;
    Vec tempVec;
    Vec vint;

    MatCreate(PETSC_COMM_WORLD, &tempMat);
    MatSetSizes(tempMat, PETSC_DECIDE, PETSC_DECIDE, 
    nNodes * 2, nNodes * 2);
    MatSetFromOptions(tempMat);
    MatSetUp(tempMat);

    VecCreate(PETSC_COMM_WORLD, &tempVec);
    VecCreate(PETSC_COMM_WORLD, &vint);
    VecSetSizes(tempVec, PETSC_DECIDE, nNodes * 2);
    VecSetSizes(vint, PETSC_DECIDE, nNodes * 2);
    VecSetFromOptions(tempVec);
    VecSetFromOptions(vint);

    //VERY SLOW **FIX**
    for(int i = 0; i < nNodes * 2; i++){
        for(int j = 0; j < nNodes * 2; j++){
            PetscInt ind = j;
            PetscScalar vecVal;
            PetscScalar matVal;
            VecGetValues(velocityVec, 1, &ind, &vecVal);
            MatGetValue(globalConvMat, i, j, &matVal);
            MatSetValue(tempMat, i, j, matVal * vecVal, INSERT_VALUES);
        }
    }

    MatAssemblyBegin(tempMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(tempMat, MAT_FINAL_ASSEMBLY);

    MatAXPY(tempMat, 1.0, globalViscMat, DIFFERENT_NONZERO_PATTERN);
    MatAXPY(tempMat, 1.0, globalMassMat, DIFFERENT_NONZERO_PATTERN);

    MatMult(globalMassMat, velocityVec, tempVec);

    //Impose Dirichlet Conditions
    applyDirichletConditions(&tempMat, &tempVec, false);

    KSP solver;
    KSPCreate(PETSC_COMM_WORLD, &solver);
    KSPSetOperators(solver, tempMat, tempMat);
    KSPSetType(solver, KSPGMRES);
    KSPSetFromOptions(solver);

    KSPSolve(solver, tempVec, vint);

    VecDestroy(&tempVec);
    MatDestroy(&tempMat);

    VecCopy(velocityVec, vint);
    VecDestroy(&vint);
}

void Solver::computeSecondStep(){
    Vec tempVec;
    Vec solVec;
    VecCreate(PETSC_COMM_WORLD, &tempVec);
    VecCreate(PETSC_COMM_WORLD, &solVec);
    VecSetSizes(tempVec, PETSC_DECIDE, nNodes * 2);
    VecSetSizes(solVec, PETSC_DECIDE, nNodes * 2 + msh->nLinear);

    VecSetFromOptions(tempVec);
    VecSetFromOptions(solVec);

    MatMult(globalMassMat, velocityVec, tempVec);

    VecZeroEntries(solVec);

    for(int i = 0; i < nNodes * 2; i++){
        PetscInt ind = i;
        PetscScalar val;
        VecGetValues(tempVec, 1, &ind, &val);
        VecSetValue(solVec, i, val, INSERT_VALUES);
    }

    VecAssemblyBegin(solVec);
    VecAssemblyEnd(solVec);

    KSP solver;
    KSPCreate(PETSC_COMM_WORLD, &solver);
    KSPSetOperators(solver, globalFullMat, globalFullMat);
    KSPSetType(solver, KSPGMRES);
    KSPSetFromOptions(solver);

    KSPSolve(solver, solVec, nodalVec);

    for(int i = 0; i < nNodes * 2; i++){
        PetscInt ind = i;
        PetscScalar val;
        VecGetValues(nodalVec, 1, &ind, &val);
        VecSetValue(velocityVec, i, val, INSERT_VALUES);
    }

    VecAssemblyBegin(velocityVec);
    VecAssemblyEnd(velocityVec);

    applyDirichletConditions(&globalFullMat, &nodalVec, true);
    applyDirichletConditions(&globalFullMat, &velocityVec, true);

    VecDestroy(&solVec);
    VecDestroy(&tempVec);
}

void Solver::computeTimeStep(){
    for(int x = 0; x < 1000; x++){
        this->computeFirstStep();
        this->computeSecondStep();
    }

    VecView(nodalVec, PETSC_VIEWER_STDOUT_WORLD);

}