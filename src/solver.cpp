#include "solver.hpp"

Solver::Solver(Mesh* msh, double dt, double viscosity){
    this->msh = msh;
    nNodes = msh->nodes.size();
    this->dt = 1/dt;
    this->viscosity = viscosity;

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

    MatAssemblyBegin(globalFullMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalFullMat, MAT_FINAL_ASSEMBLY);

    VecCreate(PETSC_COMM_WORLD, &velocityVec);
    VecSetSizes(velocityVec, PETSC_DECIDE, nNodes * 2);
    VecSetFromOptions(velocityVec);

    VecCreate(PETSC_COMM_WORLD, &nodalVec);
    VecSetSizes(nodalVec, PETSC_DECIDE, nNodes * 2 + msh->nLinear);
    VecSetFromOptions(nodalVec);
}

void Solver::applyDirichletConditions(Mat *m, Vec *v, bool expl){
    PetscInt* rows = msh->dirichletIds.data();

    for(int i = 0; i < nNodes; i++){
        Node n = msh->nodes[msh->nodeIds[i]];
        if(n.boundary || n.inlet) {
            VecSetValue(*v, i, n.velocity[0], INSERT_VALUES);
            VecSetValue(*v, i + nNodes, n.velocity[1], INSERT_VALUES);
        }
    }

    if(!expl){
        MatZeroRows(*m, msh->dirichletIds.size(), rows, 1.0, *v, *v);
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
            if(n.inlet){
                VecSetValue(*vec, x, n.velocity[0], INSERT_VALUES);
                VecSetValue(*vec, x + nNodes, n.velocity[1], INSERT_VALUES);
            }
            else{
                VecSetValue(*vec, x, 0.0, INSERT_VALUES);
                VecSetValue(*vec, x + nNodes, 0.0, INSERT_VALUES);
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
    cout << "ASSEMBLY BEGIN\n";

    localToGlobalMat(1);
    MatAssemblyBegin(globalMassMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalMassMat, MAT_FINAL_ASSEMBLY);
    MatScale(globalMassMat, dt);

    localToGlobalMat(2);
    MatAssemblyBegin(globalViscMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalViscMat, MAT_FINAL_ASSEMBLY);
    MatScale(globalViscMat, viscosity);

    localToGlobalMat(3);
    MatAssemblyBegin(globalConvMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalConvMat, MAT_FINAL_ASSEMBLY);

    localToGlobalMat(4);
    MatAssemblyBegin(globalFullMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalFullMat, MAT_FINAL_ASSEMBLY);
    Vec sol;
    VecCreate(PETSC_COMM_WORLD, &sol);
    VecSetSizes(sol, PETSC_DECIDE, nNodes * 2 + msh->nLinear);
    VecSetFromOptions(sol);
    applyDirichletConditions(&globalFullMat, &sol, false);
    VecDestroy(&sol);

    localToGlobalVec(false);
    VecAssemblyBegin(velocityVec);
    VecAssemblyEnd(velocityVec);

    localToGlobalVec(true);
    VecAssemblyBegin(nodalVec);
    VecAssemblyEnd(nodalVec);
    cout << "ASSEMBLY END\n";
}

void Solver::computeFirstStep(){
    Mat tempMat;
    Mat stabMat;
    Vec tempVec;
    Vec vint;

    MatCreate(PETSC_COMM_WORLD, &tempMat);
    MatSetSizes(tempMat, PETSC_DECIDE, PETSC_DECIDE, 
    nNodes * 2, nNodes * 2);
    MatSetFromOptions(tempMat);
    MatSetUp(tempMat);

    MatCreate(PETSC_COMM_WORLD, &stabMat);
    MatSetSizes(stabMat, PETSC_DECIDE, PETSC_DECIDE, 
    nNodes * 2, nNodes * 2);
    MatSetFromOptions(stabMat);
    MatSetUp(stabMat);

    VecCreate(PETSC_COMM_WORLD, &tempVec);
    VecCreate(PETSC_COMM_WORLD, &vint);
    VecSetSizes(tempVec, PETSC_DECIDE, nNodes * 2);
    VecSetSizes(vint, PETSC_DECIDE, nNodes * 2);
    VecSetFromOptions(tempVec);
    VecSetFromOptions(vint);

    MatConvert(globalConvMat, MATSAME, MAT_INITIAL_MATRIX, &tempMat);

    //SUPG computation
    //MatConvert(globalConvMat, MATSAME, MAT_INITIAL_MATRIX, &stabMat);
    //applyStabilisation(&stabMat);
    //MatAXPY(tempMat, 1.0, stabMat, SAME_NONZERO_PATTERN);
    //MatView(stabMat, PETSC_VIEWER_STDOUT_WORLD);
    MatDestroy(&stabMat);
    
    MatDiagonalScale(tempMat, NULL, velocityVec);
    MatAssemblyBegin(tempMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(tempMat, MAT_FINAL_ASSEMBLY);

    //Add viscous and mass matrix to system
    MatAXPY(tempMat, 1.0, globalViscMat, DIFFERENT_NONZERO_PATTERN);
    MatAXPY(tempMat, 1.0, globalMassMat, DIFFERENT_NONZERO_PATTERN);

    MatMult(globalMassMat, velocityVec, tempVec);

    //Impose Dirichlet Conditions
    applyDirichletConditions(&tempMat, &tempVec, false);

    //Solve system
    KSP solver;
    KSPCreate(PETSC_COMM_WORLD, &solver);
    KSPSetOperators(solver, tempMat, tempMat);
    KSPSetType(solver, KSPGMRES);
    KSPSetFromOptions(solver);
    KSPSolve(solver, tempVec, vint);
    VecCopy(vint, velocityVec);

    //Cleanup
    VecDestroy(&tempVec);
    MatDestroy(&tempMat);
    VecDestroy(&vint);
    KSPDestroy(&solver);
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

    applyDirichletConditions(&globalFullMat, &solVec, true);

    KSP solver;
    KSPCreate(PETSC_COMM_WORLD, &solver);
    KSPSetOperators(solver, globalFullMat, globalFullMat);
    KSPSetType(solver, KSPGMRES);
    KSPSetFromOptions(solver);

    VecZeroEntries(nodalVec);

    KSPSolve(solver, solVec, nodalVec);

    applyDirichletConditions(&globalFullMat, &nodalVec, true);
    PetscScalar max = 0;
    for(int i = 0; i < nNodes * 2; i++){
        PetscInt ind = i;
        PetscScalar val;
        VecGetValues(nodalVec, 1, &ind, &val);
        max = val > max ? val : max;
        VecSetValue(velocityVec, i, val, INSERT_VALUES);
    }
    cout << "MAX (instability metric): " << max << "\n";
    VecAssemblyBegin(velocityVec);
    VecAssemblyEnd(velocityVec);

    VecDestroy(&solVec);
    VecDestroy(&tempVec);
    KSPDestroy(&solver);
}

vector<vector<double>> Solver::computeTimeStep(int steps){
    vector<vector<double>> fluid = vector<vector<double>>(3);
    for(int x = 0; x < steps; x++){
        this->computeFirstStep();
        this->computeSecondStep();
        cout << "Step: " << x << "\n";
    }

    for(int i = 0; i < nNodes; i++){
        PetscInt iu = i;
        PetscInt iv = i + nNodes;
        PetscScalar u;
        PetscScalar v;

        VecGetValues(nodalVec, 1, &iu, &u);
        VecGetValues(nodalVec, 1, &iv, &v);

        fluid[0].push_back(u);
        fluid[1].push_back(v);
        fluid[2].push_back(-1);
    }

    for(int i = 0; i < nNodes; i++){
        PetscInt ip = i + nNodes * 2;
        PetscScalar p;
        VecGetValues(nodalVec, 1, &ip, &p);
        if(msh->nodes[msh->nodeIds[i]].pid != -1){
            fluid[2][i] = p;
        }
    }
    return fluid;
}