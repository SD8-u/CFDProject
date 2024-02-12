#include "solver.hpp"

Solver::Solver(Mesh* msh, double dt, double viscosity){
    Vec tempVec;
    this->msh = msh;

    globalBuild = new GlobalBuilder(2, dt, viscosity, msh);
    globalBuild->assembleMatrices();
    globalBuild->assembleVectors();

    VecCreate(PETSC_COMM_WORLD, &tempVec);
    VecSetSizes(tempVec, PETSC_DECIDE, msh->nNodes * 2 + msh->nLinear);
    VecSetFromOptions(tempVec);

    applyDirichletConditions(&globalBuild->globalFullMat, &tempVec, false);
    VecDestroy(&tempVec);

    KSPCreate(PETSC_COMM_WORLD, &stp1Solver);
    KSPCreate(PETSC_COMM_WORLD, &stp2Solver);
    KSPSetType(stp2Solver, KSPGMRES);
    KSPSetOperators(stp2Solver, globalBuild->globalFullMat, globalBuild->globalFullMat);
    KSPSetFromOptions(stp2Solver);
    //cout << "WHY\n";
    //MatView(globalBuild->globalViscMat, PETSC_VIEWER_STDOUT_WORLD);
    //cout << "4\n";
}

void Solver::applyDirichletConditions(Mat *m, Vec *v, bool expl){
    PetscInt* rows = msh->dirichletIds.data();

    for(int i = 0; i < msh->nNodes; i++){
        Node n = msh->nodes[msh->nodeIds[i]];
        if(n.boundary || n.inlet) {
            VecSetValue(*v, i, n.velocity[0], INSERT_VALUES);
            VecSetValue(*v, i + msh->nNodes, n.velocity[1], INSERT_VALUES);
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

void Solver::computeFirstStep(){
    Mat tempMat;
    Vec tempVec;
    Vec vint;

    MatCreate(PETSC_COMM_WORLD, &tempMat);
    MatSetSizes(tempMat, PETSC_DECIDE, PETSC_DECIDE, 
    msh->nNodes * 2, msh->nNodes * 2);
    MatSetFromOptions(tempMat);
    MatSetUp(tempMat);

    VecCreate(PETSC_COMM_WORLD, &tempVec);
    VecSetSizes(tempVec, PETSC_DECIDE, msh->nNodes * 2);
    VecSetFromOptions(tempVec);
    VecCreate(PETSC_COMM_WORLD, &vint);
    VecSetSizes(vint, PETSC_DECIDE, msh->nNodes * 2);
    VecSetFromOptions(vint);

    MatConvert(globalBuild->globalConvMat, MATSAME, MAT_INITIAL_MATRIX, &tempMat);
    
    MatDiagonalScale(tempMat, NULL, globalBuild->velocityVec);
    MatAssemblyBegin(tempMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(tempMat, MAT_FINAL_ASSEMBLY);

    //Add viscous and mass matrix to system
    MatAXPY(tempMat, 1.0, globalBuild->globalViscMat, DIFFERENT_NONZERO_PATTERN);
    MatAXPY(tempMat, 1.0, globalBuild->globalMassMat, DIFFERENT_NONZERO_PATTERN);

    MatMult(globalBuild->globalMassMat, globalBuild->velocityVec, tempVec);

    //Impose Dirichlet Conditions
    applyDirichletConditions(&tempMat, &tempVec, false);

    //Solve system
    KSPCreate(PETSC_COMM_WORLD, &stp1Solver);
    KSPSetType(stp1Solver, KSPGMRES);
    KSPSetOperators(stp1Solver, tempMat, tempMat);
    KSPSetFromOptions(stp1Solver);
    KSPSolve(stp1Solver, tempVec, vint);
    VecCopy(vint, globalBuild->velocityVec);

    //Cleanup
    VecDestroy(&tempVec);
    VecDestroy(&vint);
    MatDestroy(&tempMat);
}

void Solver::computeSecondStep(){
    Vec tempVec;
    Vec solVec;

    VecCreate(PETSC_COMM_WORLD, &tempVec);
    VecCreate(PETSC_COMM_WORLD, &solVec);
    VecSetSizes(tempVec, PETSC_DECIDE, msh->nNodes * 2);
    VecSetSizes(solVec, PETSC_DECIDE, msh->nNodes * 2 + msh->nLinear);

    VecSetFromOptions(tempVec);
    VecSetFromOptions(solVec);

    MatMult(globalBuild->globalMassMat, globalBuild->velocityVec, tempVec);

    VecZeroEntries(solVec);

    for(int i = 0; i < msh->nNodes * 2; i++){
        PetscInt ind = i;
        PetscScalar val;
        VecGetValues(tempVec, 1, &ind, &val);
        VecSetValue(solVec, i, val, INSERT_VALUES);
    }

    VecAssemblyBegin(solVec);
    VecAssemblyEnd(solVec);

    applyDirichletConditions(&globalBuild->globalFullMat, &solVec, true);
    VecZeroEntries(globalBuild->nodalVec);
    KSPSolve(stp2Solver, solVec, globalBuild->nodalVec);
    applyDirichletConditions(&globalBuild->globalFullMat, &globalBuild->nodalVec, true);
    globalBuild->updateVelocity();

    VecDestroy(&solVec);
    VecDestroy(&tempVec);
}

vector<vector<double>> Solver::computeTimeStep(int steps){
    vector<vector<double>> fluid = vector<vector<double>>(3);
    for(int x = 0; x < steps; x++){
        this->computeFirstStep();
        this->computeSecondStep();
        cout << "Step: " << x << "\n";
    }

    for(int i = 0; i < msh->nNodes; i++){
        PetscInt iu = i;
        PetscInt iv = i + msh->nNodes;
        PetscScalar u;
        PetscScalar v;

        VecGetValues(globalBuild->nodalVec, 1, &iu, &u);
        VecGetValues(globalBuild->nodalVec, 1, &iv, &v);

        fluid[0].push_back(u);
        fluid[1].push_back(v);
        fluid[2].push_back(-1);
    }

    for(int i = 0; i < msh->nNodes; i++){
        PetscInt ip = i + msh->nNodes * 2;
        PetscScalar p;
        VecGetValues(globalBuild->nodalVec, 1, &ip, &p);
        if(msh->nodes[msh->nodeIds[i]].pid != -1){
            fluid[2][i] = p;
        }
    }
    return fluid;
}