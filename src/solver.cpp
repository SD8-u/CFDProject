#include "solver.hpp"

Solver::Solver(Mesh* msh, double dt, double viscosity){
    this->msh = msh;
    nNodes = msh->nodes.size();
    this->dt = 1/dt;
    this->viscosity = viscosity;

    MatCreate(PETSC_COMM_WORLD, &globalMassMat);
    MatCreate(PETSC_COMM_WORLD, &globalMassMatF);
    MatCreate(PETSC_COMM_WORLD, &globalMassMatI);
    MatCreate(PETSC_COMM_WORLD, &globalGradMat);
    //MatCreate(PETSC_COMM_WORLD, &globalGradMat);
    MatCreate(PETSC_COMM_WORLD, &globalViscMat);
    MatCreate(PETSC_COMM_WORLD, &globalConvMat);
    MatCreate(PETSC_COMM_WORLD, &globalFullMat);

    MatSetSizes(globalMassMat, PETSC_DECIDE, PETSC_DECIDE, 
    nNodes * 2, nNodes * 2);
    MatSetSizes(globalMassMatF, PETSC_DECIDE, PETSC_DECIDE, 
    nNodes * 2, nNodes * 2);
    MatSetSizes(globalMassMatI, PETSC_DECIDE, PETSC_DECIDE, 
    nNodes * 2, nNodes * 2);
    MatSetSizes(globalGradMat, PETSC_DECIDE, PETSC_DECIDE, 
    nNodes * 2, msh->nLinear);
    //MatSetSizes(globalGradMatT, PETSC_DECIDE, PETSC_DECIDE, 
    //nNodes * 2, msh->nLinear);
    MatSetSizes(globalConvMat, PETSC_DECIDE, PETSC_DECIDE, 
    nNodes * 2, nNodes * 2);
    MatSetSizes(globalViscMat, PETSC_DECIDE, PETSC_DECIDE, 
    nNodes * 2, nNodes * 2);
    MatSetSizes(globalFullMat, PETSC_DECIDE, PETSC_DECIDE, 
    msh->nLinear, msh->nLinear);

    MatSetFromOptions(globalMassMat);
    MatSetFromOptions(globalMassMatF);
    MatSetFromOptions(globalMassMatI);
    MatSetFromOptions(globalGradMat);
    //MatSetFromOptions(globalGradMat);
    MatSetFromOptions(globalViscMat);
    MatSetFromOptions(globalConvMat);
    MatSetFromOptions(globalFullMat);

    MatSetUp(globalMassMat);
    MatSetUp(globalMassMatF);
    MatSetUp(globalMassMatI);
    MatSetUp(globalGradMat);
    //MatSetUp(globalGradMat);
    MatSetUp(globalViscMat);
    MatSetUp(globalConvMat);
    MatSetUp(globalFullMat);

    MatZeroEntries(globalMassMat);
    MatZeroEntries(globalMassMatF);
    MatZeroEntries(globalMassMatI);
    MatZeroEntries(globalGradMat);
    //MatZeroEntries(globalGradMat);
    MatZeroEntries(globalViscMat);
    MatZeroEntries(globalConvMat);
    MatZeroEntries(globalFullMat);

    MatAssemblyBegin(globalFullMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalFullMat, MAT_FINAL_ASSEMBLY);

    VecCreate(PETSC_COMM_WORLD, &velocityVec);
    VecSetSizes(velocityVec, PETSC_DECIDE, nNodes * 2);
    VecSetFromOptions(velocityVec);

    VecCreate(PETSC_COMM_WORLD, &pressureVec);
    VecSetSizes(pressureVec, PETSC_DECIDE, msh->nLinear);
    VecSetFromOptions(pressureVec);

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

//Add SUPG to convection term
void Solver::applyStabilisation(Mat* convMat){
    for(size_t elementTag : msh->elementTags[0]){
        //Compute element level SUPG parameter
        double x[3], y[3];
        double avgU, avgV;
        int n = 0;
        for(size_t node : msh->elements[elementTag]){
            if(n < 3){
                x[n] = msh->nodes[node].x;
                y[n++] = msh->nodes[node].y;
            }
            PetscInt iU = msh->nodes[node].id;
            PetscInt iV = msh->nodes[node].id + nNodes;
            PetscScalar u, v;
            VecGetValues(velocityVec, 1, &iU, &u);
            VecGetValues(velocityVec, 1, &iV, &v);
            avgU += u;
            avgV += v;
        }
        avgU /= 6; avgV /= 6;
        double elemSize = sqrt(0.5 * (x[0] * (y[1] - y[2]) + x[1] * 
        (y[2] - y[0]) + x[2] * (y[0] - y[1])));
        double convCentrVel = sqrt(avgU * avgU + avgV * avgV);

        double supgParam = (elemSize/(2 * convCentrVel)) * 
        1/(sqrt(1 + ((6 * viscosity)/(elemSize * convCentrVel))));

        //Scale elements by param in global convection matrix to form SUPG stable term
        Vec supgScaler;
        VecCreate(PETSC_COMM_WORLD, &supgScaler);
        VecSetSizes(supgScaler, PETSC_DECIDE, nNodes * 2);
        VecSetFromOptions(supgScaler);
        VecZeroEntries(supgScaler);

        for(size_t node : msh->elements[elementTag]){
            VecSetValue(supgScaler, msh->nodes[node].id, supgParam, INSERT_VALUES);
            VecSetValue(supgScaler, msh->nodes[node].id + nNodes, supgParam, INSERT_VALUES);
        }

        VecAssemblyBegin(supgScaler);
        VecAssemblyEnd(supgScaler);
        MatDiagonalScale(*convMat, supgScaler, supgScaler);
        VecDestroy(&supgScaler);
    }
}

void Solver::localToGlobalMat(int type, bool inverse=false, bool full=true){
    for(size_t elementTag : msh->elementTags[0]){

        Mat localMat;
        Mat* globalMat;

        int row = 12;
        int col = 12;
        switch(type){
            case 1:
                localMat = computeMassMatrix(elementTag, inverse, full);
                if(full){globalMat = &globalMassMat;}
                else if(inverse) {globalMat = &globalMassMatI;}
                else {globalMat = &globalMassMatF;}
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
                row = 15; col = 15;
                break;
            case 5:
                localMat = computeGradientMatrix(elementTag);
                globalMat = &globalGradMat;
                col = 3;
                break;    
        }
        for(int i = 0; i < row; i++){
            int x, y;
            x = msh->nodes[msh->elements[elementTag][i % 6]].id;
            if(i > 5){
                x += nNodes;
            }
            if(i > 11){
                x = msh->nodes[msh->elements[elementTag][i % 6]].pid;
                x += nNodes * 2;
            }
            for(int j = 0; j < col; j++){
                y = msh->nodes[msh->elements[elementTag][j % 6]].id;
                if(j > 5){
                    y += nNodes;
                }
                if(j > 11 || type == 5){
                    y = msh->nodes[msh->elements[elementTag][j % 6]].pid;
                }
                if(j > 11 && type != 5){
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
    localToGlobalMat(1, false, true);
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

    cout << "A\n";
    localToGlobalMat(1, false, false);
    MatAssemblyBegin(globalMassMatF, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalMassMatF, MAT_FINAL_ASSEMBLY);
    cout << "B\n";
    localToGlobalMat(5);
    cout << "B1\n";
    MatAssemblyBegin(globalGradMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalGradMat, MAT_FINAL_ASSEMBLY);
    MatTranspose(globalGradMat, MAT_INITIAL_MATRIX, &globalGradMatT);
    cout << "C\n";
    localToGlobalMat(1, true, false);
    MatAssemblyBegin(globalMassMatI, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalMassMatI, MAT_FINAL_ASSEMBLY);
    cout << "D\n";
    MatMatMatMult(globalGradMatT, globalMassMatI, globalGradMat, 
    MAT_INITIAL_MATRIX, PETSC_DEFAULT, &globalFullMat);

    MatAssemblyBegin(globalFullMat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalFullMat, MAT_FINAL_ASSEMBLY);
    cout << "E\n";
    localToGlobalVec(false);
    VecAssemblyBegin(velocityVec);
    VecAssemblyEnd(velocityVec);
    cout << "F\n";
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
    //MatDestroy(&stabMat);

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
    KSPSolve(solver, tempVec, velocityVec);
    //VecCopy(vint, velocityVec);

    //VecView(vint, PETSC_VIEWER_STDOUT_WORLD);
    //Cleanup
    VecDestroy(&tempVec);
    MatDestroy(&tempMat);
    VecDestroy(&vint);
    KSPDestroy(&solver);
    //VecView(velocityVec, PETSC_VIEWER_STDOUT_WORLD);
}

void Solver::computeSecondStep(){
    cout << "0\n";
    Mat tempMat1;
    Mat tempMat2;
    Vec tempVec;
    Vec solVec;
    MatCreate(PETSC_COMM_WORLD, &tempMat1);
    MatCreate(PETSC_COMM_WORLD, &tempMat2);
    MatSetSizes(tempMat1, PETSC_DECIDE, PETSC_DECIDE, 
    msh->nLinear, nNodes * 2);
    MatSetSizes(tempMat2, PETSC_DECIDE, PETSC_DECIDE, 
    nNodes * 2, msh->nLinear);
    MatSetFromOptions(tempMat1);
    MatSetFromOptions(tempMat2);
    MatSetUp(tempMat1);
    MatSetUp(tempMat2);

    cout << "1\n";
    //Solve for pressure 
    VecCreate(PETSC_COMM_WORLD, &tempVec);
    VecCreate(PETSC_COMM_WORLD, &solVec);
    VecSetSizes(tempVec, PETSC_DECIDE, nNodes * 2);
    VecSetSizes(solVec, PETSC_DECIDE, msh->nLinear);
    VecSetFromOptions(tempVec);
    VecSetFromOptions(solVec);
    VecZeroEntries(solVec);

    MatConvert(globalGradMatT, MATSAME, MAT_INITIAL_MATRIX, &tempMat1);
    MatScale(tempMat1, dt);
    MatMult(tempMat1, velocityVec, solVec);

    cout << "2\n";
    KSP solver;
    KSPCreate(PETSC_COMM_WORLD, &solver);
    KSPSetOperators(solver, globalFullMat, globalFullMat);
    KSPSetType(solver, KSPGMRES);
    KSPSetFromOptions(solver);

    KSPSolve(solver, solVec, pressureVec);

    cout << "3\n";
    //Solve for final velocity
    Vec solVec1;
    VecCreate(PETSC_COMM_WORLD, &solVec1);
    VecSetSizes(solVec1, PETSC_DECIDE, nNodes * 2);
    VecSetFromOptions(solVec1);
    MatMult(globalMassMatF, velocityVec, solVec1);
    MatConvert(globalGradMat, MATSAME, MAT_INITIAL_MATRIX, &tempMat2);
    MatScale(tempMat2, 1/dt);
    MatMult(tempMat2, pressureVec, tempVec);
    VecAXPY(solVec1, -1.0, tempVec);

    applyDirichletConditions(&globalMassMatF, &solVec1, false);

    cout << "4\n";
    KSPCreate(PETSC_COMM_WORLD, &solver);
    KSPSetOperators(solver, globalMassMatF, globalMassMatF);
    KSPSetType(solver, KSPGMRES);
    KSPSetFromOptions(solver);

    KSPSolve(solver, solVec1, velocityVec);

    cout << "5\n";
    PetscScalar max = 0;
    for(int i = 0; i < nNodes * 2; i++){
        PetscInt ind = i;
        PetscScalar val;
        VecGetValues(velocityVec, 1, &ind, &val);
        max = val > max ? val : max;
    }
    cout << "MAX (instability metric): " << max << "\n";

    //VecView(velocityVec, PETSC_VIEWER_STDOUT_WORLD);

    VecDestroy(&solVec);
    //VecDestroy(&solVec1);
    VecDestroy(&tempVec);
    MatDestroy(&tempMat1);
    MatDestroy(&tempMat2);
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

        VecGetValues(velocityVec, 1, &iu, &u);
        VecGetValues(velocityVec, 1, &iv, &v);

        fluid[0].push_back(u);
        fluid[1].push_back(v);
        fluid[2].push_back(-1);
    }

    for(int i = 0; i < msh->nLinear; i++){
        PetscInt ip = i;
        PetscScalar p;
        VecGetValues(pressureVec, 1, &ip, &p);
        if(msh->nodes[msh->nodeIds[i]].pid != -1){
            fluid[2][msh->nodes[msh->nodeIds[i]].id] = p;
        }
    }
    //MatView(globalMassMat, PETSC_VIEWER_STDOUT_WORLD);
    //MatView(globalMassMatF, PETSC_VIEWER_STDOUT_WORLD);
    return fluid;
}