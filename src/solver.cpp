#include "solver.hpp"

Solver::Solver(Mesh* msh, double dt, double viscosity){
    Vec tempVec;
    this->msh = msh;
    PC preConditioner;

    globalBuild = new GlobalBuilder(2, dt, viscosity, msh);
    globalBuild->assembleMatrices();
    globalBuild->assembleVectors();

    VecCreate(PETSC_COMM_WORLD, &tempVec);
    VecSetSizes(tempVec, PETSC_DECIDE, msh->nNodes * 2 + msh->nLinear);
    VecSetFromOptions(tempVec);

    applyDirichletConditions(&globalBuild->globalFullMat, &tempVec, false);
    VecDestroy(&tempVec);

    KSPCreate(PETSC_COMM_WORLD, &stp2Solver);
    KSPSetType(stp2Solver, KSPGMRES);
    KSPSetOperators(stp2Solver, globalBuild->globalFullMat, globalBuild->globalFullMat);
    //KSPSetInitialGuessNonzero(stp2Solver, PETSC_TRUE);
    //KSPGetPC(stp2Solver, &preConditioner);
    //PCSetType(preConditioner, PCICC);
    KSPSetFromOptions(stp2Solver);
}

Solver::~Solver(){
    KSPDestroy(&stp2Solver);
    delete(globalBuild);
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
        cout << "4\n";
        MatZeroRows(*m, msh->dirichletIds.size(), rows, 1.0, *v, *v);
    }
    MatAssemblyBegin(*m, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*m, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(*v);
    VecAssemblyEnd(*v);
}

void Solver::computeFirstStep(){
    cout << "1\n";
    Mat tempMat;
    Vec tempVec;
    PC preConditoner;

    cout << "2\n";
    VecCreate(PETSC_COMM_WORLD, &tempVec);
    VecSetSizes(tempVec, PETSC_DECIDE, msh->nNodes * 2);
    VecSetFromOptions(tempVec);

    globalBuild->assembleConvectionMatrix();
    MatConvert(globalBuild->globalMassMat, MATSAME, MAT_INITIAL_MATRIX, &tempMat);

    //MatAssemblyBegin(tempMat, MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(tempMat, MAT_FINAL_ASSEMBLY);
    cout << "3\n";
    //Subtract half viscous + convection terms from right hand side
    MatAXPY(tempMat, -1/2, globalBuild->globalViscMat, DIFFERENT_NONZERO_PATTERN);
    MatAXPY(tempMat, -1/2, globalBuild->globalConvMat, DIFFERENT_NONZERO_PATTERN);

    MatMult(tempMat, globalBuild->velocityVec, tempVec);

    //Reapply viscous and convection terms to system of equations
    MatAXPY(tempMat, 1.0, globalBuild->globalViscMat, DIFFERENT_NONZERO_PATTERN);
    MatAXPY(tempMat, 1.0, globalBuild->globalConvMat, DIFFERENT_NONZERO_PATTERN);

    //Impose Dirichlet Conditions
    applyDirichletConditions(&tempMat, &tempVec, false);

    //MatView(globalBuild->globalConvMat, PETSC_VIEWER_STDOUT_WORLD);
    //Solve system
    cout << "4\n";
    KSPCreate(PETSC_COMM_WORLD, &stp1Solver);
    KSPSetType(stp1Solver, KSPGMRES);
    KSPSetOperators(stp1Solver, tempMat, tempMat);
    cout << "5\n";
    //KSPGetPC(stp1Solver, &preConditoner);
    //PCSetType(preConditoner, PCICC);
    KSPSetFromOptions(stp1Solver);
    cout << "5a\n";
    KSPSolve(stp1Solver, tempVec, globalBuild->velocityVec);
    cout << "5b\n";
    //applyDirichletConditions(&tempMat, &globalBuild->velocityVec, true);
    cout << "5c\n";
    //Cleanup
    VecDestroy(&tempVec);
    MatDestroy(&tempMat);
    KSPDestroy(&stp1Solver);
}

void Solver::computeSecondStep(){
    Vec tempVec;
    Vec solVec;
    Mat tempMat;

    VecCreate(PETSC_COMM_WORLD, &tempVec);
    VecCreate(PETSC_COMM_WORLD, &solVec);
    VecSetSizes(tempVec, PETSC_DECIDE, msh->nNodes * 2);
    VecSetSizes(solVec, PETSC_DECIDE, msh->nNodes * 2 + msh->nLinear);

    VecSetFromOptions(tempVec);
    VecSetFromOptions(solVec);
    cout << "6\n";
    //MatDuplicate(globalBuild->globalMassMat, MAT_COPY_VALUES, &tempMat);
    //MatAXPY(tempMat, -1.0, globalBuild->globalViscMat, DIFFERENT_NONZERO_PATTERN);
    //MatMult(tempMat, globalBuild->velocityVec, tempVec);
    //MatDestroy(&tempMat);
    cout << "7\n";
    VecZeroEntries(solVec);
    VecAssemblyBegin(solVec);
    VecAssemblyEnd(solVec);
    MatMult(globalBuild->globalMassMat, globalBuild->velocityVec, solVec);
    //VecView(tempVec, PETSC_VIEWER_STDOUT_WORLD);
    //for(int i = 0; i < msh->nNodes * 2; i++){
    //    PetscInt ind = i;
    //    PetscScalar val;
    //    VecGetValues(tempVec, 1, &ind, &val);
    //    VecSetValue(solVec, i, val, INSERT_VALUES);
    //}
    //VecAssemblyBegin(solVec);
    //VecAssemblyEnd(solVec);
    cout << "8\n";
    applyDirichletConditions(&globalBuild->globalFullMat, &solVec, true);
    VecZeroEntries(globalBuild->nodalVec);
    KSPSolve(stp2Solver, solVec, globalBuild->nodalVec);

    cout << "9\n";
    //applyDirichletConditions(&globalBuild->globalFullMat, &globalBuild->nodalVec, true);
    globalBuild->updateVelocity();
    cout << "10\n";
    VecDestroy(&solVec);
    VecDestroy(&tempVec);
}

void Solver::computeTimeStep(int steps){
    for(int x = 0; x < steps; x++){
        this->computeFirstStep();
        this->computeSecondStep();
        cout << "Step: " << x << "\n";
    }
}

void Solver::interpolateValues(vector<double> &coord, 
vector<vector<double>> &solData, vector<size_t> &nodeTags){
    int nComp, nOrien;
    double pressure = 0, xv = 0, yv = 0;
    double nodePre, nodeVx, nodeVy;
    vector<double> basisFuncsVel, basisFuncsPre;

    gmsh::model::mesh::getBasisFunctions(2, coord, "Lagrange", nComp, 
    basisFuncsPre, nOrien);
    gmsh::model::mesh::getBasisFunctions(9, coord, "Lagrange", nComp, 
    basisFuncsVel, nOrien);
    
    for(int node = 0; node < nodeTags.size(); node++){
        PetscInt pi = msh->nodes[nodeTags[node]].pid + msh->nNodes * 2;
        PetscInt vxi = msh->nodes[nodeTags[node]].id;
        PetscInt vyi = vxi + msh->nNodes;
        VecGetValues(globalBuild->nodalVec, 1, &vxi, &nodeVx);
        VecGetValues(globalBuild->nodalVec, 1, &vyi, &nodeVy);
        VecGetValues(globalBuild->nodalVec, 1, &pi, &nodePre);
        if(node < 3){
            pressure += basisFuncsPre[node] * nodePre;
        }

        xv += basisFuncsVel[node] * nodeVx;
        yv += basisFuncsVel[node] * nodeVy;
    }
    
    solData[0].push_back(pressure); solData[1].push_back(xv);
    solData[2].push_back(yv);
    coord.clear();
}

vector<vector<double>> Solver::interpolateSolution(double resolution){
    size_t elementTag;
    int elementType;
    vector<size_t> nodeTags;
    double u, v, w;
    double xmax, ymax, zmax;
    double xmin, ymin, zmin;
    vector<double> coord;
    vector<vector<double>> solData = vector<vector<double>>(5);

    gmsh::model::getBoundingBox(-1, -1, xmin, ymin, zmin, xmax, ymax, zmax);
    for(double x = xmin; x < xmax; x+=resolution){
        for(double y = ymin; y < ymax; y+=resolution){
            solData[3].push_back(x); solData[4].push_back(y);

            try{
                gmsh::model::mesh::getElementByCoordinates(x, y, 0, elementTag, 
                elementType, nodeTags, u, v, w, -1, true);
            }
            catch(...) {
                solData[0].push_back(0); solData[1].push_back(0); solData[2].push_back(0);
                continue;
            }

            coord.push_back(u); coord.push_back(v); coord.push_back(w);

            interpolateValues(coord, solData, nodeTags);
        }
    }
    return solData;
}