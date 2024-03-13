#include "solver.hpp"

void Solver::updateVectors(Vec *vec1, Vec *vec2, bool vel){
    VecScatter *vecScatter = vel ? &vecScatter1 : &vecScatter2;
    VecScatterBegin(*vecScatter, *vec1, *vec2, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(*vecScatter, *vec1, *vec2, INSERT_VALUES, SCATTER_FORWARD);
    VecAssemblyBegin(*vec2);
    VecAssemblyEnd(*vec2);
}

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
    KSPSetInitialGuessNonzero(stp2Solver, PETSC_TRUE);
    KSPGetPC(stp2Solver, &preConditioner);
    PCSetType(preConditioner, PCNONE);
    KSPSetFromOptions(stp2Solver);

    PetscInt vecIndices[msh->nNodes * 2];
    for(int i = 0; i < msh->nNodes * 2; i++) vecIndices[i] = i;

    ISCreateGeneral(PETSC_COMM_WORLD, msh->nNodes * 2, vecIndices, PETSC_COPY_VALUES, &vecMapping);
    VecScatterCreate(globalBuild->nodalVec, vecMapping, globalBuild->velocityVec, vecMapping, &vecScatter1);
    VecScatterCreate(globalBuild->velocityVec, vecMapping, globalBuild->nodalVec, vecMapping, &vecScatter2);
}

Solver::~Solver(){
    KSPDestroy(&stp2Solver);
    VecScatterDestroy(&vecScatter1);
    VecScatterDestroy(&vecScatter2);
    ISDestroy(&vecMapping);
    delete(globalBuild);
}

void Solver::applyDirichletConditions(Mat *m, Vec *v, bool expl){
    PetscInt* rows = msh->dirichletIds.data();
    int high, low;
    VecGetOwnershipRange(*v, &low, &high);
    for(int i = 0; i < msh->nNodes; i++){
        Node n = msh->nodes[msh->nodeIds[i]];
        if(n.boundary || n.inlet) {
            if(i >= low && i < high){
                VecSetValue(*v, i, n.velocity[0], INSERT_VALUES);
            }
            if(i + msh->nNodes >= low && i + msh->nNodes < high) {
                VecSetValue(*v, i + msh->nNodes, n.velocity[1], INSERT_VALUES);
            }
        }
    }

    VecAssemblyBegin(*v);
    VecAssemblyEnd(*v);

    if(!expl){
        MatZeroRows(*m, msh->dirichletIds.size(), rows, 1.0, *v, *v);
    }
    MatAssemblyBegin(*m, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*m, MAT_FINAL_ASSEMBLY);
}

void Solver::computeFirstStep(){
    Mat tempMat;
    Vec tempVec;
    PC preConditioner;

    VecCreate(PETSC_COMM_WORLD, &tempVec);
    VecSetSizes(tempVec, PETSC_DECIDE, msh->nNodes * 2);
    VecSetFromOptions(tempVec);

    globalBuild->assembleConvectionMatrix();
    MatConvert(globalBuild->globalMassMat, MATSAME, MAT_INITIAL_MATRIX, &tempMat);
    //Subtract half viscous + convection terms from right hand side
    MatAXPY(tempMat, -1/2, globalBuild->globalViscMat, DIFFERENT_NONZERO_PATTERN);
    MatAXPY(tempMat, -1/2, globalBuild->globalConvMat, DIFFERENT_NONZERO_PATTERN);

    MatMult(tempMat, globalBuild->velocityVec, tempVec);
    //Reapply viscous and convection terms to system of equations
    MatAXPY(tempMat, 1.0, globalBuild->globalViscMat, DIFFERENT_NONZERO_PATTERN);
    MatAXPY(tempMat, 1.0, globalBuild->globalConvMat, DIFFERENT_NONZERO_PATTERN);

    //Impose Dirichlet Conditions
    applyDirichletConditions(&tempMat, &tempVec, false);
    //Solve system
    KSPCreate(PETSC_COMM_WORLD, &stp1Solver);
    KSPSetType(stp1Solver, KSPGMRES);
    KSPSetOperators(stp1Solver, tempMat, tempMat);
    KSPSetInitialGuessNonzero(stp1Solver, PETSC_TRUE);

    KSPGetPC(stp1Solver, &preConditioner);
    PCSetType(preConditioner, PCNONE);
    KSPSetFromOptions(stp1Solver);

    KSPSolve(stp1Solver, tempVec, globalBuild->velocityVec);
    applyDirichletConditions(&tempMat, &globalBuild->velocityVec, true);

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

    //MatDuplicate(globalBuild->globalMassMat, MAT_COPY_VALUES, &tempMat);
    //MatAXPY(tempMat, -1.0, globalBuild->globalViscMat, DIFFERENT_NONZERO_PATTERN);
    //MatMult(tempMat, globalBuild->velocityVec, tempVec);
    //MatDestroy(&tempMat);

    MatMult(globalBuild->globalMassMat, globalBuild->velocityVec, tempVec);

    updateVectors(&tempVec, &solVec, false);

    applyDirichletConditions(&globalBuild->globalFullMat, &solVec, true);
    KSPSolve(stp2Solver, solVec, globalBuild->nodalVec);
    //MatView(globalBuild->globalConvMat, PETSC_VIEWER_STDOUT_WORLD);

    applyDirichletConditions(&globalBuild->globalFullMat, &globalBuild->nodalVec, true);
    updateVectors(&globalBuild->nodalVec, &globalBuild->velocityVec, true);
    VecDestroy(&solVec);
    VecDestroy(&tempVec);
    
    PetscReal max;
    VecMax(globalBuild->velocityVec, NULL, &max);
    cout << "Max (Instability Metric): " << max << "\n";
}

void Solver::computeTimeStep(int steps){
    for(int x = 0; x < steps; x++){
        this->computeFirstStep();
        this->computeSecondStep();
        cout << "Step: " << x << "\n";
    }
}

void Solver::interpolateValues(vector<double> &coord, 
vector<vector<double>> &solData, vector<size_t> &nodeTags, Vec *solVec){
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
        VecGetValues(*solVec, 1, &vxi, &nodeVx);
        VecGetValues(*solVec, 1, &vyi, &nodeVy);
        VecGetValues(*solVec, 1, &pi, &nodePre);
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

vector<vector<double>> Solver::interpolateSolution(double resolution, int rank){
    Vec solVec;
    VecScatter vs;
    vector<vector<double>> solData = vector<vector<double>>(5);

    VecScatterCreateToZero(globalBuild->nodalVec, &vs, &solVec);
    VecScatterBegin(vs, globalBuild->nodalVec, solVec, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(vs, globalBuild->nodalVec, solVec, INSERT_VALUES, SCATTER_FORWARD);

    if(rank != 0){
        VecDestroy(&solVec);
        return solData;
    }

    size_t elementTag;
    int elementType;
    vector<size_t> nodeTags;
    double u, v, w;
    double xmax, ymax, zmax;
    double xmin, ymin, zmin;
    vector<double> coord;

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

            interpolateValues(coord, solData, nodeTags, &solVec);
        }
    }
    VecDestroy(&solVec);
    return solData;
}