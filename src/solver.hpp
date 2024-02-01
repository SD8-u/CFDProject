#include <petsc.h>

#ifndef MESH_HPP
#define MESH_HPP
#include <mesh.hpp>
#endif

class Solver {
    private:
        int nNodes;
        double dt = 1/0.01;
        Mat globalMassMat;
        Mat globalViscMat;
        Mat globalConvMat;
        Mat globalFullMat;
        Vec velocityVec;
        Vec nodalVec;
        Mesh *msh;
        void localToGlobalMat(int type);
        void localToGlobalVec(bool full);
    public:
        Solver(Mesh* msh);
        void assembleMatrices();
        void assembleVector();
        void computeFirstStep();
        void computeSecondStep();
        void computeTimeStep();
};