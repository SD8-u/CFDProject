#include <petsc.h>

#ifndef MESH_HPP
#define MESH_HPP
#include <mesh.hpp>
#endif

class Solver {
    private:
        int nNodes;
        Mat globalMassMat;
        Mat globalViscMat;
        Mat globalConvMat;
        Mat globalFullMat;
        Vec nodalVec;
        Mesh *msh;
        void localToGlobalMat(int type);
        void localToGlobalVec();
    public:
        Solver(Mesh* msh);
        void assembleMatrices();
        void assembleVector();
        void computeTimeStep();
};