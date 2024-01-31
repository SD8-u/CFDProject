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
        Mesh *msh;
        void localToGlobal(Mesh* msh, int type);
    public:
        Solver(Mesh* msh);
        void assembleMatrices();
        void computeTimeStep();
};