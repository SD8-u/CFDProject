#include <petsc.h>
#include "compute.hpp"

class Solver {
    private:
        int nNodes;
        double dt = 1/0.001;
        Mat globalMassMat;
        Mat globalViscMat;
        Mat globalConvMat;
        Mat globalFullMat;
        Vec velocityVec;
        Vec nodalVec;
        Mesh *msh;
        void localToGlobalMat(int type);
        void localToGlobalVec(bool full);
        void applyDirichletConditions(Mat *m, Vec *v, bool full);
    public:
        Solver(Mesh* msh);
        void assembleMatrices();
        void assembleVector();
        void computeFirstStep();
        void computeSecondStep();
        vector<vector<double>> computeTimeStep(int steps);
};