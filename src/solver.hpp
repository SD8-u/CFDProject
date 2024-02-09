#include <petsc.h>
#include "compute.hpp"

class Solver {
    private:
        int nNodes;
        double dt;
        double viscosity;
        Mat globalMassMat;
        Mat globalMassMatF;
        Mat globalGradMat;
        Mat globalGradMatT;
        Mat globalMassMatI;
        Mat globalViscMat;
        Mat globalConvMat;
        Mat globalFullMat;
        Vec velocityVec;
        Vec pressureVec;
        Vec nodalVec;
        Mesh *msh;
        void localToGlobalMat(int type, bool inverse, bool full);
        void localToGlobalVec(bool full);
        void applyDirichletConditions(Mat *m, Vec *v, bool full);
        void applyStabilisation(Mat* convMat);
    public:
        Solver(Mesh* msh, double dt, double visc);
        void assembleMatrices();
        void assembleVector();
        void computeFirstStep();
        void computeSecondStep();
        vector<vector<double>> computeTimeStep(int steps);
};