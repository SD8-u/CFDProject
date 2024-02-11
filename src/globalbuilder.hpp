#include "compute.hpp"
#include "localbuilder.hpp"

class GlobalBuilder {
    private:
        Mesh *msh;
        LocalBuilder localBuild;
        void localToGlobalVec();
        void localToGlobalMat();
    public:
        Mat globalMassMat;
        Mat globalViscMat;
        Mat globalConvMat;
        Mat globalFullMat;
        GlobalBuilder(int dim, Mesh *msh);
        void assembleVectors();
        void assembleMatrices();
};