#include "compute.hpp"
#include "localbuilder.hpp"

class GlobalBuilder {
    private:
        double dt;
        double viscosity;
        Mesh *msh;
        LocalBuilder *localBuild;
        void localToGlobalVec(bool full);
        void globalToLocalVec(size_t elementTag, Vec *localVec);
        void localToGlobalMat(size_t elementTag, Mat *localMat, Mat *globalMat, bool final);
    public:
        Vec nodalVec;
        Vec velocityVec;
        Mat globalMassMat;
        Mat globalViscMat;
        Mat globalConvMat;
        Mat globalFullMat;
        GlobalBuilder(int dim, double dt, double visc, Mesh *msh);
        ~GlobalBuilder();
        void assembleVectors();
        void assembleMatrices();
        void assembleConvectionMatrix();
        void updateVelocity();
};