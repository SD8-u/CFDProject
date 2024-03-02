#include "compute.hpp"

class LocalBuilder {
    private:
        const vector<double> gaussPoints = {
            0.44594849092, 0.44594849092, 0, 0.10810301817, 0.44594849092, 0, 
            0.44594849092, 0.10810301817, 0, 0.09157621351, 0.09157621351, 0, 
            0.81684757289, 0.09157621351, 0, 0.09157621351, 0.81684757289, 0
        };
        const vector<double> gaussWeights = {
            0.111690794839, 0.111690794839, 0.111690794839, 
            0.054975871827, 0.054975871827, 0.054975871827
        };

        bool conv;
        double dt;
        double viscosity;
        size_t elementTag;
        Mat localGradMat;
        vector<Mat> basisMats;
        vector<Mat> basisGradMats;
        vector<Mat> inverseJacobian;
        vector<double> basisFuncs;
        vector<double> basisFuncsPres;
        vector<double> basisFuncsGrad;
        vector<double> j;
        vector<double> jdets;

        void setUp();
        void computeGradientMatrix();
        void computeMassMatrix();
        void computeViscosityMatrix();
        void computeFinalMatrix();
        void buildBasisMatrix();
        void buildBasisGradMatrix();
        void buildInverseJacobian();
    public:
        Mat localMassMat;
        Mat localViscMat;
        Mat localConvMat;
        Mat localFullMat;

        LocalBuilder();
        LocalBuilder(double dt, double viscosity);
        ~LocalBuilder();
        void assembleMatrices(size_t elementTag);
        void computeConvectionMatrix(size_t elementTag, Vec *velocityVec);
};