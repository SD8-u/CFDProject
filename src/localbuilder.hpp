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
        void computeGradientMatrix();
    public:
        Mat localMassMat;
        Mat localViscMat;
        Mat localConvMat;
        Mat localFullMat;
        LocalBuilder(size_t elementTag);
        void computeMassMatrix();
        void computeViscosityMatrix();
        void computeConvectionMatrix();
        void computeFinalMatrix();
};