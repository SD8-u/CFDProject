#include <petsc.h>
#include "compute.hpp"
#include "globalbuilder.hpp"

class Solver {
    private:
        double dt, viscosity;
        KSP stp1Solver, stp2Solver;
        IS vecMapping;
        VecScatter vecScatter1, vecScatter2;
        Mesh *msh;
        GlobalBuilder *globalBuild;
        void applyDirichletConditions(Mat *m, Vec *v, bool full);
        void interpolateValues(vector<double> &coord, vector<vector<double>> &solData, 
        vector<size_t> &nodeTags, Vec *solVec);
        void updateVectors(Vec *vec1, Vec *vec2, bool vel);
        
    public:
        Solver(Mesh* msh, double dt, double visc);
        ~Solver();
        void computeFirstStep();
        void computeSecondStep();
        void computeTimeStep(int steps);
        vector<vector<double>> interpolateSolution(double resolution, int rank);
};