#include <petsc.h>
#include "compute.hpp"
#include "globalbuilder.hpp"

class Solver {
    private:
        double dt;
        double viscosity;
        KSP stp1Solver;
        KSP stp2Solver;
        Mesh *msh;
        GlobalBuilder *globalBuild;
        void applyDirichletConditions(Mat *m, Vec *v, bool full);
        void interpolateValues(vector<double> &coord, vector<vector<double>> &solData, vector<size_t> &nodeTags);
    public:
        Solver(Mesh* msh, double dt, double visc);
        ~Solver();
        void computeFirstStep();
        void computeSecondStep();
        void computeTimeStep(int steps);
        vector<vector<double>> interpolateSolution(double resolution);
};