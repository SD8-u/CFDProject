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
    public:
        Solver(Mesh* msh, double dt, double visc);
        void computeFirstStep();
        void computeSecondStep();
        void computeTimeStep(int steps);
        vector<vector<double>> interpolateSolution(double resolution);
};