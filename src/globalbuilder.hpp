#ifndef SRC_GLOBALBUILDER_HPP_
#define SRC_GLOBALBUILDER_HPP_

#include <vector>

#include "localbuilder.hpp"
#include "mesh.hpp"

class GlobalBuilder {
 private:
  int domainSize;
  int domainStart;
  double dt;
  double viscosity;

  Mesh *msh;
  LocalBuilder *localBuild;

  void getDomain();
  void constructVectorMapping();
  void initialiseMats();
  void initialiseVecs();

  void localToGlobalVec(bool full);
  void globalToLocalVec(size_t elementTag, Vec *localVec);
  void localToGlobalMat(size_t elementTag, Mat *localMat, Mat *globalMat,
                        bool final);

 public:
  Vec fullVec;
  Vec currVelVec;
  Mat globalMassMat;
  Mat globalViscMat;
  Mat globalConvMat;
  Mat globalFullMat;

  GlobalBuilder(int dim, double dt, double visc, Mesh *msh);
  ~GlobalBuilder();

  void assembleVectors();
  void assembleMatrices();
  void assembleConvectionMatrix();
};

#endif  // SRC_GLOBALBUILDER_HPP_
