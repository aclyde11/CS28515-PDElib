//
// Created by Austin Clyde on 4/24/18.
//

#ifndef CS28515PROJ1_LINEARPARABOLICPROBLEM_H
#define CS28515PROJ1_LINEARPARABOLICPROBLEM_H

#include <string>
#include "numvec.h"
#include "tridiag.h"
#include "simtime.h"
#include "parabolicsolver.h"
#include "utility.h"
enum BoundryCondition { vonNeumann, dirchlet };

class LinearParabolicProblem {

 public:
  NumVec U, dU;
  std::string file_name;
  std::function<double(double)> k, c;
  std::function<double(double, NumVec)> Fux;
  int nx;
  double x_0, x_n, dx;
  simTime timeControl;
  BoundryCondition bc;
  bool debug = true;

  //constructor
  LinearParabolicProblem(std::string file_name,
                         const std::function<double(double)> &init,
                         const std::function<double(double)> &k,
                         const std::function<double(double)> &c,
                         const std::function<double(double, NumVec)> &Fux,
                         int nx,
                         double x_0,
                         double x_n,
                         const simTime &timeControl,
                         BoundryCondition bc);

  void run();
  NumVec step(double t, double dt);

  void advance();
};

#endif //CS28515PROJ1_LINEARPARABOLICPROBLEM_H
