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

/**
 * Holds data for parabolic parabolic problem c(x)u_t - (k(x) u_x)_x = F(x,u)
 * @param file_name  output data for plotting
 * @param init function to initalize U at t=0
 * @param Fux function f(x,u) where it is caleld f(x, U[x])
 * @param nx number of mesh points
 * @param x_0 left spatial end point
 * @param x_n right spatial end point
 * @param timeControl simTime object
 * @param bc boundry condition
 */
class ParabolicPdeProblem {

 public:
  NumVec U, dU;
  std::string file_name;
  std::function<double(double)> k, c;
  std::function<double(double, double)> Fux;
  int nx;
  double x_0, x_n, dx;
  simTime timeControl;
  BoundryCondition bc;
  bool debug = true;

  //constructor
  ParabolicPdeProblem(std::string file_name,
                      const std::function<double(double)> &init,
                      const std::function<double(double)> &k,
                      const std::function<double(double)> &c,
                      const std::function<double(double, double)> &Fux,
                      int nx,
                      double x_0,
                      double x_n,
                      const simTime &timeControl,
                      BoundryCondition bc);
  /**
   * Runs simulation based on simTime parameters and outputs file
   */
  void run();

  /**
 * @return dU for dt step out of time t
 */
  NumVec step(double t, double dt);

  /**
 * Advances from current step to next step, using time control from simTime
 * Inplace updates dU and U
 */
  void advance();
};

#endif //CS28515PROJ1_LINEARPARABOLICPROBLEM_H
