#include <iostream>
#include "tridiag.h"
#include "parabolicsolver.h"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <vector>
#include <functional>
#include <cmath>
#include <ctime>
#include <string>
#include "LinearParabolicProblem.h"
#include "simtime.h"
#define PI 3.141592653

int main() {
  std::function<double(double)> f = [](double x) { return 0; };
  std::function<double(double)> one = [](double x) { return 1; };
  std::function<double(double)> zero = [](double x) { return 0; };
  std::function<double(double)> negOne = [](double x) { return -1; };
  std::function<double(double)> id = [](double x) { return x; };

  // c(x)u_t - k(x)u_xx + f(x,u)=0
  /*
  TriDiag M = generateMassMatrixMidpoint(zero, 5, 0.25, 0); //c(x)
  std::cout << "mass: " << std::endl << M << std::endl;

  TriDiag S = generateStiffnessMatrixMidpoint(one, 5, 0.25, 0); //k(x)
  std::cout << "stiff: " << std::endl << S << std::endl;
  */

  double x_0 = 0;
  double x_nx = 1;
  double L = (x_nx - x_0);
  int nx = 15;
  int nt = 1000;
  double tmax = 1.0;

  std::function<double(double)> init;

  std::function<double(double, NumVec)> zeroF = [](double x, NumVec U) {
    return 0;
  };

  std::function<double(double, NumVec)> F = [](double x, NumVec U) {
    double g = (0 <= x && x < 1) ? 1 : -3;
    return (1 - x) * 3 * g;
  };

  /*init = [L](double x) -> double {
    return sin(PI * x / L);
  };*/

  init = [](double x) -> double {
    return 0.1;
  };

  //k, c

  std::clock_t start;
  NumVec U(5);
  start = std::clock();
  // solveMassStiffStepDouble(one, one, x_0, x_nx, nx, nt, tmax, init, false, true);
  std::cout << "Time: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC / 1000.0) << " ms" << std::endl;

  std::cout << linearizeF(U, F, 0.25, 0.05, 0.001);
  std::cout << "\n\n\n\n" << "starting new type\n";

  simTime tc;
  LinearParabolicProblem pb("test2.txt", init, one, one, F, nx, 0, 1, tc, vonNeumann);

  start = std::clock();
  pb.run();
  std::cout << "\nTime: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC / 1000.0) << " ms" << std::endl;
  return 0;
}