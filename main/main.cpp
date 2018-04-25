#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <vector>
#include <functional>
#include <cmath>
#include <ctime>
#include <string>

#include "tridiag.h"
#include "ParabolicPdeProblem.h"
#include "simtime.h"

int main() {
  //Set up for genetic problem, xu_t - ku_xx = F(x,u)
  double x_0 = 0;
  double x_nx = 3;
  int nx = 15; //number of mesh points

  std::function<double(double)> one = [](double x) { return 1; };

  std::function<double(double, double)> F = [](double x, double U) {
    double g = (0 <= x && x < 1) ? 1 : -4;
    return U * (1 - U) * 3 * g;
  };

  std::function<double(double)> initf = [](double x) -> double {
    return 0.8;
  };

  simTime tc;
  ParabolicPdeProblem pb("test2.txt", initf, one, one, F, nx, x_0, x_nx, tc, vonNeumann);

  std::clock_t start;
  start = std::clock();
  pb.run();
  std::cout << "\nTime: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC / 1000.0) << " ms" << std::endl;
  return 0;
}