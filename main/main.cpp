#include <iostream>
#include "tridiag.h"
#include "pdesolver.h"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <vector>
#include <functional>
#include <cmath>
#include <ctime>
#include <string>

#define PI 3.141592653


int main() {
  double x_0 = 0;
  double x_nx = 1;
  double L = (x_nx - x_0);
  int nx = 15;
  int nt = 10000;
  double tmax = 1;
  double alpha = 0.1;
  std::function<double(double)> init;
  init = [L](double x) -> double { //TODO: use alpha or not?
    return sin(PI * x / L);
  };

  std::clock_t start;
  start = std::clock();

  solveHeatEquation1d(x_0, x_nx, nx, nt, tmax, alpha, init);
  std::cout << "Time: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
  return 0;
}