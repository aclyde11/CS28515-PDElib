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

int main(int argc, char *argv[]) {

  std::string file = "test1.txt";
  if (cmdOptionExists(argv, argc + argv, "-f"))
    file = getCmdOption(argv, argc + argv, "-f");

  double init = 0.8;
  if (cmdOptionExists(argv, argv + argc, "-init_value"))
    init = std::stod(getCmdOption(argv, argv + argc, "-init_value"));

  int mesh_points = 15;
  if (cmdOptionExists(argv, argv + argc, "-n"))
    mesh_points = std::stoi(getCmdOption(argv, argv + argc, "-n"));

  double x_0 = 0;
  if (cmdOptionExists(argv, argv + argc, "-x_0"))
    x_0 = std::stod(getCmdOption(argv, argv + argc, "-x_0"));

  double x_nx = 3;
  if (cmdOptionExists(argv, argv + argc, "-x_n"))
    x_nx = std::stod(getCmdOption(argv, argv + argc, "-x_n"));

  double tmax = 5.0;
  if (cmdOptionExists(argv, argv + argc, "-tmax"))
    tmax = std::stod(getCmdOption(argv, argv + argc, "-tmax"));

  std::function<double(double)> one = [](double x) { return 1; };
  std::function<double(double, double)> F = [](double x, double U) {
    double g = (0 <= x && x < 1) ? 1 : -4;
    return U * (1 - U) * 3 * g;
  };
  std::function<double(double)> initf = [init](double x) -> double {
    return init;
  };

  simTime tc;
  tc.endTime = tmax;
  std::cout << tc;
  ParabolicPdeProblem pb(file, initf, one, one, F, mesh_points, x_0, x_nx, tc, vonNeumann);

  std::clock_t start;
  start = std::clock();
  pb.run();
  std::cout << "\nTime: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC / 1000.0) << " ms" << std::endl;
  return 0;
}