//
// Created by Austin Clyde on 5/23/18.
//

//
// Created by Austin Clyde on 5/18/18.
//

#include <iostream>
#include <vector>
#include <cmath>

#include "test.h"
#include "utility.h"

#define PI 3.141592653589

int main(int argc, char *argv[]) {
    std::string file = "test1.txt";
    if (cmdOptionExists(argv, argc + argv, "-f"))
        file = getCmdOption(argv, argc + argv, "-f");

    int mesh_points = 100;
    if (cmdOptionExists(argv, argv + argc, "-n"))
        mesh_points = std::stoi(getCmdOption(argv, argv + argc, "-n"));

    double c_0 = 1;
    if (cmdOptionExists(argv, argv + argc, "-c0"))
        c_0 = std::stod(getCmdOption(argv, argv + argc, "-c0"));

    double k_0 = 1;
    if (cmdOptionExists(argv, argv + argc, "-k0"))
        c_0 = std::stod(getCmdOption(argv, argv + argc, "-k0"));

    double x_0 = 0;
    if (cmdOptionExists(argv, argv + argc, "-x_0"))
        x_0 = std::stod(getCmdOption(argv, argv + argc, "-x_0"));

    double x_nx = 2 * PI;
    if (cmdOptionExists(argv, argv + argc, "-x_n"))
        x_nx = std::stod(getCmdOption(argv, argv + argc, "-x_n"));

    int periods = 2;
    if (cmdOptionExists(argv, argv + argc, "-x_n"))
        periods = std::stoi(getCmdOption(argv, argv + argc, "-x_n"));

    double dt = 0.01;
    if (cmdOptionExists(argv, argv + argc, "-dt"))
        dt = std::stod(getCmdOption(argv, argv + argc, "-dt"));

    double tmax = 5;
    if (cmdOptionExists(argv, argv + argc, "-tmax"))
        tmax = std::stod(getCmdOption(argv, argv + argc, "-tmax"));

    std::function<double(double)> cx = [c_0](double x) {
        return (x >= 0 && x <= PI) ? 1 : c_0;
    };

    std::function<double(double)> kx = [k_0](double x) {
        return (x >= 0 && x <= PI) ? 1 : k_0;
    };

    std::function<double(double)> initf = [periods, dt](double x) -> double {
        return sin(periods * x) - cos(periods * x) * dt;
    };

    std::function<double(double)> initdU = [periods, dt](double x) -> double {
        return -1 * periods * cos(periods * x) * dt;
    };

    test_eigen();

    return 0;
}