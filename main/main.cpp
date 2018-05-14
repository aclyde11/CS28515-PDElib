#include <iostream>
#include <vector>
#include <cmath>

#include "tridiag.h"
#include "WaveEquationProblem.h"
#include "ParabolicPdeProblem.h"

#define PI 3.14592653589

void run_proj2(int argc, char *argv[]) {
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

    double tmax = 3.0;
    if (cmdOptionExists(argv, argv + argc, "-tmax"))
        tmax = std::stod(getCmdOption(argv, argv + argc, "-tmax"));

    std::function<double(double)> cx = [c_0](double x) {
        return (x >= 0 && x <= PI) ? 1 : c_0;
    };

    std::function<double(double)> kx = [k_0](double x) {
        return (x >= 0 && x <= PI) ? 1 : k_0;
    };

    std::function<double(double)> initf = [](double x) -> double {
        return sin(x);
    };

    simTime tc;
    tc.dt = 0.001;
    tc.endTime = tmax;
    WaveEquationProblem pb(file, initf, kx, cx, mesh_points, x_0, x_nx, tc);

    std::clock_t start;
    start = std::clock();
    pb.run();
    std::cout << "\nCalc Time: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC / 1000.0) << " ms" << std::endl;
    std::cout << tc << std::endl;
}

void run_proj1(int argc, char *argv[]) {
    std::string file = "test1.txt";
    if (cmdOptionExists(argv, argc + argv, "-f"))
        file = getCmdOption(argv, argc + argv, "-f");

    double init = 0.8;
    if (cmdOptionExists(argv, argv + argc, "-init_value"))
        init = std::stod(getCmdOption(argv, argv + argc, "-init_value"));

    int mesh_points = 1500;
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
    ParabolicPdeProblem pb(file, initf, one, one, F, mesh_points, x_0, x_nx, tc, vonNeumann);

    std::clock_t start;
    start = std::clock();
    pb.run();
    std::cout << "\nCalc Time: " << (std::clock() - start) / (double) (CLOCKS_PER_SEC / 1000.0) << " ms" << std::endl;
    std::cout << tc << std::endl;
}

int main(int argc, char *argv[]) {
    run_proj2(argc, argv);
    return 0;
}