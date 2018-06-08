//
// Created by Austin Clyde on 5/28/18.
//

#ifndef CS28515PROJ1_PROBLEM2D_H
#define CS28515PROJ1_PROBLEM2D_H

#include "VariMesh.h"
#include "Eigen/Dense"
#include "Eigen/Core"
#include "utility.h"
#include <string>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>

#define ITER_MAX 1000
#define EPS 1e-3

class FEMLaplacian {
public:
    VariMesh mesh;
    Eigen::VectorXd b;
    Eigen::VectorXd U;
    Eigen::VectorXd D;
    std::function<Eigen::VectorXd(Eigen::VectorXd, VariMesh)> Amul;
    std::string file;

    FEMLaplacian(VariMesh mesh, std::function<double(int, int, VariMesh)> u_init, std::string file,
                 std::function<Eigen::VectorXd(Eigen::VectorXd, VariMesh)> Amul,
                 std::function<Eigen::VectorXd(VariMesh)> GenDiag);

    void run();

    Eigen::VectorXd conjugateGradientPreconditioningNM(Eigen::VectorXd x0);

    void run(Eigen::VectorXd init);

    Eigen::VectorXd conjugateGradientPreconditioningNMEx0(Eigen::VectorXd x0);

    void write_out_mesh(std::string file);
};

Eigen::VectorXd multiplyStiffExample0(Eigen::VectorXd v, VariMesh mesh);

Eigen::VectorXd multiplyStiffExample1(Eigen::VectorXd v, VariMesh mesh);

Eigen::VectorXd conjugateGradient(Eigen::MatrixXd C, Eigen::VectorXd g, Eigen::VectorXd y0);

Eigen::VectorXd generateStiffDEx1(VariMesh mesh);

Eigen::VectorXd generateStiffDEx0(VariMesh mesh);

Eigen::VectorXd conjugateGradientPreconditioning(Eigen::MatrixXd C, Eigen::VectorXd g, Eigen::VectorXd y0);


void write_coords(Eigen::VectorXd X, VariMesh mesh, std::string file);

#endif //CS28515PROJ1_PROBLEM2D_H
