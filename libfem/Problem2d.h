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

class Problem2d {
public:
    VariMesh mesh;
    Eigen::VectorXd A;


    Problem2d();

};

Eigen::VectorXd multiplyStiff(Eigen::VectorXd v, VariMesh mesh);

Eigen::VectorXd conjugateGradient(Eigen::MatrixXd C, Eigen::VectorXd g, Eigen::VectorXd y0);

Eigen::VectorXd generateStiffD(VariMesh mesh);

Eigen::VectorXd conjugateGradientPreconditioning(Eigen::MatrixXd C, Eigen::VectorXd g, Eigen::VectorXd y0);

Eigen::VectorXd conjugateGradientPreconditioningNM(Eigen::VectorXd b, Eigen::VectorXd x0,
                                                   std::function<Eigen::VectorXd(Eigen::VectorXd, VariMesh)> Amul,
                                                   Eigen::VectorXd D, VariMesh mesh);

void write_coords(Eigen::VectorXd X, VariMesh mesh, std::string file);

#endif //CS28515PROJ1_PROBLEM2D_H
