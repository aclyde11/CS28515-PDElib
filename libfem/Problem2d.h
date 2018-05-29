//
// Created by Austin Clyde on 5/28/18.
//

#ifndef CS28515PROJ1_PROBLEM2D_H
#define CS28515PROJ1_PROBLEM2D_H

#include "VariMesh.h"
#include "Eigen/Dense"
#include "Eigen/Core"


class Problem2d {
public:
    VariMesh mesh;
    Eigen::VectorXd U;
    Eigen::VectorXd dU;

    Problem2d();

};

Eigen::VectorXd multiplyStiff(Eigen::VectorXd v, VariMesh mesh);

Eigen::VectorXd conjugateGradient(Eigen::MatrixXd C, Eigen::VectorXd g, Eigen::VectorXd y0);

Eigen::VectorXd conjugateGradientPp(Eigen::MatrixXd C, Eigen::VectorXd g, Eigen::VectorXd y0);


#endif //CS28515PROJ1_PROBLEM2D_H
