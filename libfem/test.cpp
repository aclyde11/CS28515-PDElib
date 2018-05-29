//
// Created by Austin Clyde on 5/23/18.
//
#include <iostream>
#include "Eigen/Dense"

#include "test.h"
#include "VariMesh.h"
#include "Problem2d.h"


using Eigen::MatrixXd;
void test_eigen() {


    Eigen::MatrixXd A(2, 2);
    A << 2, -1, -1, 2;

    Eigen::VectorXd b(2);
    b << 1, 0;

    Eigen::VectorXd c(2);
    c << 0, 0;

    Eigen::VectorXd result = conjugateGradientPp(A, b, c);
    std::cout << result;
}