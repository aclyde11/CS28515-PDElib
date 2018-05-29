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

    Eigen::MatrixXd result = conjugateGradient(A, b, c);
    std::cout << "Withput PC = \n " << result << std::endl;
    result = conjugateGradientPp(A, b, c);
    std::cout << "with PC = \n " << result << std::endl;

    A.resize(3, 3);
    A << 4, -1, 1,
            -1, 4, -2,
            1, -2, 4;

    b.resize(3);
    b << 12, -1, 5;

    c.resize(3);
    c << 0, 0, 0;
    result = conjugateGradient(A, b, c);
    std::cout << "Withput PC = \n " << result << std::endl;
    result = conjugateGradientPp(A, b, c);
    std::cout << "with PC = \n " << result << std::endl;
}