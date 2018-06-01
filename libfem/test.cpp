//
// Created by Austin Clyde on 5/23/18.
//
#include <iostream>
#include "Eigen/Core"

#include "test.h"
#include "VariMesh.h"
#include "Problem2d.h"
#include <cmath>


using Eigen::MatrixXd;
void test_eigen() {
    int x_nodes = 30;
    int y_nodes = 30;
    VariMesh mesh(x_nodes, y_nodes, M_PI / x_nodes, M_PI / y_nodes);
    Eigen::MatrixXd init_m(30, 30);
    for (int i = 0; i < mesh.x_nodes; i++) {
        for (int j = 0; j < mesh.y_nodes; j++) {
            init_m(i, j) = 5 * sin(mesh.node_loc(i, j).x) * sin(2 * mesh.node_loc(i, j).y);
        }
    }
    Eigen::VectorXd F = Eigen::Map<Eigen::VectorXd>(init_m.data(), init_m.size());
    Eigen::VectorXd guess(F.size());
    Eigen::VectorXd ans = conjugateGradientPreconditioningNM(F, guess, multiplyStiff, generateStiffD(mesh), mesh);
    write_coords(F, mesh, "coords.txt");
}

