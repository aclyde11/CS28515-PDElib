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

    VariMesh mesh(x_nodes, y_nodes, M_PI / (x_nodes - 1), M_PI / (y_nodes - 1));

    Eigen::MatrixXd init_m(30, 30);
    for (int i = 0; i < mesh.x_nodes; i++) {
        for (int j = 0; j < mesh.y_nodes; j++) {
            init_m(i, j) = 5 * sin(mesh.node_loc(i, j).x) * sin(2 * mesh.node_loc(i, j).y) * mesh.dx[0] * mesh.dy[0];
        }
    }
    for (int i = 0; i < mesh.x_nodes; i++) {
        init_m(i, y_nodes - 1) = 0;
        init_m(i, 0) = 0;
    }
    for (int i = 0; i < mesh.y_nodes; i++) {
        init_m(0, i) = 0;
        init_m(x_nodes - 1, i) = 0;
    }

    Eigen::VectorXd F = Eigen::Map<Eigen::VectorXd>(init_m.data(), init_m.size());
    Eigen::VectorXd guess(F.size());
    Eigen::VectorXd ans = conjugateGradientPreconditioningNM(F, guess, multiplyStiff, generateStiffD(mesh), mesh);
    write_coords(F, mesh, "coords_f.txt");
    write_coords(ans, mesh, "coords_ans.txt");
}

