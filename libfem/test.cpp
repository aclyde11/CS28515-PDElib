//
// Created by Austin Clyde on 5/23/18.
//
#include <iostream>
#include "Eigen/Core"

#include "test.h"
#include "VariMesh.h"
#include "FEMLaplacian.h"
#include <cmath>


using Eigen::MatrixXd;

double example0(int i, int j, VariMesh mesh) {
    if (i == 0 || i == mesh.x_nodes - 1 || j == 0 || j == mesh.y_nodes - 1)
        return 0;
    return 5 * sin(mesh.node_loc(i, j).x) * sin(2 * mesh.node_loc(i, j).y) * mesh.dx[0] * mesh.dy[0];
}

double example1(int i, int j, VariMesh mesh) {
    if (i == 0)
        return 1;
    else if (i == mesh.x_nodes - 1)
        return 0;
    else if (j == mesh.y_nodes - 1)
        return 1 - mesh.node_loc(i, j).x;
    else if (mesh.node_loc(i, j).x <= 0.5)
        return 5;
    return 1;
}

void test_eigen() {
    int x_nodes = 30;
    int y_nodes = 50;

    //Example 1
    VariMesh mesh(x_nodes, y_nodes, 1.0 / (x_nodes - 1), 1.0 / (y_nodes - 1));
    FEMLaplacian fm(mesh, example1, "coords_ans.txt", multiplyStiffExample0, generateStiffD);
    fm.write_out_mesh("coords_f.txt");
    fm.run();

    /*
    //Example 0
    VariMesh mesh(x_nodes, y_nodes, M_PI / (x_nodes - 1), M_PI / (y_nodes - 1));
    FEMLaplacian fm(mesh, example1, "coords_ans.txt", multiplyStiffExample0, generateStiffD);
    fm.run;
    */
}

