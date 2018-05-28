//
// Created by Austin Clyde on 5/28/18.
//

#include "Mesh.h"

Mesh::Mesh(int x, int y, double x_0, double y_0, double dx, double dy) :
        x_nodes(x), y_nodes(y), x_0(x_0), y_0(y_0) {
    nodes = x * y;

    this->dx.resize(x - 1);
    this->dy.resize(y - 1);
    for (int i = 0; i < x_nodes - 1; i++)
        this->dx(i) = dx;
    for (int j = 0; j < y_nodes - 1; j++)
        this->dy(j) = dy;
}

Mesh::Mesh(int x, int y, double dx, double dy) :
        x_nodes(x), y_nodes(y), x_0(0), y_0(0) {
    nodes = x * y;

    this->dx.resize(x - 1);
    this->dy.resize(y - 1);
    for (int i = 0; i < x_nodes - 1; i++)
        this->dx(i) = dx;
    for (int j = 0; j < y_nodes - 1; j++)
        this->dy(j) = dy;

}

Mesh::Mesh(int x, int y, Eigen::VectorXd dx, Eigen::VectorXd dy) :
        x_nodes(x), y_nodes(y), x_0(0), y_0(0), dx(dx), dy(dy) {
    nodes = x * y;
}

Mesh::Mesh(int x, int y, double x_0, double y_0, Eigen::VectorXd dx, Eigen::VectorXd dy) :
        x_nodes(x), y_nodes(y), x_0(x_0), y_0(y_0), dx(dx), dy(dy) {
    nodes = x * y;
}

Loc Mesh::node_d(int x, int y) {
    if (x >= x_nodes || x < 0 || y < 0 || y >= y_nodes)
        error("Nodes_d out of bounds!");
    else if (abs(x - y) > 1)
        error("Nodes too far");
    return Loc(dx(x), dy(y));
}

Loc Mesh::node_loc(int x, int y) {
    if (x >= x_nodes || x < 0 || y < 0 || y >= y_nodes)
        error("Nodes_d out of bounds!");
    return Loc(dx.segment(0, x).sum() + x_0, dy.segment(0, y).sum() + y_0);
}


