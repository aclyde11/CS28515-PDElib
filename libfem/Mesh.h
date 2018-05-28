//
// Created by Austin Clyde on 5/28/18.
//

#ifndef CS28515PROJ1_MESH_H
#define CS28515PROJ1_MESH_H

#include "utility.h"
#include "Eigen/Dense"

class Loc {
public:
    double x;
    double y;

    Loc() : x(0), y(0) { ; };

    Loc(double x, double y) : x(x), y(y) { ; };

    std::ostream &operator<<(std::ostream &os) {
        os << "(" << x << ", " << y << ")";
        return os;
    };
};

class Mesh {

public:
    Eigen::VectorXd dx, dy;
    double x_0, y_0;
    int x_nodes, y_nodes, nodes;

    Mesh(int x, int y, double dx, double dy);

    Mesh(int x, int y, Eigen::VectorXd dx, Eigen::VectorXd dy);

    Mesh(int x, int y, double x_0, double y_0, double dx, double dy);

    Mesh(int x, int y, double x_0, double y_0, Eigen::VectorXd dx, Eigen::VectorXd dy);


    Loc node_d(int x, int y);

    Loc node_loc(int x, int y);
};


#endif //CS28515PROJ1_MESH_H
