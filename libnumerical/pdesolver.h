//
// Created by Austin Clyde on 4/6/18.
//

#ifndef CS28515PROJ1_PDESOLVER_H
#define CS28515PROJ1_PDESOLVER_H

#include "numvec.h"
#include "tridiag.h"

TriDiag generateStiffnessMatrix(std::function<double(double)> k, int N, double mesh_width);
double stiffnessMatrixEntry(std::function<double(double)> k, int phi_i, int phi_j, double mesh_width);

TriDiag generateMassMatrix(std::function<double(double)> d, int N, double mesh_width);
double massMatrixEntry(std::function<double(double)> d, int phi_i, int phi_j, double mesh_width);



#endif //CS28515PROJ1_PDESOLVER_H
