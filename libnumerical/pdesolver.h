//
// Created by Austin Clyde on 4/6/18.
//

#ifndef CS28515PROJ1_PDESOLVER_H
#define CS28515PROJ1_PDESOLVER_H

#include "numvec.h"
#include "tridiag.h"
#include <math.h>
#include <fstream>
#include <string>
#include <iterator>

TriDiag generateStiffnessMatrix(std::function<double(double)> k, int N, double mesh_width);
double stiffnessMatrixEntry(std::function<double(double)> k, int phi_i, int phi_j, double mesh_width);

TriDiag generateMassMatrix(std::function<double(double)> d, int N, double mesh_width);
double massMatrixEntry(std::function<double(double)> d, int phi_i, int phi_j, double mesh_width);

void solveHeatEquation1d(double x_0,
                         double x_nx,
                         int nx,
                         int nt,
                         double tmax,
                         double alpha,
                         std::function<double(double)> init);
void writeParams(std::string name, std::vector<std::string> params);
void writeUpdateStep(std::string name, NumVec a);

#endif //CS28515PROJ1_PDESOLVER_H
