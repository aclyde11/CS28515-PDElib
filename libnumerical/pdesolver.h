//
// Created by Austin Clyde on 4/6/18.
//

#ifndef CS28515PROJ1_PDESOLVER_H
#define CS28515PROJ1_PDESOLVER_H

#include "numvec.h"
#include "tridiag.h"

TriDiag generateStiffnessMatrix(std::function<double(double)> f, int N);
TriDiag generateMassMatrix(std::function<double(double)> f, int N);

double midpointmethod(std::function<double(double)> f, double x_i, double x_n);


#endif //CS28515PROJ1_PDESOLVER_H
