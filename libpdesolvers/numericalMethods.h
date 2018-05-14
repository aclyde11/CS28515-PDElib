//
// Created by Austin Clyde on 4/6/18.
//

#ifndef CS28515PROJ1_numericalMethods_H
#define CS28515PROJ1_numericalMethods_H

#include "numvec.h"
#include "tridiag.h"
#include <cmath>
#include <fstream>
#include <string>
#include <iterator>

/*
 * approximates 1d interval \int_a^b f(x) dx
 */
double simpson_integration(std::function<double(double)> f, double a, double b, int n_intervals);

/*
 * Generates stiffness matrix
 */
TriDiag generateStiffnessMatrixMidpoint(const std::function<double(double)> &k, int N, double dx, double x_0);

/*
 * Generates Mass matrix
 */
TriDiag generateMassMatrixMidpoint(const std::function<double(double)> &d, int N, double dx, double x_0);

#endif //CS28515PROJ1_numericalMethods_H
