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

#endif //CS28515PROJ1_numericalMethods_H
