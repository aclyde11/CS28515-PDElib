//
// Created by Austin Clyde on 4/6/18.
//

#ifndef CS28515PROJ1_PDESOLVER_H
#define CS28515PROJ1_PDESOLVER_H

#include "numvec.h"
#include "tridiag.h"
#include <cmath>
#include <fstream>
#include <string>
#include <iterator>

/*
 * Generates stiffness matrix for cu_t - ku_xx = F using midpoint method
 */
TriDiag generateStiffnessMatrixMidpoint(std::function<double(double)> k, int N, double dx, double x_0);

/*
 * Generates Mass matrix for cu_t - ku_xx = F using midpoint method
 */
TriDiag generateMassMatrixMidpoint(std::function<double(double)> d, int N, double dx, double x_0);

/*
 * This function produces a linear F(x,u) vector from the equation u_t - u_{xx} = F(x,u)
 */
TriDiag linearizeDF(NumVec U, NumVec dU, NumVec Fvec, std::function<double(double, double)> F, double dx);

/*
 * Produces the tridiagional matrix for dF, see linearizeF
 */
NumVec linearizeF(NumVec U, std::function<double(double, double)> F, double dx);

/*
 * approximates 1d interval \int_a^b f(x) dx
 */
double simpson_integration(std::function<double(double)> f, double a, double b, int n_intervals);

/*
 * Writes PDE parameters out to file for python plotting
 */
void writeParams(std::string name, std::vector<std::string> params);

/*
 * Writes vector U at t to data file
 */
void writeUpdateStep(std::string name, NumVec a, double t);

#endif //CS28515PROJ1_PDESOLVER_H
