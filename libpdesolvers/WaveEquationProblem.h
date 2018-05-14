//
// Created by Austin Clyde on 4/29/18.
//

#ifndef CS28515PROJ1_WAVEEQUATIONPROBLEM_H
#define CS28515PROJ1_WAVEEQUATIONPROBLEM_H


//
// Created by Austin Clyde on 4/24/18.
//

#ifndef CS28515PROJ1_LINEARPARABOLICPROBLEM_H
#define CS28515PROJ1_LINEARPARABOLICPROBLEM_H

#include <string>
#include "numvec.h"
#include "tridiag.h"
#include "simtime.h"
#include "numericalMethods.h"
#include "utility.h"

#include <Eigen/Dense>

class WaveEquationProblem {

public:
    NumVec U, dU;
    std::string file_name;
    int nx;
    double x_0, x_n, dx;
    simTime timeControl;
    bool debug = true;
    std::function<double(double)> k, c;
    TriDiag M, S;

    WaveEquationProblem(std::string file_name,
                        const std::function<double(double)> &init,
                        const std::function<double(double)> &k,
                        const std::function<double(double)> &c,
                        int nx,
                        double x_0,
                        double x_n,
                        const simTime &timeControl
    );

    void run();

    NumVec step(double t, double dt);

    void advance();
};

NumVec periodic_solve(const TriDiag &A, const NumVec &R);

/*
 * Generates stiffness matrix for cu_t - ku_xx = F using midpoint method
 */
TriDiag generateStiffnessMatrixMidpoint(const std::function<double(double)> &k, int N, double dx, double x_0);

/*
 * Generates Mass matrix for cu_t - ku_xx = F using midpoint method
 */
TriDiag generateMassMatrixMidpoint(const std::function<double(double)> &d, int N, double dx, double x_0);



#endif //CS28515PROJ1_LINEARPARABOLICPROBLEM_H



#endif //CS28515PROJ1_WAVEEQUATIONPROBLEM_H
