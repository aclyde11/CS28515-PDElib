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

class WaveEquationProblem {

public:
    NumVec U, dU;
    std::string file_name;
    int nx;
    double x_0, x_n, dx;
    simTime timeControl;
    bool debug = true;

    void run();

    NumVec step(double t, double dt);

    void advance();
};


#endif //CS28515PROJ1_LINEARPARABOLICPROBLEM_H



#endif //CS28515PROJ1_WAVEEQUATIONPROBLEM_H
