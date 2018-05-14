//
// Created by Austin Clyde on 4/29/18.
//

#include "WaveEquationProblem.h"

NumVec WaveEquationProblem::step(double t, double dt) {
    return NumVec();
}

void WaveEquationProblem::run() {

}

void WaveEquationProblem::advance() {

}

NumVec periodic_solve(const TriDiag &A, const NumVec &R) {
    //update A by perodic conditions

    //set RHS(0)=RHS(N)=0, Solve X0

    //set RHS(0)=RHS(N)=1 solve X1

    //X = X0 + alpha X1

    return NumVec();
}
