//
// Created by Austin Clyde on 4/4/18.
//

#ifndef CS28515PROJ1_TRIDIAG_H
#define CS28515PROJ1_TRIDIAG_H

#include <cstdlib>
#include <iostream>
#include<vector>

#include "utility.h"
#include "numvec.h"


struct TriDiag{
  unsigned int dim;
  NumVec a;
  NumVec d;
  NumVec b;
};

TriDiag operator+(const TriDiag& A, const TriDiag& B); // A+B
TriDiag operator*( double a, const TriDiag& B); // a*B
TriDiag operator-(const TriDiag& A, const TriDiag& B); // A-B

std::ostream& operator<<(std::ostream& os, const TriDiag& A);

#endif //CS28515PROJ1_TRIDIAG_H
