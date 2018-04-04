//
// Created by Austin Clyde on 4/4/18.
//

#ifndef CS28515PROJ1_NUMVEC_H
#define CS28515PROJ1_NUMVEC_H

#include <cstdlib>
#include <iostream>
#include<vector>

#include "utility.h"

typedef std::vector<double> NumVec;


NumVec operator+(const NumVec& A, const NumVec& B); // A+B
NumVec operator*( double a, const NumVec& B); // a*B
NumVec operator-(const NumVec& A, const NumVec& B); // A-B
double operator,(const NumVec& A, const NumVec& B); // (A,B)

std::ostream& operator<<(std::ostream& os, const NumVec& A);

#endif //CS28515PROJ1_NUMVEC_H
