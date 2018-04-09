//
// Created by Todd Dupont
//

#ifndef CS28515PROJ1_NUMVEC_H
#define CS28515PROJ1_NUMVEC_H

#include <cstdlib>
#include <functional>
#include <iostream>
#include <vector>

#include "utility.h"

typedef std::vector<double> NumVec;

NumVec operator+(const NumVec& A, const NumVec& B); // A+B
NumVec operator*( double a, const NumVec& B); // a*B
NumVec operator-(const NumVec& A, const NumVec& B); // A-B
double operator,(const NumVec& A, const NumVec& B); // (A,B)

NumVec randomVector(int length, int start, int end);
NumVec randomRealVector(int length, double start, double end);

std::ostream& operator<<(std::ostream& os, const NumVec& A);

#endif //CS28515PROJ1_NUMVEC_H
