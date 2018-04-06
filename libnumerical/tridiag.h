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

class TriDiag {

 public:
  int dim;
  NumVec a;
  NumVec d;
  NumVec b;

  TriDiag(int n);
  TriDiag(const TriDiag &A);
  double &operator()(int i, int j);
  double operator()(int i, int j) const;
};

void map(TriDiag &A, std::function<double(double)> func);
TriDiag operator+(const TriDiag &A, const TriDiag &B); // A+B
TriDiag operator*(double a, const TriDiag &B); // a*B#
NumVec operator*(const TriDiag &A, const NumVec &v);
TriDiag operator-(const TriDiag &A, const TriDiag &B); // A-B

NumVec solveTriDiagMatrix(const TriDiag &A, const NumVec &d);

std::ostream &operator<<(std::ostream &os, const TriDiag A);

#endif //CS28515PROJ1_TRIDIAG_H
