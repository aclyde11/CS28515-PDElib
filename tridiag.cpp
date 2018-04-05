//
// Created by Austin Clyde on 4/4/18.
//


#include "tridiag.h"

TriDiag::TriDiag(int n) {
  dim = n;
  std::fill_n(d.begin(), n, 0.0);
  std::fill_n(a.begin(), n-1, 0.0);
  std::fill_n(b.begin(), n-1, 0.0);
}

TriDiag::TriDiag(const TriDiag& A) {
  dim = A.dim;
  d = A.d;
  a = A.a;
  b = A.b;
}

double &TriDiag::operator()(int i, int j) {
  if (i < 0 || j < 0 || i >= dim || j >= dim)
    error("TriDiag ,: indicies out of matrix bounds!");
  else if (abs(i - j) >= 2)
    error("TriDiag ,: we only know tridigonal entries! nothing else");

  if (i == j)
    return d[i];
  else if (i+1 == j)
    return a[i];
  else
    return b[j];
}

double TriDiag::operator()(int i, int j) const {
  if (i < 0 || j < 0 || i >= dim || j >= dim)
    error("TriDiag ,: indicies out of matrix bounds!");

  if (i == j)
    return d[i];
  else if (i+1 == j)
    return a[i];
  else if (j+1 == i)
    return b[j];
  else
    return 0.0;
}

void map(TriDiag& A, std::function<double(double)> func) {
  std::transform (A.d.begin(), A.d.end(), A.d.begin(), func);
  std::transform (A.a.begin(), A.a.end(), A.a.begin(), func);
  std::transform (A.b.begin(), A.b.end(), A.b.begin(), func);
}

TriDiag binary(const TriDiag& A, const TriDiag& B, std::function<double(double, double)> func) {
  if( A.dim != B.dim)
    error("TriDiag +: incompatible sizes");

  TriDiag C(A.dim);
  std::transform (A.d.begin(), A.d.end(), B.d.begin(), C.d.begin(), func);
  std::transform (A.a.begin(), A.a.end(), B.a.begin(), C.a.begin(), func);
  std::transform (A.b.begin(), A.b.end(), B.b.begin(), C.b.begin(), func);
  return C;
}

TriDiag operator*( double a, const TriDiag& B) {
  TriDiag C(B);
  map(C, [a](double b)->double { return a * b; } );
  return C;
}

TriDiag operator+(const TriDiag& A, const TriDiag& B) {
  return binary(A, B, std::plus<double>());
}

TriDiag operator-(const TriDiag& A, const TriDiag& B) {
  return binary(A, B, std::minus<double>());
}

std::ostream& operator<<(std::ostream& os, const TriDiag A){
  for (int i = 0; i < A.dim; i++){
    for (int j = 0; j < A.dim; j++) {
      std::cout << A(i,j) << "\t";
    }
    std::cout << std::endl;
  }
  os << "\n";
  return os;
}