//
// Created by Austin Clyde on 4/4/18.
//



#include "tridiag.h"

TriDiag::TriDiag() {
  dim = 1;
  d.push_back(0.0);
  a.push_back(0.0);
  b.push_back(0.0);
}

TriDiag::TriDiag(int n) {
  dim = n;
  d.resize(n);
  std::fill_n(d.begin(), n, 0.0);
  b = d;
  a = d;
}

TriDiag::TriDiag(const TriDiag &A) {
  dim = A.dim;
  d = A.d;
  a = A.a;
  b = A.b;
}

TriDiag::TriDiag(int n, double p) {
  dim = n;
  d.resize(n);
  std::fill_n(d.begin(), n, p);
  b = d;
  a = d;
}

void TriDiag::randomize() {
  randomize(1, 100);
}

void TriDiag::randomize(int start, int end) {
  a = randomVector(dim - 1, start, end);
  b = randomVector(dim - 1, start, end);
  d = randomVector(dim, start, end);
}

void TriDiag::randomizeReals(double start, double end) {
  a = randomRealVector(dim - 1, start, end);
  b = randomRealVector(dim - 1, start, end);
  d = randomRealVector(dim, start, end);
}

bool TriDiag::operator==(const TriDiag& B) {
  if (dim != B.dim)
    return false;
  return (d == B.d) && (a == B.a) && (b == B.b);
}

double &TriDiag::operator()(int i, int j) {
  if (i < 0 || j < 0 || i >= dim || j >= dim)
    error("TriDiag ,: indices out of matrix bounds!");
  else if (abs(i - j) >= 2)
    error("TriDiag ,: we only know tridiagonal entries! nothing else");

  if (i == j)
    return d[i];
  else if (i + 1 == j)
    return a[i];
  else
    return b[j];
}

double TriDiag::operator()(int i, int j) const {
  if (i < 0 || j < 0 || i >= dim || j >= dim)
    error("TriDiag ,: indicies out of matrix bounds!");

  if (i == j)
    return d[i];
  else if (i + 1 == j)
    return a[i];
  else if (j + 1 == i)
    return b[j];
  else
    return 0.0;
}

void map(TriDiag &A, std::function<double(double)> func) {
  std::transform(A.d.begin(), A.d.end(), A.d.begin(), func);
  std::transform(A.a.begin(), A.a.end(), A.a.begin(), func);
  std::transform(A.b.begin(), A.b.end(), A.b.begin(), func);
}

TriDiag binary(const TriDiag &A, const TriDiag &B, std::function<double(double, double)> func) {
  if (A.dim != B.dim)
    error("TriDiag +: incompatible sizes");

  TriDiag C(A.dim);
  std::transform(A.d.begin(), A.d.end(), B.d.begin(), C.d.begin(), func);
  std::transform(A.a.begin(), A.a.end(), B.a.begin(), C.a.begin(), func);
  std::transform(A.b.begin(), A.b.end(), B.b.begin(), C.b.begin(), func);
  return C;
}

TriDiag operator*(double a, const TriDiag &B) {
  TriDiag C(B);
  map(C, [a](double b) -> double { return a * b; });
  return C;
}

//TODO: Make better for tridiag, not general multplciation;
NumVec operator*(const TriDiag &A, const NumVec &v) {
  if (A.dim != v.size()) {
    error("TriDiag * vector: bad sizes");
  }
  NumVec y(v);
  double d;
  for (int i = 0; i < A.dim; i++) {
    d = 0.0;
    for (int j = 0; j < A.dim; j++) {
      d += A(i, j) * v[j];
    }
    y[i] = d;
  }
  return y;
}

TriDiag operator+(const TriDiag &A, const TriDiag &B) {
  return binary(A, B, std::plus<double>());
}

TriDiag operator-(const TriDiag &A, const TriDiag &B) {
  return binary(A, B, std::minus<double>());
}

double TriDiag::det() {
  return det_f(dim-1);
}

double TriDiag::det_f(int i) {
  if(i == -2)
    return 0;
  else if (i == -1)
    return 1;

  return d[i] * det_f(i-1) - a[i-1]*b[i-1]*det_f(i-2);
}

/**
 * Thomas algorithm (named after Llewellyn Thomas), is a simplified form of Gaussian elimination
 */
NumVec solveTriDiagMatrix(const TriDiag &A, const NumVec &y) {
  if (A.dim != y.size()) {
    error("SolveTriDiagMatrix, vector and matrix are different dims.");
  }

  NumVec d(y);
  NumVec x(d);
  NumVec c(A.a);
  int n = A.dim;

  //forward sweep
  c[0] = c[0] / A.d[0];
  for (int i = 1; i < n - 1; i++) {
    c[i] = c[i] / (A.d[i] - A.b[i - 1] * c[i - 1]);
  }

  d[0] = d[0] / A.d[0];
  for (int i = 1; i < n; i++) {
    d[i] = (d[i] - A.b[i - 1] * d[i - 1]) / (A.d[i] - A.b[i - 1] * c[i - 1]);
  }

  //back sub
  x[n - 1] = d[n - 1];
  for (int i = n - 2; i >= 0; i--) {
    x[i] = d[i] - c[i] * x[i + 1];
  }
  return x;
}

std::ostream &operator<<(std::ostream &os, const TriDiag A) {
  for (int i = 0; i < A.dim; i++) {
    for (int j = 0; j < A.dim; j++) {
      std::cout << A(i, j) << "\t";
    }
    std::cout << std::endl;
  }
  os << "\n";
  return os;
}