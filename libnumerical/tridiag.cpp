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

bool TriDiag::operator==(const TriDiag &B) const {
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

void TriDiag::map(std::function<double(double)> func) {
  std::transform(d.begin(), d.end(), d.begin(), func);
  std::transform(a.begin(), a.end(), a.begin(), func);
  std::transform(b.begin(), b.end(), b.begin(), func);
}

TriDiag TriDiag::binary(const TriDiag &B, std::function<double(double, double)> func) const {
  if (dim != B.dim)
    error("TriDiag +: incompatible sizes");

  TriDiag C(dim);
  std::transform(d.begin(), d.end(), B.d.begin(), C.d.begin(), func);
  std::transform(a.begin(), a.end(), B.a.begin(), C.a.begin(), func);
  std::transform(b.begin(), b.end(), B.b.begin(), C.b.begin(), func);
  return C;
}

TriDiag operator*(double a, const TriDiag &B) {
  TriDiag C(B);
  C.map([a](double b) -> double { return a * b; });
  return C;
}

//TODO: Make better for tridiag, not general multplciation;
NumVec TriDiag::operator*(const NumVec &v) const {
  if (dim != v.size()) {
    error("TriDiag * vector: bad sizes");
  }
  NumVec y(v);
  double d;
  for (int i = 0; i < dim; i++) {
    d = 0.0;
    for (int j = i - 1; j < dim && j <= i + 1; j++) {
      if (j >= 0)
        d += this->operator()(i, j) * v[j];
    }
    y[i] = d;
  }
  return y;
}

//TODO: Make better for tridiag, not general multplciation;
TriDiag TriDiag::operator*(const TriDiag &B) const {
  if (dim != B.dim)
    error("TriDiag * TriDiag: bad sizes");

  TriDiag P(dim);
  for (int i = 0; i < dim; i++) {
    for (int j = i - 1; j < dim && j <= i + 1; j++) {
      if (j >= 0)
        P(i, j) = (this->getRow(i), B.getCol(j)); //NumVec inner product
    }
  }
  return P;
}

NumVec TriDiag::getRow(int r) const {
  NumVec row(dim);
  for (int i = 0; i < dim; i++)
    row[i] = this->operator()(r, i);
  return row;
}

NumVec TriDiag::getCol(int c) const {
  NumVec col(dim);
  for (int i = 0; i < dim; i++)
    col[i] = this->operator()(i, c);
  return col;
}

TriDiag TriDiag::operator+(const TriDiag &B) const {
  return this->binary(B, std::plus<double>());
}

TriDiag TriDiag::operator-(const TriDiag &B) const {
  return this->binary(B, std::minus<double>());
}

double TriDiag::det() {
  return det_f(dim - 1);
}

double TriDiag::det_f(int i) {
  if (i == -2)
    return 0;
  else if (i == -1)
    return 1;

  return d[i] * det_f(i - 1) - a[i - 1] * b[i - 1] * det_f(i - 2);
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

std::ostream &operator<<(std::ostream &os, const TriDiag &A) {
  for (int i = 0; i < A.dim; i++) {
    for (int j = 0; j < A.dim; j++) {
      std::cout << A(i, j) << "\t";
    }
    std::cout << std::endl;
  }
  os << "\n";
  return os;
}