//
// Created by Austin Clyde on 4/6/18.
//

#include "parabolicsolver.h"

TriDiag generateStiffnessMatrixMidpoint(std::function<double(double)> k, int N, double dx, double x_0) {
  TriDiag S(N);

  for (int row = 0; row < N - 1; row++) {
    S(row, row) += k(x_0 + dx * row + dx / 2) / dx;
    S(row, row + 1) = -1 * k(x_0 + dx * row + dx / 2) / dx;
    S(row + 1, row) = -1 * k(x_0 + dx * row + dx / 2) / dx;
    S(row + 1, row + 1) = k(x_0 + dx * row + dx / 2) / dx;
  }
  return S;
}

TriDiag generateMassMatrixMidpoint(std::function<double(double)> c, int N, double dx, double x_0) {
  TriDiag M(N);
  for (int row = 0; row < N - 1; row++) {
    M(row, row) += dx * 0.25 * c(x_0 + dx * row + dx / 2);
    M(row, row + 1) = dx * 0.25 * c(x_0 + dx * row + dx / 2);
    M(row + 1, row) = dx * 0.25 * c(x_0 + dx * row + dx / 2);
    M(row + 1, row + 1) += dx * 0.25 * c(x_0 + dx * row + dx / 2);
  }
  return M;
}

NumVec linearizeF(NumVec U, std::function<double(double, double)> F, double dx) {
  NumVec f(U.size());
  int nx = U.size();
  f[0] = (dx / 2);
  F(dx * 0.5, (U[1] + U[0]) / 2);
  for (int i = 1; i < nx - 1; i++) {
    f[i] += (dx / 2) * (F(dx * i - 0.5 * dx, (U[i] + U[i - 1]) / 2) + F(dx * i + 0.5 * dx, (U[i + 1] + U[i]) / 2));
  }
  f[nx - 1] = (dx / 2) * F(dx * nx - 0.5 * dx, (U[nx - 1] + U[nx - 2]) / 2);
  std::cout << "print U: " << ((U[nx - 1] + U[nx - 2]) / 2) << "F : " << f;

  return f;
}

TriDiag linearizeDF(NumVec U, NumVec dU, NumVec Fvec, std::function<double(double, double)> F, double dx) {
  TriDiag DF(U.size());
  for (int row = 0; row < U.size() - 1; row++) {
    DF(row, row) += 0.25 * (F(row * dx + 0.5 * dx, (U[row] + U[row + 1]) / 2));
    DF(row, row + 1) = 0.25
        * (F(row * dx + 0.5 * dx, (U[row] + U[row + 1]) / 2) + F((row + 1) * dx + 0.5 * dx, (U[row] + U[row + 1]) / 2));
    DF(row + 1, row) = 0.25
        * (F(row * dx + 0.5 * dx, (U[row] + U[row + 1]) / 2) + F((row + 1) * dx + 0.5 * dx, (U[row] + U[row + 1]) / 2));
    DF(row + 1, row + 1) += 0.25 * (F((row + 1) * dx + 0.5 * dx, (U[row] + U[row + 1]) / 2)
        + F((row + 1) * dx + 0.5 * dx, (U[row] + U[row + 1]) / 2));
  }
  return DF;
}

double simpson_integration(std::function<double(double)> f, double a, double b, int n_intervals) {
  double h = (b - a) / n_intervals;
  double s = f(a) + f(b);

  for (int i = 0; i < n_intervals; i += 2)
    s += 4 * f(a + i * h);
  for (int i = 1; i < n_intervals - 1; i += 2)
    s += 2 * f(a + i * h);

  return s * h / 3;
}

void writeParams(std::string name, std::vector<std::string> params) {
  std::ofstream ofs(name);

  if (!ofs) {
    std::cout << "Error opening file for output" << std::endl;
    return;
  }
  for (int i = 0; i < params.size(); i++) {
    ofs << params[i];
    if (i != params.size() - 1)
      ofs << "\t";
  }
  ofs << std::endl;
  ofs.close();
}

void writeUpdateStep(std::string name, NumVec a, double t) {
  std::ofstream ofs(name, std::fstream::app);

  if (!ofs) {
    std::cout << "Error opening file for output" << std::endl;
    return;
  }
  ofs << t << '\t';
  ofs << a;
  ofs.close();
}
