//
// Created by Austin Clyde on 4/6/18.
//

#include "pdesolver.h"

#define PI 3.1415926535


/*
 * Generate Stiffness Matrix with linear piecewise functions
 */
TriDiag generateStiffnessMatrix(std::function<double(double)> k, int N, double mesh_width) {
  TriDiag S(N);

  for (int i = 0; i < N; i++) {
    for (int j = i - 1; j < N && j <= i + 1; j++) {
      if (j >= 0) {
        S(i, j) = stiffnessMatrixEntry(k, i, j, mesh_width);
      }
    }
  }

  return S;
}

/*
 * Uniform mesh, 1d case
 * s_{ij} = \int_{x_0}^{x_N} k(x) \phi_i' \phi_j' dx
 */
double stiffnessMatrixEntry(std::function<double(double)> k, int phi_i, int phi_j, double mesh_width) {
  return phi_i == phi_j ? 1.0 / (2.0 * mesh_width) : -1.0 / mesh_width;
}

TriDiag generateMassMatrix(std::function<double(double)> d, int N, double mesh_width) {
  TriDiag M(N);

  for (int i = 0; i < N; i++) {
    for (int j = i - 1; j < N && j <= i + 1; j++) {
      if (j >= 0) {
        M(i, j) = massMatrixEntry(d, i, j, mesh_width);
      }
    }
  }

  return M;
}

/* TODO:
 * Generate Mass Matrix with linear piecewise functions
 */
double massMatrixEntry(std::function<double(double)> d, int phi_i, int phi_j, double mesh_width) {
  return phi_i == phi_j ? 2.0 : -1.0; // \int_{x_0}^{x_N} d(x) \phi_i \phi_j dx
}

void solveHeatEquation1d(double x_0, double x_nx, int nx, int nt, double tmax, double alpha) {
  double L = x_nx - x_0;
  double dx = ((double) L) / (nx - 1);
  double dt = ((double) tmax) / (nt - 1);
  NumVec actual(nx);

  std::vector<std::string> params;
  params.push_back(std::to_string(L));
  params.push_back(std::to_string(nx));
  params.push_back(std::to_string(dx));
  params.push_back(std::to_string(tmax));
  params.push_back(std::to_string(nt));
  params.push_back(std::to_string(dt));
  params.push_back(std::to_string(alpha));
  writeParams("test.txt", params);

  TriDiag A(nx);
  for (int i = 1; i < nx - 1; i++) { //diag
    A(i, i) = (1 / dt) + (2 * alpha / (dx * dx));
  }
  A(0, 0) = 1;
  A(nx - 1, nx - 1) = 1;
  for (int i = 1; i < nx - 1; i++) { //above
    A(i, i + 1) = -1 * alpha / (dx * dx);
  }
  for (int i = 0; i < nx - 2; i++) { //below
    A(i + 1, i) = -1 * alpha / (dx * dx);
  }
  //std::cout<< A;

  NumVec U(nx);
  NumVec UOld(nx);
  for (int i = 0; i < nx - 1; i++)
    U[i] = sin(PI * i * dx / L);
  U[dx - 1] = 0;
  U[0] = 0;
  std::cout << U;
  double t = 0;
  for (int m = 1; m < nt; m++) {
    UOld = U;
    UOld = (1 / dt) * UOld;
    t = t + dt;
    U = solveTriDiagMatrix(A, UOld);
    for (int i = 1; i < nx - 1; i++) {
      actual[i] = sin(PI * i * dx / L) * exp(-1 * alpha * PI * PI * t / L);
    }
    writeUpdateStep("test.txt", U);
  }
  std::cout << "error x= " << L2norm(actual - U) << std::endl;
  std::cout << U;
  std::cout << t;
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

void writeUpdateStep(std::string name, NumVec a) {
  std::ofstream ofs(name, std::fstream::app);

  if (!ofs) {
    std::cout << "Error opening file for output" << std::endl;
    return;
  }
  ofs << a;
  ofs.close();
}
