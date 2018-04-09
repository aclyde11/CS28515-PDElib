//
// Created by Austin Clyde on 4/6/18.
//

#include "pdesolver.h"

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
  return 0; // \int_{x_0}^{x_N} d(x) \phi_i \phi_j dx
}