//
// Created by Austin Clyde on 4/6/18.
//

#include "numericalMethods.h"

double simpson_integration(std::function<double(double)> f, double a, double b, int n_intervals) {
    double h = (b - a) / n_intervals;
    double s = f(a) + f(b);

    for (int i = 0; i < n_intervals; i += 2)
        s += 4 * f(a + i * h);
    for (int i = 1; i < n_intervals - 1; i += 2)
        s += 2 * f(a + i * h);

    return s * h / 3;
}

TriDiag generateStiffnessMatrixMidpoint(const std::function<double(double)> &k, int N, double dx, double x_0) {
    TriDiag S(N);

    for (int row = 0; row < N - 1; row++) {
        S(row, row) += k(x_0 + dx * row + dx / 2) / dx;
        S(row, row + 1) = -1 * k(x_0 + dx * row + dx / 2) / dx;
        S(row + 1, row) = -1 * k(x_0 + dx * row + dx / 2) / dx;
        S(row + 1, row + 1) = k(x_0 + dx * row + dx / 2) / dx;
    }
    return S;
}

TriDiag generateMassMatrixMidpoint(const std::function<double(double)> &c, int N, double dx, double x_0) {
    TriDiag M(N);
    for (int row = 0; row < N - 1; row++) {
        M(row, row) += dx * 0.25 * c(x_0 + dx * row + dx / 2);
        M(row, row + 1) = dx * 0.25 * c(x_0 + dx * row + dx / 2);
        M(row + 1, row) = dx * 0.25 * c(x_0 + dx * row + dx / 2);
        M(row + 1, row + 1) += dx * 0.25 * c(x_0 + dx * row + dx / 2);
    }
    return M;
}
