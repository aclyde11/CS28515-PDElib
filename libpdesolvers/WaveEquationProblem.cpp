//
// Created by Austin Clyde on 4/29/18.
//

#include "WaveEquationProblem.h"

#define DEBUG 1

WaveEquationProblem::WaveEquationProblem(std::string file_name,
                                         const std::function<double(double)> &init,
                                         const std::function<double(double)> &k,
                                         const std::function<double(double)> &c,
                                         int nx,
                                         double x_0,
                                         double x_n,
                                         const simTime &timeControl
)
        : file_name(file_name),
          k(k),
          c(c),
          nx(nx),
          x_0(x_0),
          x_n(x_n),
          timeControl(timeControl) {
    dx = (x_n - x_0) / (nx - 1);

    //initializes the vector U based on the init function
    for (int i = 0; i < nx; i++) {
        U.push_back(init(i * dx));
        dU.push_back(0.0);
    }

    // mass and stiffness do not depend on U
    M = generateMassMatrixMidpoint(c, nx, dx, x_0);
    S = generateStiffnessMatrixMidpoint(k, nx, dx, x_0);
}

void WaveEquationProblem::run() {
    std::vector<std::string> params;
    params.push_back(std::to_string(x_n - x_0));
    params.push_back(std::to_string(nx));
    params.push_back(std::to_string(dx));
    params.push_back(std::to_string(timeControl.endTime));
    params.push_back(std::to_string(timeControl.endTime / timeControl.dt));
    params.push_back(std::to_string(timeControl.dt));
    params.push_back(std::to_string(k(1.0)));
    writeParams(file_name, params);
    writeUpdateStep(file_name, U, 0);

    while (timeControl.time < timeControl.endTime) {
        advance();
        if (DEBUG)
            std::cout << U;
    }
    std::cout << "Final U: " << U;
}

void WaveEquationProblem::advance() {
    dU = step(timeControl.time, timeControl.dt);
    timeControl.stepsAccepted += 1;
    timeControl.time += timeControl.dt;
    U = U + dU;
    writeUpdateStep(file_name, U, timeControl.time);
}

NumVec WaveEquationProblem::step(double t, double dt) {
    TriDiag DF, LH;
    NumVec dU_1(nx), dU_2(nx), RH(nx);
    double theta = 0.25;
    double dt_sq = dt * dt;

    LH = (1.0 / dt_sq) * (M + ((dt_sq * theta) * S));
    RH = -1.0 * (S * U) + (1.0 / dt_sq) * (M + ((dt_sq * theta) * S)) * dU;

    return periodic_solve(LH, RH); //dU new
}

NumVec periodic_solve(const TriDiag &A, NumVec R) {
    int n = A.dim;
    NumVec X(n), X1(n), X0(n);

    //Modify A
    TriDiag Amod;

    //Modify R
    R[0] += R[n - 1];
    R[n - 1] = 0;

    //Solve X0
    X0[0] = 0;
    X[n - 1] = 0;
    //solve(Amod * * X0 = R)
    //TODO

    //Solve X1
    X1[0] = 1;
    X1[n - 1] = 1;
    //solve(Amod * X1 = zeros(n)
    //TODO

    double a;
    //Solve a so X = X0 + aX1 is a solution
    a = -1.0 * (Amod(0, 0) * X0[0] + Amod(0, 1) * X0[1] + Amod(0, n - 2) * X0[n - 2]) /
        (Amod(0, 0) * X1[0] + Amod(0, 1) * X1[1] + Amod(0, n - 2) * X1[n - 2]);

    return X0 + a * X1;
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