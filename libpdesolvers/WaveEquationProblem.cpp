//
// Created by Austin Clyde on 4/29/18.
//

#include "WaveEquationProblem.h"

#define DEBUG 0

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
    U[nx - 1] = 0;

    if (DEBUG)
        std::cout << "Init U: " << U;

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

NumVec periodic_solve(TriDiag &A, NumVec R) {
    int n = A.dim;
    NumVec X(n), X0(n), X1(n);
    A(0, 0) = 1;
    A(n - 1, n - 1) = 1;

    NumVec R0(R), R1(n);
    R0[0] = 0;
    R0[n - 1] = 0;
    X0 = solveTriDiagMatrix(A, R0);

    R1[0] = 1;
    R1[n - 1] = 1;
    X1 = solveTriDiagMatrix(A, R1);

    double alpha = -1 * ((A.d[0] + A.d[n - 1]) * X0[0] + A.a[0] * X0[1] + A.b[n - 2] * X0[n - 1]) /
                   ((A.d[0] + A.d[n - 1]) * X1[0] + A.a[0] * X1[1] + A.b[n - 2] * X1[n - 1]);

    X = X0 + alpha * X1;
    return X;
}
