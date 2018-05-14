//
// Created by Austin Clyde on 4/24/18.
//

#include "ParabolicPdeProblem.h"

#define DEBUG false

ParabolicPdeProblem::ParabolicPdeProblem(std::string file_name,
                                         const std::function<double(double)> &init,
                                         const std::function<double(double)> &k,
                                         const std::function<double(double)> &c,
                                         const std::function<double(double, double)> &Fux,
                                         int nx,
                                         double x_0,
                                         double x_n,
                                         const simTime &timeControl,
                                         BoundryCondition bc)
        : file_name(file_name),
          k(k),
          c(c),
          Fux(Fux),
          nx(nx),
          x_0(x_0),
          x_n(x_n),
          timeControl(timeControl),
          bc(bc) {
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

void ParabolicPdeProblem::run() {
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

void ParabolicPdeProblem::advance() {
    double diff, newdt;
    NumVec dU_1, dU_2;
    timeControl.stepsSinceRejection = 0;
    do {
        timeControl.stepsRejected += 1;
        timeControl.stepsSinceRejection += 1;
        dU_1 = step(timeControl.time, timeControl.dt);
        dU_2 = step(timeControl.time, timeControl.dt / 2);
        dU_2 = dU_2 + step(timeControl.time + timeControl.dt / 2, timeControl.dt / 2);

        diff = L2norm(dU_1 - dU_2);
        if (diff <= 0.5 * timeControl.tol) {
            timeControl.dtOld = timeControl.dt;
            newdt = timeControl.dt * timeControl.agrow;
            timeControl.dt = (newdt > timeControl.dtmax) ? timeControl.dtmax : newdt;
        } else if (diff > timeControl.tol) {
            timeControl.dtOld = timeControl.dt;
            newdt = timeControl.dt * timeControl.ashrink;
            timeControl.dt = (newdt < timeControl.dtmin) ? timeControl.dtmin : newdt;
        }

        // Checks if error is "stuck" and proceeds with a half step after trying for a bit
        if (timeControl.stepsSinceRejection >= 100 && timeControl.dt == timeControl.dtmin) {
            std::cout << ("WARNING! taking half step\n");
            timeControl.dt = timeControl.dt / 2;
            diff = timeControl.tol;
            dU_1 = step(timeControl.time, timeControl.dt / 2);
            dU_2 = dU_1;
        }

    } while (diff > timeControl.tol);
    timeControl.stepsAccepted += 1;
    timeControl.time += timeControl.dt;
    dU = 0.5 * (dU_1 + dU_2);
    U = U + dU;
    writeUpdateStep(file_name, U, timeControl.time);
}

NumVec ParabolicPdeProblem::step(double t, double dt) {
    TriDiag DF, LH;
    NumVec F;
    NumVec RH(nx);

    F = linerizeF(U, Fux, dx);
    DF = linerizeDF(U, dU, F, Fux, dx);
    F = F + DF * dU;

    LH = (1.0 / dt) * M + S - DF;
    RH = (-1 * S) * U + F;

    if (bc == dirchlet) {
        LH(0, 0) = 1;
        LH(0, 1) = 0;
        LH(nx - 1, nx - 1) = 1;
        LH(nx - 1, nx - 2) = 0;
        RH[0] = 0;
        RH[nx - 1] = 0;
    }

    return solveTriDiagMatrix(LH, RH);
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

NumVec linerizeF(const NumVec &U, const std::function<double(double, double)> &F, double dx) {
    NumVec f(U.size());
    int nx = U.size();
    f[0] = (dx / 2);
    F(dx * 0.5, (U[1] + U[0]) / 2);
    for (int i = 1; i < nx - 1; i++) {
        f[i] += (dx / 2) * (F(dx * i - 0.5 * dx, (U[i] + U[i - 1]) / 2) + F(dx * i + 0.5 * dx, (U[i + 1] + U[i]) / 2));
    }
    f[nx - 1] = (dx / 2) * F(dx * nx - 0.5 * dx, (U[nx - 1] + U[nx - 2]) / 2);
    return f;
}

TriDiag linerizeDF(const NumVec &U,
                   const NumVec &dU,
                   const NumVec &Fvec,
                   const std::function<double(double, double)> &F,
                   double dx) {
    TriDiag DF(U.size());
    for (int row = 0; row < U.size() - 1; row++) {
        DF(row, row) += 0.25 * (F(row * dx + 0.5 * dx, (U[row] + U[row + 1]) / 2));
        DF(row, row + 1) = 0.25
                           * (F(row * dx + 0.5 * dx, (U[row] + U[row + 1]) / 2) +
                              F((row + 1) * dx + 0.5 * dx, (U[row] + U[row + 1]) / 2));
        DF(row + 1, row) = 0.25
                           * (F(row * dx + 0.5 * dx, (U[row] + U[row + 1]) / 2) +
                              F((row + 1) * dx + 0.5 * dx, (U[row] + U[row + 1]) / 2));
        DF(row + 1, row + 1) += 0.25 * (F((row + 1) * dx + 0.5 * dx, (U[row] + U[row + 1]) / 2)
                                        + F((row + 1) * dx + 0.5 * dx, (U[row] + U[row + 1]) / 2));
    }
    return DF;
}