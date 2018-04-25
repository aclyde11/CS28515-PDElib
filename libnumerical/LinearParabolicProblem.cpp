//
// Created by Austin Clyde on 4/24/18.
//

#include "LinearParabolicProblem.h"

LinearParabolicProblem::LinearParabolicProblem(std::string file_name,
                                               const std::function<double(double)> &init,
                                               const std::function<double(double)> &k,
                                               const std::function<double(double)> &c,
                                               const std::function<double(double, NumVec)> &Fux,
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

  //init U
  for (int i = 0; i < nx; i++) {
    U.push_back(init(i * dx));
    dU.push_back(0.0);
  }
}

void LinearParabolicProblem::run() {
  std::vector<std::string> params;
  params.push_back(std::to_string(x_n - x_0));
  params.push_back(std::to_string(nx));
  params.push_back(std::to_string(dx));
  params.push_back(std::to_string(timeControl.endTime));
  params.push_back(std::to_string(timeControl.endTime / timeControl.dt));
  params.push_back(std::to_string(timeControl.dt));
  params.push_back(std::to_string(k(1.0)));
  writeParams(file_name, params);

  while (timeControl.time < timeControl.endTime) {
    advance();
    std::cout << U;
  }
  std::cout << timeControl;
}

void LinearParabolicProblem::advance() {
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

    if (timeControl.stepsSinceRejection >= 50 && timeControl.dt == timeControl.dtmin) {
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

NumVec LinearParabolicProblem::step(double t, double dt) {

  TriDiag DF, LH;
  TriDiag M = generateMassMatrixMidpoint(c, nx, dx, x_0);
  TriDiag S = generateStiffnessMatrixMidpoint(k, nx, dx, x_0);
  NumVec F;
  NumVec dU_1(nx), dU_2(nx), RH(nx);

  F = linearizeF(U, Fux, dx, t, dt);
  //DF = linearizeDF(U, F, dx, dx, t, dt);
  LH = (1.0 / dt) * M + S; // -DF
  RH = (-1 * S) * U + F; // + F

  if (bc == dirchlet) {
    LH(0, 0) = 1;
    LH(0, 1) = 0;
    LH(nx - 1, nx - 1) = 1;
    LH(nx - 1, nx - 2) = 0;
    RH[0] = 0;
    RH[nx - 1] = 0;
  } else if (bc == vonNeumann) {
    RH[0] += k(x_0) * dU[0];
    RH[nx - 1] += -1 * k(x_n) * dU[nx - 1];
  }

  return solveTriDiagMatrix(LH, RH);
}
