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

void solveHeatEquation1d(double x_0,
                         double x_nx,
                         int nx,
                         int nt,
                         double tmax,
                         double alpha,
                         std::function<double(double)> init) {
  double L = x_nx - x_0;
  double dx = L / (nx - 1);
  double dt = tmax / (nt - 1);

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

  NumVec U(nx);
  NumVec UOld(nx);
  NumVec actual(nx);

  //init conditions
  for (int i = 1; i < nx - 1; i++)
    U[i] = init(i * dx);
  U[nx - 1] = 0;
  U[0] = 0;

  double t = 0;
  for (int m = 1; m < nt; m++) {
    UOld = U;
    UOld = (1 / dt) * UOld;
    t = t + dt;
    U = solveTriDiagMatrix(A, UOld);
    writeUpdateStep("test.txt", U, t);
  }
  std::cout << "error x= " << L2norm(actual - U) << std::endl;
}

void solveMassStiffStepDouble(std::function<double(double)> k, std::function<double(double)> c, double x_0,
                              double x_nx,
                              int nx,
                              int nt,
                              double tmax,
                              std::function<double(double)> init, bool debug, bool dirch) {
  double L = x_nx - x_0;
  double dx = L / (nx - 1);
  double dt = tmax / (nt);
  double t = 0;

  NumVec U(nx);
  NumVec dU(nx), dU_1(nx), dU_2(nx);

  //init conditions
  for (int i = 0; i < nx; i++)
    U[i] = init(i * dx);

  std::vector<std::string> params;
  params.push_back(std::to_string(L));
  params.push_back(std::to_string(nx));
  params.push_back(std::to_string(dx));
  params.push_back(std::to_string(tmax));
  params.push_back(std::to_string(nt));
  params.push_back(std::to_string(dt));
  params.push_back(std::to_string(k(1.0)));
  writeParams("test.txt", params);


  NumVec RH(nx);
  for (int m = 1; m <= nt; m++) {
    dU_1 = step(k, c, m * dt, x_0, x_nx, L, dx, dt, nx, U, dU, debug, dirch);
    dU_2 = step(k, c, m * dt, x_0, x_nx, L, dx, dt / 2, nx, U, dU, debug, dirch);
    dU_2 = dU_2 + step(k, c, m * dt + dt / 2, x_0, x_nx, L, dx, dt / 2, nx, U, dU_2, debug, dirch);
    U = U + dU_1;
    writeUpdateStep("test.txt", U, t);
    std::cout << "error " << L2norm(dU_1 - dU_2) << std::endl;
  }
  std::cout << "U at t = " << t << ": " << U;
}

NumVec step(std::function<double(double)> k, std::function<double(double)> c,
    //  std::function<double(double, std::function<double(double, double)>)> F,
            double t,
            double x_0,
            double x_nx,
            double L,
            double dx,
            double dt,
            int nx,
            NumVec U, NumVec dU, bool debug, bool dirch) {

  TriDiag DF, LH;
  TriDiag M = generateMassMatrixMidpoint(c, nx, dx, x_0);
  TriDiag S = generateStiffnessMatrixMidpoint(k, nx, dx, x_0);
  NumVec F;
  NumVec dU_1(nx), dU_2(nx), RH(nx);

  if (debug) {
    std::cout << "L = " << L << " nx = " << nx << " dx = " << dx << " dt = " << dt << std::endl;
    std::cout << "initial vector for  U: " << U << std::endl;
    std::cout << "Mass: \n" << M << "Stiff: \n" << S << std::endl;
  }

  t += dt;

//  F = linearizeF(U, dx, dx, t, dt);
//  DF = linearizeDF(U, F, dx, dx, t, dt);
  LH = (1.0 / dt) * M + S - DF;
  RH = ((-1 * S) * U + F);

  if (dirch) {
    LH(0, 0) = 1;
    LH(0, 1) = 0;
    LH(nx - 1, nx - 1) = 1;
    LH(nx - 1, nx - 2) = 0;
    RH[0] = 0;
    RH[nx - 1] = 0;
  } else { // von
    RH[0] += k(0) * dU[0];
    RH[nx - 1] += -1 * k(L) * dU[nx - 1];
  }

  dU_1 = solveTriDiagMatrix(LH, RH);

  if (debug) {
    std::cout << "t = " << t << std::endl;
    std::cout << "left side of equation\n" << LH;
    std::cout << "right side of equation: " << RH;
    std::cout << "dU step " << t << ": " << dU;
    std::cout << "U step " << t << ": " << U << std::endl << std::endl;
  }
  return dU_1;
}

void solveMassStiff(std::function<double(double)> k, std::function<double(double)> c, double x_0,
                    double x_nx,
                    int nx,
                    int nt,
                    double tmax,
                    std::function<double(double)> init) {
  double L = x_nx - x_0;
  double dx = L / (nx - 1);
  double dt = tmax / (nt);
  double t = 0;
  bool debug = false;
  bool dirch = false;
  TriDiag M = generateMassMatrixMidpoint(c, nx, dx, x_0);
  TriDiag S = generateStiffnessMatrixMidpoint(k, nx, dx, x_0);
  TriDiag DF;
  NumVec F;
  NumVec U(nx);
  NumVec dU(nx);

  //init conditions
  for (int i = 0; i < nx; i++)
    U[i] = init(i * dx);

  std::vector<std::string> params;
  params.push_back(std::to_string(L));
  params.push_back(std::to_string(nx));
  params.push_back(std::to_string(dx));
  params.push_back(std::to_string(tmax));
  params.push_back(std::to_string(nt));
  params.push_back(std::to_string(dt));
  params.push_back(std::to_string(k(1.0)));
  writeParams("test.txt", params);

  if (debug) {
    std::cout << "L = " << L << " nx = " << nx << " dx = " << dx << " dt = " << dt << std::endl;
    std::cout << "initial vector for  U: " << U << std::endl;
    std::cout << "Mass: \n" << M << "Stiff: \n" << S << std::endl;
  }
  TriDiag LH = (1.0 / dt) * M + S;

  NumVec RH(nx);
  for (int m = 1; m <= nt; m++) {
    writeUpdateStep("test.txt", U, t);
    t += dt;

//    F = linearizeF(U, m * dx, dx, t, dt);
    //  DF = linearizeDF(U, F, m * dx, dx, t, dt);
    LH = (1.0 / dt) * M + S - DF;
    RH = ((-1 * S) * U + F);

    if (dirch) {
      LH(0, 0) = 1;
      LH(0, 1) = 0;
      LH(nx - 1, nx - 1) = 1;
      LH(nx - 1, nx - 2) = 0;
      RH[0] = 0;
      RH[nx - 1] = 0;
    } else { // von
      RH[0] += k(0) * dU[0];
      RH[nx - 1] += -1 * k(L) * dU[nx - 1];
    }

    dU = solveTriDiagMatrix(LH, RH);
    U = U + dU;
    if (debug) {
      std::cout << "t = " << t << std::endl;
      std::cout << "left side of equation\n" << LH;
      std::cout << "right side of equation: " << RH;
      std::cout << "dU step " << m << ": " << dU;
      std::cout << "U step " << m << ": " << U << std::endl << std::endl;
    }
  }
  std::cout << "U at t = " << t << ": " << U;
}

NumVec linearizeF(NumVec U, std::function<double(double, NumVec)> F, double dx, double t_i, double dt) {
  NumVec f(U.size());
  f[0] = -1 * F(dx * 0.5, U);
  for (int i = 1; i < U.size() - 1; i++) {
    f[i] += F((i - 1) * dx + dx * 0.5, U) - F((i) * dx + 0.5 * dx, U);
  }
  f[f.size() - 1] = F((U.size() - 2) * dx + 0.5 * dx, U);
  return -1 * f;
}

TriDiag linearizeDF(NumVec U, NumVec F, double x_i, double dx, double t_i, double dt) {
  TriDiag DF(U.size());

  return DF;
}

void solveHeatEquation1dStepDoubling(double x_0,
                                     double x_nx,
                                     int nx,
                                     int nt,
                                     double tmax,
                                     double alpha,
                                     std::function<double(double)> init) {
  double L = x_nx - x_0;
  double dx = L / (nx - 1);
  double dt = tmax / (nt - 1);
  double acc = pow(10, -3);

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

  NumVec U(nx);
  NumVec UOld(nx);
  NumVec actual(nx);

  //init conditions
  for (int i = 1; i < nx - 1; i++)
    U[i] = init(i * dx);
  U[nx - 1] = 0;
  U[0] = 0;
  std::cout << A;
  double t = 0;
  for (int m = 1; m < nt; m++) {
    UOld = U;
    UOld = (1 / dt) * UOld;
    t = t + dt;
    U = solveTriDiagMatrix(A, UOld);
    writeUpdateStep("test.txt", U, t);
  }
  std::cout << "old: " << U;
  std::cout << "error x= " << L2norm(actual - U) << std::endl;
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
