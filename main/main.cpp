#include <iostream>
#include "tridiag.h"
#include "pdesolver.h"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <vector>
#include <functional>
#include <cmath>

#include <string>
#define PI 3.141592653
/*

TriDiag ftcs() {
  int L = 1;
  double tmax = 1;
  int nx = 10;
  int nt = 100;
  double alpha = 0.1;
  double dx = ((double) L)/(nx-1);
  double dt = ((double) tmax)/(nt-1);
  double r = alpha * dt / (dx * dx);
  double r2 = 1 - 2*r;
  std::cout << dx << ", " << dt << ", " << r <<std::endl;;
  NumVec U(nx);
  NumVec UOld(nx);
  UOld[0] = 0;
  UOld[nx-1] = 0;
  for(int i = 0; i < nx-1; i++)
    U[i] = sin(PI * i*dx / L);

  TriDiag A(nx);
  for(int i = 0; i < nx-2; i++) {
    A(i+1,i) = r;
  }
  for(int i = 0; i < nx; i++) {
    A(i,i) = r2;
  }
  for(int i = 1; i < nx-1; i++)
    A(i, i+1) = r;
  A(0,0) = 1;
  A(nx-1,nx-1) = 1;
  std::cout << A;
  std::cout << "init UOld" << std::endl << U << std::endl;
  NumVec actual(nx);
  double t = 0;
  for(int m = 1; m < nt; m++) {
    UOld = U;
    t = t+dt;
    U = A * UOld;
    for(int i = 1; i < nx-1; i++) {
      actual[i] = sin(PI * i * dx ) * exp(-1*alpha * PI * PI *t);
    }
    //std::cout << "t = " << t << std::endl << U << "actual" << std::endl << actual;
    //std::cout << "error = " << L2norm(actual - U) << std::endl << std::endl;
    //writeVector(U);
    //std::cout << "t = " << t << ", x = " ;
  }
  std::cout << "error L2= " << L2norm(actual - U) << std::endl;


} */

int main() {
  double x_0 = 0;
  double x_nx = 1;
  double L = (x_nx - x_0);
  int nx = 100;
  int nt = 1000;
  double tmax = 1;
  double alpha = 0.1;
  std::function<double(double)> init;
  init = [alpha, L](double x) -> double {
    return sin(PI * x / L);
  };
  solveHeatEquation1d(x_0, x_nx, nx, nt, tmax, alpha, init);
  return 0;
}