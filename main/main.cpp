#include <iostream>
#include "tridiag.h"
#include "pdesolver.h"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <vector>
#include <functional>
#include <math.h>

#include <string>

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


}

void btcs() {
  int L = 1;
  double tmax = 2;
  int nx = 10;
  int nt = 200;
  double alpha = 0.2;
  double dx = ((double) L)/(nx-1);
  double dt = ((double) tmax)/(nt-1);
  NumVec actual(nx);

  TriDiag A(nx);
  for (int i = 1; i < nx -1; i++) { //diag
    A(i,i) = (1/dt) + (2*alpha / (dx * dx));
  }
  A(0,0 ) =1;
  A(nx-1,nx-1) =1;
  for(int i = 1; i < nx - 1; i++) { //above
    A(i,i+1) = -1 * alpha / (dx * dx);
  }
  for(int i = 0; i < nx-2; i++) { //below
    A(i+1,i) = -1 * alpha / (dx * dx);
  }
  std::cout<< A;

  NumVec U(nx);
  NumVec UOld(nx);
  for(int i = 0; i < nx-1; i++)
    U[i] = sin(PI * i*dx / L);
  U[dx-1] =0;
  U[0] = 0;
  std::cout << U;
  double t = 0;
  for(int m = 1; m < nt; m++) {
    UOld = U;
    UOld = (1/dt) * UOld;
    t = t + dt;
    U = solveTriDiagMatrix(A, UOld);
    for(int i = 1; i < nx-1; i++) {
      actual[i] = sin(PI * i * dx ) * exp(-1*alpha * PI * PI *t);
    }
    //std::cout << "error = " << L2norm(actual - U) << std::endl;
    //std::cout << U;
  }
  std::cout << "error x= " << L2norm(actual - U) << std::endl;
  std::cout << U;
  writeVector(U);
}

void actual() {

}
*/
int main() {
  //ftcs();
  // btcs();
  double x_0 = 0;
  double x_nx = 1;
  int nx = 1000;
  int nt = 4000;
  double tmax = 4;
  double alpha = 0.2;
  solveHeatEquation1d(x_0, x_nx, nx, nt, tmax, alpha);
  return 0;
}