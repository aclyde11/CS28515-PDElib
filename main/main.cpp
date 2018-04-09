#include <iostream>
#include "tridiag.h"

int main() {
  /*
  TriDiag A(5);
  for(int i = 0; i<4; i++) {
    A(i,i) = i+1;
    A(i+1,i) = i*0.25;
    A(i,i+1) = i*2;
  }
  A(4,4) = 5;
  std::cout << A << std::endl;


  NumVec d(5);
  for(int i = 0; i < 5; i++) {
    d[i] = i;
  }
  std::cout << d << std::endl;

  NumVec x = solveTriDiagMatrix(A, d);
  std::cout << x << std::endl;
  std::cout << (A*x) << std::endl;

  TriDiag test(3);
  test.randomize(1, 10);
  std::cout << test << std::endl;

  std::cout << test.det() << std::endl;
   */

  TriDiag A(5);
  TriDiag B(5);
  A.randomize(1, 100);
  B.randomize(1, 100);
  NumVec v = randomVector(5, 0, 100);
  std::cout << A << std::endl;
  std::cout << v << std::endl;
  std::cout << A*v << std::endl;
  std::cout << A << std::endl;
  std::cout << v << std::endl;

  std::cout << A*B << std::endl;
}