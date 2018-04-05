#include <iostream>
#include "tridiag.h"

int main() {
  TriDiag A(3);
  A(1,1) = 2;
  std::cout << (10.0 * A)<< std::endl;

  TriDiag B(3);
  B(1,1) = 1;
  B(0,0) = 1;
  std::cout << A+B << std::endl;
}