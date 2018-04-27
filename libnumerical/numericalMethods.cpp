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


