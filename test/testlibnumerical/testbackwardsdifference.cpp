//
// Created by Austin Clyde on 4/16/18.
//

#include "tridiagtest.h"

TEST_F(BackwardsDifferenceTest, basic_heat_eq
) {
std::vector<double> a = {1, 2, 3, 4};
std::vector<double> d = {1, 2, 3, 4, 5};
std::vector<double> b = {1, 2, 3, 4};
std::vector<double> r = {1, 2, 3, 4, 5};
std::vector<double> ans = {-8.0 / 49, 57.0 / 49, -4.0 / 49, 15.0 / 49, 37.0 / 49};
TriDiag M(5);
M.
a = a;
M.
b = b;
M.
d = d;
assert(L2norm(solveTriDiagMatrix(M, r) - ans) < 0.001);
}