#ifndef CS28515PROJ1_TRIDIAGTEST_H
#define CS28515PROJ1_TRIDIAGTEST_H

#include <stdio.h>
#include <iostream>
#include "gtest/gtest.h"
#include "tridiag.h"


class TriDiagTest : public ::testing::Test {
 public:
  TriDiag m1;
  TriDiag m2;
  TriDiag result;
  TriDiag testMatrix;
};

#endif