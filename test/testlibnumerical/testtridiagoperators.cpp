#include "tridiagtest.h"

TEST_F(TriDiagTest, MatrixEqual) {
  m1 = TriDiag(3);
  m2 = TriDiag(3);
  m1(1,1) = 2;
  m2(1,1) = 2;
  EXPECT_TRUE(m1 == m1);
}

TEST_F(TriDiagTest, MatrixEqualId) {
  m1 = TriDiag(3);
  EXPECT_TRUE(m1 == m1);
}

TEST_F(TriDiagTest, MatrixNotEqualDim) {
  m1 = TriDiag(3);
  m2 = TriDiag(4);
  m1(1,1) = 2;
  m2(1,1) = 2;
  EXPECT_FALSE(m1 == m2);
}

TEST_F(TriDiagTest, MatrixNotEqualVal) {
  m1 = TriDiag(3);
  m2 = TriDiag(3);
  m1(2,2) = -1;
  m1(1,1) = -1;
  EXPECT_FALSE(m1 == m2);
}

TEST_F(TriDiagTest, MatrixAdd) {
  m1 = TriDiag(3, 1);
  m2 = TriDiag(3, 1);
  result = TriDiag(3,2);
  EXPECT_TRUE((m1+m2) == result);
}