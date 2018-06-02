//
// Created by Austin Clyde on 5/23/18.
//

#ifndef CS28515PROJ1_TEST_H
#define CS28515PROJ1_TEST_H

#include <iostream>
#include "Eigen/Core"

#include "test.h"
#include "VariMesh.h"
#include "FEMLaplacian.h"
#include <cmath>

void example0();

void example1();

double example0_init(int i, int j, VariMesh mesh);

double example1_init(int i, int j, VariMesh mesh);

#endif //CS28515PROJ1_TEST_H
