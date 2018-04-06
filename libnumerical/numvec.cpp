//
// Created by Todd Dupont
//
#include <algorithm>
#include <random>
#include "numvec.h"

NumVec operator+(const NumVec& A, const NumVec& B) { // A+B
  if( A.size() != B.size() )
    error("NumVec +: incompatible sizes");
  NumVec C(A);
  for(int i=0; i<A.size(); i++)
    C[i] += B[i];
  return C;
}

NumVec operator*(double a, const NumVec& B) { // a*B
  NumVec C(B);
  for(int i=0; i<B.size(); i++ )
    C[i] *= a;
  return C;
}

NumVec operator-(const NumVec& A, const NumVec& B) { // A-B
  if( A.size() != B.size() )
    error("NumVec -: incompatible sizes");
  NumVec C(A);
  for(int i=0; i<A.size(); i++)
    C[i] -= B[i];
  return C;
}

double operator,(const NumVec& A, const NumVec& B) { // (A,B)
  if( A.size() != B.size() )
    error("NumVec ,: incompatible sizes");
  double sum = 0.0;
  for(int i=0; i< A.size(); i++ )
    sum += A[i]*B[i];
  return sum;
}

std::ostream& operator<<(std::ostream& os, const NumVec& A) {
   for(int i=0; i<A.size(); i++){
     os << A[i] << "\t";
   }
 os << "\n";
 return os;
}

NumVec randomVector(int length, int start, int end) {
  std::random_device rnd_device;
  std::mt19937 mersenne_engine(rnd_device());
  std::uniform_int_distribution<int> distr(start, end);
  auto gen = std::bind(distr, mersenne_engine);
  NumVec r(length);
  generate(r.begin(), r.end(), gen);
  return r;
}

NumVec randomRealVector(int length, double start, double end) {
  std::random_device rnd_device;
  std::mt19937 mersenne_engine(rnd_device());
  std::uniform_real_distribution<double> distr(start, end);
  auto gen = std::bind(distr, mersenne_engine);
  NumVec r(length);
  generate(r.begin(), r.end(), gen);
  return r;
}