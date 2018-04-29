//
// Created by Todd Dupont
//

#ifndef CS28515PROJ1_SIMTIME_H
#define CS28515PROJ1_SIMTIME_H

#include <ostream>
class simTime { // class for time step control
 public:
  double time;
  double dt;
  double dtOld;
  double tol;
  double agrow;
  double ashrink;
  double dtmin;
  double dtmax;
  double endTime;
  int stepsSinceRejection;
  int stepsRejected;
  int stepsAccepted;
  bool newSimulation;

  simTime() : time(0.0), dt(0.001), dtOld(0.001), tol(1.0e-4), agrow(1.2),
              ashrink(1.0 / 1.2), dtmin(1.0e-6), dtmax(1e-1), endTime(10.0),
              stepsSinceRejection(0), stepsRejected(0), stepsAccepted(true), newSimulation(true) {}

  friend std::ostream &operator<<(std::ostream &os, const simTime &time1);
};

#endif //CS28515PROJ1_SIMTIME_H
