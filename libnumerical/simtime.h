//
// Created by Austin Clyde on 4/16/18.
//

#ifndef CS28515PROJ1_SIMTIME_H
#define CS28515PROJ1_SIMTIME_H

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
  simTime() : time(0.0), dt(0.0), dtOld(0.0), tol(1.0e-3), agrow(1.25),
              ashrink(1.0 / 1.25), dtmin(1.0e-6), dtmax(1.0), endTime(1.0),
              stepsSinceRejection(0), stepsRejected(0), stepsAccepted(true), newSimulation(true) {};
};

#endif //CS28515PROJ1_SIMTIME_H
