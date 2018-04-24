//
// Created by Austin Clyde on 4/16/18.
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
  simTime() : time(0.0), dt(0.001), dtOld(0.001), tol(1.0e-3), agrow(1.25),
              ashrink(1.0 / 1.25), dtmin(1.0e-7), dtmax(1.0), endTime(1.0),
              stepsSinceRejection(0), stepsRejected(0), stepsAccepted(true), newSimulation(true) {}
  friend std::ostream &operator<<(std::ostream &os, const simTime &time1) {
    os << "time: " << time1.time << " dt: " << time1.dt << " dtOld: " << time1.dtOld << " tol: " << time1.tol
       << " agrow: " << time1.agrow << " ashrink: " << time1.ashrink << " dtmin: " << time1.dtmin << " dtmax: "
       << time1.dtmax << " endTime: " << time1.endTime << " stepsSinceRejection: " << time1.stepsSinceRejection
       << " stepsRejected: " << time1.stepsRejected << " stepsAccepted: " << time1.stepsAccepted << " newSimulation: "
       << time1.newSimulation;
    return os;
  };
};

#endif //CS28515PROJ1_SIMTIME_H
