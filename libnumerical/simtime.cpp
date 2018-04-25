//
// Created by Austin Clyde on 4/16/18.
//

#include "simtime.h"

std::ostream &operator<<(std::ostream &os, const simTime &time1) {
  os << "time: " << time1.time << " dt: " << time1.dt << " dtOld: " << time1.dtOld << " tol: " << time1.tol
     << " agrow: " << time1.agrow << " ashrink: " << time1.ashrink << " dtmin: " << time1.dtmin << " dtmax: "
     << time1.dtmax << " endTime: " << time1.endTime << " stepsSinceRejection: " << time1.stepsSinceRejection
     << " stepsRejected: " << time1.stepsRejected << " stepsAccepted: " << time1.stepsAccepted << " newSimulation: "
     << time1.newSimulation;
  return os;
}
