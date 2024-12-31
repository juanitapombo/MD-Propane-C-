#ifndef FORCES2_H
#define FORCES2_H

#include <iostream>
#include <vector>
#include <cmath>

void forces2(int NATOM, double DENS, double LJCUT, double SIGMA, double EPSILON,
             std::vector<double>& RX, std::vector<double>& RY, std::vector<double>& RZ,
             double& LSIMBOX, double& LJCUTSQ, std::vector<double>& FX, std::vector<double>& FY, std::vector<double>& FZ,
             double& EPOT, double& PRESS); 

#endif
