#ifndef VVERLET_H
#define VVERLET_H

#include <iostream>
#include <vector>
#include <cmath>

void vverlet(std::vector<double>& RX, std::vector<double>& RY, std::vector<double>& RZ,
             std::vector<double>& VX, std::vector<double>& VY, std::vector<double>& VZ,
             std::vector<double>& FX, std::vector<double>& FY, std::vector<double>& FZ,
             int NATOM, double DENS, double LJCUT, double SIGMA, double EPSILON,
             double MASS, double DT, double kT, bool RESCALE, double& LSIMBOX,
             double& EPOT, double& PRESS, double& EKIN, double& ETOT); 


#endif

