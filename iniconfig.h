#ifndef INICONFIG_H
#define INICONFIG_H

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <cstdio>

void iniconfig(std::vector<double>& RX, std::vector<double>& RY, std::vector<double>& RZ,
               std::vector<double>& VX, std::vector<double>& VY, std::vector<double>& VZ,
               std::vector<double>& FX, std::vector<double>& FY, std::vector<double>& FZ,
               int NATOM, double DENS, double kT, double SIGMA, double EPSILON,
               double MASS, double LJCUT); 

#endif
