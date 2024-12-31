#include <iostream>
#include <vector>
#include <cmath>

#include "forces2.h"
#include "vverlet.h"

int main() {
    int NATOM = 3;
    double DENS = 0.05;
    double LJCUT = 2.0;
    double SIGMA = 1.0;
    double EPSILON = 1.0;
    double kT = 100; 
    double MASS = 12; 
    bool RESCALE = true;
    double DT = 0.005;

    std::vector<double> RX;
    RX.push_back(1.0);
    RX.push_back(5.0);
    RX.push_back(3.0);

    std::vector<double> RY;
    RY.push_back(6.0);
    RY.push_back(1.0);
    RY.push_back(3.0);

    std::vector<double> RZ;
    RZ.push_back(2.0);
    RZ.push_back(5.0);
    RZ.push_back(3.0);

    std::vector<double> FX;
    FX.push_back(0.0927);
    FX.push_back(-0.0927);
    FX.push_back(0);

    std::vector<double> FY;
    FY.push_back(-1.1815);
    FY.push_back(1.1815);
    FY.push_back(0);

    std::vector<double> FZ;
    FZ.push_back(-0.9961);
    FZ.push_back(0.9961);
    FZ.push_back(0);

    std::vector<double> VX;
    VX.push_back(0.10);
    VX.push_back(0.20);
    VX.push_back(0.22);

    std::vector<double> VY;
    VY.push_back(0.40);
    VY.push_back(0.27);
    VY.push_back(0.24);

    std::vector<double> VZ;
    VZ.push_back(0.35);
    VZ.push_back(0.15);
    VZ.push_back(0.31);

    double LSIMBOX, EPOT, ETOT, EKIN, PRESS;  
   
    vverlet(RX, RY, RZ, VX, VY, VZ, FX, FY, FZ, NATOM, DENS, LJCUT, SIGMA, EPSILON, MASS, DT, kT, RESCALE, LSIMBOX, EPOT, PRESS, EKIN, ETOT);

    // Print the variables
    std::cout << "NATOM: " << NATOM << std::endl;
    std::cout << "DENS: " << DENS << std::endl;
    std::cout << "LJCUT: " << LJCUT << std::endl;
    std::cout << "LSIMBOX: " << LSIMBOX << std::endl;
    std::cout << "EKIN: " << EKIN << std::endl;
    std::cout << "EPOT: " << EPOT << std::endl;
    std::cout << "ETOT: " << ETOT << std::endl;
    std::cout << "PRESS: " << PRESS << std::endl;

    // Print FX, FY, and FZ
    std::cout << "FX: ";
    for (int i = 0; i < NATOM; ++i) {
        std::cout << FX[i] << " ";
    }
    std::cout << std::endl;
    
    std::cout << "FY: ";
    for (int i = 0; i < NATOM; ++i) {
        std::cout << FY[i] << " ";
    }
    std::cout << std::endl;
    
    std::cout << "FZ: ";
    for (int i = 0; i < NATOM; ++i) {
        std::cout << FZ[i] << " ";
    }
    std::cout << std::endl;

    // Print VX, VY, and VZ
    std::cout << "VX: ";
    for (int i = 0; i < NATOM; ++i) {
        std::cout << VX[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "VY: ";
    for (int i = 0; i < NATOM; ++i) {
        std::cout << VY[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "VZ: ";
    for (int i = 0; i < NATOM; ++i) {
        std::cout << VZ[i] << " ";
    }
    std::cout << std::endl;

    // Print RX, RY, and RZ
    std::cout << "RX: ";
    for (int i = 0; i < NATOM; ++i) {
        std::cout << RX[i] << " ";
    }
    std::cout << std::endl;
    
    std::cout << "RY: ";
    for (int i = 0; i < NATOM; ++i) {
        std::cout << RY[i] << " ";
    }
    std::cout << std::endl;
    
    std::cout << "RZ: ";
    for (int i = 0; i < NATOM; ++i) {
        std::cout << RZ[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}


