#include <iostream>
#include <vector>
#include <cmath>

#include "forces2.h"

int main() {
    int NATOM = 3;
    double DENS = 0.05;
    double LJCUT = 2.0;
    double SIGMA = 1.0;
    double EPSILON = 1.0;

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

    double LSIMBOX, LJCUTSQ, EPOT, PRESS;
    std::vector<double> FX(NATOM, 0.0);
    std::vector<double> FY(NATOM, 0.0);
    std::vector<double> FZ(NATOM, 0.0);

    forces2(NATOM, DENS, LJCUT, SIGMA, EPSILON, RX, RY, RZ, LSIMBOX, LJCUTSQ, FX, FY, FZ, EPOT, PRESS);

    // Print the variables
    std::cout << "NATOM: " << NATOM << std::endl;
    std::cout << "DENS: " << DENS << std::endl;
    std::cout << "LJCUT: " << LJCUT << std::endl;
    std::cout << "LSIMBOX: " << LSIMBOX << std::endl;
    std::cout << "LJCUTSQ: " << LJCUTSQ << std::endl;
    std::cout << "EPOT: " << EPOT << std::endl;
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

    return 0;
}


