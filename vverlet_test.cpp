#include <iostream>
#include <vector>
#include <cmath>

// Define the forces function signature
void forces2(int NATOM, double DENS, double LJCUT, double SIGMA, double EPSILON,
             std::vector<double>& RX, std::vector<double>& RY, std::vector<double>& RZ,
             double& LSIMBOX, double& LJCUTSQ, std::vector<double>& FX, std::vector<double>& FY, std::vector<double>& FZ,
             double& EPOT, double& PRESS);

void test_vverlet() {
    // define testing parameters
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

    // test the functions here
    double LSIMBOX = pow(NATOM / DENS, 1.0 / 3.0);
    double DT2 = DT/2.0;

    for (int i = 0; i < NATOM; ++i) {
        printf("FX[%d]: %f, FY[%d]: %f, FZ[%d]: %f\n", i, FX[i], i, FY[i], i, FZ[i]);
    }

    printf("LSIMBOX: %f, DT2: %f\n", LSIMBOX, DT2);

    // ADVANCE POSITIONS BY DT AND VELOCITIES BY DT/2
    for (int I = 0; I < NATOM; ++I) {
        VX[I] = VX[I] + DT2 * FX[I] / MASS;
        VY[I] = VY[I] + DT2 * FY[I] / MASS;
        VZ[I] = VZ[I] + DT2 * FZ[I] / MASS;
        RX[I] = RX[I] + DT * VX[I];
        RY[I] = RY[I] + DT * VY[I];
        RZ[I] = RZ[I] + DT * VZ[I];
        printf("VX[%d]: %f, VY[%d]: %f, VZ[%d]: %f, RX[%d]: %f, RY[%d]: %f, RZ[%d]: %f\n",
               I, VX[I], I, VY[I], I, VZ[I], I, RX[I], I, RY[I], I, RZ[I]);
    }

    double EPOT, PRESS, LJCUTSQ;
    // CALCULATE FORCES
    forces2(NATOM, DENS, LJCUT, SIGMA, EPSILON, RX, RY, RZ, LSIMBOX, LJCUTSQ, FX, FY, FZ, EPOT, PRESS);

    for (int i = 0; i < NATOM; ++i) {
        printf("FX[%d]: %f, FY[%d]: %f, FZ[%d]: %f\n", i, FX[i], i, FY[i], i, FZ[i]);
    }

    printf("EPOT: %f, PRESS: %f\n", EPOT, PRESS);
 
    // ADVANCE VELOCITIES BY THE REMAINING DT/2, WITH NEW FORCES
    for (int I = 0; I < NATOM; ++I) {
        VX[I] = VX[I] + DT2 * FX[I] / MASS;
        VY[I] = VY[I] + DT2 * FY[I] / MASS;
        VZ[I] = VZ[I] + DT2 * FZ[I] / MASS;

        printf("VX[%d]: %f, VY[%d]: %f, VZ[%d]: %f\n",I, VX[I], I, VY[I], I, VZ[I]); 
    }

    // APPLY PERIODIC BOUNDARY CONDITIONS
    for (int I = 0; I < NATOM; ++I) {
        RX[I] = RX[I] - LSIMBOX * std::round(RX[I]/LSIMBOX);
        RY[I] = RY[I] - LSIMBOX * std::round(RY[I]/LSIMBOX);
        RZ[I] = RZ[I] - LSIMBOX * std::round(RZ[I]/LSIMBOX);
        
        printf("RX[%d]: %f, RY[%d]: %f, RZ[%d]: %f\n",I, RX[I], I, RY[I], I, RZ[I]);

    }

    // RESCALE THE VELOCITIES TO THE OBTAIN THE SETPOINT TEMPERATURE. DO
    // THIS ONLY DURING THE EQUILIBRATION
    if (RESCALE) {
        double VELSQ = 0.0;
        for (int I = 0; I < NATOM; ++I) {
            VELSQ = VELSQ + MASS * (pow(VX[I], 2) + pow(VY[I], 2) + pow(VZ[I], 2));
            printf("VELSQ[%d]: %f\n",I, VELSQ);
        }
        double FACTOR = sqrt( 3 * NATOM * kT / VELSQ );
        for (int I = 0; I < NATOM; ++I) {
            VX[I] = VX[I] * FACTOR;
            VY[I] = VY[I] * FACTOR;
            VZ[I] = VZ[I] * FACTOR;
	    printf("VX[%d]: %f, VY[%d]: %f, VZ[%d]: %f\n",I, VX[I], I, VY[I], I, VZ[I]);
        }
    }

   // CALCULATE CURRENT KINETIC ENERGY
   double EKIN = 0.0;
   for (int I = 0; I < NATOM; ++I) {
       EKIN = EKIN + pow(VX[I], 2) + pow(VY[I], 2) + pow(VZ[I], 2);
       printf("EKIN[%d]: %f\n",I, EKIN);
   }
   EKIN = EKIN * MASS / 2;
   printf("EKIN: %f\n", EKIN); 

   // CALCULATE TOTAL ENERGY AS KINETIC ENERGY + POTENTIAL ENERGY
   double ETOT = EKIN + EPOT;
   printf("EKIN: %f, EPOT: %f, ETOT: %f\n", EKIN, EPOT, ETOT); 

   // ADD CURRENT KINETIC (IDEAL GAS) CONTRIBUTION TO CURRENT PRESSURE
   printf("PRESS: %f\n", PRESS); 
   PRESS = PRESS + 2.0 * EKIN * DENS / 3.0 / NATOM;
   printf("PRESS: %f\n", PRESS); 

}

#include "forces2.h"

int main() {
    test_vverlet();
    return 0;
}

