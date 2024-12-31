#include <iostream>
#include <vector>
#include <cmath>

// Define the forces function signature
void forces2(int NATOM, double DENS, double LJCUT, double SIGMA, double EPSILON,
             std::vector<double>& RX, std::vector<double>& RY, std::vector<double>& RZ,
             double& LSIMBOX, double& LJCUTSQ, std::vector<double>& FX, std::vector<double>& FY, std::vector<double>& FZ,
             double& EPOT, double& PRESS); 

void vverlet(std::vector<double>& RX, std::vector<double>& RY, std::vector<double>& RZ,
	     std::vector<double>& VX, std::vector<double>& VY, std::vector<double>& VZ,
 	     std::vector<double>& FX, std::vector<double>& FY, std::vector<double>& FZ,
	     int NATOM, double DENS, double LJCUT, double SIGMA, double EPSILON,
	     double MASS, double DT, double kT, bool RESCALE, double& LSIMBOX, 
             double& EPOT, double& PRESS, double& EKIN, double& ETOT) {

    // CALCULATE SIMULATION BOX LENGTH AND SQUARE OF LJCUT
    LSIMBOX = pow(NATOM / DENS, 1.0 / 3.0);
    
    // CALCULATE HALF TIMESTEP
    double DT2 = DT/2.0;

    // ADVANCE POSITIONS FROM T TO T+DT AND VELOCITIES FROM T TO T+DT/2 
    // ACCORDING TO THE VELOCITY VERLET INTEGRATOR:
    //  V(T+DT/2) = V(T) + 1/2*DT*A(T)
    //  R(T+DT)   = R(T) + DT*V(T) + 1/2*(DT**2)*A(T)
    //  V(T+DT)   = V(T+DT/2) + 1/2*DT*A(T+DT)


    // ADVANCE POSITIONS BY DT AND VELOCITIES BY DT/2
    for (int I = 0; I < NATOM; ++I) {
        VX[I] = VX[I] + DT2 * FX[I] / MASS;
        VY[I] = VY[I] + DT2 * FY[I] / MASS;
        VZ[I] = VZ[I] + DT2 * FZ[I] / MASS;
        RX[I] = RX[I] + DT * VX[I];
        RY[I] = RY[I] + DT * VY[I];
        RZ[I] = RZ[I] + DT * VZ[I];
    }

    // CALCULATE FORCES
    double LJCUTSQ;
    forces2(NATOM, DENS, LJCUT, SIGMA, EPSILON, RX, RY, RZ, LSIMBOX, LJCUTSQ, FX, FY, FZ, EPOT, PRESS);
 
    // ADVANCE VELOCITIES BY THE REMAINING DT/2, WITH NEW FORCES
    for (int I = 0; I < NATOM; ++I) {
	VX[I] = VX[I] + DT2 * FX[I] / MASS;
	VY[I] = VY[I] + DT2 * FY[I] / MASS;
	VZ[I] = VZ[I] + DT2 * FZ[I] / MASS;
    }

    // APPLY PERIODIC BOUNDARY CONDITIONS
    for (int I = 0; I < NATOM; ++I) {
        RX[I] = RX[I] - LSIMBOX * std::round(RX[I]/LSIMBOX);
        RY[I] = RY[I] - LSIMBOX * std::round(RY[I]/LSIMBOX);
        RZ[I] = RZ[I] - LSIMBOX * std::round(RZ[I]/LSIMBOX);
    }

    // RESCALE THE VELOCITIES TO THE OBTAIN THE SETPOINT TEMPERATURE. DO 
    // THIS ONLY DURING THE EQUILIBRATION
    if (RESCALE) {
        double VELSQ = 0.0; 
        for (int I = 0; I < NATOM; ++I) {
            VELSQ = VELSQ + MASS * (pow(VX[I], 2) + pow(VY[I], 2) + pow(VZ[I], 2));
        }
	double FACTOR = pow( 3 * NATOM * kT / VELSQ, 1/2); 
        for (int I = 0; I < NATOM; ++I) {
            VX[I] = VX[I] * FACTOR; 
            VY[I] = VY[I] * FACTOR; 
            VZ[I] = VZ[I] * FACTOR; 
        } 
    }

   // CALCULATE CURRENT KINETIC ENERGY 
   EKIN = 0.0;  
   for (int I = 0; I < NATOM; ++I) {
       EKIN = EKIN + pow(VX[I], 2) + pow(VY[I], 2) + pow(VZ[I], 2);
   }
   EKIN = EKIN * MASS / 2; 

   // CALCULATE TOTAL ENERGY AS KINETIC ENERGY + POTENTIAL ENERGY
   ETOT = EKIN + EPOT;

   // ADD CURRENT KINETIC (IDEAL GAS) CONTRIBUTION TO CURRENT PRESSURE
   PRESS = PRESS + 2.0 * EKIN * DENS / 3.0 / NATOM;

   }

