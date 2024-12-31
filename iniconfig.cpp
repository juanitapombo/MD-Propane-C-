#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <cstdio>

// Define the forces function signature
void forces2(int NATOM, double DENS, double LJCUT, double SIGMA, double EPSILON,
             std::vector<double>& RX, std::vector<double>& RY, std::vector<double>& RZ,
             double& LSIMBOX, double& LJCUTSQ, std::vector<double>& FX, std::vector<double>& FY, std::vector<double>& FZ,
             double& EPOT, double& PRESS);

// FUNCTION TO INITIALIZE POSITIONS AND VELOCITIES 

void iniconfig(std::vector<double>& RX, std::vector<double>& RY, std::vector<double>& RZ,
               std::vector<double>& VX, std::vector<double>& VY, std::vector<double>& VZ,
               std::vector<double>& FX, std::vector<double>& FY, std::vector<double>& FZ,
               int NATOM, double DENS, double kT, double SIGMA, double EPSILON, 
               double MASS, double LJCUT) { 

    // INITIALIZE POSITIONS!  

    // PREALLOCATE AND INITIALIZE VECTORS TO ZERO 
    RX.assign(NATOM, 0.0);
    RY.assign(NATOM, 0.0);
    RZ.assign(NATOM, 0.0);
    VX.assign(NATOM, 0.0);
    VY.assign(NATOM, 0.0);
    VZ.assign(NATOM, 0.0);

    // LENGTH OF SIMULATION BOX
    double LSIMBOX = pow(NATOM / DENS, 1.0 / 3.0);
    double LSIMBOX2 = LSIMBOX/2.0;

    // ASSIGN DISTANCE CUTOFF TO PREVENT OVERLAPS
    double SIGMACUT2 = 0.9 * SIGMA; 
    SIGMACUT2 = SIGMACUT2 * SIGMACUT2; 

    // SEED THE RANDOM NUMBER GENERATOR 
    std::random_device rd;
    std::mt19937 gen(rd());

    // DEFINE THE RANGE FOR THE RANDOM NUMBERS 
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // ASSIGN TEMPORARY POSITION OF FIRST ATOM 
    double POSX = dis(gen);
    double POSY = dis(gen);
    double POSZ = dis(gen);

    // TRANSFORM TEMPORARY POS (RANDOM NUMBER BETWEEN 0 AND 1)
    //  TO A LENGTH BETWEEN -1/2 * LSIMBOX and +1/2 * LSIMBOX: (-L/2,+L/2)
    RX[0] = LSIMBOX * POSX - LSIMBOX2; 
    RY[0] = LSIMBOX * POSY - LSIMBOX2; 
    RZ[0] = LSIMBOX * POSZ - LSIMBOX2; 

    // MAKE SURE TWO ATOMS FO NOT OVERLAP BY MORE THAN 10% OF SIGMA 
    for (int I = 1; I < NATOM; ++I) { 
        printf("INSERTION OF MOL %4i SUCCESSFULL \n", I);
        bool REPEAT = true;

        while (REPEAT) {
           REPEAT = false; 
           // ASSIGN THE POSITION OF THE I'TH ATOM WITHIN (-L/2,+L/2)
           POSX = LSIMBOX * dis(gen) - LSIMBOX2;
           POSY = LSIMBOX * dis(gen) - LSIMBOX2;
           POSZ = LSIMBOX * dis(gen) - LSIMBOX2;
    
           // CALCULATE DISTANCES FROM PREVIOUSLY INSERTED ATOMS 
           for (int K = 0; K < I; ++K) {
               double DRX = POSX - RX[K]; 
               double DRY = POSY - RY[K]; 
               double DRZ = POSZ - RZ[K];
               // MINIMUM IMAGE CONVENTION 
               DRX = DRX - LSIMBOX * std::round(DRX/LSIMBOX); 
               DRY = DRY - LSIMBOX * std::round(DRY/LSIMBOX); 
               DRZ = DRZ - LSIMBOX * std::round(DRZ/LSIMBOX); 
               double DELRSQ = (pow(DRX, 2) + pow(DRY, 2) + pow(DRZ, 2)); 
              
               // TRY A DIFFERENT RANDOM POSITION FOR THE CURRENT ATOM IF IT IS 
               // TO CLOSE TO ANY OF THE PREVIOUSLY INSERTED ATOMS 
               if (DELRSQ < SIGMACUT2) {
                   REPEAT = true;  
               }
           }
        }
 
        RX[I] = POSX; 
        RY[I] = POSY; 
        RZ[I] = POSZ; 
    }  

    // INITIALIZE VELOCITIES!

    // CALCULATE FORCES 

    double LJCUTSQ, EPOT, PRESS;
    forces2(NATOM, DENS, LJCUT, SIGMA, EPSILON, RX, RY, RZ, LSIMBOX, LJCUTSQ, FX, FY, FZ, EPOT, PRESS);

    // RESCALE VELOCITIES TO (-1,1) & CALCULATE SUM OF X, Y, Z VELOCITY 
    double SUMX = 0.0;   
    double SUMY = 0.0;   
    double SUMZ = 0.0;  
    double VELSQ = 0.0; 
    for (int I = 0; I < NATOM; ++I) { 
        VX[I] = 2.0 * dis(gen) - 1.0; 
        VY[I] = 2.0 * dis(gen) - 1.0; 
        VZ[I] = 2.0 * dis(gen) - 1.0; 
        double SUMX = SUMX + VX[I];  
        double SUMY = SUMY + VY[I];  
        double SUMZ = SUMZ + VZ[I];  
    } 

    // ZERO THE TOTAL LINEAR MOMENTUM BY REMOVING/ADDING MOMENTUM
    double NATOMI = 1.0/NATOM;  
    for (int I = 0; I < NATOM; ++I) {
        VX[I] = VX[I] - SUMX * NATOMI;
        VY[I] = VY[I] - SUMY * NATOMI;
        VZ[I] = VZ[I] - SUMZ * NATOMI;
        VELSQ = VELSQ + MASS*(pow(VX[I], 2) + pow(VY[I], 2) + pow(VZ[I], 2)); 
    }
     
    // SCALE VELOCITIES TO SET-POINT TEMPERATURE 
    double FACTOR = pow( 3 * NATOM * kT / VELSQ, 1/2);
    for (int I = 0; I < NATOM; ++I) {
        VX[I] = VX[I] * FACTOR; 
        VY[I] = VY[I] * FACTOR;
        VZ[I] = VZ[I] * FACTOR;    
    }
 
    }
