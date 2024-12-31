#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <cstdio>
#include <fstream>

#include "iniconfig.h"
#include "forces2.h"
#include "vverlet.h"

int main() {

    // PARAMETERS TO BE SET BY STUDENTS 
    int NATOM = 200;                // Number of atoms 
    double kB = 1.0;                // Boltzmann constant [amu, ps, A]
    double SIGMA = 3.0;             // LJ size parameter [A]
    double EPSILON = 200.0 * kB;    // LJ energy parameter [amu, A^2, ps^-2]
    double MASS = 10.0;             // Atom mass [amu]
    double DENS = 0.005;            // System density [atoms/A^3]
    double TEM = 100.0;             // Temperature [K]
    int EQSTEPS = 10000;            // Equilibration time steps 
    int PRODSTEPS = 50000;          // Production time steps 

    // DO NOT CHANGE THE PARAMETERS BELOW UNLESS YOU KNOW WHAT YOU ARE DOING 
    double DT = 0.002;              // Integration step [ps] 
    double LJCUT = 2.5 * SIGMA;     // LJ potential cutoff 
    double kT = kB * TEM;           // kB*T
    int IE = 1;                     // Writing intervals - energy 
    int ITAPE = 100;                // Writing intervals - positions and velocities 
    int IANIMATE = 100;             // Animation intervals 
    
    // OPEN FILES TO STORE OUTPUT DATA 
    std::ofstream tapePfile("TAPE_POS.dat");    // STORES POSITIONS EVERY ITAPE STEPS 
    std::ofstream tapeVfile("TAPE_VEL.dat");    // STORES VELOCITIES EVERY ITAPE STEPS 
    std::ofstream moviefile("MOVIE.xyz");       // POSITIONS IN MOVIE FORMAT EVERY IANIMATE 
    std::ofstream energyfile("ENERGY.dat");     // STORES ENERGIES EVERY IE STEPS 
    std::ofstream pressurefile("PRESSURE.dat"); // STORES PRESSURES EVERY IE STEPS 

    // CALCULATE LONG RANGE CORRECTIONS FOR ENERGY AND PRESSURE A PRIORI
    double ELRC = (1.0/3.0) * pow(SIGMA/LJCUT, 9) - pow(SIGMA/LJCUT, 3);
    ELRC = ELRC * (8.0/3.0) * NATOM * M_PI * EPSILON * DENS * pow(SIGMA, 3);
    double PLRC = (2.0/3.0) * pow(SIGMA/LJCUT, 9) - pow(SIGMA/LJCUT, 3);
    PLRC = PLRC * (16.0/3.0) * M_PI * EPSILON * DENS * DENS * pow(SIGMA, 3);
 
    // INITIALIZE THE POSITIONS AND VELOCITIES OF ATOMS 
    printf("STARTING INITIALIZATION \n"); 
    printf("\n");

    // THE POSITIONS ARE CHOSEN SUCH THAT THERE ARE NO MAJOR OVERLAPS OF ATOMS 
    // CALCULATE SIMULATION BOX LENGTH ACCORDING TO DENSITY AND NATOM    
    double LSIMBOX = pow(NATOM / DENS, 1.0 / 3.0);
    double LSIMBOX2 = LSIMBOX/2.0;
    if (LJCUT > LSIMBOX2) {
        std::cout << "******** WARNING: LJCUT > L/2 ********" << std::endl;
        std::cout << "LJCUT = " << LJCUT << ", L/2 = " << LSIMBOX2 << std::endl;
        std::cout << "******** STOPPING SIMULATION ********" << std::endl;
        exit(1); // Exit the program with status 1
    } 

    // FUNCTION FOR INITIALIZING POSITIONS AND VELOCITIES

    std::vector<double> FX(NATOM, 0.0);
    std::vector<double> FY(NATOM, 0.0);
    std::vector<double> FZ(NATOM, 0.0);
    std::vector<double> VX(NATOM, 0.0);
    std::vector<double> VY(NATOM, 0.0);
    std::vector<double> VZ(NATOM, 0.0);
    std::vector<double> RX(NATOM, 0.0);
    std::vector<double> RY(NATOM, 0.0);
    std::vector<double> RZ(NATOM, 0.0);

    iniconfig(RX, RY, RZ, VX, VY, VZ, FX, FY, FZ, NATOM, DENS, kT, SIGMA, EPSILON, MASS, LJCUT); 

    // WRITE OUT AN XYZ FILE
    std::ofstream xyzfile("init.xyz");
    // Write the number of atoms
    xyzfile << NATOM << "\n";
    xyzfile << "Comment line\n"; // You can modify this line as needed
    // Write the atom type and coordinates for each atom
    for (int i = 0; i < NATOM; ++i) {
        xyzfile << "H " << RX[i] << " " << RY[i] << " " << RZ[i] << "\n";
    }
    xyzfile.close();


    // EQUILIBRATION STAGE 
 
    double EKIN, EPOT, ETOT, PRESS;

    printf("\n");
    printf("STARTING EQUILIBRATION \n");
    printf("\n");

    bool RESCALE = true;                           // RESCALE VELOCITIES TO KEEP TEMPERATURE FIXED 

    for (int T = 0; T < EQSTEPS; ++T) {
        // WRITE TIMESTEP TO SCREEN 
        if (T % 1000 == 0) {
            printf("EQ STEP %6i \n", T);
        }
        
        // UPDATE POSITIONS AND VELOCITIES OF THE ATOMS 
        vverlet(RX, RY, RZ, VX, VY, VZ, FX, FY, FZ, NATOM, DENS, LJCUT, SIGMA, EPSILON, MASS, DT, kT, RESCALE);
        
        // ADD LONG-RANGE CORRECTIONS TO ENERGY AND PRESSURE 
        ETOT = ETOT + ELRC; 
        EPOT = EPOT + ELRC;
        PRESS = PRESS + PLRC; 

        // WRITE OUT ENERGIES AND PRESSURES EVERY IE STEPS 
        if (T % IE == 0) {
            energyfile << T << " " << EKIN << " " << EPOT << " " << ETOT << "\n";
            pressurefile << T << " " << PRESS << "\n";
        }
    }

  
    // PRODUCTION STEP 

    printf("\n");
    printf("STARTING PRODUCTION \n");
    printf("\n");

    RESCALE = false;                            // STOP RESCALING VELOCITIES (NVE ENSEMBLE)

    for (int T = 0; T < PRODSTEPS; ++T) {
        // WRITE TIMESTEP TO SCREEN
        if (T % 1000 == 0) {
            printf("PROD STEP %6i \n", T);
        }

        // UPDATE POSITIONS AND VELOCITIES OF THE ATOMS
        vverlet(RX, RY, RZ, VX, VY, VZ, FX, FY, FZ, NATOM, DENS, LJCUT, SIGMA, EPSILON, MASS, DT, kT, RESCALE);     

        // ADD LONG-RANGE CORRECTIONS TO ENERGY AND PRESSURE
        ETOT = ETOT + ELRC;
        EPOT = EPOT + ELRC;
        PRESS = PRESS + PLRC;

	// WRITE DOWN POSITIONS OF ATOMS IN THE TAPE FILE EVERY ITAPE STEPS
        if (T % ITAPE == 0) {
            for (int I = 0; I < NATOM; ++I) {
                tapePfile << " " << RX[I] << " " << RY[I] << " " << RZ[I] << "\n";
                tapeVfile << " " << VX[I] << " " << VY[I] << " " << VZ[I] << "\n";
            }
        }

        // WRITE DOWN THE POSITIONS OF THE ATOMS IN AN "XMOL" FORMAT TO MAKE
        // THE MOVIE EVERY IANIMATE STEPS
        if (T % IANIMATE == 0) {
            moviefile << NATOM << "\n";
            moviefile << T << "Comment line\n";
            // Write the atom type and coordinates for each atom
            for (int i = 0; i < NATOM; ++i) { 
                moviefile << "H " << RX[i] << " " << RY[i] << " " << RZ[i] << "\n";
            }
            
        }

        // WRITE OUT ENERGIES AND PRESSURES EVERY IE STEPS
        if (T % IE == 0) {
            energyfile << T + EQSTEPS << " " << EKIN << " " << EPOT << " " << ETOT << "\n";
            pressurefile << T + EQSTEPS << " " << PRESS << "\n";
        }

    }

    // CLOSE FILES 
    tapePfile.close();
    tapeVfile.close();
    moviefile.close();
    energyfile.close();
    pressurefile.close();
    
    return 0; 

}    
