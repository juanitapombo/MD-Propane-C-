#include <iostream>
#include <vector>
#include <cmath>

void forces2(int NATOM, double DENS, double LJCUT, double SIGMA, double EPSILON,
             std::vector<double>& RX, std::vector<double>& RY, std::vector<double>& RZ,
             double& LSIMBOX, double& LJCUTSQ, std::vector<double>& FX, std::vector<double>& FY, std::vector<double>& FZ,
             double& EPOT, double& PRESS) {
    LSIMBOX = pow(NATOM / DENS, 1.0 / 3.0);
    LJCUTSQ = LJCUT * LJCUT;

    // SET FORCES ON ALL ATOMS, POTENTIAL ENERGY, AND PRESSURE TO ZERO
    FX.assign(NATOM, 0.0);
    FY.assign(NATOM, 0.0);
    FZ.assign(NATOM, 0.0);
    EPOT = 0.0;
    PRESS = 0.0;

    for (int I = 0; I < NATOM - 1; ++I) {
        double FXI = FX[I];
        double FYI = FY[I];
        double FZI = FZ[I];
        double RXI = RX[I];
        double RYI = RY[I];
        double RZI = RZ[I];
        for (int J = I+1; J < NATOM; ++J) {
            double RXIJ = RXI - RX[J];
            double RYIJ = RYI - RY[J];
            double RZIJ = RZI - RZ[J];
            // APPLY MINIMUM IMAGE CONVENTION
            RXIJ = RXIJ - LSIMBOX * std::round(RXIJ/LSIMBOX);
            RYIJ = RYIJ - LSIMBOX * std::round(RYIJ/LSIMBOX);
            RZIJ = RZIJ - LSIMBOX * std::round(RZIJ/LSIMBOX);
            double RIJSQ = RXIJ*RXIJ + RYIJ*RYIJ + RZIJ*RZIJ;
            if (RIJSQ < LJCUTSQ) {
                double R6I = std::pow(SIGMA, 6)/(RIJSQ * RIJSQ * RIJSQ);
                double R12I = R6I * R6I;
                EPOT = EPOT + (-R6I + R12I);
                double FIJ = -(-12.0 * R12I + 6.0 * R6I) / RIJSQ;
                PRESS = PRESS + FIJ*RIJSQ;
                double FXIJ  = RXIJ * FIJ;
                double FYIJ  = RYIJ * FIJ;
                double FZIJ  = RZIJ * FIJ;
                FXI   = FXI + FXIJ;
                FYI   = FYI + FYIJ;
                FZI   = FZI + FZIJ;
                FX[J] = FX[J] - FXIJ;
                FY[J] = FY[J] - FYIJ;
                FZ[J] = FZ[J] - FZIJ;
             }
        }
        FX[I] = FXI;
        FY[I] = FYI;
        FZ[I] = FZI;
    }

    // MULTIPLY BACK THE MISSING 4*EPS FACTOR IN THE COMPUTED LJ ENERGY
    EPOT = 4.0 * EPSILON * EPOT;

    // MULTIPLY BACK THE MISSING 4*EPS FACTOR IN THE COMPUTED RIJ*FIJ
    // AND DIVIDE BY 3*VOLUME
    PRESS = 4.0 * EPSILON * DENS * PRESS / 3.0 / NATOM;

    for (int I = 0; I < NATOM - 1; ++I) {
        FX[I] = 4.0 * EPSILON * FX[I];
        FY[I] = 4.0 * EPSILON * FY[I];
        FZ[I] = 4.0 * EPSILON * FZ[I];
    }
}

