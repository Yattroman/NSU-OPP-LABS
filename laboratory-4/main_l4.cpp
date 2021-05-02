#include <iostream>
#include <cstring>
#include "mpi.h"

#define EPSILON 1e-10
#define A 1e5
#define START_VALUE_N 8

#define X0 -1
#define XNX 1

#define Y0 -1
#define YNY 1

#define Z0 -1
#define ZNY 1

#define DX 2
#define DY 2
#define DZ 2

using namespace std;

void calculatePhiArguments(int i, int j, int k, int &x, int &y, int &z, const int &Hx, const int &Hy, const int &Hz){
    x = X0 + i*Hx;
    y = Y0 + j*Hy;
    z = Z0 + k*Hz;
}

double ro(double phiValue){
    return 6 - A*phiValue;
}

double phi(double x, double y, double z){
    return x*x + y*y + z*z;
}

void fillStartPhiValues(double ** phiSolid, int Nx, int Ny, int Nz, const int& Hx, const int& Hy, const int& Hz){
    memset(phiSolid[0], 0, Nx*Ny*Nz);

    int x, y, z;

    int marginalValues[2] = {0, 1};

    for (int s = 0; s < 2; ++s) {
        for (int i1 = 0; i1 < Ny; ++i1) {
            for (int i2 = 0; i2 < Nx; ++i2) {
                calculatePhiArguments(i2, i1, marginalValues[s], x, y, z, Hx, Hy, Hz);
                phiSolid[0][i1*Ny + i2] = phi(x, y, z);
            }
        }

        for (int i1 = 0; i1 < Nz; ++i1) {
            for (int i2 = 0; i2 < Nx; ++i2) {
                calculatePhiArguments(i2, marginalValues[s], i1, x, y, z, Hx, Hy, Hz);
                phiSolid[0][i1*Nz + i2] = phi(x, y, z);
            }
        }

        for (int i1 = 0; i1 < Nz; ++i1) {
            for (int i2 = 0; i2 < Ny; ++i2) {
                calculatePhiArguments(marginalValues[s], i2, i1, x, y, z, Hx, Hy, Hz);
                phiSolid[0][i1*Nz + i2] = phi(x, y, z);
            }
        }
    }
}

char calculateCompleteCondition(){

}

void calculateMarginalPhiValues(){

}

void calculateNextPhiValue(double ** phiPart, int NxPart, int NyPart, int NzPart, int procSize){
    double medianValue[3]; // z - med, y - startValue, x - startValue
}

int main(){
    double * phiPart[2]; // 0 ~ phiPart[0], 1 ~ phiPart[1]
    double * phiMarginalValues[2];

    int procCount = 1;
    int procRank = 1;

    int Nx = START_VALUE_N;
    int Ny = START_VALUE_N;
    int Nz = START_VALUE_N;

    int NxPart = Nx / procCount;
    int NyPart = Ny / procCount;
    int NzPart = Nz / procCount;

    double Hx = DX / (Nx - 1);
    double Hy = DY / (Ny - 1);
    double Hz = DZ / (Nz - 1);

    double * phiSolid = new double[Nx*Ny*Nz];

    for(auto & i : phiPart) {
        i = new double[NxPart*NyPart*NzPart];
    }

    fillStartPhiValues(phiSolid, Nx, Ny, Nz, Hx, Hy, Hz);

}
