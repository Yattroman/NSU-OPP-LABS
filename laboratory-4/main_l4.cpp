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

void calculatePhiArguments(int i, int j, int k, double &x, double &y, double &z, const double &Hx, const double &Hy, const double &Hz){
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

void fillStartPhiValues(double * phiSolid, int Nx, int Ny, int Nz, const double& Hx, const double& Hy, const double& Hz){
    memset(phiSolid, 0, Nx*Ny*Nz);

    double x, y, z;

    int marginalValues[2] = {0, 1};

    for (int s = 0; s < 2; ++s) {
        // phi(xi, yj, 1), phi(xi, yj, -1)
        for (int i1 = 0; i1 < Ny; ++i1) {
            for (int i2 = 0; i2 < Nx; ++i2) {
                calculatePhiArguments(i2, i1, Nz-1, x, y, z, Hx, Hy, Hz);
                phiSolid[i1*Nx + i2 + marginalValues[s]*(Nz-1)*Nx*Ny] = phi(x, y, z);
            }
        }

        // phi(xi, 1, zk), phi(xi, -1, zk)
        for (int i1 = 0; i1 < Nz; ++i1) {
            for (int i2 = 0; i2 < Nx; ++i2) {
                calculatePhiArguments(i2, Ny-1, i1, x, y, z, Hx, Hy, Hz);
                phiSolid[i1*Nx*Ny + i2 + marginalValues[s]*(Ny-1)*Nx] = phi(x, y, z);
            }
        }

        // phi(1, yj, zk), phi(-1, yj, zk)
        for (int i1 = 0; i1 < Nz; ++i1) {
            for (int i2 = 0; i2 < Ny; ++i2) {
                calculatePhiArguments(Nx-1, i2, i1, x, y, z, Hx, Hy, Hz);
                phiSolid[i1*Nx*Ny + i2*Nx + marginalValues[s]*(Nx-1)] = phi(x, y, z);
            }
        }
    }
}

char calculateCompleteCondition(double ** phiValuesPart,  const double& Hx, const double& Hy, const double& Hz){

}

void calculateMPlusOnePhiValue(double ** phiValuesPart, int NxPart, int NyPart, int NzPart, const double& Hx, const double& Hy, const double& Hz){
    int NzAddition = NzPart % 2;
    double medianValueZ = NzPart / 2 + NzAddition;
    double a[4];

    double denominator = 1 / ( 2/Hx*Hx + 2/Hy*Hy + 2/Hz*Hz + A );

    if(NzAddition == 1) {
        // NzPart mod 2 = 1
        for(int k = medianValueZ; k < NzPart; ++k){
            for (int i = 0; i < NyPart; ++i) {
                for (int j = 0; j < NxPart; ++j) {
                    a[0] = 1;
                    a[1] = 1;
                    a[2] = 1;
                    a[3] = 1;
                    phiValuesPart[1][j + i*NxPart + k*NxPart*NyPart] = (a[0] + a[1] + a[2] + a[3]) / denominator;
                }
            }
        }
    } else {
        // NzPart mod 2 = 0

    }
}

int main(){
    double * phiValuesPart[2]; // 0 ~ phiPart[0], 1 ~ phiPart[1]

    int procCount = 1;
    int procRank = 1;

    int Nx = START_VALUE_N;
    int Ny = START_VALUE_N;
    int Nz = START_VALUE_N;

    int NxPart = Nx / procCount;

    double Hx = DX / (Nx - 1);
    double Hy = DY / (Ny - 1);
    double Hz = DZ / (Nz - 1);

    double * phiSolid = new double[Nx*Ny*Nz];

    for(auto & i : phiValuesPart) {
        i = new double[NxPart*Ny*Nz];
    }

    fillStartPhiValues(phiSolid, Nx, Ny, Nz, Hx, Hy, Hz);

    /*for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nx; ++j) {
            cout << phiSolid[i*Ny + j + (Nz-2)*Nx*Nz] << " ";
        }
        cout << endl;
    }*/

}
