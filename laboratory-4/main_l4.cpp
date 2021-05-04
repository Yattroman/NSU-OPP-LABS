#include <iostream>
#include <cstring>
#include "mpi.h"

#define EPSILON 1e-10
#define MAX_DIFF_DELTA 1e2

#define A 1e5
#define START_VALUE_N 11

#define X0 -1.0
#define XNX 1.0

#define Y0 -1.0
#define YNY 1.0

#define Z0 -1.0
#define ZNZ 1.0

#define DX XNX - X0
#define DY YNY - Y0
#define DZ ZNZ - Z0

using namespace std;

int isVld(int value, int limit, int addition){
    if(value + addition < 0 || value + addition >= limit){
        return 0;
    }

    return 1;
}

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
                calculatePhiArguments(i2, i1, marginalValues[s]*(Nz-1), x, y, z, Hx, Hy, Hz);
                phiSolid[i2 + i1*Nx + marginalValues[s]*(Nz-1)*Nx*Ny] = phi(x, y, z);
            }
        }

        // phi(xi, 1, zk), phi(xi, -1, zk)
        for (int i1 = 0; i1 < Nz; ++i1) {
            for (int i2 = 0; i2 < Nx; ++i2) {
                calculatePhiArguments(i2, marginalValues[s]*(Ny-1), i1, x, y, z, Hx, Hy, Hz);
                phiSolid[i2 + marginalValues[s]*(Ny-1)*Nx + i1*Nx*Ny] = phi(x, y, z);
            }
        }

        // phi(1, yj, zk), phi(-1, yj, zk)
        for (int i1 = 0; i1 < Nz; ++i1) {
            for (int i2 = 0; i2 < Ny; ++i2) {
                calculatePhiArguments(marginalValues[s]*(Nx-1), i2, i1, x, y, z, Hx, Hy, Hz);
                phiSolid[marginalValues[s]*(Nx-1) + i2*Nx + i1*Nx*Ny] = phi(x, y, z);
            }
        }
    }
}

void calculateMPlusOnePhiValue(double ** phiValuesPart, int NxPart, int NyPart, int NzPart, const double& Hx, const double& Hy, const double& Hz, double &maxDiff, double &delta){
    int NzAddition = NzPart % 2;
    int medValZ = NzPart / 2;
    double a[4];
    double b[4];
    double precisePhiValue[2];
    int K[2];

    double denominator = 1 / ( 2/(Hx*Hx) + 2/(Hy*Hy) + 2/(Hz*Hz) + A );

    double x, y, z;

    for(int k = 0; k < medValZ + NzAddition; ++k){
        for (int j = 0; j < NyPart; ++j) {
            for (int i = 0; i < NxPart; ++i) {
                K[0] = medValZ+k; // NzPart MOD 2 == 1 -> K[0] = medValZ + k | NzPart MOD 2 == 0 -> K[0] = medValZ + k
                // phi[M]_{i+1, j, k} + phi[M]_{i-1, j, k}
                a[0] = isVld(i, NxPart, +1)*phiValuesPart[0][(i+1) + j*NxPart + K[0]*NxPart*NyPart]
                       + isVld(i, NxPart, -1)*phiValuesPart[0][(i-1) + j*NxPart + K[0]*NxPart*NyPart];

                // phi[M]_{i, j+1, k} + phi[M]_{i, j-1, k}
                a[1] = isVld(j, NyPart, +1)*phiValuesPart[0][i + (j+1)*NxPart + K[0]*NxPart*NyPart]
                       + isVld(j, NyPart, -1)*phiValuesPart[0][i + (j-1)*NxPart + K[0]*NxPart*NyPart];

                // phi[M]_{i, j, k+1} + phi[M]_{i, j, k-1}
                a[2] = isVld(K[0], NzPart, +1)*phiValuesPart[0][i + j*NxPart + (K[0]+1)*NxPart*NyPart]
                       + isVld(K[0], NzPart, -1)*phiValuesPart[0][i + j*NxPart + (K[0]-1)*NxPart*NyPart];

                // ro_{i, j, k}
                calculatePhiArguments(i, j, K[0], x, y, z, Hx, Hy, Hz);
                precisePhiValue[0] = phi(x, y, z);
                a[3] = ro(precisePhiValue[0]);

                K[1] = (medValZ-1+NzAddition-k); // NzPart MOD 2 == 1 -> K[0] = medValZ - k | NzPart MOD 2 == 0 -> K[0] = medValZ - 1 - k
                // phi[M]_{i+1, j, k} + phi[M]_{i-1, j, k}
                b[0] = isVld(i, NxPart, +1)*phiValuesPart[0][(i+1) + j*NxPart + K[1]*NxPart*NyPart]
                       + isVld(i, NxPart, -1)*phiValuesPart[0][(i-1) + j*NxPart + K[1]*NxPart*NyPart];

                // phi[M]_{i, j+1, k} + phi[M]_{i, j-1, k}
                b[1] = isVld(j, NyPart, +1)*phiValuesPart[0][i + (j+1)*NxPart + K[1]*NxPart*NyPart]
                       + isVld(j, NyPart, -1)*phiValuesPart[0][i + (j-1)*NxPart + K[1]*NxPart*NyPart];

                // phi[M]_{i, j, k+1} + phi[M]_{i, j, k-1}
                b[2] = isVld(K[1], NzPart, +1)*phiValuesPart[0][i + j*NxPart + (K[1]+1)*NxPart*NyPart]
                       + isVld(K[1], NzPart, -1)*phiValuesPart[0][i + j*NxPart + (K[1]-1)*NxPart*NyPart];

                // ro_{i, j, k}
                calculatePhiArguments(i, j, K[1], x, y, z, Hx, Hy, Hz);
                precisePhiValue[1] = phi(x, y, z);
                b[3] = ro(precisePhiValue[1]);

                // phi[M+1]_{i, j, k}
                phiValuesPart[1][i + j*NxPart + K[0]*NxPart*NyPart] = (a[0]/(Hx*Hx) + a[1]/(Hy*Hy) + a[2]/(Hz*Hz) + a[3]) / denominator;

                // phi[M+1]_{i, j, k}
                phiValuesPart[1][i + j*NxPart + K[1]*NxPart*NyPart] = (b[0]/(Hx*Hx) + b[1]/(Hy*Hy) + b[2]/(Hz*Hz) + b[3]) / denominator;

                if(abs(phiValuesPart[1][i + j*NxPart + K[0]*NxPart*NyPart] - phiValuesPart[0][i + j*NxPart + K[0]*NxPart*NyPart]) > maxDiff){
                    maxDiff = abs(phiValuesPart[1][i + j*NxPart + K[0]*NxPart*NyPart] - phiValuesPart[0][i + j*NxPart + K[0]*NxPart*NyPart]);
                }

                if(abs(phiValuesPart[1][i + j*NxPart + K[1]*NxPart*NyPart] - phiValuesPart[0][i + j*NxPart + K[1]*NxPart*NyPart]) > maxDiff){
                    maxDiff = abs(phiValuesPart[1][i + j*NxPart + K[1]*NxPart*NyPart] - phiValuesPart[0][i + j*NxPart + K[1]*NxPart*NyPart]);
                }

                if(abs(phiValuesPart[0][i + j*NxPart + K[0]*NxPart*NyPart] - precisePhiValue[0]) > delta){
                    delta = abs(phiValuesPart[0][i + j*NxPart + K[0]*NxPart*NyPart] - precisePhiValue[0]);
                }

                if(abs(phiValuesPart[0][i + j*NxPart + K[1]*NxPart*NyPart] - precisePhiValue[1]) > delta){
                    delta = abs(phiValuesPart[0][i + j*NxPart + K[1]*NxPart*NyPart] - precisePhiValue[1]);
                }
            }
        }
    }

    memcpy(phiValuesPart[0], phiValuesPart[1], NxPart*NyPart*NzPart);
}

int main(){
    double * phiValuesPart[2]; // 0 ~ phiPart[0], 1 ~ phiPart[1]

    int procCount = 1;
    int procRank = 1;

    int Nx = START_VALUE_N;
    int Ny = START_VALUE_N;
    int Nz = START_VALUE_N;

    int NzPart = Nz / procCount;

    double Hx = DX / (Nx - 1);
    double Hy = DY / (Ny - 1);
    double Hz = DZ / (Nz - 1);

    double maxDiff = numeric_limits<double>::min();
    double delta = numeric_limits<double>::min();

    for(auto & i : phiValuesPart) {
        i = new double[NzPart*Ny*Nx];
    }

    fillStartPhiValues(phiValuesPart[0], Nx, Ny, Nz, Hx, Hy, Hz);

    do{
        calculateMPlusOnePhiValue(phiValuesPart, Nx, Ny, NzPart, Hx, Hy, Hz, maxDiff, delta);

        if(delta > EPSILON*MAX_DIFF_DELTA){
            cout << "DELTA is bigger than EPSILON*1e2.\n";
            break;
        }
    } while (maxDiff > EPSILON);

    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nx; ++j) {
            cout << phiValuesPart[1][i*Ny + j + (0)*Nx*Nz] << " ";
        }
        cout << endl;
    }

}
