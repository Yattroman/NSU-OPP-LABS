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

double calculateMPlusOnePhiValue(double ** phiValuesPart, int NxPart, int NyPart, int NzPart, const double& Hx, const double& Hy, const double& Hz){
    int NzAddition = NzPart % 2;
    int medValZ = NzPart / 2 + NzAddition;
    double a[4];
    double b[4];

    double maxDiff = numeric_limits<double>::max();

    double denominator = 1 / ( 2/Hx*Hx + 2/Hy*Hy + 2/Hz*Hz + A );

    double x, y, z;

    if(NzAddition == 1) {
        // NzPart mod 2 = 1
        for(int k = 0; k < medValZ; ++k){
            for (int j = 0; j < NyPart; ++j) {
                for (int i = 0; i < NxPart; ++i) {
                    // phi[M]_{i+1, j, k} + phi[M]_{i-1, j, k}
                    a[0] = isVld(i, NxPart, +1)*phiValuesPart[0][(i+1) + j*NxPart + (medValZ-1+k)*NxPart*NyPart]
                            + isVld(i, NxPart, -1)*phiValuesPart[0][(i-1) + j*NxPart + (medValZ-1+k)*NxPart*NyPart];

                    // phi[M]_{i, j+1, k} + phi[M]_{i, j-1, k}
                    a[1] = isVld(j, NyPart, +1)*phiValuesPart[0][i + (j+1)*NxPart + (medValZ-1+k)*NxPart*NyPart]
                           + isVld(j, NyPart, -1)*phiValuesPart[0][i + (j-1)*NxPart + (medValZ-1+k)*NxPart*NyPart];

                    // phi[M]_{i, j, k+1} + phi[M]_{i, j, k-1}
                    a[2] = isVld(medValZ-1+k, NzPart, +1)*phiValuesPart[0][i + j*NxPart + (medValZ-1+k+1)*NxPart*NyPart]
                           + isVld(medValZ-1+k, NzPart, -1)*phiValuesPart[0][i + j*NxPart + (medValZ-1+k-1)*NxPart*NyPart];

                    // ro_{i, j, k}
                    calculatePhiArguments(i, j, medValZ-1+k, x, y, z, Hx, Hy, Hz);
                    a[3] = ro(phi(x, y, z));

                    // phi[M]_{i+1, j, k} + phi[M]_{i-1, j, k}
                    b[0] = isVld(i, NxPart, +1)*phiValuesPart[0][(i+1) + j*NxPart + (medValZ-1-k)*NxPart*NyPart]
                           + isVld(i, NxPart, -1)*phiValuesPart[0][(i-1) + j*NxPart + (medValZ-1-k)*NxPart*NyPart];

                    // phi[M]_{i, j+1, k} + phi[M]_{i, j-1, k}
                    b[1] = isVld(j, NyPart, +1)*phiValuesPart[0][i + (j+1)*NxPart + (medValZ-1-k)*NxPart*NyPart]
                           + isVld(j, NyPart, -1)*phiValuesPart[0][i + (j-1)*NxPart + (medValZ-1-k)*NxPart*NyPart];

                    // phi[M]_{i, j, k+1} + phi[M]_{i, j, k-1}
                    b[2] = isVld(medValZ-1-k, NzPart, +1)*phiValuesPart[0][i + j*NxPart + (medValZ-1-k+1)*NxPart*NyPart]
                           + isVld(medValZ-1-k, NzPart, -1)*phiValuesPart[0][i + j*NxPart + (medValZ-1-k-1)*NxPart*NyPart];

                    // ro_{i, j, k}
                    calculatePhiArguments(i, j, medValZ-1-k, x, y, z, Hx, Hy, Hz);
                    b[3] = ro(phi(x, y, z));

                    // phi[M+1]_{i, j, k}
                    phiValuesPart[1][j + i*NxPart + (medValZ-1+k)*NxPart*NyPart] = (a[0]/(Hx*Hx) + a[1]/(Hy*Hy) + a[2]/(Hz*Hz) + a[3]) / denominator;

                    // phi[M+1]_{i, j, k}
                    phiValuesPart[1][j + i*NxPart + (medValZ-1-k)*NxPart*NyPart] = (b[0]/(Hx*Hx) + b[1]/(Hy*Hy) + b[2]/(Hz*Hz) + b[3]) / denominator;

                    if(abs(phiValuesPart[1][i + j*NxPart + (medValZ-1+k)*NxPart*NyPart] - phiValuesPart[0][i + j*NxPart + (medValZ-1+k)*NxPart*NyPart]) > maxDiff){
                        maxDiff = abs(phiValuesPart[1][i + j*NxPart + (medValZ-1+k)*NxPart*NyPart] - phiValuesPart[0][i + j*NxPart + (medValZ-1+k)*NxPart*NyPart]);
                    }

                    if(abs(phiValuesPart[1][i + j*NxPart + (medValZ-1-k)*NxPart*NyPart] - phiValuesPart[0][i + j*NxPart + (medValZ-1-k)*NxPart*NyPart]) > maxDiff){
                        maxDiff = abs(phiValuesPart[1][i + j*NxPart + (medValZ-1-k)*NxPart*NyPart] - phiValuesPart[0][i + j*NxPart + (medValZ-1-k)*NxPart*NyPart]);
                    }
                }
            }
        }
    } else {
        // NzPart mod 2 = 0
        for(int k = 0; k < medValZ; ++k){
            for (int j = 0; j < NyPart; ++j) {
                for (int i = 0; i < NxPart; ++i) {
                    // phi[M]_{i+1, j, k} + phi[M]_{i-1, j, k}
                    a[0] = isVld(i, NxPart, +1)*phiValuesPart[0][(i+1) + j*NxPart + (medValZ+k)*NxPart*NyPart]
                           + isVld(i, NxPart, -1)*phiValuesPart[0][(i-1) + j*NxPart + (medValZ+k)*NxPart*NyPart];

                    // phi[M]_{i, j+1, k} + phi[M]_{i, j-1, k}
                    a[1] = isVld(j, NyPart, +1)*phiValuesPart[0][i + (j+1)*NxPart + (medValZ+k)*NxPart*NyPart]
                           + isVld(j, NyPart, -1)*phiValuesPart[0][i + (j-1)*NxPart + (medValZ+k)*NxPart*NyPart];

                    // phi[M]_{i, j, k+1} + phi[M]_{i, j, k-1}
                    a[2] = isVld(medValZ+k, NzPart, +1)*phiValuesPart[0][i + j*NxPart + (medValZ+k+1)*NxPart*NyPart]
                           + isVld(medValZ+k, NzPart, -1)*phiValuesPart[0][i + j*NxPart + (medValZ+k-1)*NxPart*NyPart];

                    // ro_{i, j, k}
                    calculatePhiArguments(i, j, medValZ+k, x, y, z, Hx, Hy, Hz);
                    a[3] = ro(phi(x, y, z));

                    // phi[M]_{i+1, j, k} + phi[M]_{i-1, j, k}
                    b[0] = isVld(i, NxPart, +1)*phiValuesPart[0][(i+1) + j*NxPart + (medValZ-k-1)*NxPart*NyPart]
                           + isVld(i, NxPart, -1)*phiValuesPart[0][(i-1) + j*NxPart + (medValZ-k-1)*NxPart*NyPart];

                    // phi[M]_{i, j+1, k} + phi[M]_{i, j-1, k}
                    b[1] = isVld(j, NyPart, +1)*phiValuesPart[0][i + (j+1)*NxPart + (medValZ-k-1)*NxPart*NyPart]
                           + isVld(j, NyPart, -1)*phiValuesPart[0][i + (j-1)*NxPart + (medValZ-k-1)*NxPart*NyPart];

                    // phi[M]_{i, j, k+1} + phi[M]_{i, j, k-1}
                    b[2] = isVld(medValZ-k-1, NzPart, +1)*phiValuesPart[0][i + j*NxPart + (medValZ-k-1+1)*NxPart*NyPart]
                           + isVld(medValZ-k-1, NzPart, -1)*phiValuesPart[0][i + j*NxPart + (medValZ-k-1-1)*NxPart*NyPart];

                    // ro_{i, j, k}
                    calculatePhiArguments(i, j, medValZ-k-1, x, y, z, Hx, Hy, Hz);
                    b[3] = ro(phi(x, y, z));

                    // phi[M+1]_{i, j, k}
                    phiValuesPart[1][j + i*NxPart + (medValZ+k)*NxPart*NyPart] = (a[0]/(Hx*Hx) + a[1]/(Hy*Hy) + a[2]/(Hz*Hz) + a[3]) / denominator;

                    // phi[M+1]_{i, j, k}
                    phiValuesPart[1][j + i*NxPart + (medValZ-k-1)*NxPart*NyPart] = (b[0]/(Hx*Hx) + b[1]/(Hy*Hy) + b[2]/(Hz*Hz) + b[3]) / denominator;

                    if(abs(phiValuesPart[1][i + j*NxPart + (medValZ+k)*NxPart*NyPart] - phiValuesPart[0][i + j*NxPart + (medValZ+k)*NxPart*NyPart]) > maxDiff){
                        maxDiff = abs(phiValuesPart[1][i + j*NxPart + (medValZ+k)*NxPart*NyPart] - phiValuesPart[0][i + j*NxPart + (medValZ+k)*NxPart*NyPart]);
                    }

                    if(abs(phiValuesPart[1][i + j*NxPart + (medValZ-k-1)*NxPart*NyPart] - phiValuesPart[0][i + j*NxPart + (medValZ-k-1)*NxPart*NyPart]) > maxDiff){
                        maxDiff = abs(phiValuesPart[1][i + j*NxPart + (medValZ-k-1)*NxPart*NyPart] - phiValuesPart[0][i + j*NxPart + (medValZ-k-1)*NxPart*NyPart]);
                    }
                }
            }
        }
    }

    return maxDiff;
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

    double maxDiff;

    for(auto & i : phiValuesPart) {
        i = new double[NzPart*Ny*Nx];
    }

    fillStartPhiValues(phiValuesPart[0], Nx, Ny, Nz, Hx, Hy, Hz);

    do{
        maxDiff = calculateMPlusOnePhiValue(phiValuesPart, Nx, Ny, NzPart, Hx, Hy, Hz);
    } while (maxDiff < EPSILON);

    /*for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nx; ++j) {
            cout << phiSolid[i*Ny + j + (Nz-2)*Nx*Nz] << " ";
        }
        cout << endl;
    }*/

}
