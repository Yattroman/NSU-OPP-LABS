#include <iostream>
#include <cstring>
#include "mpi.h"

#define EPSILON 1e-8
#define MAX_DIFF_DELTA 1e2

#define A 1e5
#define START_VALUE_N 300.0

#define X0 -1.0

#define Y0 -1.0

#define Z0 -1.0

#define DX 2.0
#define DY 2.0
#define DZ 2.0

#define TRUE 1

using namespace std;

void cutStartPhiValuesMatrix(long double * phiValuesPart, long double * phiValuesSolid, const int& NxPart, const int& NyPart, const int& NzPart){
    // Nx MOD procCount == 0 | Ny MOD procCount == 0 | Ny MOD procCount == 0
//    MPI_Scatter(phiValuesSolid, NxPart*NyPart*NzPart, MPI_long double, phiValuesPart, NxPart*NyPart*NzPart, MPI_long double, 0, MPI_COMM_WORLD);
}

int isVld(int value, int limit, int addition){
    if(value + addition < 0 || value + addition >= limit){
        return 0;
    }

    return 1;
}

void calculatePhiArguments(int i, int j, int k, long double &x, long double &y, long double &z, const long double &Hx, const long double &Hy, const long double &Hz){
    x = X0 + i*Hx;
    y = Y0 + j*Hy;
    z = Z0 + k*Hz;
}

long double rho(long double phiValue){
    return (6.0 - A*phiValue);
}

long double phi(long double x, long double y, long double z){
    return (x*x + y*y + z*z);
}

void fillStartPhiValues(long double ** phiValuesPart, int Nx, int Ny, int Nz, const long double& Hx, const long double& Hy, const long double& Hz){
    memset(phiValuesPart[0], 0, Nx*Ny*Nz);

    long double x, y, z;

    for (int s = 0; s < 2; ++s) {
        // phi(xi, yj, 1), phi(xi, yj, -1)
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                calculatePhiArguments(i, j, s*(Nz-1), x, y, z, Hx, Hy, Hz);
                phiValuesPart[0][i + j*Nx + s*(Nz-1)*Nx*Ny] = phi(x, y, z);
            }
        }

        // phi(xi, 1, zk), phi(xi, -1, zk)
        for (int k = 0; k < Nz; ++k) {
            for (int i = 0; i < Nx; ++i) {
                calculatePhiArguments(i, s*(Ny-1), k, x, y, z, Hx, Hy, Hz);
                phiValuesPart[0][i +s*(Ny-1)*Nx + k*Nx*Ny] = phi(x, y, z);
            }
        }

        // phi(1, yj, zk), phi(-1, yj, zk)
        for (int k = 0; k < Nz; ++k) {
            for (int j = 0; j < Ny; ++j) {
                calculatePhiArguments(s*(Nx-1), j, k, x, y, z, Hx, Hy, Hz);
                phiValuesPart[0][s*(Nx-1) + j*Nx + k*Nx*Ny] = phi(x, y, z);
            }
        }
    }

    memcpy(phiValuesPart[1], phiValuesPart[0], Nx*Ny*Nz*sizeof(long double));
}

void calculateMPlusOnePhiValue(long double ** phiValuesPart, int NxPart, int NyPart, int NzPart, const long double& Hx, const long double& Hy, const long double& Hz, long double &maxDiff, long double &delta){
    int NzAddition = NzPart % 2;
    int medValZ = NzPart / 2;
    int K[2];

    long double a[4];
    long double b[4];
    long double precisePhiValue[2];

    maxDiff = 0;
    delta = 0;

    long double denominator = ( 2.0/(Hx*Hx) + 2.0/(Hy*Hy) + 2.0/(Hz*Hz) + A );

    long double x, y, z;

    for(int k = 0; k < medValZ + NzAddition - 1; ++k){
        for (int j = 1; j < NyPart - 1; ++j) {
            for (int i = 1; i < NxPart - 1; ++i) {
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

                // rho_{i, j, k}
                calculatePhiArguments(i, j, K[0], x, y, z, Hx, Hy, Hz);
                precisePhiValue[0] = phi(x, y, z);
                a[3] = -rho(precisePhiValue[0]);

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

                // rho_{i, j, k}
                calculatePhiArguments(i, j, K[1], x, y, z, Hx, Hy, Hz);
                precisePhiValue[1] = phi(x, y, z);
                b[3] = -rho(precisePhiValue[1]);

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

    memcpy(phiValuesPart[0], phiValuesPart[1], NxPart*NyPart*NzPart*sizeof(long double));
}

int main(int argc, char* argv[]){
    long double * phiValuesPart[2];

//    MPI_Init(&argc,&argv);

    int procCount = 1;
    int procRank = 1;

    int Nx = START_VALUE_N;
    int Ny = START_VALUE_N;
    int Nz = START_VALUE_N;

    int NxPart = Nx / 1;
    int NyPart = Ny / 1;
    int NzPart = Nz / procCount;

    long double Hx = DX / ((long double) Nx - 1);
    long double Hy = DY / ((long double) Ny - 1);
    long double Hz = DZ / ((long double) Nz - 1);

    long double maxDiff;
    long double delta;

    for(auto & i : phiValuesPart) {
        i = new long double[NzPart*NyPart*NxPart];
    }

    fillStartPhiValues(phiValuesPart, Nx, Ny, Nz, Hx, Hy, Hz);

    do{
        calculateMPlusOnePhiValue(phiValuesPart, NxPart, NyPart, NzPart, Hx, Hy, Hz, maxDiff, delta);

        if(maxDiff < EPSILON){
            break;
        }

    } while (TRUE);

    if(delta > EPSILON*MAX_DIFF_DELTA){
        cout << "DELTA is bigger than EPSILON*1e2.\n";
    } else {
        cout << "Function is found";
    }

//    cout << phiValuesPart[1][2 + 3*NxPart + 4*NxPart*NyPart] << endl;
//    cout << Hx << endl;
    cout << delta;

/*    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            cout << phiValuesPart[0][i + j*Nx + (1)*Nx*Ny] << " ";
        }
        cout << endl;
    }*/

}
