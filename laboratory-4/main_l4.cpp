#include <iostream>
#include <cstring>
#include <cmath>
#include "mpi.h"

#define EPSILON 1e-8
#define MAX_DIFF_deltaLocal 1e2

#define A 1e5
#define START_VALUE_N 12.0

#define X0 -1.0

#define Y0 -1.0

#define Z0 -1.0

#define DX 2.0
#define DY 2.0
#define DZ 2.0

#define TRUE 1

using namespace std;

void waitMessages(MPI_Request * reqr, MPI_Request * reqs, int procRankSend, int procCount){
    if(procRankSend != 0) {
        MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        MPI_Wait(&reqr[1], MPI_STATUS_IGNORE);
    }
    if(procRankSend != procCount - 1) {
        MPI_Wait(&reqs[1], MPI_STATUS_IGNORE);
        MPI_Wait(&reqr[0], MPI_STATUS_IGNORE);
    }
}

void getBoundaries(long double ** boundaries, int procRankSend, int count, MPI_Request * reqr, int procCount){
    if(procRankSend != 0) {
        MPI_Irecv(boundaries[0], count, MPI_LONG_DOUBLE, procRankSend - 1, 0, MPI_COMM_WORLD,&reqr[1]);
    }

    if(procRankSend != procCount - 1) {
        MPI_Irecv(boundaries[1], count, MPI_LONG_DOUBLE, procRankSend + 1, 0, MPI_COMM_WORLD, &reqr[0]);
    }
}

void sendBoundaries(long double ** phiValuesPart, int const * boundariesOffset, int procRankSend, int count, MPI_Request * reqs, int procCount){
    if(procRankSend != 0) {
        MPI_Isend(&phiValuesPart[0][boundariesOffset[0]], count, MPI_LONG_DOUBLE, procRankSend - 1, 0, MPI_COMM_WORLD, &reqs[0]);
    }

    if(procRankSend != procCount - 1) {
        MPI_Isend(&phiValuesPart[0][boundariesOffset[1]], count, MPI_LONG_DOUBLE, procRankSend + 1, 0, MPI_COMM_WORLD, &reqs[1]);
    }
}

void cutStartPhiValuesMatrix(long double ** phiValuesPart, long double * phiValuesSolid, const int& NxPart, const int& NyPart, const int& NzPart){
    // Nx MOD procCount == 0 | Ny MOD procCount == 0 | Ny MOD procCount == 0
    MPI_Scatter(phiValuesSolid, NxPart*NyPart*NzPart, MPI_LONG_DOUBLE, phiValuesPart[0], NxPart*NyPart*NzPart, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);

    // Start Phi: phi[M] values and phi[M+1] values
    memcpy(phiValuesPart[1], phiValuesPart[0], NxPart*NyPart*NzPart*sizeof(long double));
}

int isVld(int value, int limit, int addition){
    return (value + addition < 0 || value + addition >= limit) ? 0 : 1;
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

void fillStartPhiValues(long double * phiSolid, int Nx, int Ny, int Nz, const long double& Hx, const long double& Hy, const long double& Hz){
    memset(phiSolid, 0, Nx*Ny*Nz);

    long double x, y, z;

    for (int s = 0; s < 2; ++s) {
        // phi(xi, yj, 1), phi(xi, yj, -1)
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                calculatePhiArguments(i, j, s*(Nz-1), x, y, z, Hx, Hy, Hz);
                phiSolid[i + j*Nx + s*(Nz-1)*Nx*Ny] = phi(x, y, z);
            }
        }

        // phi(xi, 1, zk), phi(xi, -1, zk)
        for (int k = 0; k < Nz; ++k) {
            for (int i = 0; i < Nx; ++i) {
                calculatePhiArguments(i, s*(Ny-1), k, x, y, z, Hx, Hy, Hz);
                phiSolid[i +s*(Ny-1)*Nx + k*Nx*Ny] = phi(x, y, z);
            }
        }

        // phi(1, yj, zk), phi(-1, yj, zk)
        for (int k = 0; k < Nz; ++k) {
            for (int j = 0; j < Ny; ++j) {
                calculatePhiArguments(s*(Nx-1), j, k, x, y, z, Hx, Hy, Hz);
                phiSolid[s*(Nx-1) + j*Nx + k*Nx*Ny] = phi(x, y, z);
            }
        }
    }

}

void calculateMaxDiffLocalAndDeltaLocal(long double ** phiValuesPart, long double precisePhiValue, int NxPart, int NyPart,
                             long double& maxDiffLocal, long double& deltaLocal, int i, int j, int k){

    if(abs(phiValuesPart[1][i + j*NxPart + k*NxPart*NyPart] - phiValuesPart[0][i + j*NxPart + k*NxPart*NyPart]) > maxDiffLocal){
        maxDiffLocal = abs(phiValuesPart[1][i + j*NxPart + k*NxPart*NyPart] - phiValuesPart[0][i + j*NxPart + k*NxPart*NyPart]);
    }
    if(abs(phiValuesPart[0][i + j*NxPart + k*NxPart*NyPart] - precisePhiValue) > deltaLocal){
        deltaLocal = abs(phiValuesPart[0][i + j*NxPart + k*NxPart*NyPart] - precisePhiValue);
    }

}

void calculateBoundaries(long double ** phiValuesPart, int NxPart, int NyPart, int NzPart, const long double& Hx, const long double& Hy, const long double& Hz,
                         long double &maxDiffLocal, long double &deltaLocal, long double ** boundaries, int procRank, int procCount, int * H, int * L){

    long double a[4], b[4];
    long double precisePhiValue[2];
    long double x, y, z;
    long double denominator = ( 2.0/(Hx*Hx) + 2.0/(Hy*Hy) + 2.0/(Hz*Hz) + A );

    for (int j = 1; j < NyPart - 1; ++j) {
        for (int i = 1; i < NxPart - 1; ++i) {
            if( procRank != 0 ){
                a[0] = phiValuesPart[0][(i+1) + j*NxPart + H[0]*NxPart*NyPart] + phiValuesPart[0][(i-1) + j*NxPart + H[0]*NxPart*NyPart];
                a[1] = phiValuesPart[0][i + (j+1)*NxPart + H[0]*NxPart*NyPart] + phiValuesPart[0][i + (j-1)*NxPart + H[0]*NxPart*NyPart];
                a[2] = phiValuesPart[0][i + j*NxPart + (H[0]+1)*NxPart*NyPart] + boundaries[0][i + j*NxPart];

                calculatePhiArguments(i, j, L[0], x, y, z, Hx, Hy, Hz);
                precisePhiValue[0] = phi(x, y, z);
                a[3] = -rho(precisePhiValue[0]);

                phiValuesPart[1][i + j*NxPart + H[0]*NxPart*NyPart] = (a[0]/(Hx*Hx) + a[1]/(Hy*Hy) + a[2]/(Hz*Hz) + a[3]) / denominator;

                calculateMaxDiffLocalAndDeltaLocal(phiValuesPart, precisePhiValue[0], NxPart, NyPart, maxDiffLocal, deltaLocal, i, j, H[0]);
            }

            if( procRank != procCount-1 ){
                b[0] = phiValuesPart[0][(i+1) + j*NxPart + H[1]*NxPart*NyPart] + phiValuesPart[0][(i-1) + j*NxPart + H[1]*NxPart*NyPart];
                b[1] = phiValuesPart[0][i + (j+1)*NxPart + H[1]*NxPart*NyPart] + phiValuesPart[0][i + (j-1)*NxPart + H[1]*NxPart*NyPart];
                b[2] = boundaries[1][i + j*NxPart] + phiValuesPart[0][i + j*NxPart + (H[1]-1)*NxPart*NyPart];

                calculatePhiArguments(i, j, L[1], x, y, z, Hx, Hy, Hz);
                precisePhiValue[1] = phi(x, y, z);
                b[3] = -rho(precisePhiValue[1]);

                phiValuesPart[1][i + j*NxPart + H[1]*NxPart*NyPart] = (b[0]/(Hx*Hx) + b[1]/(Hy*Hy) + b[2]/(Hz*Hz) + b[3]) / denominator;

                calculateMaxDiffLocalAndDeltaLocal(phiValuesPart, precisePhiValue[1], NxPart, NyPart, maxDiffLocal, deltaLocal, i, j, H[1]);
            }

        }
    }

}

void calculateMPlusOnePhiValue(long double ** phiValuesPart, int NxPart, int NyPart, int NzPart, const long double& Hx, const long double& Hy, const long double& Hz,
                               long double &maxDiffLocal, long double &deltaLocal, long double ** boundaries, int procRank, int procCount, int Nz){
    int NzAddition = NzPart % 2;
    int medValZ = NzPart / 2;
    int K[2];

    MPI_Request reqs[2] = {0, 0};
    MPI_Request reqr[2] = {0, 0};

    long double a[4];
    long double b[4];
    long double precisePhiValue[2];

    maxDiffLocal = 0;
    deltaLocal = 0;

    long double denominator = ( 2.0/(Hx*Hx) + 2.0/(Hy*Hy) + 2.0/(Hz*Hz) + A );

    long double x, y, z;

    // Boundaries 0 - Upper, 1 - Lower
    int H[2] = {0, NzPart-1};
    int L[2] = {procRank*NzPart, (procRank+1)*NzPart - 1};
    int boundariesOffsetS[2] = { 0, (NzPart-1)*NxPart*NyPart };

    sendBoundaries(phiValuesPart, boundariesOffsetS, procRank, NxPart*NyPart, reqs, procCount);
    getBoundaries(boundaries, procRank, NxPart*NyPart, reqr, procCount);

    for(int k = 0; k < medValZ + NzAddition - 1; ++k){
        for (int j = 1; j < NyPart - 1; ++j) {
            for (int i = 1; i < NxPart - 1; ++i) {
                K[0] = medValZ+k; // NzPart MOD 2 == 1 -> K[0] = medValZ + k | NzPart MOD 2 == 0 -> K[0] = medValZ + k

                // phi[M]_{i+1, j, k} + phi[M]_{i-1, j, k}
                a[0] = phiValuesPart[0][(i+1) + j*NxPart + K[0]*NxPart*NyPart] + phiValuesPart[0][(i-1) + j*NxPart + K[0]*NxPart*NyPart];

                // phi[M]_{i, j+1, k} + phi[M]_{i, j-1, k}
                a[1] = phiValuesPart[0][i + (j+1)*NxPart + K[0]*NxPart*NyPart] + phiValuesPart[0][i + (j-1)*NxPart + K[0]*NxPart*NyPart];

                // phi[M]_{i, j, k+1} + phi[M]_{i, j, k-1}
                a[2] = phiValuesPart[0][i + j*NxPart + (K[0]+1)*NxPart*NyPart] + phiValuesPart[0][i + j*NxPart + (K[0]-1)*NxPart*NyPart];

                // rho_{i, j, k}
                calculatePhiArguments(i, j, K[0], x, y, z, Hx, Hy, Hz);
                precisePhiValue[0] = phi(x, y, z);
                a[3] = -rho(precisePhiValue[0]);

                // phi[M+1]_{i, j, k}
                phiValuesPart[1][i + j*NxPart + K[0]*NxPart*NyPart] = (a[0]/(Hx*Hx) + a[1]/(Hy*Hy) + a[2]/(Hz*Hz) + a[3]) / denominator;

                if(isnan(phiValuesPart[1][i + j*NxPart + K[0]*NxPart*NyPart]) && procRank == 2){
                    cout << i << " " << j << " " << K[0] << "\n";
                }

                calculateMaxDiffLocalAndDeltaLocal(phiValuesPart, precisePhiValue[0], NxPart, NyPart, maxDiffLocal, deltaLocal, i, j, K[0]);

                K[1] = medValZ-1+NzAddition-k; // NzPart MOD 2 == 1 -> K[0] = medValZ - k | NzPart MOD 2 == 0 -> K[0] = medValZ - 1 - k
                // phi[M]_{i+1, j, k} + phi[M]_{i-1, j, k}
                b[0] = phiValuesPart[0][(i+1) + j*NxPart + K[1]*NxPart*NyPart] + phiValuesPart[0][(i-1) + j*NxPart + K[1]*NxPart*NyPart];

                // phi[M]_{i, j+1, k} + phi[M]_{i, j-1, k}
                b[1] = phiValuesPart[0][i + (j+1)*NxPart + K[1]*NxPart*NyPart] + phiValuesPart[0][i + (j-1)*NxPart + K[1]*NxPart*NyPart];

                // phi[M]_{i, j, k+1} + phi[M]_{i, j, k-1}
                b[2] = phiValuesPart[0][i + j*NxPart + (K[1]+1)*NxPart*NyPart] + phiValuesPart[0][i + j*NxPart + (K[1]-1)*NxPart*NyPart];

                // rho_{i, j, k}
                calculatePhiArguments(i, j, K[1], x, y, z, Hx, Hy, Hz);
                precisePhiValue[1] = phi(x, y, z);
                b[3] = -rho(precisePhiValue[1]);

                // phi[M+1]_{i, j, k}
                phiValuesPart[1][i + j*NxPart + K[1]*NxPart*NyPart] = (b[0]/(Hx*Hx) + b[1]/(Hy*Hy) + b[2]/(Hz*Hz) + b[3]) / denominator;

                // deb
                if(isnan(phiValuesPart[1][i + j*NxPart + K[1]*NxPart*NyPart]) && procRank == 2){
                    cout << i << " " << j << " " << K[1] << "\n";
                }

                calculateMaxDiffLocalAndDeltaLocal(phiValuesPart, precisePhiValue[1], NxPart, NyPart, maxDiffLocal, deltaLocal, i, j, K[1]);
            }
        }
    }

    waitMessages(reqr, reqs, procRank, procCount);

    calculateBoundaries(phiValuesPart, NxPart, NyPart, NzPart, Hx, Hy, Hz, maxDiffLocal, deltaLocal, boundaries, procRank, procCount, H, L);

    /*if(procRank == 0){
        for (int j = 0; j < NyPart; ++j) {
            for (int i = 0; i < NxPart; ++i) {
                cout << boundaries[1][i + j*NxPart] << " ";
            }
            cout << endl;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(procRank == 1){
        cout << endl;

        for (int j = 0; j < NyPart; ++j) {
            for (int i = 0; i < NxPart; ++i) {
                cout << phiValuesPart[0][i + j*NxPart + boundariesOffsetS[0]] << " ";
            }
            cout << endl;

        }
    }*/

    memcpy(phiValuesPart[0], phiValuesPart[1], NxPart*NyPart*NzPart*sizeof(long double));
}

int main(int argc, char* argv[]){
    long double * phiValuesPart[2];
    long double * boundary[2]; // 0 - upper boundary, 1 - lower boundary
    long double * phiValuesSolid;
    long double maxDiff;
    long double delta;

    MPI_Init(&argc,&argv);

    int procCount;
    int procRank;

    MPI_Comm_size(MPI_COMM_WORLD, &procCount);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

    int Nx = START_VALUE_N;
    int Ny = START_VALUE_N;
    int Nz = START_VALUE_N;

    if( Nz%procCount != 0 ){
        cout << "Nz MOD processes count != 0. Change Ni or processes count\n";
        return EXIT_SUCCESS;
    }

    int NxPart = Nx / 1;
    int NyPart = Ny / 1;
    int NzPart = Nz / procCount;

    long double Hx = DX / ((long double) Nx - 1.0);
    long double Hy = DY / ((long double) Ny - 1.0);
    long double Hz = DZ / ((long double) Nz - 1.0);

    long double maxDiffLocal;
    long double deltaLocal;


    for(auto & i : phiValuesPart) {
        i = new long double[NzPart*NyPart*NxPart];
    }

    for (int i = 0; i < 2; ++i) {
        boundary[i] = new long double[NyPart*NxPart];
        phiValuesPart[i] = new long double[NzPart*NyPart*NxPart];
    }

    phiValuesSolid = new long double[Nx*Ny*Nz];

    if(procRank == 0){
        fillStartPhiValues(phiValuesSolid, Nx, Ny, Nz, Hx, Hy, Hz);
    }

    cutStartPhiValuesMatrix(phiValuesPart, phiValuesSolid, NxPart, NyPart, NzPart);

    do{
        calculateMPlusOnePhiValue(phiValuesPart, NxPart, NyPart, NzPart, Hx, Hy, Hz, maxDiffLocal, deltaLocal, boundary, procRank, procCount, Nz);

        MPI_Allreduce(&maxDiffLocal, &maxDiff, 1, MPI_LONG_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&deltaLocal, &delta, 1, MPI_LONG_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        if(procRank == 0){
//            cout << 1 << " ";
//            cout << delta << endl;
//            cout << maxDiff << endl;
        }

        if(maxDiff < EPSILON){
            break;
        }

    } while (TRUE);

    if(procRank == 0){

        if(delta > EPSILON*MAX_DIFF_deltaLocal){
            cout << "delta is bigger than EPSILON*1e2.\n";
        } else {
            cout << "Function is found\n";
        }

        cout << delta << endl;

    }

}
