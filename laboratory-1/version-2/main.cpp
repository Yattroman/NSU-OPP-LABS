#include <mpi.h>
#include <iostream>
#include <unistd.h>
#include <cstring>
#include "mv/mvOperations_2nd.h"
#include "mv/mvInit_2nd.h"

#define EPSILON 1e-5

void distributeData(double* vectorX, double* vectorB, int* sendcounts, int* recvcounts, int* displs, int procRank, double* vecBPart, int N){
    MPI_Bcast(vectorX, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Allgatherv(vecBPart, sendcounts[procRank], MPI_DOUBLE, vectorB, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
}

void fillDisplsAndRecvcountsTables(int* displs, int* recvcounts, int* sendcounts, int rowNum, int lastRowAdding, int procSize){
    for (int i = 1; i < procSize; ++i) {
        displs[i] = i*rowNum+lastRowAdding;
        recvcounts[i] = rowNum;
        sendcounts[i] = rowNum;
    }
    displs[0] = 0;
    sendcounts[0] = rowNum+lastRowAdding;
    recvcounts[0] = rowNum+lastRowAdding;
}

int main(int argc, char** argv)
{
    int N = atoi(argv[1]);

    int procSize, procRank;

    int repeats = 0;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

    int rowNum = N/procSize;
    int lastRowAdding = N%procSize;

    int growStatus = 0;

    double* mProcRows;
    double* vecBPart;

    double* r[] = {NULL, NULL};
    double* z[] = {NULL, NULL};
    double* x[] = {NULL, NULL};

    double* vecU = (double*) calloc(N, sizeof(double));
    double* vecB = (double*) calloc(N, sizeof(double));

    for (int i = 0; i < 2; ++i) {
        r[i] = (double*) calloc(N, sizeof(double));
        z[i] = (double*) calloc(N, sizeof(double));
        x[i] = (double*) calloc(N, sizeof(double));
    }

    double* AzkPart = (double*) calloc(N, sizeof(double));

    double alpha[] = {0, 0};
    double beta[] = {0, 0};

    int* displs = (int*) calloc(procSize, sizeof(int));
    int* recvcounts = (int*) calloc(procSize, sizeof(int));
    int* sendcounts = (int*) calloc(procSize, sizeof(int));

    initVectorU(N, vecU);

    fillDisplsAndRecvcountsTables(displs, recvcounts, sendcounts, rowNum, lastRowAdding, procSize);

    int rowNumMod = (procRank == 0) ? rowNum+lastRowAdding : rowNum;

    mProcRows = (double*) calloc(rowNumMod*N, sizeof(double));
    vecBPart = (double*) calloc(rowNumMod, sizeof(double));

    initMatrixProcRows(rowNum, N, mProcRows, procRank, lastRowAdding);
    mulMatrixAndVector(rowNumMod, N, mProcRows, vecU, vecBPart);

    distributeData(x[0], vecB, sendcounts, recvcounts, displs, procRank, vecBPart, N);

    std::memcpy( r[0], vecB, N*sizeof(double) );                    // r0 = b - Ax0, где x0 - нулевой вектор
    std::memcpy( z[0], r[0], N*sizeof(double) );                // z0 = r0

    double* temp[4];

    for (int ui = 0; ui < 4; ++ui) {
        if(ui == 0)
            temp[ui] = (double*) calloc(rowNumMod, sizeof(double));
        else
            temp[ui] = (double*) calloc(N, sizeof(double));
    }

    while (1){

        mulMatrixAndVector(rowNumMod, N, mProcRows, z[0], temp[0]);                      // Az(k)

        MPI_Allgatherv(temp[0], sendcounts[procRank], MPI_DOUBLE, AzkPart, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);

        alpha[1] = scalarVectorAndVector(N, r[0], r[0])
                   / scalarVectorAndVector(N, AzkPart, z[0]);               // alpha(k+1) = (r(k), r(k)) / (Az(k), z(k))

        mulVectorAndScalar(N, alpha[1], z[0], temp[1]);
        sumVectorAndVector(N, x[0], temp[1], x[1]);            // x(k+1) = x(k) + alpha(k+1)z(k)

        mulVectorAndScalar(N, alpha[1], AzkPart, temp[2]);
        subVectorAndVector(N, r[0], temp[2], r[1]);            // r(k+1) = r(k) - alpha(k+1)Az(k)

        beta[1] = scalarVectorAndVector(N, r[1], r[1])
                  / scalarVectorAndVector(N, r[0], r[0]);         // b(k+1) = (r(k+1), r(k+1)) / (r(k), r(k))

        mulVectorAndScalar(N, beta[1], z[0], temp[3]);
        sumVectorAndVector(N, r[1], temp[3], z[1]);            // z(k+1) = r(k+1) + beta(k+1)z(k)

        repeats++;

        if( (vectorLength(N, r[0]) / vectorLength(N, vecB) ) < EPSILON){    // |r(k)| / |b| < EPSILON
            break;
        }

        if(growStatus > 10){
            break;
        } else if( vectorLength(N, r[0]) < vectorLength(N, r[1]) ){
            growStatus++;
        } else if( vectorLength(N, r[0]) > vectorLength(N, r[1]) ){
            growStatus = 0;
        }

        std::memcpy(x[0], x[1], N*sizeof(double));
        std::memcpy(r[0], r[1], N*sizeof(double));
        std::memcpy(z[0], z[1], N*sizeof(double));

    }

    if(procRank == 0 && growStatus <= 10){
        printVector(vecU, N, procRank);
        printVector(x[1], N, procRank);
        std::cout << "Repeats in total: " << repeats << "\n";
    } else if(procRank == 0 && growStatus > 10) {
        std::cout << "There are no roots!\n";
    }

    for (size_t ui = 0; ui < 4; ++ui) {
        free(temp[ui]);
    }

    for (int i = 0; i < 2; ++i) {
        free(x[i]);
        free(z[i]);
        free(r[i]);
    }

    MPI_Finalize();

    return 0;
}