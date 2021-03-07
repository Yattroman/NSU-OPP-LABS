#include <mpi.h>
#include <iostream>
#include <unistd.h>
#include <cstring>
#include "mv/mvOperations_3rd.h"
#include "mv/mvInit_3rd.h"

#define EPSILON 0.0001

void distributeData(){ }

void fillDisplsAndRecvcountsTables(int* displs, int* recvcounts, int* sendcounts, int rowNum, int lastRowAdding, int procSize){
    for (int i = 1; i < procSize; ++i) {
        displs[i] = (i+lastRowAdding)*rowNum;
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

    int  procSize, procRank;

    int repeats = 0;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

    int rowNum = N/procSize;
    int lastRowAdding = N%procSize;

    double* mProcRows;
    double* vecBPart;
    double* vecXPart;

    double* vecU = (double*) calloc(N, sizeof(double));

    double* rPart[] = {NULL, NULL};
    double* zPart[] = {NULL, NULL};
    double* xPart[] = {NULL, NULL};

    double alpha[] = {0, 0};
    double beta[] = {0, 0};

    int* displs = (int*) calloc(procSize, sizeof(int));
    int* recvcounts = (int*) calloc(procSize, sizeof(int));
    int* sendcounts = (int*) calloc(procSize, sizeof(int));

    initVectorU(N, vecU);

    fillDisplsAndRecvcountsTables(displs, recvcounts, sendcounts, rowNum, lastRowAdding, procSize);

    int rowNumMod = (procRank == 0) ? rowNum+lastRowAdding : rowNum;

    mProcRows = (double*) calloc(rowNumMod*N, sizeof(double));
    vecBPart = (double*) calloc(rowNum+lastRowAdding, sizeof(double));
    vecXPart = (double*) calloc(rowNum+lastRowAdding, sizeof(double));

    initMatrixProcRows(rowNum, N, mProcRows, procRank, lastRowAdding);
    vecBPart = mulMatrixAndVector(rowNum, lastRowAdding, rowNumMod, N, mProcRows, vecU, recvcounts); // init vector B part

    //printVector(vecBPart, rowNumMod, procRank);

    rPart[0] = (double*) calloc(rowNum, sizeof(double));
    xPart[0] = (double*) calloc(rowNum, sizeof(double));
    zPart[0] = (double*) calloc(rowNum, sizeof(double));

    std::memcpy(rPart[0], vecBPart, sizeof(double)*rowNumMod); // r0 = b - Ax0, где x0 - нулевой вектор
    std::memcpy(zPart[0], rPart[0], sizeof(double)*rowNumMod);  // z0 = r0

    //printProcRows(mProcRows, rowNumMod, N);

    //std::cout << scalarVectorAndVector(rowNumMod, rPart[0], rPart[0]);

    while (1){
        double* temp[] = {NULL, NULL, NULL, NULL};

        temp[0] = mulMatrixAndVector(rowNum, lastRowAdding, rowNumMod, N, mProcRows, zPart[0], recvcounts);
        alpha[1] = scalarVectorAndVector(rowNumMod, rPart[0], rPart[0])
                    / scalarVectorAndVector(rowNumMod, temp[0], zPart[0]);

        temp[1] = mulVectorAndScalar(rowNumMod, alpha[1], zPart[0]);
        xPart[1] = sumVectorAndVector(rowNumMod, xPart[0], temp[1]);

        temp[2] = mulVectorAndScalar(rowNumMod, alpha[1], temp[0]);
        rPart[1] = subVectorAndVector(rowNumMod, rPart[0], temp[2]);

        beta[1] = scalarVectorAndVector(rowNumMod, rPart[1], rPart[1])
                    / scalarVectorAndVector(rowNumMod, rPart[0], rPart[0]);

        temp[3] = mulVectorAndScalar(rowNumMod, beta[1], zPart[0]);
        zPart[1] = sumVectorAndVector(rowNumMod, rPart[1], temp[3]);

        for (int ui = 0; ui < 4; ++ui) {
            free(temp[ui]);
        }

        if(procRank == 0){
            repeats++;
        }

        if( (vectorLength(rowNumMod, rPart[0]) / vectorLength(rowNumMod, vecBPart) ) < EPSILON){    // |r(k+1)| / |b| < EPSILON
            free(xPart[0]);
            free(zPart[0]);
            free(rPart[0]);
            break;
        }

        free(xPart[0]);
        free(zPart[0]);
        free(rPart[0]);
        xPart[0] = xPart[1];
        zPart[0] = zPart[1];
        rPart[0] = rPart[1];
    }

    if(procRank == 0){
        printVector(vecU, N, procRank);

        std::cout << "Repeats in total: " << repeats << "\n";
    }
    printVector(xPart[1], N, procRank);

    MPI_Finalize();

    return 0;
}