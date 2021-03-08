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

    int procSize, procRank;

    int repeats = 0;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

    int rowNum = N/procSize;
    int lastRowAdding = N%procSize;

    double* mProcRows;
    double* vecBPart;
    double* vecUPart;

    double* vecXRes;

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

    vecUPart = (double*) calloc(rowNum+lastRowAdding, sizeof(double));
    mProcRows = (double*) calloc(rowNumMod*N, sizeof(double));
    vecXRes = (double*) calloc(N, sizeof(double));

    initMatrixProcRows(rowNum, N, mProcRows, procRank, lastRowAdding);
    MPI_Scatterv(vecU, sendcounts, displs, MPI_DOUBLE, vecUPart, recvcounts[procRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    vecBPart = mulMatrixAndVector(rowNum, lastRowAdding, rowNumMod, N, mProcRows, vecUPart, recvcounts); // init vector B part

    //printVector(vecBPart, rowNumMod, procRank);

    rPart[0] = (double*) calloc(rowNum, sizeof(double));
    xPart[0] = (double*) calloc(rowNum, sizeof(double));
    zPart[0] = (double*) calloc(rowNum, sizeof(double));

    std::memcpy(rPart[0], vecBPart, sizeof(double)*rowNumMod); // r0 = b - Ax0, где x0 - нулевой вектор
    std::memcpy(zPart[0], rPart[0], sizeof(double)*rowNumMod);  // z0 = r0

    // printProcRows(mProcRows, rowNumMod, N);

   /* double * temp = subVectorAndVector(rowNumMod, vecBPart, vecBPart);
    double * tempres = (double*) calloc(N, sizeof(double));
    MPI_Allgatherv(temp, sendcounts[procRank], MPI_DOUBLE, tempres, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
    printVector(tempres, N, procRank);*/

    // std::cout << scalarVectorAndVector(rowNumMod, rPart[0], rPart[0]);
    double* temp[] = {NULL, NULL, NULL, NULL};

    while (1){

        temp[0] = mulMatrixAndVector(rowNum, lastRowAdding, rowNumMod, N, mProcRows, zPart[0], recvcounts);
        alpha[1] = scalarVectorAndVector(rowNumMod, rPart[0], rPart[0])
                    / scalarVectorAndVector(rowNumMod, temp[0], zPart[0]);          // "alpha(k+1) = (r(k), r(k)) / (Az(k), z(k))"

        temp[1] = mulVectorAndScalar(rowNumMod, alpha[1], zPart[0]);
        xPart[1] = sumVectorAndVector(rowNumMod, xPart[0], temp[1]);                // "x(k+1) = x(k) + alpha(k+1)z(k)"

        temp[2] = mulVectorAndScalar(rowNumMod, alpha[1], temp[0]);
        rPart[1] = subVectorAndVector(rowNumMod, rPart[0], temp[2]);                // "r(k+1) = r(k) - alpha(k+1)Az(k)"

        beta[1] = scalarVectorAndVector(rowNumMod, rPart[1], rPart[1])
                    / scalarVectorAndVector(rowNumMod, rPart[0], rPart[0]);         // "b(k+1) = (r(k+1), r(k+1)) / (r(k), r(k))"

        temp[3] = mulVectorAndScalar(rowNumMod, beta[1], zPart[0]);
        zPart[1] = sumVectorAndVector(rowNumMod, rPart[1], temp[3]);                // "z(k+1) = r(k+1) + beta(k+1)z(k)"

        for (int ui = 0; ui < 4; ++ui) {
            free(temp[ui]);
        }

        if(procRank == 0){
            repeats++;
        }

        if( (vectorLength(rowNumMod, rPart[0]) / vectorLength(rowNumMod, vecBPart) ) < EPSILON){    // |r(k)| / |b| < EPSILON
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

    MPI_Allgatherv(xPart[0], sendcounts[procRank], MPI_DOUBLE, vecXRes, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);

    if(procRank == 0){
        printVector(vecU, N, procRank);
        printVector(vecXRes, N, procRank);

        std::cout << "Repeats in total: " << repeats << "\n";
    }

    MPI_Finalize();

    return 0;
}