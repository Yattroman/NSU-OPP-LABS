#include <mpi.h>
#include <iostream>
#include <cstring>
#include "mv/mvOperations_3rd.h"
#include "mv/mvInit_3rd.h"

#define EPSILON 1e-6

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
    int growStatus = 0;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

    int rowNum = N/procSize;
    int lastRowAdding = N%procSize;

    double* mProcRows;
    double* vecBPart;
    double* vecUPart;

    double* vecXRes = (double*) calloc(N, sizeof(double));
    double* vecU = (double*) calloc(N, sizeof(double));

    double* rPart[] = {NULL, NULL};
    double* zPart[] = {NULL, NULL};
    double* xPart[] = {NULL, NULL};

    double alpha[] = {0, 0};
    double beta[] = {0, 0};

    int* displs = (int*) calloc(procSize, sizeof(int));
    int* recvcounts = (int*) calloc(procSize, sizeof(int));
    int* sendcounts = (int*) calloc(procSize, sizeof(int));

    int rowNumMod = (procRank == 0) ? rowNum+lastRowAdding : rowNum;

    fillDisplsAndRecvcountsTables(displs, recvcounts, sendcounts, rowNum, lastRowAdding, procSize);

    initVectorU(N, vecU);
    vecUPart = (double*) calloc(N, sizeof(double));
    mProcRows = (double*) calloc(rowNumMod*N, sizeof(double));

    initMatrixProcRows(rowNum, N, mProcRows, procRank, lastRowAdding);
    MPI_Scatterv(vecU, sendcounts, displs, MPI_DOUBLE, vecUPart, recvcounts[procRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    vecBPart = (double*) calloc(N, sizeof(double));
    mulMatrixAndVector(rowNumMod, N, mProcRows, vecUPart, recvcounts, vecBPart); // init vector B part

    for (int i = 0; i < 2; ++i) {
        rPart[i] = (double*) calloc(rowNumMod, sizeof(double));
        xPart[i] = (double*) calloc(rowNumMod, sizeof(double));
        zPart[i] = (double*) calloc(rowNumMod, sizeof(double));
    }

    std::memcpy(rPart[0], vecBPart, sizeof(double)*rowNumMod); // r0 = b - Ax0, где x0 - нулевой вектор
    std::memcpy(zPart[0], rPart[0], sizeof(double)*rowNumMod);  // z0 = r0

    double* temp[4];

    for (int ui = 0; ui < 4; ++ui) {
        temp[ui] = (double*) calloc(N, sizeof(double));
    }

    while (1){

        mulMatrixAndVector(rowNumMod, N, mProcRows, zPart[0], recvcounts, temp[0]);

        alpha[1] = scalarVectorAndVector(rowNumMod, rPart[0], rPart[0])
                    / scalarVectorAndVector(rowNumMod, temp[0], zPart[0]);          // "alpha(k+1) = (r(k), r(k)) / (Az(k), z(k))"

        mulVectorAndScalar(rowNumMod, alpha[1], zPart[0], temp[1]);
        sumVectorAndVector(rowNumMod, xPart[0], temp[1], xPart[1]);                // "x(k+1) = x(k) + alpha(k+1)z(k)"

        mulVectorAndScalar(rowNumMod, alpha[1], temp[0], temp[2]);
        subVectorAndVector(rowNumMod, rPart[0], temp[2], rPart[1]);                // "r(k+1) = r(k) - alpha(k+1)Az(k)"

        beta[1] = scalarVectorAndVector(rowNumMod, rPart[1], rPart[1])
                    / scalarVectorAndVector(rowNumMod, rPart[0], rPart[0]);         // "b(k+1) = (r(k+1), r(k+1)) / (r(k), r(k))"

        mulVectorAndScalar(rowNumMod, beta[1], zPart[0], temp[3]);
        sumVectorAndVector(rowNumMod, rPart[1], temp[3], zPart[1]);                // "z(k+1) = r(k+1) + beta(k+1)z(k)"

        if(procRank == 0){
            repeats++;
            //std::cout << repeats;
        }

        if( (vectorLength(rowNumMod, rPart[0]) / vectorLength(rowNumMod, vecBPart) ) < EPSILON){    // |r(k)| / |b| < EPSILON
            break;
        }

        /*MPI_Allgatherv(xPart[1], sendcounts[procRank], MPI_DOUBLE, vecXRes, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
        printVector(vecXRes, N, procRank);*/

        if(growStatus > 10){
            break;
        } else if( vectorLength(rowNumMod, rPart[0]) < vectorLength(rowNumMod, rPart[1]) ){
            growStatus++;
        } else if( vectorLength(rowNumMod, rPart[0]) > vectorLength(rowNumMod, rPart[1]) ){
            growStatus = 0;
        }

        std::memcpy(xPart[0], xPart[1], sizeof(double)*rowNumMod);
        std::memcpy(rPart[0], rPart[1], sizeof(double)*rowNumMod);
        std::memcpy(zPart[0], zPart[1], sizeof(double)*rowNumMod);

        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Allgatherv(xPart[1], sendcounts[procRank], MPI_DOUBLE, vecXRes, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);

    if(procRank == 0 && growStatus <= 10){
        printVector(vecU, N, procRank);
        printVector(vecXRes, N, procRank);
        std::cout << "Repeats in total: " << repeats << "\n";
    } else if(procRank == 0 && growStatus > 10) {
        std::cout << "There are no roots!\n";
    }

    MPI_Finalize();

    for (int ui = 0; ui < 4; ++ui) {
        free(temp[ui]);
    }

    for (int i = 0; i < 2; ++i) {
        free(xPart[i]);
        free(zPart[i]);
        free(rPart[i]);
    }

    free(vecBPart);
    free(vecXRes);
    free(vecUPart);
    free(mProcRows);
    free(vecU);

    free(displs);
    free(sendcounts);
    free(recvcounts);

    return EXIT_SUCCESS;
}