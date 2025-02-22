#include <mpi.h>
#include <iostream>
#include <cstring>
#include "mv/mvOperations_l1_3rd.h"
#include "mv/mvInit_l1_3rd.h"

#define EPSILON 1e-10

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
    int procSize, procRank;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

    int N = atoi(argv[1]);

    int repeats = 0;
    int growStatus = 0;
    double startTime = MPI_Wtime();

    int rowNum = N/procSize;
    int lastRowAdding = N%procSize;
    int rowNumMod = (procRank == 0) ? rowNum+lastRowAdding : rowNum;

    double* mProcRows = (double*) calloc(N*rowNumMod, sizeof(double));
    double* vecBPart = (double*) calloc(rowNumMod, sizeof(double));
    // double* vecUPart = (double*) calloc(rowNumMod+lastRowAdding, sizeof(double)); FOR TEST

    double* vecXRes = (double*) calloc(N, sizeof(double));
    // double* vecU = (double*) calloc(N, sizeof(double)); FOR TEST

    double* rPart[] = {NULL, NULL};
    double* zPart[] = {NULL, NULL};
    double* xPart[] = {NULL, NULL};

    double alpha[] = {0, 0};
    double beta[] = {0, 0};

    int* displs = (int*) calloc(procSize, sizeof(int));
    int* recvcounts = (int*) calloc(procSize, sizeof(int));
    int* sendcounts = (int*) calloc(procSize, sizeof(int));

    fillDisplsAndRecvcountsTables(displs, recvcounts, sendcounts, rowNum, lastRowAdding, procSize);

    // initVectorU(N, vecU); FOR TEST
    initMatrixProcRows(rowNum, N, mProcRows, procRank, lastRowAdding, displs);
    // MPI_Scatterv(vecU, sendcounts, displs, MPI_DOUBLE, vecUPart, recvcounts[procRank], MPI_DOUBLE, 0, MPI_COMM_WORLD); FOR TEST
    // mulMatrixAndVector(rowNumMod, N, mProcRows, vecUPart, recvcounts, vecBPart); FOR TEST
    initVectorBPart(rowNumMod, vecBPart, displs, procRank);

    for (int i = 0; i < 2; ++i) {
        rPart[i] = (double*) calloc(rowNumMod, sizeof(double));
        xPart[i] = (double*) calloc(rowNumMod, sizeof(double));
        zPart[i] = (double*) calloc(rowNumMod, sizeof(double));
    }

    std::memcpy(rPart[0], vecBPart, sizeof(double)*rowNumMod); // r0 = b - Ax0, где x0 - нулевой вектор
    std::memcpy(zPart[0], rPart[0], sizeof(double)*rowNumMod);  // z0 = r0

    double* temp[4];

    for (int ui = 0; ui < 4; ++ui) {
        temp[ui] = (double*) malloc(rowNumMod*sizeof(double));
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

        repeats++;

        if( (vectorLength(rowNumMod, rPart[0]) / vectorLength(rowNumMod, vecBPart) ) < EPSILON){    // |r(k)| / |b| < EPSILON
            break;
        }

        if(growStatus > 10){
            break;
        } else if( vectorLength(rowNumMod, rPart[0]) < vectorLength(rowNumMod, rPart[1]) ){
            growStatus++;
        } else if( vectorLength(rowNumMod, rPart[0]) > vectorLength(rowNumMod, rPart[1]) ){
            growStatus = 0;
        }

        std::memmove(xPart[0], xPart[1], sizeof(double)*rowNumMod);
        std::memmove(rPart[0], rPart[1], sizeof(double)*rowNumMod);
        std::memmove(zPart[0], zPart[1], sizeof(double)*rowNumMod);

    }

    MPI_Allgatherv(xPart[0], sendcounts[procRank], MPI_DOUBLE, vecXRes, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);

    if(procRank == 0 && growStatus <= 10){
        // printVector(vecU, N, procRank); FOR TEST
        printVector(vecXRes, N, procRank);
        std::cout << "Repeats in total: " << repeats << "\n";
    } else if(procRank == 0 && growStatus > 10) {
        std::cout << "There are no roots!\n";
    }

    double endTime = MPI_Wtime();
    double minimalStartTime;
    double maximumEndTime;
    MPI_Reduce( &endTime, &maximumEndTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Reduce( &startTime, &minimalStartTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    if ( procRank == 0 ) {
        printf( "Total time spent in seconds id %f\n", maximumEndTime - minimalStartTime );
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
    // free(vecUPart); FOR TEST
    free(mProcRows);
    // free(vecU); FOR TEST

    free(displs);
    free(sendcounts);
    free(recvcounts);

    return EXIT_SUCCESS;
}