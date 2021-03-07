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
    vecBPart = mulMatrixAndVector(rowNumMod, N, mProcRows, vecU, recvcounts); // init vector B part

    printVector(vecBPart, rowNumMod, procRank);

    rPart[0] = (double*) calloc(rowNum, sizeof(double));
    xPart[0] = (double*) calloc(rowNum, sizeof(double));
    zPart[0] = (double*) calloc(rowNum, sizeof(double));

    std::memcpy(rPart[0], vecBPart, sizeof(double)*rowNumMod); // r0 = b - Ax0, где x0 - нулевой вектор
    std::memcpy(zPart[0], rPart[0], sizeof(double)*rowNumMod);  // z0 = r0

    //std::cout << scalarVectorAndVector(rowNumMod, rPart[0], rPart[0]);

    /*while (1){
        double* temp[] = {NULL, NULL, NULL, NULL};

        //temp[0] = mulMatrixAndVector(rowNumMod, mProcRows, zPart[0])
        //alpha[1] = scalarVectorAndVector(rowNumMod, rPart[0], rPart[0]) /
                                            scalarVectorAndVector(rowNumMod, temp[0], zPart[0]);
    }*/

    MPI_Finalize();

    return 0;
}