#include <mpi.h>
#include <iostream>
#include <unistd.h>
#include <cstring>
#include "mv/mvOperations_3rd.h"
#include "mv/mvInit_3rd.h"

#define EPSILON 0.0001

void distributeData(){ }

void fillDisplsAndRecvcountsTables(int* displs, int* recvcounts, int* sendcounts, int rowNum, int lastRowAdding, int N, int procSize){
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

    double* Azk = (double*) calloc(N, sizeof(double));

    double* rPart[] = {NULL, NULL};
    double* zPart[] = {NULL, NULL};
    double* xPart[] = {NULL, NULL};

    double alpha[] = {0, 0};
    double beta[] = {0, 0};

    int* displs = (int*) calloc(procSize, sizeof(int));
    int* recvcounts = (int*) calloc(procSize, sizeof(int));
    int* sendcounts = (int*) calloc(procSize, sizeof(int));

    initVectorU(N, vecU);

    fillDisplsAndRecvcountsTables(displs, recvcounts, sendcounts, rowNum, lastRowAdding, N, procSize);

    if(procRank == 0){
        mProcRows = (double*) calloc((rowNum+lastRowAdding)*N, sizeof(double));
        vecBPart = (double*) calloc(rowNum+lastRowAdding, sizeof(double));
        vecXPart = (double*) calloc(rowNum+lastRowAdding, sizeof(double));

        initMatrixProcRows(rowNum, N, mProcRows, procRank, lastRowAdding);
        vecBPart = mulMatrixAndVector(rowNum+lastRowAdding, N, mProcRows, vecU);
    } else {
        mProcRows = (double*) calloc(rowNum*N, sizeof(double));
        vecBPart = (double*) calloc(rowNum, sizeof(double));
        vecXPart = (double*) calloc(rowNum, sizeof(double));

        initMatrixProcRows(rowNum, N, mProcRows, procRank, lastRowAdding);
        vecBPart = mulMatrixAndVector(rowNum, N, mProcRows, vecU);
    }

    distributeData();

                  // r0 = b - Ax0, где x0 - нулевой вектор
                   // z0 = r0

    int rowNumMod = (procRank == 0) ? rowNum+lastRowAdding : rowNum;

    while (1){
        double* temp[] = {NULL, NULL, NULL, NULL};

    }

    MPI_Finalize();

    return 0;
}