#include <mpi.h>
#include <iostream>
#include <unistd.h>
#include <cstring>
#include "mv/mvOperations_2nd.h"
#include "mv/mvInit_2nd.h"

#define EPSILON 0.0001

void testPartialResults(double* pProcResult, int procSize, int procRank, int rowNum) {
    for (int i=0; i<procSize; i++) {
        if (procRank == i) {
            printf("ProcRank = %d \n", procRank);
            printf("Part of result vector: \n");
            printVector(pProcResult, rowNum);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void distributeData(double* vectorX, double* vectorB, int N){
    MPI_Bcast(vectorX, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(vectorB, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void fillDisplsAndRecvcountsTables(int* displs, int* recvcounts, int* sendcounts, int rowNum, int lastRowAdding, int N, int procSize){
    for (int i = 0; i < procSize; ++i) {
        displs[i] = N*rowNum*i;
        recvcounts[i] = rowNum;
        sendcounts[i] = rowNum;
    }

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

    double* vecU = new double[N];
    double* x0 = new double[N];
    double* vecB = new double[N+N];
    double* r0; // r0 = b - Ax0, где x0 - нулевой вектор
    double* z0; // z0 = r0

    double* r[] = {r0, NULL};
    double* z[] = {z0, NULL};
    double* x[] = {x0, NULL};
    double alpha[] = {0, 0};
    double beta[] = {0, 0};

    int* displs = new int[procSize];
    int* recvcounts = new int[procSize];
    int* sendcounts = new int[procSize];

    initVectorU(N, vecU);

    fillDisplsAndRecvcountsTables(displs, recvcounts, sendcounts, rowNum, lastRowAdding, N, procSize);

    if(procRank == 0){
        mProcRows = new double[(rowNum+lastRowAdding)*N];
        vecBPart = new double[rowNum+lastRowAdding];

        initVectorX(N, x0);

        initMatrixProcRows(rowNum, N, mProcRows, procRank, lastRowAdding);
        vecBPart = mulMatrixAndVector(rowNum+lastRowAdding, N, mProcRows, vecU);

        printVector(vecBPart, rowNum+lastRowAdding);
    } else {
        mProcRows = new double[rowNum*N];
        vecBPart = new double[rowNum];

        initMatrixProcRows(rowNum, N, mProcRows, procRank, lastRowAdding);
        vecBPart = mulMatrixAndVector(rowNum, N, mProcRows, vecU);
        printVector(vecBPart, rowNum);
    }

    MPI_Gatherv(vecBPart, sendcounts[procRank], MPI_DOUBLE, vecB, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //initVectorBPart()

    //distributeData(x0, vecB, N);

    //r0 = vecB; // r0 = b - Ax0, где x0 - нулевой вектор
    //memcpy(z0, r0, N); // z0 = r0

    //testPartialResults(mulMatrixAndVector(rowNum, N, mA, vecB), procSize, procRank, rowNum);

    //printProcRows(mProcRows, rowNum, N);
    if(procRank == 0){
        printVector(vecB, N);
    }

    /*

    while (1){
        double* temp[] = {NULL, NULL, NULL, NULL};

        temp[0] = mulMatrixAndVector(rowNum, N, mProcRows, z[0]);                      // Az(k)
        alpha[1] = scalarVectorAndVector(N, r[0], r[0])
                   / scalarVectorAndVector(N, temp[0], z[0]);      // alpha(k+1) = (r(k), r(k)) / (Az(k), z(k))

        temp[1] = mulVectorAndScalar(N, alpha[1], z[0]);
        x[1] = sumVectorAndVector(N, x[0], temp[1]);            // x(k+1) = x(k) + alpha(k+1)z(k)

        temp[2] = mulVectorAndScalar(N, alpha[1], temp[0]);
        r[1] = subVectorAndVector(N, r[0], temp[2]);            // r(k+1) = r(k) - alpha(k+1)Az(k)

        beta[1] = scalarVectorAndVector(N, r[1], r[1])
                  / scalarVectorAndVector(N, r[0], r[0]);         // b(k+1) = (r(k+1), r(k+1)) / (r(k), r(k))

        temp[3] = mulVectorAndScalar(N, beta[1], z[0]);
        z[1] = sumVectorAndVector(N, r[1], temp[3]);            // z(k+1) = r(k+1) + beta(k+1)z(k)

        alpha[0] = alpha[1];
        beta[0] = beta[1];

        delete[] x[0];
        delete[] z[0];
        delete[] r[0];

        x[0] = x[1];
        z[0] = z[1];
        r[0] = r[1];

        for (size_t ui = 0; ui < 4; ++ui) {
            delete[] temp[ui];
        }

        ++repeats;
        if( (vectorLength(N, r[1]) / vectorLength(N, vecB) ) < EPSILON){    // |r(k+1)| / |b| < EPSILON
            break;
        }
    }

//    printVector(vecU, N);
//    printVector(x[1], N);

    if(procRank == 0){
        delete[] vecU;
    }

     */

    MPI_Finalize();

    return 0;
}