#include <mpi.h>
#include <iostream>
#include <unistd.h>
#include <cstring>
#include "mv/mvOperations_2nd.h"
#include "mv/mvInit_2nd.h"

#define EPSILON 0.0001

void distributeData(double* vectorX, double* vectorB, int N){
    MPI_Bcast(vectorX, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(vectorB, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

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

    double* vecU = (double*) calloc(N, sizeof(double));
    double* x0 = (double*) malloc(N*sizeof(double));
    double* vecB = (double*) calloc(N, sizeof(double));
    double* r0; // r0 = b - Ax0, где x0 - нулевой вектор
    double* z0; // z0 = r0

    double* r[] = {r0, NULL};
    double* z[] = {z0, NULL};
    double* x[] = {x0, NULL};
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

        initVectorX(N, x0);

        initMatrixProcRows(rowNum, N, mProcRows, procRank, lastRowAdding);
        vecBPart = mulMatrixAndVector(rowNum+lastRowAdding, N, mProcRows, vecU);

        printVector(vecBPart, rowNum+lastRowAdding);
    } else {
        mProcRows = (double*) calloc(rowNum*N, sizeof(double));
        vecBPart = (double*) calloc(rowNum, sizeof(double));

        initMatrixProcRows(rowNum, N, mProcRows, procRank, lastRowAdding);
        vecBPart = mulMatrixAndVector(rowNum, N, mProcRows, vecU);
        printVector(vecBPart, rowNum);
    }

    std::cout << "Rank: " << procRank << " - " << displs[procRank] << "\n";

    //distributeData(x0, vecB, N);

    //r0 = vecB; // r0 = b - Ax0, где x0 - нулевой вектор
    //memcpy(z0, r0, N); // z0 = r0

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

        free(x[0]);
        free(z[0]);
        free(r[0]);

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
     */

    MPI_Finalize();

    return 0;
}