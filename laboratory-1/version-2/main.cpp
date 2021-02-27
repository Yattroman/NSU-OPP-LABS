#include <mpi.h>
#include <iostream>
#include <unistd.h>
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

void distributeData(double* matrixA, double* vectorX, double* vectorB, double* matrixProcRows, int N, int rowNum, int procSize, int procRank){
    MPI_Bcast(vectorX, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(vectorB, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int* displs = new int[procSize];
    int* scounts = new int[procSize];

    for (size_t i=0; i< procSize; ++i) {
        displs[i] = i*N*rowNum;
        scounts[i] = N*rowNum;
    }

    MPI_Scatterv(matrixA, scounts, displs, MPI_DOUBLE, matrixProcRows, rowNum*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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
    double* mProcRows = new double[rowNum*N];

    double* mA;

    double* vecU = new double[N];
    double* x0 = new double[N];
    double* vecB = new double[N];
    double* r0; // r0 = b - Ax0, где x0 - нулевой вектор
    double* z0; // z0 = r0

    double* r[] = {r0, NULL};
    double* z[] = {z0, NULL};
    double* x[] = {x0, NULL};
    double alpha[] = {0, 0};
    double beta[] = {0, 0};

    if(procRank == 0){
        mA = new double[N*N];

        initVectorX(N, x0);
        initMatrixA(N, mA);
        initVectorU(N, vecU);
        initVectorB(N, mA, vecU, vecB);
        r0 = vecB; // r0 = b - Ax0, где x0 - нулевой вектор
        z0 = vecB; // z0 = r0
    }

    distributeData(mA, x0, vecB, mProcRows, N, rowNum, procSize, procRank);

    //std::cout << "test\n" << procRank;
    //printVector(vecB, N);
    //printVector(x0, N);

    testPartialResults(mulMatrixAndVector(rowNum, N, mA, x0), procSize, procRank, rowNum);
    MPI_Barrier(MPI_COMM_WORLD);

    /* while (1){
        double* temp[] = {NULL, NULL, NULL, NULL};

        temp[0] = mulMatrixAndVector(rowNum, N, mA, z[0]);                      // Az(k)
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
     */

//    printVector(vecU, N);
//    printVector(x[1], N);

    if(procRank == 0){
        delete[] mA;
        delete[] vecU;
    }

    MPI_Finalize();

    return 0;
}