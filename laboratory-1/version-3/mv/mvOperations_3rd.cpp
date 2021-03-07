#include <cstring>
#include <mpi.h>
#include "mvOperations_3rd.h"

void printVector(double* vector, int N, int procRank){
    printf("proc rank: %d. ", procRank);
    for(size_t i = 0; i < N; ++i){
        printf("%f ", vector[i]);
    }
    printf("\n");
}

void printProcRows(double* matrix, int M, int N){
    for(size_t i = 0; i < M; ++i){
        for(size_t j = 0; j < N; ++j){
            std::cout << matrix[i*N + j] << ' ';
        }
        std::cout << '\n';
    }
}

double* subVectorAndVector(int rowNumMod, const double* vectorLPart, const double* vectorRPart){
    double* res = (double*) calloc(rowNumMod, sizeof(double));

    for(size_t j = 0; j < rowNumMod; ++j){
        res[j] = vectorLPart[j] - vectorRPart[j];
    }

    return res;
}

double* sumVectorAndVector(int rowNumMod, const double* vectorLPart, const double* vectorRPart){
    double* res = (double*) calloc(rowNumMod, sizeof(double));

    for(size_t j = 0; j < rowNumMod; ++j){
        res[j] = vectorLPart[j] + vectorRPart[j];
    }

    return res;
}

double* mulMatrixAndVector(int rowNumMod, int N, const double* matrixPart, const double* vectorPart, int* recvcounts){
    double* res = (double*) calloc(rowNumMod, sizeof(double));

    double temp[N];

    for (int i = 0; i < rowNumMod; ++i) {
        for (int j = 0; j < N; ++j) {
            temp[j] += matrixPart[i*rowNumMod + j]*vectorPart[i];
        }
    }

    MPI_Reduce_scatter(temp, res, recvcounts, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return res;
}

double scalarVectorAndVector(int rowNumMod, const double* vectorLPart, const double* vectorRPart) { // OK.
    // res = tempRes1 + tempRes2 + tempRes3 + ... + tempRes(procRank) [between processes]
    double res = 0;

    // tempRes = a1*b1 + a2*b2 + a3*b3 + ...+ a(rowNumMod)*b(rowNumMod) [between vectorL and vector R]
    double tempRes = 0;

    for (int i = 0; i < rowNumMod; ++i) {
        tempRes += vectorLPart[i] * vectorRPart[i];
    }

    MPI_Allreduce(&tempRes, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return res;
}

double vectorLength(int rowNumMod, const double* vectorPart){ // OK.
    double res = sqrt(scalarVectorAndVector(rowNumMod, vectorPart, vectorPart));
    return res;
}

double* mulVectorAndScalar(int rowNumMod, double scalar, const double* vector){
    double* res = (double*) calloc(rowNumMod, sizeof(double));

    for(size_t i = 0; i < rowNumMod; ++i) {
        res[i] = scalar * vector[i];
    }

    return res;
}