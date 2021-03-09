#include <cstring>
#include <mpi.h>
#include "mvOperations_3rd.h"

void printVector(const double* vector, int N, int procRank){
    printf("proc rank: %d. ", procRank);
    for(size_t i = 0; i < N; ++i){
        printf("%f ", vector[i]);
    }
    printf("\n");
}

void printProcRows(const double* matrix, int M, int N){
    for(size_t i = 0; i < M; ++i){
        for(size_t j = 0; j < N; ++j){
            std::cout << matrix[i*N + j] << ' ';
        }
        std::cout << '\n';
    }
}

void subVectorAndVector(int rowNumMod, const double* vectorLPart, const double* vectorRPart, double * res){
    for(size_t j = 0; j < rowNumMod; ++j){
        res[j] = vectorLPart[j] - vectorRPart[j];
    }
}

void sumVectorAndVector(int rowNumMod, const double* vectorLPart, const double* vectorRPart, double * res){
    for(size_t j = 0; j < rowNumMod; ++j){
        res[j] = vectorLPart[j] + vectorRPart[j];
    }
}

void mulMatrixAndVector(int rowNum, int lastRowAdding, int rowNumMod, int N, const double* matrixPart, const double* vectorPart, int* recvcounts, double * res){ // OK.
    double temp[N];

    for (int i = 0; i < rowNumMod; ++i) {
        for (int j = 0; j < N; ++j) {
            temp[j] += matrixPart[i*N+j]*vectorPart[i];
        }
    }

    MPI_Reduce_scatter(temp, res, recvcounts, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

double scalarVectorAndVector(int rowNumMod, const double* vectorLPart, const double* vectorRPart) { // OK.
    // res = tempRes1 + tempRes2 + tempRes3 + ... + tempRes(procRank) [between processes]
    double res = 0;

    // tempRes = a1*b1 + a2*b2 + a3*b3 + ...+ a(rowNumMod)*b(rowNumMod) [between vectorL and vectorR]
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

void mulVectorAndScalar(int rowNumMod, double scalar, const double* vectorPart, double * res){
    for(size_t i = 0; i < rowNumMod; ++i) {
        res[i] = scalar * vectorPart[i];
    }
}