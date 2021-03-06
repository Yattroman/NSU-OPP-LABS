#include <cstring>
#include <mpi.h>
#include "mvOperations_3rd.h"

void printMatrix(double* matrix, int M, int N){
    for(size_t i = 0; i < M; ++i){
        for(size_t j = 0; j < N; ++j){
            printf("%f ", matrix[i*N + j]);
        }
        printf("\n");
    }
}

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

double* subVectorAndVector(int N, double* vectorL, double* vectorR){
    double* res = (double*) calloc(N, sizeof(double));

    for(size_t j = 0; j < N; ++j){
        res[j] = vectorL[j] - vectorR[j];
    }

    return res;
}

double* sumVectorAndVector(int N, double* vectorL, double* vectorR){
    double* res = (double*) calloc(N, sizeof(double));

    for(size_t j = 0; j < N; ++j){
        res[j] = vectorL[j] + vectorR[j];
    }

    return res;
}

double* mulMatrixAndVector(int M, int N, double* matrix, double* vector){
    double* res = (double*) calloc(M, sizeof(double));

    for(size_t i = 0; i < M; ++i) {
        res[i] = 0;
        for (size_t j = 0; j < N; ++j) {
            res[i] += matrix[i*N + j] * vector[j];
        }
    }

    return res;
}

double scalarVectorAndVector(int rowNumMod, const double* vectorL, const double* vectorR){ // Memory OK
    // res = c1 + c2 + c3 + ...
    double res = 0;

    // temp1 = a1*b1=c1 | a2*b2=c2 | a3*b3=c3 | ...
    double* temp1 = (double*) calloc(rowNumMod, sizeof(double));

    // temp2 = с1 | с2 | с3 | ...
    double* temp2 = (double*) calloc(rowNumMod, sizeof(double));

    for(size_t i = 0; i < rowNumMod; ++i) {
        temp1[i] += vectorL[i] * vectorR[i];
    }

    MPI_Reduce(temp1, temp2, rowNumMod, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    for(size_t i = 0; i < rowNumMod; ++i) {
        res += temp1[i] * temp2[i];
    }

    MPI_Reduce(temp1, temp2, rowNumMod, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


    return res;
}

double vectorLength(int N, const double* vector){ // Memory OK
    double res = 0;

    for(size_t i = 0; i < N; ++i) {
        res += vector[i] * vector[i];
    }

    res = sqrt(res);

    return res;
}

double* mulVectorAndScalar(int N, double scalar, double* vector){
    double* res = (double*) calloc(N, sizeof(double));

    for(size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            res[j] = scalar * vector[j];
        }
    }

    return res;
}