#include <cstring>
#include "mvOperations_2nd.h"

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

double scalarVectorAndVector(int N, const double* vectorL, const double* vectorR){ // Memory OK
    double res = 0;

    for(size_t i = 0; i < N; ++i) {
        res += vectorL[i] * vectorR[i];
    }

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