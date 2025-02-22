#include <cstring>
#include "mvOperations_l2_2nd.h"

/*
void printMatrix(double* matrix, int N){
    for(size_t i = 0; i < N; ++i){
        for(size_t j = 0; j < N; ++j){
            printf("%f ", matrix[j + i*N]);
        }
        printf("\n");
    }
}
*/

void printVector(double* vector, int N){
    for(size_t i = 0; i < N; ++i){
        printf("%f ", vector[i]);
    }
    printf("\n");
}

void subVectorAndVector(int N, const double* vectorL, const double* vectorR, double * res){
    #pragma omp parallel for schedule(runtime)
    for(size_t j = 0; j < N; ++j){
        res[j] = vectorL[j] - vectorR[j];
    }
}

void sumVectorAndVector(int N, const double* vectorL, const double* vectorR, double * res){
    #pragma omp parallel for schedule(runtime)
    for(size_t j = 0; j < N; ++j){
        res[j] = vectorL[j] + vectorR[j];
    }
}

void mulMatrixAndVector(int N, const double* matrix, const double* vector, double * res){
    #pragma omp parallel for schedule(runtime)
    for(size_t i = 0; i < N; ++i) {
        res[i] = 0;
        #pragma omp parallel for schedule(runtime)
        for (size_t j = 0; j < N; ++j) {
            res[i] += matrix[i*N + j] * vector[j];
        }
    }
}

double scalarVectorAndVector(int N, const double* vectorL, const double* vectorR){ // Memory OK
    double res = 0;

    #pragma omp parallel for reduction (+: res) schedule(runtime)
    for(size_t i = 0; i < N; ++i) {
        res += vectorL[i] * vectorR[i];
    }

    return res;
}

double vectorLength(int N, const double* vector){ // Memory OK
    return sqrt(scalarVectorAndVector(N, vector, vector));
}

void mulVectorAndScalar(int N, double scalar, const double* vector, double * res){
    #pragma omp parallel for schedule(runtime)
    for(size_t i = 0; i < N; ++i) {
        res[i] = scalar * vector[i];
    }
}