#include <cstring>
#include "mvOperations_l2_1st.h"

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
    for(size_t j = 0; j < N; ++j){
        res[j] = vectorL[j] - vectorR[j];
    }
}

void sumVectorAndVector(int N, const double* vectorL, const double* vectorR, double * res){
    for(size_t j = 0; j < N; ++j){
        res[j] = vectorL[j] + vectorR[j];
    }
}

void mulMatrixAndVector(int N, const double* matrix, const double* vector, double * res){
    for(size_t i = 0; i < N; ++i) {
        res[i] = 0;
        for (size_t j = 0; j < N; ++j) {
            res[i] += matrix[i*N + j] * vector[j];
        }
    }
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

void mulVectorAndScalar(int N, double scalar, const double* vector, double * res){
    for(size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            res[j] = scalar * vector[j];
        }
    }
}