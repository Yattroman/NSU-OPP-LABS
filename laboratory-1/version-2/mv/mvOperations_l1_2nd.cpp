#include <cstring>
#include "mvOperations_l1_2nd.h"

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

void mulMatrixAndVector(int rowNumMod, int N, const double* matrix, const double* vector, double * res){

    for(size_t i = 0; i < rowNumMod; ++i) {
        res[i] = 0;
        for (size_t j = 0; j < N; ++j) {
            res[i] += matrix[i*N + j] * vector[j];
        }
    }

}

double scalarVectorAndVector(int N, const double* vectorL, const double* vectorR){
    double res = 0;

    for(size_t i = 0; i < N; ++i) {
        res += vectorL[i] * vectorR[i];
    }

    return res;
}

double vectorLength(int N, const double* vector){
    double res = sqrt(scalarVectorAndVector(N, vector, vector));
    return res;
}

void mulVectorAndScalar(int N, double scalar, const double* vector, double * res){

    for(size_t i = 0; i < N; ++i) {
        res[i] = scalar * vector[i];
    }

}