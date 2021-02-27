#include <cstring>
#include "mvInit_2nd.h"

void initMatrixA(int N, double* matrixA){
    for(size_t i = 0; i < N; ++i){
        for(size_t j = 0; j < N; ++j){
            if(i == j) {
                matrixA[i*N + j] = 2.0;
            } else {
                matrixA[i*N + j] = 1.0;
            }
        }
    }
}

void initVectorU(int N, double* vectorU){
    for(size_t i = 0; i < N; ++i){
        //vectorU[i] = sin(2*M_PI*(i+1)/N );
        vectorU[i] = 1;
    }
}

void initVectorB(int N, double* matrixA, double* vectorU, double* vectorB){
    memset(vectorB, 0, N);
    for(size_t i = 0; i < N; ++i){
        for(size_t j = 0; j < N; ++j){
            vectorB[j] += matrixA[j + i*N] * vectorU[j];
        }
    }
}

void initVectorX(int N, double* vectorX){
    memset(vectorX, 0, N);

    for(size_t i = 0; i < N; ++i){
        //vectorU[i] = sin(2*M_PI*(i+1)/N );
        vectorX[i] = 1;
    }
}