#include "mvInit.h"

double* initMatrixA(int N){
    double* mA = (double*) malloc(sizeof(double)*N*N);

    for(size_t i = 0; i < N; ++i){
        for(size_t j = 0; j < N; ++j){
            if(i == j){
                mA[i*N + j] = 2.0;
            } else {
                mA[i*N + j] = 1.0;
            }
        }
    }

    return mA;
}

double* initVectorU(int N){
    double* vecU = (double*) malloc(sizeof(double)*N);

    for(size_t i = 0; i < N; ++i){
        vecU[i] = sin(2*M_PI*i/N );
    }

    return vecU;
}

double* initVectorB(int N, double* mA, double* vecU){
    double* vecB = (double*) calloc(N,sizeof(double)*N);

    for(size_t i = 0; i < N; ++i){
        for(size_t j = 0; j < N; ++j){
            vecB[j] += mA[j + i*N] * vecU[j];
        }
    }

    free(vecU);
    return vecB;
}

double* initVectorX(int N){
    double* vecX = (double*) calloc(N, sizeof(double)*N);
    return vecX;
}