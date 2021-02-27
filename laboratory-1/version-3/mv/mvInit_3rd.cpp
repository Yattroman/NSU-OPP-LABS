#include <cstring>
#include "mvInit_3rd.h"

double* initMatrixA(int N){
    double* mA = new double[N*N];

    for(size_t i = 0; i < N; ++i){
        for(size_t j = 0; j < N; ++j){
            if(i == j) {
                mA[i*N + j] = 300.0;
            } else {
                mA[i*N + j] = 1.0;
            }
        }
    }

    return mA;
}

double* initVectorU(int N){
    double* vecU = new double[N];

    for(size_t i = 0; i < N; ++i){
        vecU[i] = sin(2*M_PI*(i+1)/N );
    }

    return vecU;
}

double* initVectorB(int N, double* mA, double* vecU){
    double* vecB = new double[N];
    memset(vecB, 0, N);

    for(size_t i = 0; i < N; ++i){
        for(size_t j = 0; j < N; ++j){
            vecB[j] += mA[j + i*N] * vecU[j];
        }
    }

    //delete[] vecU;
    return vecB;
}

double* initVectorX(int N){
    double* vecX = new double[N];
    memset(vecX, 0, N);

    return vecX;
}