#include <cstring>
#include "mvInit_1st.h"

double* initMatrixA(int N){
    double* mA = (double*) calloc(N*N, sizeof(double));

    for(size_t i = 0; i < N; ++i){
        for(size_t j = 0; j < N; ++j){
            if(i == j) {
                mA[i*N + j] = 2.0;
            } else {
                mA[i*N + j] = 1.0;
            }
        }
    }

    return mA;
}

double* initVectorU(int N){
    double* vecU = (double*) calloc(N, sizeof(double));

    for(size_t i = 0; i < N; ++i){
        vecU[i] = sin(2*M_PI*i/N );
//        vecU[i] = i;
    }

    return vecU;
}

double* initVectorB(int N, double* mA, double* vecU){
    double* vecB = (double*) calloc(N, sizeof(double));

    for(size_t i = 0; i < N; ++i){
        vecB[i] = 0;
        for(size_t j = 0; j < N; ++j){
            vecB[i] += mA[i*N + j] * vecU[j];
        }
    }

    return vecB;
}

double* initVectorX(int N){
    double* vecX = (double*) calloc(N, sizeof(double));
    return vecX;
}