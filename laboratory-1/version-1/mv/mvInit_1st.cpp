#include <cstring>
#include "mvInit_1st.h"

double* initMatrixA(int N){
    double* mA = new double[N*N];

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
    double* vecU = new double[N];

    for(size_t i = 0; i < N; ++i){
        //vecU[i] = sin(2*M_PI*(i+1)/N );
        vecU[i] = 100;
    }

    return vecU;
}

double* initVectorB(int N, double* mA, double* vecU){
    double* vecB = new double[N];

    for(size_t i = 0; i < N; ++i){
        vecB[i] = 0;
        for(size_t j = 0; j < N; ++j){
            vecB[i] += mA[i*N + j] * vecU[j];
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