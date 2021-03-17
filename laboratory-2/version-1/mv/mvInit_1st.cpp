#include <cstring>
#include "mvInit_1st.h"

double* initMatrixA(int N){
    double* mA = (double*) calloc(N*N, sizeof(double));

    for(size_t i = 0; i < N; ++i){
        srand(i);
        for(size_t j = 0; j < N; ++j){
            if(i == j) {
                mA[i*N + j] = rand() % 300 + 1001;
            } else {
                mA[i*N + j] = 999.0;
            }
        }
    }

    return mA;
}

double* initVectorB(int N){
    double* vecB = (double*) calloc(N, sizeof(double));

    for(size_t i = 0; i < N; ++i){
        srand(i);
        vecB[i] = rand() % 300;
    }

    return vecB;
}