#include <cstring>
#include "mvInit_l2_3rd.h"

void initMatrixA(int N, double* matrixA){
    for(size_t i = 0; i < N; ++i){
        srand(i);
        for(size_t j = 0; j < N; ++j){
            if(i == j) {
                matrixA[i*N + j] = rand() % 300 + 1001;
            } else {
                matrixA[i*N + j] = 999.0;
            }
        }
    }
}

void initVectorB(int N, double* vectorB){
    for(size_t i = 0; i < N; ++i){
        srand(i);
        vectorB[i] = rand() % 300;
    }
}