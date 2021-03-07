#include <cstring>
#include "mvInit_3rd.h"

void initVectorU(int N, double* vectorU){
    for(size_t i = 0; i < N; ++i){
        //vectorU[i] = sin(2*M_PI*(i+1)/N );
        vectorU[i] = 1;
    }
}

void initMatrixProcRows(int M, int N, double* matrixProcRows, int procRank, int lastRowAdding){
    if(procRank == 0)
        M += lastRowAdding;

    for(size_t i = 0; i < M; ++i){
        for(size_t j = 0; j < N; ++j){
            matrixProcRows[N*i + j] = 1;
            if(procRank == 0){
                matrixProcRows[N*i + i] = 2;
            } else {
                matrixProcRows[N*i + M*procRank + lastRowAdding + i] = 2;
            }
        }
    }
}