#include <cstring>
#include "mvInit_3rd.h"

void initVectorU(int N, double* vectorU){
    for(size_t i = 0; i < N; ++i){
        vectorU[i] = sin(2*M_PI*i/N );
    }
}

void initMatrixProcRows(int rowNum, int N, double* matrixProcRows, int procRank, int lastRowAdding, int * displs){
    if(procRank == 0)
        rowNum += lastRowAdding;

    for(size_t i = 0; i < rowNum; ++i){
        srand(displs[procRank]+ i);
        for(size_t j = 0; j < N; ++j){
            if( (displs[procRank] + i) !=  j){
                matrixProcRows[N*i + j] = 999.0;
            } else {
                if(procRank == 0){
                    matrixProcRows[N*i + i] = rand() % 300 + 1001;
                } else {
                    matrixProcRows[N*i + rowNum*procRank + lastRowAdding + i] = rand() % 300 + 1001;
                }
            }

        }
    }
}

void initVectorBPart(int rowNumMod, double* vectorBPart, int* displs, int procRank){
    for(size_t i = 0; i < rowNumMod; ++i){
        srand(displs[procRank] + i);
        vectorBPart[i] = rand() % 300;
    }
}