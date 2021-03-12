#include <cstring>
#include "mvInit_2nd.h"

/*void initVectorU(int N, double* vectorU){
    for(size_t i = 0; i < N; ++i){
        vectorU[i] = sin(2*M_PI*i/N );
    }
}
 FOR TEST */

void initMatrixProcRows(int rowNumMod, int N, double* matrixProcRows, int procRank, int lastRowAdding){
    if(procRank == 0)
        rowNumMod += lastRowAdding;

    for(size_t i = 0; i < rowNumMod; ++i){
        srand(procRank*rowNumMod + i);
        for(size_t j = 0; j < N; ++j){
            if( (procRank*rowNumMod + i) !=  j){
                matrixProcRows[N*i + j] = 999.0;
            } else {
                if(procRank == 0){
                    matrixProcRows[N*i + i] = rand() % 300 + 20001;
                } else {
                    matrixProcRows[N*i + rowNumMod*procRank + lastRowAdding + i] = rand() % 300 + 20001;
                }
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