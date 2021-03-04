#ifndef OPP_NSU_LABS_MVINIT_2ND_H
#define OPP_NSU_LABS_MVINIT_2ND_H

#include <cstdlib>
#include <cmath>

void initMatrixA(int N, double* matrixA);
void initVectorU(int N, double* vectorU);
void initVectorBPart(int N, double* matrixA, double* vectorU, double* vectorB);
void initVectorX(int N, double* vectorX);
void initMatrixProcRows(int M, int N, double* matrixProcRows, int procRank, int lastRowAdding);

#endif //OPP_NSU_LABS_MVINIT_2ND_H
