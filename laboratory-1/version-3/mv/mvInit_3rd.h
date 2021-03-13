#ifndef OPP_NSU_LABS_MVINIT_H
#define OPP_NSU_LABS_MVINIT_H

#include <cstdlib>
#include <cmath>

void initVectorU(int N, double* vectorU);
void initMatrixProcRows(int rowNumMod, int N, double* matrixProcRows, int procRank, int lastRowAdding, int* displs);
void initVectorBPart(int rowNumMod, double* vectorBPart, int* displs, int procRank);

#endif //OPP_NSU_LABS_MVINIT_H
