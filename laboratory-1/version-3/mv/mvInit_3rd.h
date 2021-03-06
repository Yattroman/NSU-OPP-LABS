#ifndef OPP_NSU_LABS_MVINIT_H
#define OPP_NSU_LABS_MVINIT_H

#include <cstdlib>
#include <cmath>

void initVectorU(int N, double* vectorU);
void initMatrixProcRows(int M, int N, double* matrixProcRows, int procRank, int lastRowAdding);

#endif //OPP_NSU_LABS_MVINIT_H
