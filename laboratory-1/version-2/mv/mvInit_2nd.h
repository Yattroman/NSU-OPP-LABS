#ifndef OPP_NSU_LABS_MVINIT_2ND_H
#define OPP_NSU_LABS_MVINIT_2ND_H

#include <cstdlib>
#include <cmath>

/*void initVectorU(int N, double* vectorU); FOR TEST*/
void initVectorB(int N, double * vectorB);
void initMatrixProcRows(int rowNumMod, int N, double* matrixProcRows, int procRank, int lastRowAdding);

#endif //OPP_NSU_LABS_MVINIT_2ND_H
