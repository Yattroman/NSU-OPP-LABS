#ifndef OPP_NSU_LABS_MVINIT_2ND_H
#define OPP_NSU_LABS_MVINIT_H

#include <cstdlib>
#include <cmath>

double* initMatrixA(int N);
double* initVectorU(int N);
double* initVectorB(int N, double* mA, double* vU);
double* initVectorX(int N);

#endif //OPP_NSU_LABS_MVINIT_2ND_H
