#ifndef OPP_NSU_LABS_MVOPERATIONS_2ND_H
#define OPP_NSU_LABS_MVOPERATIONS_H

#include <cstdlib>
#include <iostream>
#include <cmath>

void printMatrix(double* matrix, int N);
void printVector(double * vector, int N);

void mulMatrixAndVector(int N, const double* matrix, const double* vector, double * res);

void subVectorAndVector(int N, const double* vectorL, const double* vectorR, double * res);
void sumVectorAndVector(int N, const double* vectorL, const double* vectorR, double * res);

void mulVectorAndScalar(int N, double scalar, const double* vector, double * res);

#endif
