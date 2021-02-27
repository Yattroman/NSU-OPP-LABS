#ifndef OPP_NSU_LABS_MVOPERATIONS_2ND_H
#define OPP_NSU_LABS_MVOPERATIONS_H

#include <cstdlib>
#include <iostream>
#include <cmath>

void printMatrix(double* matrix, int N);
void printVector(double * vector, int N);

double* mulMatrixAndVector(int N, double* matrix, double* vector);

double* subVectorAndVector(int N, double* vectorL, double* vectorR);
double* sumVectorAndVector(int N, double* vectorL, double* vectorR);
double scalarVectorAndVector(int N, const double* vectorL, const double* vectorR);

double vectorLength(int N, const double* vector);

double* mulVectorAndScalar(int N, double scalar, double* vector);

#endif
