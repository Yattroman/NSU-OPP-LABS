#ifndef OPP_NSU_LABS_MVOPERATIONS_3RD_H
#define OPP_NSU_LABS_MVOPERATIONS_3RD_H

#include <cstdlib>
#include <iostream>
#include <cmath>

void printMatrix(double* matrix, int M, int N);
void printProcRows(double* matrix, int M, int N);
void printVector(double * vector, int N, int procRank);

double* mulMatrixAndVector(int M, int N, double* matrix, double* vector);

double* subVectorAndVector(int N, double* vectorL, double* vectorR);
double* sumVectorAndVector(int N, double* vectorL, double* vectorR);
double scalarVectorAndVector(int rowNumMod, const double* vectorL, const double* vectorR);

double vectorLength(int N, const double* vector);

double* mulVectorAndScalar(int N, double scalar, double* vector);

#endif
