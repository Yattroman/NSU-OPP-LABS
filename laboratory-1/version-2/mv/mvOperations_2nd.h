#ifndef OPP_NSU_LABS_MVOPERATIONS_2ND_H
#define OPP_NSU_LABS_MVOPERATIONS_2ND_H

#include <cstdlib>
#include <iostream>
#include <cmath>

void printMatrix(double* matrix, int M, int N);
void printProcRows(double* matrix, int M, int N);
void printVector(double * vector, int N, int procRank);

void mulMatrixAndVector(int rowNumMod, int N, const double* matrix, const double* vector, double * res);

void subVectorAndVector(int N, const double* vectorL, const double* vectorR, double * res);
void sumVectorAndVector(int N, const double* vectorL, const double* vectorR, double * res);
double scalarVectorAndVector(int N, const double* vectorL, const double* vectorR);

double vectorLength(int N, const double* vector);

void mulVectorAndScalar(int N, double scalar, const double* vector, double * res);

#endif
