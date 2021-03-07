#ifndef OPP_NSU_LABS_MVOPERATIONS_3RD_H
#define OPP_NSU_LABS_MVOPERATIONS_3RD_H

#include <cstdlib>
#include <iostream>
#include <cmath>

void printProcRows(const double* matrix, int M, int N);
void printVector(const double * vector, int N, int procRank);

double* mulMatrixAndVector(int rowNum, int lastRowAdding, int rowNumMod, int N, const double* matrixPart, const double* vectorPart, int* recvcounts);

double* subVectorAndVector(int rowNumMod, const double* vectorLPart, const double* vectorRPart);
double* sumVectorAndVector(int rowNumMod, const double* vectorLPart, const double* vectorRPart);
double scalarVectorAndVector(int rowNumMod, const double* vectorLPart, const double* vectorRPart);

double vectorLength(int rowNumMod, const double* vectorPart);

double* mulVectorAndScalar(int N, double scalar, const double* vector);

#endif
