#ifndef OPP_NSU_LABS_MVOPERATIONS_3RD_H
#define OPP_NSU_LABS_MVOPERATIONS_3RD_H

#include <cstdlib>
#include <iostream>
#include <cmath>

void printProcRows(const double* matrix, int M, int N);
void printVector(const double * vector, int N, int procRank);

void mulMatrixAndVector(int rowNum, int lastRowAdding, int rowNumMod, int N, const double* matrixPart, const double* vectorPart, int* recvcounts, double * res);

void subVectorAndVector(int rowNumMod, const double* vectorLPart, const double* vectorRPart, double * res);
void sumVectorAndVector(int rowNumMod, const double* vectorLPart, const double* vectorRPart, double * res);
double scalarVectorAndVector(int rowNumMod, const double* vectorLPart, const double* vectorRPart);

double vectorLength(int rowNumMod, const double* vectorPart);

void mulVectorAndScalar(int N, double scalar, const double* vector, double * res);

#endif
