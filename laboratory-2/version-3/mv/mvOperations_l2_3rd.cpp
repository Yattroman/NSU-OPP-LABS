#include <cstring>
#include "mvOperations_l2_3rd.h"

/*
void printMatrix(double* matrix, int N){
    for(size_t i = 0; i < N; ++i){
        for(size_t j = 0; j < N; ++j){
            printf("%f ", matrix[j + i*N]);
        }
        printf("\n");
    }
}
*/

void printVector(double *vector, int N) {
    for (size_t i = 0; i < N; ++i) {
        printf("%f ", vector[i]);
    }
    printf("\n");
}

void subVectorAndVector(int N, const double *vectorL, const double *vectorR, double *res) {
#pragma omp for
    for (size_t j = 0; j < N; ++j) {
        res[j] = vectorL[j] - vectorR[j];
    }
}

void sumVectorAndVector(int N, const double *vectorL, const double *vectorR, double *res) {
#pragma omp for
    for (size_t j = 0; j < N; ++j) {
        res[j] = vectorL[j] + vectorR[j];
    }
}

void mulMatrixAndVector(int N, const double *matrix, const double *vector, double *res) {
    double tempRes;

    for (size_t i = 0; i < N; ++i) {
        tempRes = 0;
#pragma omp for
        for (size_t j = 0; j < N; ++j) {
            tempRes += matrix[i * N + j] * vector[j];
        }
#pragma omp reduction(+:res[i])
        {
            res[i] += tempRes;
        }
    }
}

/*double scalarVectorAndVector(int N, const double *vectorL, const double *vectorR, double &resHolder) {
    double tempRes = resHolder;

#pragma omp for reduction(+: tempRes)
    for (size_t i = 0; i < N; ++i) {
        tempRes += vectorL[i] * vectorR[i];
    }

#pragma omp critical
    std::cout << "resHolder: '" << resHolder << "' ";

//    #pragma omp critical
//    memcpy(&tempRes, &resHolder, sizeof(double) );
//    tempRes = resHolder;

//#pragma omp critical
//    std::cout << "tempRes: '" << tempRes << "' ";
//    std::cout << "resHolder: '" << resHolder << "' ";

#pragma omp master
    resHolder = 0;

    return tempRes;
}*/

/*double vectorLength(int N, const double *vector, double &resHolder) {
    return sqrt(scalarVectorAndVector(N, vector, vector, resHolder));
}*/

void mulVectorAndScalar(int N, double scalar, const double *vector, double *res) {
#pragma omp for
    for (size_t i = 0; i < N; ++i) {
        res[i] = scalar * vector[i];
    }
}