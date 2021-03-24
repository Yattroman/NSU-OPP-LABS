#include <iostream>
#include <cstdlib>
#include <cstring>
#include <omp.h>
#include "mv/mvOperations_l2_3rd.h"
#include "mv/mvInit_l2_3rd.h"

#define EPSILON 1e-10
#define ACCURACY 1e-10

int main(int argc, char *argv[]) {
    int N = atoi(argv[1]);
    int repeats = 0;

    double resHolder1 = 0;
    double resHolder2 = 0;

    int growStatus = 0;

    double *mA = (double *) calloc(N * N, sizeof(double));
    initMatrixA(N, mA);

    double *vecB = (double *) calloc(N, sizeof(double));
    initVectorB(N, vecB);

    double *r[2];
    double *z[2];
    double *x[2];

    for (int i = 0; i < 2; ++i) {
        r[i] = (double *) calloc(N, sizeof(double));
        z[i] = (double *) calloc(N, sizeof(double));
        x[i] = (double *) calloc(N, sizeof(double));
    }

    double alpha[2];
    double beta[2];

    std::memcpy(r[0], vecB, N * sizeof(double)); // r0 = b - Ax0, где x0 - нулевой вектор
    std::memcpy(z[0], r[0], N * sizeof(double)); // z0 = r0

    double *temp[4];

    for (size_t ui = 0; ui < 4; ++ui) {
        temp[ui] = (double *) calloc(N, sizeof(double));
    }

#pragma omp parallel shared(temp, r, z, x, resHolder1, resHolder2, mA, vecB, alpha, beta)
    {
//        #pragma omp single
//        printVector(vecB, N);

//        std::cout << scalarVectorAndVector(N, vecB, vecB, resHolder) << std::endl;
//        std::cout << scalarVectorAndVector(N, vecB, vecB, resHolder) << std::endl;

#pragma omp for reduction(+: resHolder1, resHolder2)
        for (size_t i = 0; i < N; ++i) {
            resHolder1 += vecB[i] * vecB[i];
            resHolder2 += vecB[i] * vecB[i];
        }
#pragma omp barrier
#pragma omp critical
        std::cout << "'" << resHolder1 << "'";
#pragma omp critical
        std::cout << "'" << resHolder2 << "'";
#pragma omp barrier
#pragma omp single
        {
            resHolder1 = 0;
            resHolder2 = 0;
        }

//        printVector(r[0], N);
        /*while (1) {
            mulMatrixAndVector(N, mA, z[0], temp[0]);                      // Az(k)
            alpha[1] = scalarVectorAndVector(N, r[0], r[0])
                       / scalarVectorAndVector(N, temp[0], z[0]);      // alpha(k+1) = (r(k), r(k)) / (Az(k), z(k))

            mulVectorAndScalar(N, alpha[1], z[0], temp[1]);
            sumVectorAndVector(N, x[0], temp[1], x[1]);            // x(k+1) = x(k) + alpha(k+1)z(k)

            mulVectorAndScalar(N, alpha[1], temp[0], temp[2]);
            subVectorAndVector(N, r[0], temp[2], r[1]);            // r(k+1) = r(k) - alpha(k+1)Az(k)

            beta[1] = scalarVectorAndVector(N, r[1], r[1])
                      / scalarVectorAndVector(N, r[0], r[0]);         // b(k+1) = (r(k+1), r(k+1)) / (r(k), r(k))

            mulVectorAndScalar(N, beta[1], z[0], temp[3]);
            sumVectorAndVector(N, r[1], temp[3], z[1]);            // z(k+1) = r(k+1) + beta(k+1)z(k)

    #pragma omp master
            {
                alpha[0] = alpha[1];
                beta[0] = beta[1];
            }

            if ((vectorLength(N, r[0]) / vectorLength(N, vecB)) < EPSILON) {    // |r(k)| / |b| < EPSILON
                break;
            }

            if (growStatus > 10) {
                break;
            } else if (vectorLength(N, r[0]) < vectorLength(N, r[1])) {
                growStatus++;
            } else if (vectorLength(N, r[0]) > vectorLength(N, r[1])) {
                growStatus = 0;
            }

    #pragma omp master
            {
                std::memcpy(x[0], x[1], N * sizeof(double));
                std::memcpy(r[0], r[1], N * sizeof(double));
                std::memcpy(z[0], z[1], N * sizeof(double));
                ++repeats;
            }*/
    }

    if (growStatus <= 10) {
//        printVector(x[1], N);
//printVector(vecB, N);
        std::cout << "Repeats in total: " << repeats << "\n";
    } else {
        std::cout << "There are no roots!\n";
    }


    for (size_t ui = 0; ui < 4; ++ui) {
        free(temp[ui]);
    }

    for (size_t i = 0; i < 2; ++i) {
        free(x[i]);
        free(z[i]);
        free(r[i]);
    }

    free(mA);
    free(vecB);

    return EXIT_SUCCESS;
}