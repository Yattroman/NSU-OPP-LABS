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

    double breakCondition = 0;

    double resHolder1 = 0;
    double resHolder2 = 0;
    double resHolder3 = 0;

    double value = 0;

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

#pragma omp parallel shared(temp, r, z, x, resHolder1, resHolder2, mA, vecB, alpha, beta) private(value)
    {

/*#pragma omp for reduction(+:resHolder1)
        for (size_t i = 0; i < N; ++i) {
            value = 0;
            for (size_t j = 0; j < N; ++j) {
                value += mA[i * N + j] * z[0][j];
            }
            resHolder1 += value;
            temp[0][i] = resHolder1;

            resHolder1 = 0;
        }
#pragma omp single
        printVector(temp[0], N);*/

        while (1) {
            // Az(k)
#pragma omp for reduction(+:resHolder1)
            for (size_t i = 0; i < N; ++i) {
                value = 0;
                for (size_t j = 0; j < N; ++j) {
                    value += mA[i * N + j] * z[0][j];
                }
                resHolder1 += value;
                temp[0][i] = resHolder1;

                resHolder1 = 0;
            }

            // alpha(k+1) = (r(k), r(k)) / (Az(k), z(k))
#pragma omp for reduction(+: resHolder1, resHolder2)
            for (size_t i = 0; i < N; ++i) {
                resHolder1 += r[0][i] * r[0][i];
                resHolder2 += temp[0][i] * z[0][i];
            }
#pragma omp single
            {
                alpha[1] = resHolder1 / resHolder2;

                resHolder1 = 0;
                resHolder2 = 0;
            }

            mulVectorAndScalar(N, alpha[1], z[0], temp[1]);
            sumVectorAndVector(N, x[0], temp[1], x[1]);            // x(k+1) = x(k) + alpha(k+1)z(k)

            mulVectorAndScalar(N, alpha[1], temp[0], temp[2]);
            subVectorAndVector(N, r[0], temp[2], r[1]);            // r(k+1) = r(k) - alpha(k+1)Az(k)

            // beta(k+1) = (r(k+1), r(k+1)) / (r(k), r(k))
#pragma omp for reduction(+: resHolder1, resHolder2)
            for (size_t i = 0; i < N; ++i) {
                resHolder1 += r[1][i] * r[1][i];
                resHolder2 += r[0][i] * r[0][i];
            }
#pragma omp single
            {
                beta[1] = resHolder1 / resHolder2;

                resHolder1 = 0;
                resHolder2 = 0;
            }

            mulVectorAndScalar(N, beta[1], z[0], temp[3]);
            sumVectorAndVector(N, r[1], temp[3], z[1]);            // z(k+1) = r(k+1) + beta(k+1)z(k)

#pragma omp for reduction(+: resHolder1, resHolder2, resHolder3)
            for (size_t i = 0; i < N; ++i) {
                resHolder1 += r[0][i] * r[0][i]; // |r(k)|
                resHolder2 += vecB[i] * vecB[i]; // |b|
                resHolder3 += r[1][i] * r[1][i]; // |r(k+1)|
            }
#pragma omp single
            breakCondition = resHolder1 / resHolder2;

            if (breakCondition < EPSILON) {    // |r(k)| / |b| < EPSILON
                break;
            }

            if (growStatus > 10) {
                break;
            } else if (resHolder1 < resHolder3) {
                growStatus++;
            } else if (resHolder1 > resHolder3) {
                growStatus = 0;
            }
#pragma omp single
            {
                std::memcpy(x[0], x[1], N * sizeof(double));
                std::memcpy(r[0], r[1], N * sizeof(double));
                std::memcpy(z[0], z[1], N * sizeof(double));

                resHolder1 = 0;
                resHolder2 = 0;
                resHolder3 = 0;

                ++repeats;
            }
        }
    }

        if (growStatus <= 10) {
          printVector(x[1], N);
          printVector(vecB, N);
          printMatrix(mA, N);
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