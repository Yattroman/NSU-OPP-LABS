#include <iostream>
#include <cstdlib>
#include <cstring>
#include "mv/mvOperations_1st.h"
#include "mv/mvInit_1st.h"

#define EPSILON 1e-5

int main(int argc, char* argv[]){
    int N = atoi(argv[1]);
    int repeats = 0;

    int growStatus = 0;

    double* mA = initMatrixA(N);
    double* vecU = initVectorU(N);
    double* vecB = initVectorB(N, mA, vecU);

    double* x0 = initVectorX(N);
    double* r0 = (double*) calloc(N, sizeof(double));
    double* z0 = (double*) calloc(N, sizeof(double));

    std::memcpy(r0, vecB, N*sizeof(double)); // r0 = b - Ax0, где x0 - нулевой вектор
    std::memcpy(z0, r0, N*sizeof(double)); // z0 = r0

    double* r[] = {r0, NULL};
    double* z[] = {z0, NULL};
    double* x[] = {x0, NULL};

    r[1] = (double*) calloc(N, sizeof(double));
    z[1] = (double*) calloc(N, sizeof(double));
    x[1] = (double*) calloc(N, sizeof(double));

    double alpha[] = {0, 0};
    double beta[] = {0, 0};

    double* temp[4];

    for (size_t ui = 0; ui < 4; ++ui) {
        temp[ui] = (double*) calloc(N, sizeof(double));
    }

    while (1){
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

        alpha[0] = alpha[1];
        beta[0] = beta[1];

        if( (vectorLength(N, r[0]) / vectorLength(N, vecB) ) < EPSILON){    // |r(k)| / |b| < EPSILON
            break;
        }

        if(growStatus > 10){
            break;
        } else if( vectorLength(N, r[0]) < vectorLength(N, r[1]) ){
            growStatus++;
        } else if( vectorLength(N, r[0]) > vectorLength(N,  r[1]) ){
            growStatus = 0;
        }

        std::memcpy(x[0], x[1], N*sizeof(double));
        std::memcpy(r[0], r[1], N*sizeof(double));
        std::memcpy(z[0], z[1], N*sizeof(double));

        ++repeats;
    }

    if(growStatus <= 10){
        printVector(vecU, N);
        printVector(x[1], N);
        std::cout << "Repeats in total: " << repeats << "\n";
    } else {
        std::cout << "There are no roots!\n";
    }


    for (size_t ui = 0; ui < 4; ++ui) {
        free(temp[ui]);
    }

    for (int i = 0; i < 2; ++i) {
        free(x[i]);
        free(z[i]);
        free(r[i]);
    }

    free(mA);
    free(vecU);

    return 0;
}
