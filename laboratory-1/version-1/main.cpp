#include <stdio.h>
#include "mv/mvOperations.h"
#include "mv/mvInit.h"

#define EPSILON 0.00001

int main(int argc, char* argv[]){
    int N = 3;

    double* x0 = initVectorX(N);
    double* mA = initMatrixA(N);
    double* vecU = initVectorU(N);
    double* vecB = initVectorB(N, mA, vecU);
    double* r0 = subVectorAndVector(N, vecB,mulMatrixAndVector(N, mA, x0));
    double* z0 = subVectorAndVector(N, vecB,mulMatrixAndVector(N, mA, x0));

    double* r[] = {r0, NULL};
    double* z[] = {z0, NULL};
    double* x[] = {x0, NULL};
    double alpha[] = {0, 0};
    double beta[] = {0, 0};

    double* temp[] = {NULL, NULL, NULL, NULL, NULL};

    while (1){
        temp[0] = mulMatrixAndVector(N, mA, z[0]);
        alpha[1] = scalarVectorAndVector(N, r[0], r[0])
                   / scalarVectorAndVector(N, temp[0], z[0]);      // alpha(k+1) = (r(k), r(k)) / (Az(k), z(k))

        temp[1] = mulVectorAndScalar(N, alpha[1], z[0]);
        x[1] = sumVectorAndVector(N, x[0], temp[1]);            // x(k+1) = x(k) + alpha(k+1)z(k)

        temp[2] = mulMatrixAndVector(N,mA, z[0]);
        temp[3] = mulVectorAndScalar(N, alpha[1], temp[2]);
        r[1] = subVectorAndVector(N, r[0], temp[3]);            // r(k+1) = r(k) - alpha(k+1)Az(k)

        beta[1] = scalarVectorAndVector(N, r[1], r[1])
                  / scalarVectorAndVector(N, r[0], r[0]);         // b(k+1) = (r(k+1), r(k+1)) / (r(k), r(k))

        temp[4] = mulVectorAndScalar(N, beta[1], z[0]);
        z[1] = sumVectorAndVector(N, r[1], temp[4]);            // z(k+1) = r(k+1) - beta(k+1)z(k)

        if(vectorLength(N, r[1]) / vectorLength(N, vecB) < EPSILON){    // |r(k+1)| / |b| < EPSILON
            break;
        }

        for(size_t ti = 0; ti < 5; ++ti){
            free(temp[ti]);
        }

        free(x[0]);
        free(z[0]);
        free(r[0]);

        alpha[0] = alpha[1];
        beta[0] = beta[1];
        x[0] = x[1];
        z[0] = z[1];
        r[0] = r[1];
    }

    printVector(vecU, N);
    printVector(x[1], N);

    free(x[1]);
    free(z[1]);
    free(r[1]);

    return 0;
}
