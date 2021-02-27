#include <iostream>
#include "mv/mvOperations.h"
#include "mv/mvInit.h"

#define EPSILON 0.0001

int main(int argc, char* argv[]){
    int N = atoi(argv[1]);
    int repeats = 0;

    double* x0 = initVectorX(N);
    double* mA = initMatrixA(N);
    double* vecU = initVectorU(N);
    double* vecB = initVectorB(N, mA, vecU);
    double* r0 = vecB; // r0 = b - Ax0, где x0 - нулевой вектор
    double* z0 = vecB; // z0 = r0

    double* r[] = {r0, NULL};
    double* z[] = {z0, NULL};
    double* x[] = {x0, NULL};
    double alpha[] = {0, 0};
    double beta[] = {0, 0};

    while (1){
        double* temp[] = {NULL, NULL, NULL, NULL};

        temp[0] = mulMatrixAndVector(N, mA, z[0]);                      // Az(k)
        alpha[1] = scalarVectorAndVector(N, r[0], r[0])
                   / scalarVectorAndVector(N, temp[0], z[0]);      // alpha(k+1) = (r(k), r(k)) / (Az(k), z(k))

        temp[1] = mulVectorAndScalar(N, alpha[1], z[0]);
        x[1] = sumVectorAndVector(N, x[0], temp[1]);            // x(k+1) = x(k) + alpha(k+1)z(k)

        temp[2] = mulVectorAndScalar(N, alpha[1], temp[0]);
        r[1] = subVectorAndVector(N, r[0], temp[2]);            // r(k+1) = r(k) - alpha(k+1)Az(k)

        beta[1] = scalarVectorAndVector(N, r[1], r[1])
                  / scalarVectorAndVector(N, r[0], r[0]);         // b(k+1) = (r(k+1), r(k+1)) / (r(k), r(k))

        temp[3] = mulVectorAndScalar(N, beta[1], z[0]);
        z[1] = sumVectorAndVector(N, r[1], temp[3]);            // z(k+1) = r(k+1) + beta(k+1)z(k)

        alpha[0] = alpha[1];
        beta[0] = beta[1];

        delete[] x[0];
        delete[] z[0];
        //delete[] r[0];

        x[0] = x[1];
        z[0] = z[1];
        r[0] = r[1];

//        for (size_t ui = 0; ui < 4; ++ui) {
//            delete[] temp[ui];
//        }

        ++repeats;
        if( (vectorLength(N, r[1]) / vectorLength(N, vecB) ) < EPSILON){    // |r(k+1)| / |b| < EPSILON
            break;
        }
    }

    std::cout << '\n';

    printVector(vecU, N);
    printVector(x[1], N);

    printVector(mulMatrixAndVector(N, mA, vecU), N);

    std::cout << repeats;

    return 0;
}
