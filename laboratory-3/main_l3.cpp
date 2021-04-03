#include <mpi.h>
#include <cstdlib>
#include <iostream>

#define P1 3
#define P2 2
#define DIMENSION 2

using namespace std;

void fillXandYComms(MPI_Comm* yComms, MPI_Comm* xComms, MPI_Comm comm2d){
    int remainDimsY[2] = {true, false};
    int remainDimsX[2] = {false, true};

    MPI_Cart_sub(comm2d, remainDimsX, xComms);
    MPI_Cart_sub(comm2d, remainDimsY, yComms);
}

void printMatrix(double* matrix, int N1, int N2){
    for(int i = 0; i < N1; ++i){
        for(int j = 0; j < N2; ++j){
            cout << matrix[i*N2 + j] << " ";
        }
        cout << endl;
    }
}

void initMatrixA(double* matrixA, int N1, int N2){
    for (int i = 0; i < N1; ++i) {
        for (int j = 0; j < N2; ++j) {
            matrixA[i*N2 + j] = i*N2 + j + 1;
        }
    }
}

void initMatrixB(double* matrixB, int N2, int N3){
    for (int i = 0; i < N2; ++i) {
        for (int j = 0; j < N3; ++j) {
            matrixB[i*N3 + j] = 2;
        }
    }
}

int main(int argc, char* argv[]){
    int procSize;

    MPI_Init(&argc,&argv);
    MPI_Comm comm2d;
    MPI_Comm_size(MPI_COMM_WORLD,&procSize);

    int dims[DIMENSION] = {0, 0};
    int periods[DIMENSION] = {0, 0};
    int procCoords[DIMENSION];
    int reorder = 1;

    MPI_Dims_create(procSize, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &comm2d);
    MPI_Cart_get(comm2d, 2, dims, periods, procCoords);

//    MPI_Comm yComms[procCoords[0]];
//    MPI_Comm xComms[procCoords[1]];

    int N1 = atoi(argv[1]);
    int N2 = atoi(argv[2]);
    int N3 = atoi(argv[3]);

    double* matrixA;
    double* matrixB;
    double* matrixC;

    double * matrixAPart = new double[N1*N2/procCoords[0]];
    double * matrixBPart = new double[N2*N3/procCoords[1]];
    double * matrixCPart = new double[112]; // !!!

    cout << "Y: " << procCoords[0] << " " << "X: " << procCoords[1] << '\n';

//    if( procCoords[0] == 0 && procCoords[1] == 0 ){
//        matrixA = new double[N1*N2];
//        matrixB = new double[N2*N3];
//
//        initMatrixA(matrixA, N1, N2);
//        initMatrixB(matrixB, N2, N3);

//        printMatrix(matrixA, N1, N2);
//        cout << endl;
//        printMatrix(matrixB, N2, N3);
//    }

//    fillXandYComms(yComms, xComms, comm2d);
//    MPI_Scatter(matrixA, N1*N2 / procCoords[0], MPI_DOUBLE, matrixAPart, N1*N2 / procCoords[0], MPI_DOUBLE, 0, yComms[0]);

//    printMatrix(matrixAPart, N1, N2/procCoords[0]);
//    cout << endl;
//    printMatrix(matrixAPart, N2/procCoords[1], N3);
                                   
    MPI_Finalize();

    return EXIT_SUCCESS;
}

