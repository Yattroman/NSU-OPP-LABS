#include <mpi.h>
#include <cstdlib>
#include <iostream>

#define P1 3
#define P2 2
#define DIMENSION 2

using namespace std;

int isN1andN3SimilarP1andP2Respectively(int& N1, int& N3, int& p1, int& p2){
    return (N1%p1 + N3%p2 != 0) ? 0 : 1;
}

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
            matrixB[i*N3 + j] = i*N3 + j + 1;
        }
    }
}

int main(int argc, char* argv[]){
    int procSize;

    MPI_Init(&argc,&argv);

    MPI_Comm comm2d;
    MPI_Comm_size(MPI_COMM_WORLD, &procSize);

    int dims[DIMENSION] = {0, 0};
    int periods[DIMENSION] = {0, 0};
    int procCoords[DIMENSION];
    int reorder = 1;

    MPI_Dims_create(procSize, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &comm2d);
    MPI_Cart_get(comm2d, 2, dims, periods, procCoords);

    MPI_Comm yComms[dims[1]];
    MPI_Comm xComms[dims[0]];

    int N1 = atoi(argv[1]);
    int N2 = atoi(argv[2]);
    int N3 = atoi(argv[3]);

    MPI_Datatype bMatColumn, bMatColumnType;
//    MPI_Type_vector(); // N1 / p1

    MPI_Type_vector(N3, 1, N3, MPI_DOUBLE, &bMatColumn);
    MPI_Type_commit(&bMatColumn);
    MPI_Type_create_resized(bMatColumn, 0, sizeof(double), &bMatColumnType);
    MPI_Type_commit(&bMatColumnType);

    double* matrixA;
    double* matrixB;
    double* matrixC;

    double * matrixAPart = new double[N1*N2/dims[0]];
    double * matrixBPart = new double[N2*N3/dims[1]];
//    double * matrixCPart = new double[112]; // !!!

    if( procCoords[0] == 0 && procCoords[1] == 0 ){
        matrixA = new double[N1*N2];
        matrixB = new double[N2*N3];

        initMatrixA(matrixA, N1, N2);
        initMatrixB(matrixB, N2, N3);

//        printMatrix(matrixA, N1, N2);
        printMatrix(matrixB, N2, N3);
        cout << endl;
    }

    fillXandYComms(yComms, xComms, comm2d);
    MPI_Scatter(matrixA, N1*N2/dims[0], MPI_DOUBLE, matrixAPart, N1*N2/dims[0], MPI_DOUBLE, 0, yComms[0]);
    MPI_Scatter(matrixB, N3/dims[1], bMatColumnType, matrixBPart, N3/dims[1], bMatColumnType, 0, xComms[0]);

//    for (int i = 0; i < dims[1]; ++i) {
//        MPI_Bcast(matrixAPart, N1*N2/dims[0], MPI_DOUBLE, 0, xComms[i]);
//    }

    if(procCoords[0] == 0) {
        cout << "Y: " << procCoords[0] << " " << "X: " << procCoords[1] << '\n';
        printMatrix(matrixAPart, N1 / dims[0], N2);
        printMatrix(matrixBPart, N2, N3/dims[1]);
        cout << endl;
    }
//    cout << N1/dims[0] << " " << N2;

    MPI_Finalize();

    return EXIT_SUCCESS;
}