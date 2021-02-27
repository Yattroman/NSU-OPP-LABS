#include <mpi.h>
#include <iostream>
#include <unistd.h>

int main(int argc, char** argv)
{
    int  size, rank;
    char host[32];

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    gethostname(host,32);
    printf("Hello, MPI world! I'm number %d from %d and I run on host %s.\n", rank, size, host);

    MPI_Finalize();
    return 0;
}