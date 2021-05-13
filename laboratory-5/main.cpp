#include <cstdlib>
#include <cmath>
#include <iostream>
#include "mpi.h"

#define MAX_LISTS_COUNT 5
#define TASKS_NUMBER 100
#define L 100
#define TRUE 1

using namespace std;

typedef struct TaskList{
    int repeatNumber;
} TaskList;

void printResults(int& pRank, int& tasksCompleted, double& globalRes, double& timeSpent){
    cout << "Process rank: " << pRank << endl;
    cout << "Tasks completed: " << tasksCompleted << endl;
    cout << "Global result: " << globalRes << endl;
    cout << "Time spent: " << timeSpent << endl;

    cout << endl;
}

void fillTasksList(int pRank, int pSize, TaskList * TL, int iCounter){
    for (int i = 0; i < TASKS_NUMBER; ++i) {
        TL[i].repeatNumber = abs(TASKS_NUMBER/2 - i) * abs(pRank - iCounter % pSize) * L;
    }
}

void doTasks(int& iCounter, TaskList * TL, double& globalResult, int& pRank, int& pSize){
    double startTime, endTime, pDiff, timeSpent;
    int tasksCompleted = 0;

    while(iCounter < MAX_LISTS_COUNT){
        startTime = MPI_Wtime();

        fillTasksList(pRank, pSize, TL, iCounter);

        for (int i = 0; i < TASKS_NUMBER; ++i) {
            for (int j = 0; j < TL[i].repeatNumber; ++j) {
                globalResult += sin(i);
            }
            tasksCompleted++;
        }

        endTime = MPI_Wtime();
        pDiff = endTime - startTime;
        MPI_Reduce(&pDiff, &timeSpent, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        // Синхронизация процессов
        MPI_Barrier(MPI_COMM_WORLD);
        // Вывод результатов с этой итерации
        printResults(pRank, tasksCompleted, globalResult, timeSpent);

        iCounter++;
    }
}

int main(int argc, char ** argv){
    MPI_Init(&argc, &argv);

    int pRank, pSize;

    MPI_Comm_size(MPI_COMM_WORLD, &pSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &pRank);

    int iCounter = 0; // Iterations Counter
    double globalResult = 0;

    TaskList TL[TASKS_NUMBER];

    doTasks(iCounter, TL, globalResult, pRank, pSize);

    MPI_Finalize();
}