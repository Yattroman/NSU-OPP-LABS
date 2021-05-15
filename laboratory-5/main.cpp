#include <cstdlib>
#include <cmath>
#include <iostream>
#include "mpi.h"

#define THREADS_NUM 2

#define MAX_LISTS_COUNT 2
#define UPLOADING_START_VALUE 9
#define TASKS_NUMBER 10
#define L 10

using namespace std;

double globalResult = 0;

pthread_t pthreads[THREADS_NUM];

typedef struct TaskList{
    int repeatNumber;
} TaskList;

void printResults(int& pRank, int& tasksCompleted, double& timeSpent){
    cout << "Process rank: " << pRank << endl;
    cout << "Tasks completed: " << tasksCompleted << endl;
    cout << "Global result: " << globalResult << endl;
    cout << "Time spent: " << timeSpent << endl;

    cout << endl;
}

void fillTasksList(int pRank, int pSize, TaskList * TL, int iCounter){
    for (int i = 0; i < TASKS_NUMBER; ++i) {
        TL[i].repeatNumber = abs(TASKS_NUMBER/2 - i) * abs(pRank - iCounter % pSize) * L;
    }
}

void * executeTasks(void * args){
    // Executor deals
    auto TL = (TaskList*) args;

    pthread_mutex_t globalResMutex;
    pthread_mutex_init(&globalResMutex, nullptr);

    double localResult = 0;

    for (int i = 0; i < TASKS_NUMBER; ++i) {
        if(TL[i].repeatNumber < UPLOADING_START_VALUE){
            // Analyzer вычисляет, у какого процесса по сравнению с этим больше всего разница в незавершенных заданиях
            // Loader начинает подгрузку заданий с подходящего процесса на фоне вычислений
            // Ставим флаг подрузки, чтобы снова не заходить в эту секцию
            // Проблема может быть в круговой постоянной пересылке данных, если процессы выполняют задания +- одновременно
        }
        for (int j = 0; j < TL[i].repeatNumber; ++j) {
//            localResult += sin(i);
            localResult += 1;
        }
    }

    pthread_mutex_lock(&globalResMutex);
    globalResult += localResult;
    pthread_mutex_unlock(&globalResMutex);

    pthread_exit(nullptr);
}

void doTasks(int& iCounter, TaskList * TL, int& pRank, int& pSize){
    double startTime, endTime, pDiff, timeSpent;

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    int tasksCompleted = 0;

    while(iCounter < MAX_LISTS_COUNT){
        startTime = MPI_Wtime();

        fillTasksList(pRank, pSize, TL, iCounter);

        pthread_create(&pthreads[0], &attr, executeTasks, (void*) TL);
        pthread_join(pthreads[0], nullptr);

        double commonResult;
        MPI_Allreduce(&globalResult, &commonResult, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // Update global result on this iteration
        globalResult = commonResult;

        endTime = MPI_Wtime();

        // Calculate iteration time
        pDiff = endTime - startTime;
        MPI_Reduce(&pDiff, &timeSpent, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        // Processes sync
        MPI_Barrier(MPI_COMM_WORLD);

        // Print this iteration's results
        printResults(pRank, tasksCompleted, timeSpent);

        iCounter++;
    }

    pthread_attr_destroy(&attr);
}

int main(int argc, char ** argv){
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    int pRank, pSize;

    MPI_Comm_size(MPI_COMM_WORLD, &pSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &pRank);

    int iCounter = 0; // Iterations Counter

    TaskList TL[TASKS_NUMBER];

//    double temp = 0;
//    fillTasksList(pRank, pSize, TL, iCounter, temp);

//    cout << temp << endl;

    doTasks(iCounter, TL, pRank, pSize);

    MPI_Finalize();
}