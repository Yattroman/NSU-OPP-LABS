#include <cstdlib>
#include <cmath>
#include <iostream>
#include <cstring>
#include <pthread.h>
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

typedef struct AnalyzeInfo{
    int uncompletedTasks;
    int procRankRequest;
    int procSize;

    int& maxDiff;
    int& maxDiffRank;
} AnalyzeInfo;

typedef struct LoadInfo{
    TaskList * addTL;
    TaskList * mainTL;
    int taskSent;
} LoadInfo;

void atError(){
    MPI_Finalize();
    exit(EXIT_FAILURE);
}

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

void * waitingNewTasks(void * args){

}

void * sendNewTasks(void * args){
    // Analyze which process has maximal load
    auto analyzeInf = (AnalyzeInfo*) args;

    int uncompletedTasks = analyzeInf->uncompletedTasks;

    int procRankRequest = analyzeInf->procRankRequest;
    int uncompletedTasksRequest = analyzeInf->uncompletedTasks;

    int procSize = analyzeInf->procSize;

    int localDiff;
    int differs[procSize];
    int maxDiff = 0;
    int maxDiffRank;

    memset(differs, -1, procSize*sizeof(double));

    MPI_Bcast(&uncompletedTasksRequest, 1, MPI_INT, procRankRequest, MPI_COMM_WORLD);

    localDiff = uncompletedTasks - uncompletedTasksRequest;
    MPI_Allgather(&localDiff, 1, MPI_INT, differs, 1, MPI_INT, MPI_COMM_WORLD);

    for (int i = 0; i < procSize; ++i) {
        if (differs[i] > maxDiff) {
            maxDiff = differs[i];
            maxDiffRank = i;
        }
    }

    if(maxDiff > 3){
        analyzeInf->maxDiff = maxDiff;
        analyzeInf->maxDiffRank = maxDiffRank;
    } else {
        analyzeInf->maxDiff = -1;
        analyzeInf->maxDiffRank = -1;
    }

    pthread_exit(nullptr);
}

void * executeTasks(void * args){
    // Executor deals
    auto TL = (TaskList*) args;

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    double localResult = 0;

    pthread_create(&pthreads[1], &attr, waitingNewTasks, (void*) TL);

    for (int i = 0; i < TASKS_NUMBER; ++i) {
        for (int j = 0; j < TL[i].repeatNumber; ++j) {
//            localResult += sin(i);
            localResult += 1;
        }
    }

    pthread_join(pthreads[1], nullptr);

    // Analyzer вычисляет, у какого процесса по сравнению с этим больше всего разница в незавершенных заданиях
    // Loader начинает подгрузку заданий с подходящего процесса на фоне вычислений

    globalResult += localResult;

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

        // Update global result on this iteration
        double commonResult;
        MPI_Allreduce(&globalResult, &commonResult, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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

        MPI_Barrier(MPI_COMM_WORLD);
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

    doTasks(iCounter, TL, pRank, pSize);

    MPI_Finalize();
}