#include <cstdlib>
#include <cmath>
#include <iostream>
#include <cstring>
#include <pthread.h>
#include <vector>
#include "mpi.h"

#define THREADS_NUM 2
#define REQUEST_TAG 900
#define RESPONSE_TAG 90
#define TASK_TAG 9
#define STOP_FL -1
#define MAX_LISTS_COUNT 3
#define TASKS_NUMBER 10
#define L 10
#define TRUE 1
#define UPLOAD_LIMIT 2
#define NO_AVAILABLE_TASKS 2

using namespace std;

typedef struct Task{
    int repeatNumber;

    explicit Task(int number) : repeatNumber(number) {};
    Task() : repeatNumber(0) {};
} Task;

int pSize;
int pRank;

int iCounter;
unsigned long iTask;

int tasksSent;
int tasksReceived;

double globalResult = 0;

pthread_t receiverThread;
pthread_t executorThread;

pthread_mutex_t taskListMutex;
pthread_mutex_t receiverMutex;

vector<Task> TaskList;

/*typedef struct AnalyzeInfo{
    int uncompletedTasks;
    int procRankRequest;
    int procSize;

    int& maxDiff;
    int& maxDiffRank;
} AnalyzeInfo;*/

void atError(){
    MPI_Finalize();
    exit(EXIT_FAILURE);
}

void printResults(double& timeSpent){
    cout << "Process rank: " << pRank << endl;
    cout << "Tasks completed: " << iTask << endl;
    cout << "Global result: " << globalResult << endl;
    cout << "Time spent: " << timeSpent << endl;

    cout << endl;
}

void fillTasksList(){
    TaskList.clear();

    for (int i = 0; i < TASKS_NUMBER; ++i) {
        TaskList.emplace_back(abs(TASKS_NUMBER/2 - i) * abs(pRank - iCounter % pSize) * L);
    }
}

int analyzeProcessesLoad(int uncompletedTasks, int& maxDiff){
    // Analyze which process has maximal load
    int uncompletedTasksRequest = uncompletedTasks;

    int localDiff;
    int differents[pSize];

    int maxDiffRank = -1;

    memset(differents, -1, pSize*sizeof(int));

    MPI_Bcast(&uncompletedTasksRequest, 1, MPI_INT, pRank, MPI_COMM_WORLD);

    localDiff = uncompletedTasks - uncompletedTasksRequest;
    MPI_Allgather(&localDiff, 1, MPI_INT, differents, 1, MPI_INT, MPI_COMM_WORLD);

    for (int i = 0; i < pSize; ++i) {
        if (differents[i] > maxDiff && maxDiffRank > 3) {
            maxDiff = differents[i];
            maxDiffRank = i;
        }

    }

    return maxDiffRank;
}

void executeTasks(Task task){
    double localResult = 0;

    for (int i = 0; i < TASKS_NUMBER; ++i) {
        for (int j = 0; j < task.repeatNumber; ++j) {
            localResult += 1;
        }
    }

    globalResult += localResult;
}

int receiveTasks(int reqProc){
    if(reqProc == pRank || reqProc == -1){
        return NO_AVAILABLE_TASKS;
    }
    int availabilityFl = 1;

    MPI_Send(&availabilityFl, 1, MPI_INT, reqProc, REQUEST_TAG, MPI_COMM_WORLD);
    MPI_Recv(&availabilityFl, 1, MPI_INT, reqProc, REQUEST_TAG, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

    if(availabilityFl == 0){
        return NO_AVAILABLE_TASKS;
    }

    // Receive task
    Task receivedTask;
    MPI_Recv(&receivedTask, sizeof(receivedTask), MPI_BYTE, reqProc, TASK_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    pthread_mutex_lock(&taskListMutex);
        TaskList.push_back(receivedTask);
    pthread_mutex_unlock(&taskListMutex);

    tasksReceived++;

    return availabilityFl;
}

void * receiveRequests(void * args){
    int availabilityFl;
    MPI_Status status;

    while(iCounter < MAX_LISTS_COUNT){

        MPI_Recv(&availabilityFl, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST_TAG, MPI_COMM_WORLD, &status);

        if(availabilityFl == STOP_FL){
            return nullptr;
        }

        if(TaskList.size() - iTask > UPLOAD_LIMIT){
            availabilityFl = 0;
        }

        MPI_Send(&availabilityFl, 1, MPI_INT, status.MPI_SOURCE, RESPONSE_TAG, MPI_COMM_WORLD);

        if(availabilityFl == 0){
            continue;
        }

        pthread_mutex_lock(&receiverMutex);
            Task taskToSend = TaskList.back();
            TaskList.pop_back();
        pthread_mutex_unlock(&receiverMutex);

        MPI_Send(&taskToSend, sizeof(Task), MPI_BYTE, status.MPI_SOURCE, TASK_TAG, MPI_COMM_WORLD);
        ++tasksSent;
    }

    pthread_exit(nullptr);
}

void doTasks(){
    double startTime, endTime, timeSpent;
    int requestProcess;
    int maxDiff = 0;
    int uncompletedTasks;
    for(iCounter = 0; iCounter < MAX_LISTS_COUNT; ++iCounter){
        pthread_mutex_lock(&taskListMutex);

        fillTasksList();
        startTime = MPI_Wtime();

        while(TRUE){
            for(; iTask < TaskList.size(); ++iTask){
                pthread_mutex_unlock(&taskListMutex);
                // Execute tasks by main (executor) thread
                executeTasks(TaskList[iTask]);
                pthread_mutex_lock(&taskListMutex);
            }
            pthread_mutex_unlock(&taskListMutex);

            /*requestProc = (pRank + proc) % pSize;
            if(receiveTasks(requestProc) == NO_TASKS_AVAILABLE){
                ++proc;
            }*/

            uncompletedTasks = TaskList.size() - 1 - iTask;
            requestProcess = analyzeProcessesLoad(uncompletedTasks, maxDiff);

            cout << "'" << "!" << "'" << endl;

            if(receiveTasks(requestProcess) == NO_AVAILABLE_TASKS){
                break;
            }

            pthread_mutex_lock(&taskListMutex);
        }
        pthread_mutex_unlock(&taskListMutex);

        endTime = MPI_Wtime();
        timeSpent = endTime - startTime;

        printResults(timeSpent);

        MPI_Barrier(MPI_COMM_WORLD);
    }

    int stopFl= STOP_FL;
    MPI_Send(&stopFl, 1, MPI_INT, pRank, REQUEST_TAG, MPI_COMM_WORLD);
}

void initThreads(){
    if(pthread_create(&receiverThread, nullptr, receiveRequests, nullptr)){
        std::cerr << "Error creating thread" << std::endl;
        atError();
    }

    doTasks();

    if(pthread_join(receiverThread, nullptr)){
        std::cerr << "Error joining thread" << std::endl;
        atError();
    }
}

void go(){
    MPI_Comm_size(MPI_COMM_WORLD, &pSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &pRank);

    pthread_mutex_init(&taskListMutex, nullptr);
    initThreads();
    pthread_mutex_destroy(&taskListMutex);
}

int main(int argc, char ** argv){
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    go();

    MPI_Finalize();

    return EXIT_SUCCESS;
}