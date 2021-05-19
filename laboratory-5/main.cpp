#include <iostream>
#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <vector>

#define TASKS_LIST_COUNT 4
#define DEFAULT_TASKS_NUMBER 100

#define REQUEST_TAG 987
#define RESPONSE_TAG 98
#define TASK_TAG 9

#define NO_TASKS_IN_STOCK 2
#define STOP_FL -1
#define L 77777

#define DEFAULT_TASK_VALUE 0

struct Task {
    int repeatNumber;

    explicit Task(int num) : repeatNumber(num) {}
    Task() : Task(DEFAULT_TASK_VALUE) {}
};

using namespace std;

int pSize;
int pRank;
int iTask;
int iCounter;

int inadmissibleValue;

double globalResult;
double ImbalanceTime;
double ImbalancePortion;

int tasksReceived;
int tasksSent;

vector<Task> TaskList;

pthread_t receiverThread;
pthread_t executorThread;

pthread_mutex_t ExecutorMutex;
pthread_mutex_t ReceiverMutex;

void errorHandler(){
    MPI_Finalize();
    exit(EXIT_FAILURE);
}

void printResults(double& timeSpent){
    cout << "Process rank:\t" << pRank << "; Iteration:\t" << iCounter << "; Tasks completed:\t" << iTask << "; Global result:\t" << globalResult
    << "; Time spent:\t" << timeSpent << "; Task sent:\t" << tasksSent << "; Task received:\t" << tasksReceived << endl;
}

void fillTaskList(){
    TaskList.clear();
    for (int i = 0; i < DEFAULT_TASKS_NUMBER; ++i) {
        TaskList.emplace_back(abs(DEFAULT_TASKS_NUMBER/2 - i) * abs(pRank - iCounter % pSize) * L);
    }
}

void subTask1(Task task, double& intermediateValue){
    intermediateValue = pow(sqrt(task.repeatNumber) * sqrt(task.repeatNumber), -3.4);
}

void subTask2(double& intermediateValue){
    intermediateValue = intermediateValue / 29062001.47;
}

void subTaskFinal(double& intermediateValue){
    globalResult = globalResult + 1e9 * pow(intermediateValue, 1/2);
}

void executeTasks(Task task){
    double intermediateValue;

    for(int i = 0; i  < task.repeatNumber; ++i){
        subTask1(task, intermediateValue);
        subTask2(intermediateValue);
        subTaskFinal(intermediateValue);
    }
}

void calculateImbalanceTimeAndPortion(double &timeSpent){
    double maxTime, minTime;
    MPI_Allreduce(&timeSpent, &maxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&timeSpent, &minTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    ImbalanceTime = maxTime - minTime;
    ImbalancePortion = ImbalanceTime / maxTime * 100;
}

int receiveTasks(int requestProcess){
    if(requestProcess == pRank){
        return NO_TASKS_IN_STOCK;
    }
    int inStockFl = 1;

    MPI_Send(&inStockFl,1, MPI_INT, requestProcess, REQUEST_TAG, MPI_COMM_WORLD);

    MPI_Recv(&inStockFl,1, MPI_INT, requestProcess, RESPONSE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if(inStockFl == 0){
        return NO_TASKS_IN_STOCK;
    }

    Task receivedTask;
    MPI_Recv(&receivedTask, sizeof(receivedTask), MPI_BYTE, requestProcess, TASK_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    pthread_mutex_lock(&ExecutorMutex);
        TaskList.push_back(receivedTask);
    pthread_mutex_unlock(&ExecutorMutex);

    tasksReceived++;

    return inStockFl;
}

void * doTasks(void * args){ // Executor routine
    int requestProcess;
    int addition = 0;
    double startTime, endTime, timeSpent;

    for(iCounter = 0; iCounter < TASKS_LIST_COUNT; ++iCounter, iTask = 0, addition = 0){
        tasksReceived = 0;

        pthread_mutex_lock(&ExecutorMutex);
        fillTaskList();

        startTime = MPI_Wtime();
        while(addition < pSize){
            while(iTask < TaskList.size()){
                pthread_mutex_unlock(&ExecutorMutex);
                    executeTasks(TaskList[iTask]);
                pthread_mutex_lock(&ExecutorMutex);

                ++iTask;
            }

            pthread_mutex_unlock(&ExecutorMutex);
                requestProcess = (pRank + addition) % pSize;
                if(receiveTasks(requestProcess) == NO_TASKS_IN_STOCK){
                    ++addition;
                }
            pthread_mutex_lock(&ExecutorMutex);
        }

        pthread_mutex_unlock(&ExecutorMutex);
        endTime = MPI_Wtime();

        timeSpent = endTime - startTime;

        printResults(timeSpent);

        tasksSent = 0;
        MPI_Barrier(MPI_COMM_WORLD);

        calculateImbalanceTimeAndPortion(timeSpent);


        if(pRank == 0) {
            std::cout << std::endl;
        }
        if(pRank == pSize-1){
            std::cout << "Imbalance time: " << ImbalanceTime << "; Imbalance portion: " << ImbalancePortion << std::endl;
        }
    }

    int stopFlag = STOP_FL;
    MPI_Send(&stopFlag, 1, MPI_INT, pRank, REQUEST_TAG, MPI_COMM_WORLD);

    pthread_exit(nullptr);
}

void* receiveRequests(void* args){ // Receiver routine
    int inStockFl;
    MPI_Status status;
    while(iCounter < TASKS_LIST_COUNT){
        MPI_Recv(&inStockFl, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST_TAG, MPI_COMM_WORLD, &status);
        if(inStockFl == STOP_FL){
            pthread_exit(nullptr);
        }

        pthread_mutex_lock(&ReceiverMutex);
            if(iTask + inadmissibleValue > TaskList.size()){
                inStockFl = 0;
            }
        pthread_mutex_unlock(&ReceiverMutex);

        MPI_Send(&inStockFl, 1, MPI_INT, status.MPI_SOURCE, RESPONSE_TAG, MPI_COMM_WORLD);
        if(inStockFl == 0){
            continue;
        }

        pthread_mutex_lock(&ReceiverMutex);
            Task taskToSend = TaskList.back();
            TaskList.pop_back();
        pthread_mutex_unlock(&ReceiverMutex);

        MPI_Send(&taskToSend, sizeof(Task), MPI_BYTE, status.MPI_SOURCE, TASK_TAG, MPI_COMM_WORLD);

        ++tasksSent;
    }

    pthread_exit(nullptr);
}

void startWorkWithPthreads(){
    pthread_mutex_init(&ExecutorMutex, nullptr);
    pthread_mutex_init(&ReceiverMutex, nullptr);

    int statusOne = pthread_create(&receiverThread, nullptr, receiveRequests, nullptr);
    int statusTwo = pthread_create(&executorThread, nullptr, doTasks, nullptr);

    if(statusOne && statusTwo){
        std::cerr << "Error creating threads" << std::endl;
        errorHandler();
    }

    statusOne = pthread_join(executorThread,nullptr);
    statusTwo = pthread_join(receiverThread,nullptr);

    if(statusOne && statusTwo){
        std::cerr << "Error joining thread" << std::endl;
        errorHandler();
    }

    pthread_mutex_destroy(&ExecutorMutex);
    pthread_mutex_destroy(&ReceiverMutex);
}

int main(int argc, char* argv[]){
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    MPI_Comm_size(MPI_COMM_WORLD, &pSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &pRank);

    inadmissibleValue = pSize + 1;

    startWorkWithPthreads();

    MPI_Finalize();

    return EXIT_SUCCESS;
}