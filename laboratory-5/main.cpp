#include <iostream>
#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <vector>

#define TASKS_LIST_COUNT 4
#define DEFAULT_TASKS_NUMBER 100
#define REQUEST_TAG 900
#define RESPONSE_TAG 90
#define TASK_TAG 9
#define NO_AVAILABLE_TASKS 2
#define STOP_FL -1
#define INADMISSIBLE_VALUE 5
#define L 77777

struct Task {
    int repeatNumber;

    explicit Task(int num) : repeatNumber(num) {}
    Task() : Task(0) {}
};

using namespace std;

int pSize;
int pRank;
int iTask;
int iCounter;

double globalResult;
double ImbalanceTime;
double ImbalancePortion;

int tasksReceived;
int tasksSent;

vector<Task> TaskList;

pthread_t receiverThread;
pthread_t executorThread;

pthread_mutex_t TaskListMutex;

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

void executeTasks(Task task){
    double intermediateValue;

    for(int i = 0; i  < task.repeatNumber ; ++i){
        intermediateValue = sqrt(task.repeatNumber) * sqrt(task.repeatNumber);
        intermediateValue = pow(intermediateValue, -3.4);
        intermediateValue /= 29062001.47;
        intermediateValue += 21301.47*324.342;
        intermediateValue -= 21132.21*823.702;
        globalResult += sqrt(intermediateValue) * 1e9;
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
        return NO_AVAILABLE_TASKS;
    }
    int availabilityFl = 1;

    MPI_Send(&availabilityFl,1, MPI_INT, requestProcess, REQUEST_TAG, MPI_COMM_WORLD);
    MPI_Recv(&availabilityFl,1, MPI_INT, requestProcess, RESPONSE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if(!availabilityFl){
        return NO_AVAILABLE_TASKS;
    }

    Task receivedTask;
    MPI_Recv(&receivedTask, sizeof(receivedTask), MPI_BYTE, requestProcess, TASK_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    pthread_mutex_lock(&TaskListMutex);
        TaskList.push_back(receivedTask);
    pthread_mutex_unlock(&TaskListMutex);

    tasksReceived++;

    return availabilityFl;
}

void * doTasks(void * args){
    int requestProcess;
    int addition = 0;
    double startTime, endTime, timeSpent;

    for(iCounter = 0; iCounter < TASKS_LIST_COUNT; ++iCounter, iTask = 0, addition = 0){
        tasksReceived = 0;

        pthread_mutex_lock(&TaskListMutex);
        fillTaskList();

        startTime = MPI_Wtime();
        while(addition < pSize){
            while(iTask < TaskList.size()){
                pthread_mutex_unlock(&TaskListMutex);
                    executeTasks(TaskList[iTask]);
                pthread_mutex_lock(&TaskListMutex);

                ++iTask;
            }

            pthread_mutex_unlock(&TaskListMutex);

            requestProcess = (pRank + addition) % pSize;
            if(receiveTasks(requestProcess) == NO_AVAILABLE_TASKS){
                ++addition;
            }
            pthread_mutex_lock(&TaskListMutex);
        }

        pthread_mutex_unlock(&TaskListMutex);
        endTime = MPI_Wtime();

        timeSpent = endTime - startTime;

        printResults(timeSpent);

        tasksSent = 0;
        MPI_Barrier(MPI_COMM_WORLD);

        if(!pRank) {
            std::cout << std::endl;
        }

        calculateImbalanceTimeAndPortion(timeSpent);

        if(pRank == pSize-1){
            std::cout << "Imbalance time: " << ImbalanceTime << "; Imbalance portion: " << ImbalancePortion << std::endl;
        }
    }

    int stopFlag = STOP_FL;
    MPI_Send(&stopFlag, 1, MPI_INT, pRank, REQUEST_TAG, MPI_COMM_WORLD);

    return nullptr;
}

void* receiveRequests(void* args){
    int availabilityFl;
    MPI_Status status;
    while(iCounter < TASKS_LIST_COUNT){
        MPI_Recv(&availabilityFl, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST_TAG, MPI_COMM_WORLD, &status);
        if(availabilityFl == STOP_FL){
            return nullptr;
        }

        pthread_mutex_lock(&TaskListMutex);
            if(iTask > TaskList.size() - INADMISSIBLE_VALUE){
                availabilityFl = 0;
            }
        pthread_mutex_unlock(&TaskListMutex);

        MPI_Send(&availabilityFl, 1, MPI_INT, status.MPI_SOURCE, RESPONSE_TAG, MPI_COMM_WORLD);
        if(!availabilityFl){
            continue;
        }

        pthread_mutex_lock(&TaskListMutex);
            Task taskToSend = TaskList.back();
            TaskList.pop_back();
        pthread_mutex_unlock(&TaskListMutex);

        MPI_Send(&taskToSend, sizeof(Task), MPI_BYTE, status.MPI_SOURCE, TASK_TAG, MPI_COMM_WORLD);

        ++tasksSent;
    }

    return nullptr;
}

void createThreads(){
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
}

int main(int argc, char* argv[]){
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    MPI_Comm_size(MPI_COMM_WORLD, &pSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &pRank);

    pthread_mutex_init(&TaskListMutex, nullptr);

    createThreads();

    pthread_mutex_destroy(&TaskListMutex);

    MPI_Finalize();

    return EXIT_SUCCESS;
}