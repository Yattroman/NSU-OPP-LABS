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
    cout << "Time spent: " << time << endl;
}

void fillTasksList(int pRank, int pSize, TaskList * TL, int iCounter){
    for (int i = 0; i < TASKS_NUMBER; ++i) {
        TL[i].repeatNumber = abs(TASKS_NUMBER/2 - i) * abs(pRank - iCounter % pSize) * L;
    }
}

void doTasks(int iCounter, TaskList * TL, double globalResult, int pRank, int pSize){
    while(iCounter < MAX_LISTS_COUNT){
        fillTasksList(pRank, pSize, TL, iCounter);
        for (int i = 0; i < TASKS_NUMBER; ++i) {
            for (int j = 0; j < TL[i].repeatNumber; ++j) {
                globalResult += sin(i);
            }
        }
        // Синхронизация процессов
        // Вывод результатов с этой итерации
        iCounter++;
    }
}

int main(int argc, char ** argv){
    int pRank, pSize;
    int iCounter = 0; // Iterations Counter
    double globalResult = 0;

    TaskList TL[TASKS_NUMBER];
}