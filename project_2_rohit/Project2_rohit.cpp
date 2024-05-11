#include <iostream>
#include <vector>
#include <stack>
#include <algorithm>
#include <limits>
#include <numeric>
#include <iomanip>
#include <cstring>
#include <chrono>

using namespace std;

struct Task {
  int Taskid;  
  vector<int> TasklocalExecutions; 
  vector<int> TaskpreceedingId; 
  vector<int> TasksucceedingId; 
  vector<int> TaskfinishingTimes;  
  vector<int> TaskreadyingTimes;  
  int TaskstartingTime; 
  int allotedCore; 
  double Taskweight; 
  double Taskpriority; 
};
Task task{0, {}, {}, {}, {}, {}, 0, -1, 0.0, 0.0};


const int coresNumber = 3; //number of cores
const int tasksNumber = 20;  //Number of tasks                       (change this to 20 for example 3,4,5)
vector<int> cloudTimes = {3, 1, 1}; //Cloud timing
vector<float> corePower = { 1, 2, 4, 0.5}; //Core1, core 2, core 3, cloud sending power


vector<Task> PerformTasks() {
//{Task(i), {Core1,core2,core3},{Task(i-1)},{Task(i+1)}}
//Example 1
return {
    {},
    {1, {9, 7, 5}, {}, {2, 3, 4, 5, 6}},
    {2, {8, 6, 5}, {1}, {8, 9}},
    {3, {6, 5, 4}, {1}, {7}},
    {4, {7, 5, 3}, {1}, {8, 9}},
    {5, {5, 4, 2}, {1}, {9}},
    {6, {7, 6, 4}, {1}, {8}},
    {7, {8, 5, 3}, {3}, {10}},
    {8, {6, 4, 2}, {2, 4, 6}, {10}},
    {9, {5, 3, 2}, {2, 4, 5}, {10}},
    {10, {7, 4, 2}, {7, 8, 9}, {}}
};


//Example 2
// return {
//     {},
//     {1, {9, 7, 5}, {}, {2, 3, 4}},
//     {2, {8, 6, 5}, {1}, {5}},
//     {3, {6, 5, 4}, {1}, {6, 7}},
//     {4, {7, 5, 3}, {1}, {5}},

//     {5, {5, 4, 2}, {2, 4}, {8, 9}},
//     {6, {7, 6, 4}, {3}, {8}},
//     {7, {8, 5, 3}, {3}, {9}},

//     {8, {6, 4, 2}, {5, 6}, {10}},
//     {9, {5, 3, 2}, {5, 7}, {10}},
//     {10, {7, 4, 2}, {8, 9}, {}}
// };




//Example 3
// return {
//     {},
//     {1, {9, 7, 5}, {}, {2, 3, 4, 5, 6}},

//     {2, {8, 6, 5}, {1}, {8}},
//     {3, {6, 5, 4}, {1}, {9}},
//     {4, {7, 5, 3}, {1}, {10}},
//     {5, {5, 4, 2}, {1}, {7}},
//     {6, {7, 6, 4}, {1}, {8}},

//     {7, {8, 5, 3}, {5}, {12}},
//     {8, {6, 4, 2}, {2, 6}, {13}},
//     {9, {5, 3, 2}, {3}, {14}},
//     {10, {7, 4, 2}, {4}, {11}},

//     {11, {4, 3, 2}, {10}, {16}},
//     {12, {5, 2, 1}, {7}, {17}},
//     {13, {6, 3, 2}, {8}, {15}},
//     {14, {7, 3, 2}, {9}, {16}},

//     {15, {5, 3, 2}, {13}, {18}},
//     {16, {7, 2, 1}, {11, 14}, {19}},
//     {17, {7, 5, 4}, {12}, {18}},

//     {18, {8, 6, 3}, {15,17}, {20}},
//     {19, {5, 4, 3}, {16}, {20}},

//     {20, {4, 3, 2}, {18, 19}, {}}
// };



// Example 4
// return {
//     {},
//     {1, {9, 7, 5}, {}, {4, 5, 6}},
//     {2, {8, 6, 5}, {}, {5,6}},
//     {3, {6, 5, 4}, {}, {4,6}},

//     {4, {7, 5, 3}, {1,3}, {8,10}},
//     {5, {5, 4, 2}, {1,2}, {7,9,10}},
//     {6, {7, 6, 4}, {1,2,3}, {8}},

//     {7, {8, 5, 3}, {5}, {11,12}},
//     {8, {6, 4, 2}, {4,6}, {12,13}},
//     {9, {5, 3, 2}, {5}, {11,13,14}},
//     {10, {7, 4, 2}, {4,5}, {11}},

//     {11, {4, 3, 2}, {7,9,10}, {15,16,17}},
//     {12, {5, 2, 1}, {7,8}, {17}},
//     {13, {6, 3, 2}, {8,9}, {15,17}},
//     {14, {7, 3, 2}, {9}, {16,17}},

//     {15, {5, 3, 2}, {11,13}, {18,19}},
//     {16, {7, 2, 1}, {11, 14}, {18,19}},
//     {17, {7, 5, 4}, {11,12,13,14}, {18}},

//     {18, {8, 6, 3}, {15,16,17}, {20}},
//     {19, {5, 4, 3}, {15,16}, {20}},

//     {20, {4, 3, 2}, {18, 19}, {}}
// };



// Example 5
// return {
//     {},
//     {1, {9, 7, 5}, {}, {4, 5, 6}},
//     {2, {8, 6, 5}, {}, {5,6}},
//     {3, {6, 5, 4}, {}, {4,6}},

//     {4, {7, 5, 3}, {1,3}, {8,10}},
//     {5, {5, 4, 2}, {1,2}, {7,9,10}},
//     {6, {7, 6, 4}, {1,2,3}, {8}},

//     {7, {8, 5, 3}, {5}, {11,12}},
//     {8, {6, 4, 2}, {4,6}, {12,13}},
//     {9, {5, 3, 2}, {5}, {11,13,14}},
//     {10, {7, 4, 2}, {4,5}, {11}},

//     {11, {4, 3, 2}, {7,9,10}, {15,16,17}},
//     {12, {5, 2, 1}, {7,8}, {17}},
//     {13, {6, 3, 2}, {8,9}, {15,17}},
//     {14, {7, 3, 2}, {9}, {16,17}},

//     {15, {5, 3, 2}, {11,13}, {18,19,20}},
//     {16, {7, 2, 1}, {11, 14}, {18,19,20}},
//     {17, {7, 5, 4}, {11,12,13,14}, {18,20}},

//     {18, {8, 6, 3}, {15,16,17}, {}},
//     {19, {5, 4, 3}, {15,16}, {}},

//     {20, {4, 3, 2}, {15,16,17}, {}}
// };

}













vector<Task> tasks;
vector<int> wirelessSendQueue, exitTasks;
vector<vector<int>> cores(coresNumber);
void setInitialTimes(vector<Task>& tasks);
void updateTaskPriority(Task &task);
int calculateSendResponseTime(const vector<Task> &tasks, const vector<int> &queue);
int calculateCoreResponseTime(const vector<Task> &tasks, const vector<int> &queue);
void calculateCloudTimes(const vector<Task> &tasks, Task &task, const vector<int> wirelessQueue);
void calculateCoreTimes(const vector<Task> &tasks, Task &task, const vector<int> coreQueue, int coreNumber);
double calculateTotalEnergy(const vector<Task> &tasks, const vector<int> &sendQueue, const vector<vector<int>> &coreQueue);
int calculateTotalTime(const vector<Task> &tasks);
void printTaskAssignment(const vector<Task> &tasks, const vector<int> &wirelessQueue, const vector<vector<int>> &coreQueues, double energyTotal, int timeTotal);
void removeTaskFromCore(vector<vector<int>> &scheduling, int core, int taskId);
void addTaskToCore(vector<vector<int>> &scheduling, int core, int taskId, const vector<Task> &tasks);
void resetTaskScheduling(vector<Task> &tasks);
void PerformTasksScheduling(stack<int> &toassignTaskToCores, vector<int> &temporaryVectorF, vector<int> &temporaryVectorS, const vector<Task> &tasks, const vector<vector<int>> & currentTask, int tasksNumber, int coresNumber);
void processScheduling(vector<vector<int>> &schedulingQueue, stack<int> &toassignTaskToCores, vector<Task> &copyTasks, vector<int> &temporaryVectorF, vector<int> &temporaryVectorS, vector<vector<int>> & currentTask, int tasksNumber, int coresNumber);
void calculateTotalEnergyAndTime(const vector<Task>& tasks, const vector<int>& wirelessQueue, const vector<vector<int>>& schedulingQueue, double &energyTotal, int &timeTotal);
vector<Task> rearrangeTask(int targetTask, int targetCore, double &energyTotal, int &timeTotal, vector<vector<int>> currentTask, vector<vector<vector<int>>> &schedulingQueues);
void FindExitTasks(const vector<Task>& tasks, vector<int>& exitTasks);
void assignTasksAndCalculateLoad(vector<Task>& tasks, const vector<int>& cloudTimes, int coresNumber);
vector<Task> sortTasksByTaskPriority(const vector<Task>& tasks);
void assignTaskToCore(Task& task, const vector<Task>& tasks, vector<int>& wirelessSendQueue, vector<vector<int>>& cores, int coresNumber);
void performInitialAssignment(vector<Task>& tasks, const vector<int>& cloudTimes, int coresNumber, vector<int>& wirelessSendQueue, vector<vector<int>>& cores);




int main() {
    using namespace std::chrono;
    auto startTime = high_resolution_clock::now();

    tasks = PerformTasks();
    setInitialTimes(tasks);
    FindExitTasks(tasks, exitTasks);

    double initialEnergyTotal, finalEnergyTotal;
    int initialTimeTotal, finalTimeTotal;

    performInitialAssignment(tasks, cloudTimes, coresNumber, wirelessSendQueue, cores);
    initialEnergyTotal = calculateTotalEnergy(tasks, wirelessSendQueue, cores);
    initialTimeTotal = calculateTotalTime(tasks);

    double energyTotal = initialEnergyTotal;
    int timeTotal = initialTimeTotal;

    auto initialTimeStop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(initialTimeStop - startTime);

    cout << string(155, '~') << endl;
    cout << "                                                            Example 5: Initial Scheduling Table                                                            \n";
    cout << string(155, '~') << endl;

    printTaskAssignment(tasks, wirelessSendQueue, cores, energyTotal, timeTotal);
    cout << "Initial Scheduling algorithm Execution Time: " << duration.count() << " ms\n" << endl;

    double T_max = 1.5 * timeTotal; 
    vector<Task> taskSorted = sortTasksByTaskPriority(tasks);
    auto finalTimeStart = high_resolution_clock::now();

    for (int run = 0; run < 1; ++run) { 
        for (const Task &task : taskSorted) {
            vector<double> energies;
            vector<int> times;
            vector<vector<Task>> assignments;
            vector<vector<vector<int>>> schedulingQueues;
            vector<vector<int>> currentTask(coresNumber + 1);
            currentTask[0] = wirelessSendQueue;
            for (int i = 1; i <= coresNumber; ++i) {
                currentTask[i] = cores[i - 1];
            }
            for (int targetCore = 0; targetCore <= coresNumber; ++targetCore) {
                double testingEnergyTotal = energyTotal; 
                int testingTimeTotal = timeTotal;      

                assignments.push_back(rearrangeTask(task.Taskid, targetCore, testingEnergyTotal, testingTimeTotal, currentTask, schedulingQueues));
                energies.push_back(testingEnergyTotal);
                times.push_back(testingTimeTotal);
            }
            int assignmentBest = tasks[task.Taskid].allotedCore; 
            double best_energy = energyTotal; 
            for (int targetCore = 0; targetCore <= coresNumber; ++targetCore) {

                if (times[targetCore] <= timeTotal && energies[targetCore] < best_energy) {
                    assignmentBest = targetCore;
                    best_energy = energies[targetCore];
                }
            }
            if (assignmentBest == tasks[task.Taskid].allotedCore) {
                double best_ratio = 0;

                for (int targetCore = 0; targetCore <= coresNumber; ++targetCore) {

                    if (targetCore == tasks[task.Taskid].allotedCore) {
                        continue;
                    }

                    if (times[targetCore] <= T_max && energies[targetCore] < energyTotal) {
                        double ratioOfEnergyReduction = (energyTotal - energies[targetCore]) / (times[targetCore] - timeTotal);
                        if (ratioOfEnergyReduction > best_ratio) {
                            assignmentBest = targetCore;
                            best_ratio = ratioOfEnergyReduction;
                        }
                    }
                }
            }
            if (assignmentBest != tasks[task.Taskid].allotedCore) {
                tasks = assignments[assignmentBest];
                wirelessSendQueue = schedulingQueues[assignmentBest][0];
                schedulingQueues[assignmentBest].erase(schedulingQueues[assignmentBest].begin());
                cores = schedulingQueues[assignmentBest];
                energyTotal = energies[assignmentBest];
                timeTotal = times[assignmentBest];
            }
        }
    }

    finalEnergyTotal = calculateTotalEnergy(tasks, wirelessSendQueue, cores);
    finalTimeTotal = calculateTotalTime(tasks);
    auto finalTimeStop = high_resolution_clock::now();
    auto finalDuration = duration_cast<microseconds>(finalTimeStop - finalTimeStart);


    cout << string(155, '~') << endl;
    cout << "                                                              Example 5: Final Scheduling Table                                                            \n";
    cout << string(155, '~') << endl;
    printTaskAssignment(tasks, wirelessSendQueue, cores, energyTotal, timeTotal);
    cout << "Final Scheduling algorithm Execution Time: " << duration.count() << " ms\n" << endl;
    cout << string(37, '~') << endl;
    cout << "       Example 5: Summary Table      \n";
    cout << string(37, '~') << endl;
    cout << left << setw(20) << "      " << setw(10) << "Initial" << setw(10) << "Final" << endl;
    cout << left << setw(20) << "Total Energy" << setw(10) << initialEnergyTotal << setw(10) << finalEnergyTotal << endl;
    cout << left << setw(20) << "Total Time" << setw(10) << initialTimeTotal << setw(10) << finalTimeTotal << endl;
    cout << string(37, '~') << endl;

    return 0;
}






































void setInitialTimes(vector<Task>& tasks) {
    for (Task& task : tasks) {
        task.TaskfinishingTimes.resize(4, 0); 
        task.TaskreadyingTimes.resize(4, 0); 
    }
}



void updateTaskPriority(Task &task) {

if (task.Taskpriority == 0) {

    if (task.TasksucceedingId.empty()) {
    task.Taskpriority = task.Taskweight;
    } else {
    double maxPriority = 0;

    for (int successor : task.TasksucceedingId) {
        updateTaskPriority(tasks[successor]);

        maxPriority = max(maxPriority, tasks[successor].Taskpriority);
    }

    task.Taskpriority = task.Taskweight + maxPriority;
    }
}
}



int calculateSendResponseTime(const vector<Task> &tasks, const vector<int> &queue) {
if (!queue.empty()) {
    int lastIndex = queue.back();

    return tasks[lastIndex].TaskfinishingTimes[1];
}

return 0;
}



int calculateCoreResponseTime(const vector<Task> &tasks, const vector<int> &queue) {
if (!queue.empty()) {
    int lastIndex = queue.back();

    return tasks[lastIndex].TaskfinishingTimes[0];
}

return 0;
}



void calculateCloudTimes(const vector<Task> &tasks, Task &task, const vector<int> wirelessQueue) {

int maximumFTSending = 0;
for (int j : task.TaskpreceedingId) {
    int max_FT = max(tasks[j].TaskfinishingTimes[0], tasks[j].TaskfinishingTimes[1]);
    maximumFTSending = max(maximumFTSending, max_FT);
}

task.TaskreadyingTimes[1] = maximumFTSending;
task.TaskstartingTime = max(task.TaskreadyingTimes[1], calculateSendResponseTime(tasks, wirelessQueue));
task.TaskfinishingTimes[1] = task.TaskstartingTime + cloudTimes[0];

int maximumFTCloud = 0;
for (int j : task.TaskpreceedingId) {
    maximumFTCloud = max(maximumFTCloud, tasks[j].TaskfinishingTimes[2]);
}

task.TaskreadyingTimes[2] = max(task.TaskfinishingTimes[1], maximumFTCloud);
task.TaskfinishingTimes[2] = task.TaskreadyingTimes[2] + cloudTimes[1];

task.TaskreadyingTimes[3] = task.TaskfinishingTimes[2];
task.TaskfinishingTimes[3] = task.TaskfinishingTimes[2] + cloudTimes[2];
}




void calculateCoreTimes(const vector<Task> &tasks, Task &task, const vector<int> coreQueue, int coreNumber) {
int maxPredFT = 0;  

for (int j : task.TaskpreceedingId) {

    int currentMaxFT = (tasks[j].TaskfinishingTimes[0] > tasks[j].TaskfinishingTimes[3]) ? tasks[j].TaskfinishingTimes[0] : tasks[j].TaskfinishingTimes[3];
    if (currentMaxFT > maxPredFT) {
    maxPredFT = currentMaxFT;
    }
}

task.TaskreadyingTimes[0] = maxPredFT;

task.TaskstartingTime = max(task.TaskreadyingTimes[0], calculateCoreResponseTime(tasks, coreQueue));

task.TaskfinishingTimes[0] = task.TaskstartingTime + task.TasklocalExecutions[coreNumber];
}




double calculateTotalEnergy(const vector<Task> &tasks, const vector<int> &sendQueue, const vector<vector<int>> &coreQueue) {
double energy = 0.0;

double wirelessEnergy = static_cast<double>(sendQueue.size()) * cloudTimes[0] * corePower[3];
energy += wirelessEnergy;

for (int i = 0; i < coresNumber; ++i) {
    for (int j : coreQueue[i]) {

    double executionTime = tasks[j].TaskfinishingTimes[0] - tasks[j].TaskstartingTime;

    double coreEnergy = executionTime * corePower[i];
    energy += coreEnergy;
    }
}

return energy;
}






int calculateTotalTime(const vector<Task> &tasks) {
int maxOverallFT = 0;  

for (int i : exitTasks) {

    int maxFTtask = max(tasks[i].TaskfinishingTimes[0], tasks[i].TaskfinishingTimes[3]);

    maxOverallFT = max(maxOverallFT, maxFTtask);
}

return maxOverallFT;
}




void printTaskAssignment(const vector<Task> &tasks, const vector<int> &wirelessQueue,
                const vector<vector<int>> &coreQueues, double energyTotal, int timeTotal) {
int maxTimeline = timeTotal;
int numCores = coreQueues.size();

maxTimeline = max(maxTimeline, 45);

int headerLength = 0;
for (int core = 0; core < numCores; ++core) {
    headerLength = max(headerLength, static_cast<int>(to_string(core + 1).length()));
}
headerLength = max(headerLength, static_cast<int>(strlen("WL Sending")));

cout << "TIME            :";
for (int i = 0; i <= maxTimeline; ++i) {
    cout << setw(3) << i;
}
cout << endl;

cout << string(155, '~') << endl;

for (int core = 0; core < numCores; ++core) {
    cout << "CORE " << core + 1 << setw(5) << "          :";
    for (int i = 0; i <= maxTimeline; ++i) { 
        bool taskScheduled = false;
        for (int taskId : coreQueues[core]) {
            const Task &task = tasks[taskId];
            if (task.TaskstartingTime <= i && i < task.TaskfinishingTimes[0]) {
                cout << setw(3) << taskId;
                taskScheduled = true;
                break;
            }
        }
        if (!taskScheduled) {
            cout << "  -";
        }
    }
        cout << endl;
    }

    cout << "CLOUD Sending   :";
    for (int i = 0; i <= maxTimeline; ++i) { 
        bool taskScheduled = false;
        for (int taskId : wirelessQueue) {
            const Task &task = tasks[taskId];
            if (task.TaskstartingTime <= i && i < task.TaskfinishingTimes[1]) {
                cout << setw(3) << taskId;
                taskScheduled = true;
                break;
            }
        }
        if (!taskScheduled) {
            cout << "  -";
        }
    }
    cout << endl;

    cout << "CLOUD Computing :";
    for (int i = 0; i <= maxTimeline; ++i) { 
        bool taskScheduled = false;
        for (int taskId : wirelessQueue) {
            const Task &task = tasks[taskId];
            if (task.TaskreadyingTimes[2] <= i && i < task.TaskfinishingTimes[2]) {
                cout << setw(3) << taskId;
                taskScheduled = true;
                break;
            }
        }
        if (!taskScheduled) {
            cout << "  -";
        }
    }
    cout << endl;

    cout << "CLOUD Receiving :";
    for (int i = 0; i <= maxTimeline; ++i) { 
        bool taskScheduled = false;
        for (int taskId : wirelessQueue) {
            const Task &task = tasks[taskId];
            if (task.TaskreadyingTimes[3] <= i && i < task.TaskfinishingTimes[3]) {
                cout << setw(3) << taskId;
                taskScheduled = true;
                break;
            }
        }
        if (!taskScheduled) {
            cout << "  -";
        }
    }
    cout << endl;

    cout << string(155, '~') << endl;
    cout << "Total Energy : " << energyTotal << endl;
    cout << "Total Time : " << timeTotal << endl;
}

















void removeTaskFromCore(vector<vector<int>> &scheduling, int core, int taskId) {
    auto &queue = scheduling[core];
    queue.erase(remove(queue.begin(), queue.end(), taskId), queue.end());
}



void addTaskToCore(vector<vector<int>> &scheduling, int core, int taskId, const vector<Task> &tasks) {
    auto &queue = scheduling[core];

    auto insertPos = find_if(queue.begin(), queue.end(),
                            [&](int Taskid) { return tasks[Taskid].TaskstartingTime >= tasks[taskId].TaskreadyingTimes[0]; });

    queue.insert(insertPos, taskId);
}



void resetTaskScheduling(vector<Task> &tasks) {
    for (auto &task : tasks) {
        fill(task.TaskreadyingTimes.begin(), task.TaskreadyingTimes.end(), 0);
        fill(task.TaskfinishingTimes.begin(), task.TaskfinishingTimes.end(), 0);
        task.TaskstartingTime = 0;
    }
}




void PerformTasksScheduling(stack<int> &toassignTaskToCores, vector<int> &temporaryVectorF, vector<int> &temporaryVectorS, const vector<Task> &tasks, const vector<vector<int>> & currentTask, int tasksNumber, int coresNumber) {

    temporaryVectorF.push_back(1); 
    for (int i = 1; i <= tasksNumber; ++i) {
        temporaryVectorF.push_back(static_cast<int>(tasks[i].TaskpreceedingId.size()));
    }

    temporaryVectorS.assign(tasksNumber + 1, 1);
    for (int i = 0; i <= coresNumber; ++i) {
        if (!currentTask[i].empty()) {
            temporaryVectorS[currentTask[i][0]] = 0;
        }
    }

    for (int i = 0; i <= tasksNumber; ++i) {
        if (temporaryVectorF[i] == 0 && temporaryVectorS[i] == 0) {
            toassignTaskToCores.push(i);
            temporaryVectorF[i] = temporaryVectorS[i] = -1; 
        }
    }
}




void processScheduling(vector<vector<int>> &schedulingQueue, stack<int> &toassignTaskToCores, vector<Task> &copyTasks, vector<int> &temporaryVectorF, vector<int> &temporaryVectorS, vector<vector<int>> & currentTask, int tasksNumber, int coresNumber) {

    schedulingQueue.resize(coresNumber + 1);

    while (!toassignTaskToCores.empty()) {
        Task &task = copyTasks[toassignTaskToCores.top()];
        toassignTaskToCores.pop();

        int k = task.allotedCore;
        if (k == 0) {

            calculateCloudTimes(copyTasks, task, schedulingQueue[k]);
        } else {

            calculateCoreTimes(copyTasks, task, schedulingQueue[k], k - 1);
        }
        schedulingQueue[k].push_back(task.Taskid);

        for (int i : task.TasksucceedingId) {
            --temporaryVectorF[i];
        }

        currentTask[k].erase(currentTask[k].begin());

        if (!currentTask[k].empty()) {
            temporaryVectorS[currentTask[k][0]] = 0;
        }

        for (int i = 1; i <= tasksNumber; ++i) {
            if (temporaryVectorF[i] == 0 && temporaryVectorS[i] == 0) {
                toassignTaskToCores.push(i);
                temporaryVectorF[i] = temporaryVectorS[i] = -1; 
            }
        }
    }
}






void calculateTotalEnergyAndTime(const vector<Task>& tasks, const vector<int>& wirelessQueue, 
                            const vector<vector<int>>& schedulingQueue, double &energyTotal, int &timeTotal) {

    energyTotal = calculateTotalEnergy(tasks, wirelessQueue, schedulingQueue);

    timeTotal = calculateTotalTime(tasks);
}

 
 
 
vector<Task> rearrangeTask(int targetTask, int targetCore, double &energyTotal, int &timeTotal, vector<vector<int>> currentTask,
                    vector<vector<vector<int>>> &schedulingQueues) {

    int originalCore = tasks[targetTask].allotedCore;

    if (originalCore == targetCore) { 
    schedulingQueues.push_back(currentTask);
    return tasks;
    }

    energyTotal = 0;
    timeTotal = 0;

    auto iter = currentTask[originalCore].begin();
    removeTaskFromCore(currentTask, originalCore, targetTask);

    addTaskToCore(currentTask, targetCore, targetTask, tasks);

    vector<Task> copyTasks = tasks;
    resetTaskScheduling(copyTasks);

    copyTasks[targetTask].allotedCore = targetCore;

    vector<vector<int>> schedulingQueue;
    vector<int> temporaryVectorF, temporaryVectorS;
    stack<int> toassignTaskToCores;

    PerformTasksScheduling(toassignTaskToCores, temporaryVectorF, temporaryVectorS, tasks, currentTask, tasksNumber, coresNumber);

    processScheduling(schedulingQueue, toassignTaskToCores, copyTasks, temporaryVectorF, temporaryVectorS, currentTask, tasksNumber, coresNumber);

    schedulingQueues.push_back(schedulingQueue);

    vector<int> wirelessQueue = schedulingQueue[0];
    schedulingQueue.erase(schedulingQueue.begin());

    calculateTotalEnergyAndTime(copyTasks, wirelessQueue, schedulingQueue, energyTotal, timeTotal);

    return copyTasks;
}






void FindExitTasks(const vector<Task>& tasks, vector<int>& exitTasks) {

    for (const Task& task : tasks) {

        if (task.TasksucceedingId.empty()) {

            exitTasks.push_back(task.Taskid);
        }
    }
}









void assignTasksAndCalculateLoad(vector<Task>& tasks, const vector<int>& cloudTimes, int coresNumber) {
    int totalTimeCloudExecution = cloudTimes[0] + cloudTimes[1] + cloudTimes[2]; 

    for (Task& task : tasks) {
        int totalTimeLocalExecution = numeric_limits<int>::max();

        for (int executionTime : task.TasklocalExecutions) {
            totalTimeLocalExecution = min(totalTimeLocalExecution, executionTime);
        }

        if (totalTimeCloudExecution < totalTimeLocalExecution) {
            task.allotedCore = 0; 
            task.Taskweight = totalTimeCloudExecution; 
        } else {

            int totalLocalExecutionTime = accumulate(task.TasklocalExecutions.begin(), task.TasklocalExecutions.end(), 0);

            task.Taskweight = static_cast<double>(totalLocalExecutionTime) / coresNumber;
        }
    }
}





vector<Task> sortTasksByTaskPriority(const vector<Task>& tasks) {

    vector<Task> taskSorted(tasks);

    sort(taskSorted.begin(), taskSorted.end(), [](const Task &Vector1, const Task &Vector2) {
        return Vector1.Taskpriority > Vector2.Taskpriority;
    });

    if (!taskSorted.empty()) {
        taskSorted.pop_back(); 
    }

    return taskSorted;
}







void assignTaskToCore(Task& task, const vector<Task>& tasks, vector<int>& wirelessSendQueue, vector<vector<int>>& cores, int coresNumber) {
    int allotedCore = -1;
    int MinimumFinalTime = numeric_limits<int>::max();

    vector<Task> CopyTaskCore(coresNumber, task);
    for (int j = 0; j < coresNumber; ++j) {
        calculateCoreTimes(tasks, CopyTaskCore[j], cores[j], j);
        if (CopyTaskCore[j].TaskfinishingTimes[0] < MinimumFinalTime) {
            MinimumFinalTime = CopyTaskCore[j].TaskfinishingTimes[0];
            allotedCore = j + 1;
        }
    }

    Task CopyTaskWireless = task;
    calculateCloudTimes(tasks, CopyTaskWireless, wirelessSendQueue);
    if (CopyTaskWireless.TaskfinishingTimes[3] < MinimumFinalTime) {
        MinimumFinalTime = CopyTaskWireless.TaskfinishingTimes[3];
        allotedCore = 0; 
    }

    if (allotedCore == 0) {
        task = CopyTaskWireless;
    } else {
        task = CopyTaskCore[allotedCore - 1];
    }
    task.allotedCore = allotedCore;

    if (allotedCore == 0) {
        wirelessSendQueue.push_back(task.Taskid);
    } else {
        cores[allotedCore - 1].push_back(task.Taskid);
    }
}





void performInitialAssignment(vector<Task>& tasks, const vector<int>& cloudTimes, int coresNumber, vector<int>& wirelessSendQueue, vector<vector<int>>& cores) {

    assignTasksAndCalculateLoad(tasks, cloudTimes, coresNumber);

    for (int i = 1; i <= tasksNumber; ++i) {
        updateTaskPriority(tasks[i]);
    }

    vector<Task> taskSorted = sortTasksByTaskPriority(tasks);

    for (int i = 0; i < tasksNumber; ++i) {
        assignTaskToCore(tasks[taskSorted[i].Taskid], tasks, wirelessSendQueue, cores, coresNumber);
    }
}
