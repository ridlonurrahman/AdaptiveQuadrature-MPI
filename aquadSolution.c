#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>
#include "stack.h"

#define EPSILON 1e-3
#define F(arg)  cosh(arg)*cosh(arg)*cosh(arg)*cosh(arg)
#define A 0.0
#define B 5.0

#define SLEEPTIME 1

int *tasks_per_process;

double farmer(int);

void worker(int);

int num_workers_available(int*,int);

int main(int argc, char **argv ) {
  int i, myid, numprocs;
  double area, a, b;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  if(numprocs < 2) {
    fprintf(stderr, "ERROR: Must have at least 2 processes to run\n");
    MPI_Finalize();
    exit(1);
  }

  if (myid == 0) { // Farmer
    // init counters
    tasks_per_process = (int *) malloc(sizeof(int)*(numprocs));
    for (i=0; i<numprocs; i++) {
      tasks_per_process[i]=0;
    }
  }

  if (myid == 0) { // Farmer
    area = farmer(numprocs);
  } else { //Numworkers
    worker(myid);
  }

  if(myid == 0) {
    fprintf(stdout, "Area=%lf\n", area);
    fprintf(stdout, "\nTasks Per Process\n");
    for (i=0; i<numprocs; i++) {
      fprintf(stdout, "%d\t", i);
    }
    fprintf(stdout, "\n");
    for (i=0; i<numprocs; i++) {
      fprintf(stdout, "%d\t", tasks_per_process[i]);
    }
    fprintf(stdout, "\n");
    free(tasks_per_process);
  }
  MPI_Finalize();
  return 0;
}

int num_workers_available(int* workers, int numworkers){
  int available = 0;
  int i = 0;
  for (i=0;i<numworkers;i++){
    available = available + workers[i];
  }
  return available;
}

double farmer(int numprocs) {
  MPI_Status status;
  stack* stack;
  stack = new_stack();
  double temp[] = {A,B}; push(temp, stack);
  int numworkers = numprocs - 1;
  
  int* workers = (int*) malloc(sizeof(int)*(numworkers));
  int i,j;
  for (i=0;i<numworkers;i++){
    workers[i] = 0;
  }
  int available_workers = 0;
  double result = 0;

  while(1) {
    MPI_Recv(temp, 2, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    fprintf(stdout, "Source = %d\t", status.MPI_SOURCE);
    workers[status.MPI_SOURCE - 1] = 1;

    if (status.MPI_TAG == 0) {
      push(temp,stack);
      MPI_Recv(temp, 2, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
      fprintf(stdout, "Source = %d\t", status.MPI_SOURCE);
      push(temp,stack);
    }
    else {
      result += temp[0];
      fprintf(stdout, "Temp[0] = %f\t", temp[0]);
        }
    for (j=0;j<numworkers;j++){
      if (!is_empty(stack) && workers[j]) {
        MPI_Send(pop(stack), 2, MPI_DOUBLE, j+1, 0, MPI_COMM_WORLD);
        workers[j] = 0;
        fprintf(stdout, "Worker Just Assigned = %d\t", j+1);
        tasks_per_process[j+1]++;
            }
    }
    available_workers = num_workers_available(workers, numworkers);
    fprintf(stdout, "Down avail = %d \n", num_workers_available(workers, numworkers));
    if (is_empty(stack) && available_workers == numworkers){
      break;
    }
  }
  for (i=0;i<numworkers;i++) {
    temp[0] = 0; temp[1] = 0;
    MPI_Send(temp, 2, MPI_DOUBLE, i+1, 2, MPI_COMM_WORLD);
  }
  return result;
}

void worker(int mypid) {
  // You must complete this function
  MPI_Status status;
  double temp[] = {0,0};
  MPI_Send(temp, 2, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
  while (1) {
    MPI_Recv(temp, 2, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    if (status.MPI_TAG != 2) {
      double left = temp[0];
      double right = temp[1];
      double lrarea = (F(left) + F(right)) * (right - left) / 2;
      double mid, fmid, larea, rarea;

      mid = (left + right) / 2;
      fmid = F(mid);
      larea = (F(left) + fmid) * (mid - left) / 2;
      rarea = (fmid + F(right)) * (right - mid) / 2;
      
      usleep(2000);
      
      if ( fabs((larea + rarea) - lrarea) > EPSILON) {
        temp[0] = left; temp[1] = mid;
        MPI_Send(temp, 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); // larea
        temp[0] = mid; temp[1] = right;
        MPI_Send(temp, 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); // rarea
      }
      else {
        temp[0] = larea + rarea; temp[1] = 0;
        MPI_Send(temp, 2, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
          }
    }
    else {
      break;
    }
    }
}
