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
  /* Count how many workers available
  Args:
    workers: List of workers (1 available, 0 not available)
    numworkers: Number of workers (numprocs-1)

  Return:
    available: number of available worker
  */
  int available = 0;
  int i = 0;
  for (i=0; i<numworkers; i++){
    available = available + workers[i];
  }
  return available;
}

double farmer(int numprocs) {
  /* Farmer function
  Args:
    numprocs: number of process

  Return:
    result: area
  */

  MPI_Status status;
  stack* stack;
  stack = new_stack();
  double temp[] = {A,B}; push(temp, stack);
  int numworkers = numprocs - 1;
  
  int* workers = (int*) malloc(sizeof(int)*(numworkers));
  int i;
  for (i=0; i<numworkers; i++){
    workers[i] = 0;
  }
  double result = 0;

  while (1) {
    MPI_Recv(temp, 2, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status); // receive from a worker
    workers[status.MPI_SOURCE - 1] = 1;

    // check if worker sends left/mid/right (MPI_TAG = 0) else larea+rarea (MPI_TAG = 1)
    if (status.MPI_TAG == 0) {
      push(temp,stack);
      MPI_Recv(temp, 2, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status); // receive again from a worker 2/2
      push(temp,stack);
    }
    else {
      result = result + temp[0]; // add the calculation to result
    }

    // distributing tasks for available/free workers
    for (i=0; i<numworkers; i++){
      if (!is_empty(stack) && workers[i]) {
        MPI_Send(pop(stack), 2, MPI_DOUBLE, i+1, 0, MPI_COMM_WORLD);
        tasks_per_process[i+1]++;
        workers[i] = 0;
      }
    }

    // calculation complete condition. when stack is empty and all workers are available (distribution tasks complete).
    if (is_empty(stack) && num_workers_available(workers, numworkers) == numworkers){
      break;
    }

  }

  // send finish tag (MPI_TAG = 2)
  for (i=0; i<numworkers; i++) {
    temp[0] = 0; temp[1] = 0;
    MPI_Send(temp, 2, MPI_DOUBLE, i+1, 2, MPI_COMM_WORLD);
  }

  return result; // return the result
}

void worker(int mypid) {
  /* Worker function
  Args:
    mypid
  Return:
    None
  */

  MPI_Status status;
  double temp[] = {0,0};
  MPI_Send(temp, 2, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

  while (1) {
    MPI_Recv(temp, 2, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status); // receive computation condition from farmer. If MPI_TAG==2 then break.

    if (status.MPI_TAG != 2) {
      // Aquad function. It's just the aquadsequential.c when written in MPI workers.
      double left = temp[0];
      double right = temp[1];
      double lrarea = (F(left) + F(right)) * (right - left)/2;
      double mid, fmid, larea, rarea;

      mid = (left + right)/2;
      fmid = F(mid);
      larea = (F(left) + fmid) * (mid - left)/2;
      rarea = (fmid + F(right)) * (right - mid)/2;
      
      usleep(2000); // delay so that the farmer is not hijacked by a worker.
      
      if ( fabs((larea + rarea) - lrarea) > EPSILON) {
        temp[0] = left; temp[1] = mid;
        MPI_Send(temp, 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        temp[0] = mid; temp[1] = right;
        MPI_Send(temp, 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
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
