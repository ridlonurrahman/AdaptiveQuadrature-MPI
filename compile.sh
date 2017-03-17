/usr/lib64/openmpi/bin/mpicc -o aquadSolution aquadSolution.c stack.h stack.c -lm
/usr/lib64/openmpi/bin/mpirun -c 5 ./aquadSolution
