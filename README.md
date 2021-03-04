# jacobi_mpi_parallel

## jacobi_series.c

jacobi_series.c 是简单的jacobi迭代的C语言串行实现，用来给MPI并行实现提供思路

+ compile and run:

`gcc jacobi_series.c -o jacobi_series -g -Wall`  

`./jacobi_series`

## jacobi_parallel.c 

jacobi_parallel.c 是简单的jacobi迭代的MPI并行C语言实现

+ preparation:

install MPICH

+ compile and run:

`mpicc jacobi_parallel.c -o jacobi_parallel -Wall`  

`mpirun -n [2/3/4/5] ./jacobi_parallel`



