/*
Jacobi迭代的MPI实现(C语言)

运行：
mpicc -g jacobi_parallel.c -o jacobi_parallel -Wall
mpirun -n [2/3/4/5] ./jacobi_parallel

例子：(注意Jacobi迭代的收敛条件：A选择不可约对角占优阵，或严格对角占优阵)
A = [10, 1, 1, 1;
	 1, 10, 1, 1;
	 1, 1, 10, 1;
	 1, 1, 1, 10]
b = [1, 2, 3, 4]
x_0 = [0, 0, 0, 0]
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include "mpi.h"
#define SIZE 4
#define MASTER 0
#define FROM_MASTER 0
#define FROM_WORKER 1
#define TRUE 1
#define FALSE 0
#define MAXTIMES 20
#define TOL 1e-5

typedef int BOOL;

/**
 * @brief 打印指定矩阵
*/
void Print_matix(float *buf, int nr, int nc);


int main(int argc, char *argv[]){
    int numtasks, // 进程数
		numworkers, // 工作进程数
        rank, // 进程ID
        source, // 信息来源
        dest, // 信息目的地
        tag, // 信息标签
        i, k,
		averows, // 每个工作进程平均分配的进程数
		extrarow, // 平均分配剩余的进程数
		rows, // 每个worker进程分配到需要计算的行数
		offset, // 每个worker进程分配到的函数的起始偏移量
		tmp_offset, 
		times; // 迭代次数
	BOOL first; // 标识是否是首次传递消息(首次需要传递系数矩阵A的行块和b)
	float A[SIZE][SIZE] = {
		{10.0, 1.0, 1.0, 1.0},
		{1.0, 10.0, 1.0, 1.0},
		{1.0, 1.0, 10.0, 1.0},
		{1.0, 1.0, 1.0, 10.0}
	};
	float *rp_A[SIZE]; // 存储每行首地址的指针数组(row point of A)
	for(i=0; i<SIZE; i++){
		rp_A[i] = &(A[i][0]);
	}
    float b[SIZE] = {1.0, 2.0, 3.0, 4.0};
	float x[SIZE] = {0.0, 0.0, 0.0, 0.0};
    float x_pre[SIZE] = {0.0, 0.0, 0.0, 0.0}; // 记录上一次迭代所得的x
    float *buf_A, *result; // worker进程中需要动态分配法空间
    float diff = 1.0; // 向量两次迭代x的L1误差

	MPI_Status stat;
	MPI_Datatype coltype; // 存储矩阵的一行或者一列的数据类型
	MPI_Aint extent;
	int size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	numworkers = numtasks-1;

    // 定义派生数据类型并打印相关信息
    MPI_Type_contiguous(SIZE, MPI_FLOAT, &coltype);
    MPI_Type_commit(&coltype);
	if(rank == MASTER){
		MPI_Type_extent(MPI_FLOAT, &extent);
    	printf("extent of MPI_FLOAT: %d\n", (int)extent);
		MPI_Type_extent(coltype, &extent);
    	printf("extent of coltype: %d\n", (int)extent);
		MPI_Type_size(coltype, &size); // 返回跨度减去类型中空隙后的空间的大小
		printf("size of coltype=%d\n", size);

		/* !!!巨坑:不要用sizeof作用在MPI派生数据类型上
			sizeof(coltype) ！= SIZE*sizeof(float) 
							!= extent(coltype)
			sizeof(coltype) == 4 == sizeof(float)		
		*/
		// printf("sizeof(float)=%ld\n", sizeof(float));
		// printf("sizeof(MPI_FLOAT)=%ld\n", sizeof(MPI_FLOAT));
		// printf("sizeof(coltype)=%ld\n", sizeof(coltype));
	}

	// 均匀分配数据
	averows = SIZE / numworkers;
	extrarow = SIZE % numworkers;

    times = 0;
    first = TRUE;
    // Jacobi迭代的循环
    while(times < MAXTIMES && diff > TOL){
        if(first == TRUE){
            // MASTER将数据传递给workers：
            MPI_Bcast(b, 1, coltype, MASTER, MPI_COMM_WORLD); 
            if(rank == MASTER){
                tag = FROM_MASTER;
                offset = 0;
                for(dest=1; dest<=numworkers; dest++){
                    if(dest <= extrarow)
                        rows = averows + 1;
                    else
                        rows = averows;
                    MPI_Sendrecv(MPI_BOTTOM, 0, MPI_INT, dest, 14, 
                                MPI_BOTTOM, 0, MPI_INT, dest, 14, 
                                MPI_COMM_WORLD, &stat);
                    MPI_Send(&offset, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
                    MPI_Send(&rows, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
                    MPI_Send(rp_A[offset], rows, coltype, dest, tag, MPI_COMM_WORLD);	
                    offset += rows;
                }
            }
            if(rank != MASTER){
                // worker接收MASTER传出的计算所需的基本信息
                tag = FROM_MASTER;
                MPI_Sendrecv(MPI_BOTTOM, 0, MPI_INT, MASTER, 14, 
                            MPI_BOTTOM, 0, MPI_INT, MASTER, 14, 
                            MPI_COMM_WORLD, &stat);
                MPI_Recv(&offset, 1, MPI_INT, MASTER, tag, MPI_COMM_WORLD, &stat);
                MPI_Recv(&rows, 1, MPI_INT, MASTER, tag, MPI_COMM_WORLD, &stat);
                buf_A = (float *)malloc(sizeof(float)*SIZE*rows);
                /* bug:
                    float *buf_A = (float *)malloc(sizeof(coltype)*rows);
                    报错：double free or corruption (out)
                */
                memset(buf_A, 0, sizeof(float)*SIZE*rows);
                MPI_Recv(buf_A, rows, coltype, MASTER, tag, MPI_COMM_WORLD, &stat);
                // printf("rank=%d, offset=%d\n", rank, offset);
                // printf("rank=%d, rows=%d\n", rank, rows);
                // printf("rank=%d, b=", rank);	Print_matix(b, 1, SIZE);
                // printf("rank=%d, x=", rank); 	Print_matix(x, 1, SIZE);
                // printf("rank=%d, buf_A=", rank); 	Print_matix(buf_A, rows, SIZE);
            }
            first = FALSE;
        }

        // 将每次迭代得到的x广播到每个进程
        MPI_Bcast(x, 1, coltype, MASTER, MPI_COMM_WORLD);

        // MASTER主控进程接收worker进程的计算结果，且在分配；并计算终止依据
        if(rank == MASTER){
            // 记录上一次迭代的x
            // printf("x_pre=");
            for(i=0; i<SIZE; i++){
                x_pre[i] = x[i];
                // printf("%f ", x_pre[i]);
            }
            // printf("\n");

            // MASTER接收worker回传的数据
            tag = FROM_WORKER;
            for(source=1; source<=numworkers; source++){
                MPI_Sendrecv(MPI_BOTTOM, 0, MPI_INT, source, 14, 
                            MPI_BOTTOM, 0, MPI_INT, source, 14, 
                            MPI_COMM_WORLD, &stat);
                MPI_Recv(&rows, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &stat);
                MPI_Recv(&offset, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &stat);
                MPI_Recv(&(x[offset]), rows, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &stat);
            }
            // 打印每次迭代的x
            printf("x=");
            for(i=0; i<SIZE; i++){
            	printf("%f ", x[i]);
            }
            printf("\n");

            // 此次迭代和上一次的差距(判断收敛的依据)
            diff = 0.0;
            for(i=0; i<SIZE; i++){
                if(x[i] > x_pre[i])
                    diff += x[i]-x_pre[i];
                else   
                    diff += x_pre[i]-x[i];
            }
            printf("diff=%f\n", diff);
        }

        // worker 进程计算自己的结果，并返回给MASTER进程
        if(rank != MASTER){
            // 分布计算Jacobi迭代
            tmp_offset = offset;
            result = (float *)malloc(sizeof(float)*rows);
            memset(result, 0, sizeof(float)*rows);
            for(i=0; i<rows; i++){
                for(k=0; k<SIZE; k++){
                    if(k != tmp_offset){
                        result[i] += (-1) * (buf_A[i*SIZE+k] / buf_A[i*SIZE+tmp_offset]) * x[k];
                    }
                }
                result[i] += b[tmp_offset] / buf_A[i*SIZE+tmp_offset];
                tmp_offset++;
                printf("rank=%d result[%d]=%f \n", rank, tmp_offset, result[i]);
            }

            // workers将计算结果result传回给MASTER
            tag = FROM_WORKER;
            MPI_Sendrecv(MPI_BOTTOM, 0, MPI_INT, MASTER, 14, 
                        MPI_BOTTOM, 0, MPI_INT, MASTER, 14, 
                        MPI_COMM_WORLD, &stat);
            MPI_Send(&rows, 1, MPI_INT, MASTER, tag, MPI_COMM_WORLD);
            MPI_Send(&offset, 1, MPI_INT, MASTER, tag, MPI_COMM_WORLD);
            MPI_Send(result, rows, MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD);
        }
        // 此次迭代和上一次的差距(判断收敛的依据)，广播到每个进程
        MPI_Bcast(&diff, 1, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
        times++;
    }
    
    // 打印迭代的总次数
    if(rank == MASTER){
        printf("There are %d times Jacobi iterations!\n", times);
    }

    // 释放动态分配的内存
    if(rank != MASTER){
        free(buf_A);
        free(result);
    }
    MPI_Type_free(&coltype);
    MPI_Finalize();
}


// function---------------------------------------------------------------------------
void Print_matix(float *buf, int nr, int nc){
	int i, j;
	for(i=0; i<nr; i++){
		for(j=0; j<nc; j++){
			printf("%f ", *(buf+i*nc+j));
		}
		printf("\n");
	}
}
