#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>


void generate_random_matrix(double* matrix, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            matrix[i * size + j] = (double)(rand()%10000)/100;
        }
    }
}

int verify_result(double* A, double* B, double* C, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            double expected = 0.0;
            for (int k = 0; k < size; k++) {
                expected += A[i * size + k] * B[k * size + j];
            }
            if (fabs(C[i * size + j] - expected) > 1e-5) {
                return 0;
            }
        }
    }
    return 1;
}

void gemm(double* A, double* B, double* C, int size, int rows_per_proc) {
    for (int i = 0; i < rows_per_proc; i++) {
        for (int j = 0; j < size; j++) {
            double sum = 0.0;
            for (int k = 0; k < size; k++) {
                sum += A[i * size + k] * B[k * size + j];
            }
            C[i * size + j] = sum;
        }
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // Input matrix size
    if (argc < 2) {
        printf("Input the matrix size as a parameter\n");
        MPI_Finalize();
        return 1;
    }

    int size = atoi(argv[1]);
    if (size <= 0) {
        printf("Matrix size must be a positive integer.\n");
        MPI_Finalize();
        return 1;
    }

    int rows_per_proc = size / num_procs;
    int remaining_rows = size % num_procs;

    // Allocating and generating matrixies 
    double* A = NULL;
    double* B = NULL;
    double* C = NULL;
    double* local_A = malloc(rows_per_proc * size * sizeof(double));
    double* local_C = malloc(rows_per_proc * size * sizeof(double));

    if (rank == 0) {
        A = malloc(size * size * sizeof(double));
        B = malloc(size * size * sizeof(double));
        C = malloc(size * size * sizeof(double));

        //srand(time(NULL));
        generate_random_matrix(A, size);
        generate_random_matrix(B, size);

    } else {
        B = malloc(size * size * sizeof(double));
    }

    // Start execution
    double start_time = MPI_Wtime();

    MPI_Bcast(B, size * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(A, rows_per_proc * size, MPI_DOUBLE, local_A, rows_per_proc * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    gemm(local_A, B, local_C, size, rows_per_proc);

    
    MPI_Gather(local_C, rows_per_proc * size, MPI_DOUBLE, C, rows_per_proc * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0 && remaining_rows > 0) {
        for (int i = size - remaining_rows; i < size; i++) {
            for (int j = 0; j < size; j++) {
                double sum = 0.0;
                for (int k = 0; k < size; k++) {
                    sum += A[i * size + k] * B[k * size + j];
                }
                C[i * size + j] = sum;
            }
        }
    }

    // End execution
    double end_time = MPI_Wtime();
    if (rank == 0) {
        printf("Execution time: %f seconds\n", end_time - start_time);

        // Verify result
        if(argc >= 3 && strcmp(argv[2], "-v") == 0){
            if (verify_result(A, B, C, size)) {
                printf("Matrix multiplication is correct.\n");
            } else {
                printf("Matrix multiplication is incorrect.\n");
            }
        }

        free(A);
        free(B);
        free(C);
    }
    
    free(local_A);
    free(local_C);
    MPI_Finalize();
    
    return 0;
}
