#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>


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

void gemm(double* A, double* B, double* C, int size) {
    for (int i = 0; i < size; i++) {
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
    // Input matrix size
    if (argc < 2) {
        printf("Input the matrix size as a parameter\n");
        return 1;
    }

    int size = atoi(argv[1]);
    if (size <= 0) {
        printf("Matrix size must be a positive integer.\n");
        return 1;
    }

    // Allocating and generating matrixies 
    double* A = (double*)malloc(size * size * sizeof(double));
    double* B = (double*)malloc(size * size * sizeof(double));
    double* C = (double*)malloc(size * size * sizeof(double));

    //srand(time(NULL));
    generate_random_matrix(A, size);
    generate_random_matrix(B, size);

    // Start execution
    double start_time = omp_get_wtime();

    gemm(A, B, C, size);

    // End execution
    double end_time = omp_get_wtime();
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

    return 0;
}
