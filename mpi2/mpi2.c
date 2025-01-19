#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>


void generate_random_matrix(double* matrix, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            matrix[i * size + j] = (double)(rand() % 10000) / 100.0;
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
            for (int k = 0; k < size; k++) {
                C[i * size + j] += A[i * size + k] * B[k * size + j];
            }
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
        if (rank == 0) printf("Input the matrix size as a parameter\n");
        MPI_Finalize();
        return 1;
    }

    int size = atoi(argv[1]);
    int sqrt_p = (int)sqrt(num_procs);
    if (size <= 0 || sqrt(num_procs) != sqrt_p || size % sqrt_p != 0) {
        if (rank == 0) printf("Matrix size must be divisible by sqrt of number of processes.\n");
        MPI_Finalize();
        return 1;
    }

    // Allocating and generating matrixies 
    double* A = NULL;
    double* B = NULL;
    double* C = NULL;
    int block_size;
    int bCastData[4];
    int procDim;
    if (rank == 0) {
        A = malloc(size * size * sizeof(double));
        B = malloc(size * size * sizeof(double));
        C = malloc(size * size * sizeof(double));

        //srand(time(NULL));
        generate_random_matrix(A, size);
        generate_random_matrix(B, size);

        procDim = sqrt_p;
        block_size = size / sqrt_p;
        bCastData[0] = procDim;
        bCastData[1] = block_size;
        bCastData[2] = size;
        bCastData[3] = size;
    }

    // Start execution
    double start_time = MPI_Wtime();
    
    MPI_Bcast(&bCastData, 4, MPI_INT, 0, MPI_COMM_WORLD);
    procDim = bCastData[0];
    block_size = bCastData[1];
    size = bCastData[2];
    size = bCastData[3];

    int dim[2] = {procDim, procDim};
    int period[2] = {1, 1};
    int reorder = 1;
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &cart_comm);

    double* local_A = malloc(block_size * block_size * sizeof(double));
    double* local_B = malloc(block_size * block_size * sizeof(double));

    int global_size[2] = { size, size };
    int local_size[2] = { block_size, block_size };
    int starts[2] = { 0,0 };
    MPI_Datatype type, subarrtype;
    MPI_Type_create_subarray(2, global_size, local_size, starts, MPI_ORDER_C, MPI_DOUBLE, &type);
    MPI_Type_create_resized(type, 0, block_size * sizeof(double), &subarrtype);
    MPI_Type_commit(&subarrtype);

    int* sendCounts = (int*)malloc(sizeof(int) * num_procs);
    int* displacements = (int*)malloc(sizeof(int) * num_procs);

    if (rank == 0) {
        for (int i = 0; i < num_procs; i++) {
            sendCounts[i] = 1;
        }
        int disp = 0;
        for (int i = 0; i < procDim; i++) {
            for (int j = 0; j < procDim; j++) {
                displacements[i * procDim + j] = disp;
                disp += 1;
            }
            disp += (block_size - 1) * procDim;
        }
    }

    MPI_Scatterv(A, sendCounts, displacements, subarrtype, local_A, block_size * block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(B, sendCounts, displacements, subarrtype, local_B, block_size * block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double* local_C = calloc(block_size, block_size * sizeof(double));
    int left, right, up, down, coords[2];

    MPI_Cart_coords(cart_comm, rank, 2, coords);
    MPI_Cart_shift(cart_comm, 1, coords[0], &left, &right);
    MPI_Sendrecv_replace(local_A, block_size * block_size, MPI_DOUBLE, left, 1, right, 1, cart_comm, MPI_STATUS_IGNORE);
    MPI_Cart_shift(cart_comm, 0, coords[1], &up, &down);
    MPI_Sendrecv_replace(local_B, block_size * block_size, MPI_DOUBLE, up, 1, down, 1, cart_comm, MPI_STATUS_IGNORE);

    for (int k = 0; k < procDim; k++) {
        gemm(local_A, local_B, local_C, block_size);

        MPI_Cart_shift(cart_comm, 1, 1, &left, &right);
        MPI_Cart_shift(cart_comm, 0, 1, &up, &down);
        MPI_Sendrecv_replace(local_A, block_size * block_size, MPI_DOUBLE, left, 1, right, 1, cart_comm, MPI_STATUS_IGNORE);
        MPI_Sendrecv_replace(local_B, block_size * block_size, MPI_DOUBLE, up, 1, down, 1, cart_comm, MPI_STATUS_IGNORE);
    }
    
    MPI_Gatherv(local_C, block_size * block_size, MPI_DOUBLE, C, sendCounts, displacements, subarrtype, 0, MPI_COMM_WORLD);

    free(local_C);

    // End execution
    double end_time = MPI_Wtime();
    if (rank == 0) {
        printf("Execution time: %f seconds\n", end_time - start_time);

        // Verify result
        if (argc >= 3 && strcmp(argv[2], "-v") == 0) {
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

    free(sendCounts);
    free(displacements);
    MPI_Type_free(&subarrtype);
    MPI_Finalize();

    return 0;
}
