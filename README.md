# Matrix Multiplication with OpenMP
Tested with gcc-9.1.0 and gcc-14.2.0
Use these commands on the folder `/ODMM`
First execute `chmod +x *.sh ; find . -type f -name "*.pbs" -exec dos2unix {} \;`, this command makes script for testing executable and parses to unix PBS files.
Each PBS file and C file tries to mantain consistency with code, comments and line number over each version.

### Compilation commands:
- `gcc -o blas/blas blas/blas.c -fopenmp`
- `gcc -o serial/serial serial/serial.c -fopenmp`
- `gcc -o openmp/openmp openmp/openmp.c -fopenmp -Ofast`
- `gcc -o openmp2/openmp2 openmp2/openmp2.c -fopenmp -Ofast`
- `mpicc -o mpi/mpi mpi/mpi.c -fopenmp -lm`
- `mpicc -o mpi2/mpi2 mpi2/mpi2.c -fopenmp -lm`

### To test different threads:
- On linux: `export NUM_THREADS = 4`
- On cmd: `set NUM_THREADS=4`
- On Powershell: `$env:OMP_NUM_THREADS=4`

You can add this on the PBS file to test different threads.

### Execution commands:
- `./blas/blas 512`
- `./serial/serial 512`
- `./openmp/openmp 512`
- `./openmp2/openmp2 512`
- `mpirun -np 4 ./mpi/mpi 512`
- `mpirun -np 4 ./mpi2/mpi2 512`

The first argument is the size of the matrixes.
Use the flag `-v` as second argument to verify the correct execution of the code.
E.g. `./openmp/openmp 512 -v`

Use `test.sh` script for testing different sizes of matrixies.
Use `testThreads.sh` script for testing different sizes of matrixies and different number of threads.
Use `testMPI.sh` script for testing different sizes of matrixies and different number of processors.
E.g. `./test.sh serial/serial`, `./testThreads.sh openmp/openmp`
The output consists of the command with the corresponding argument, the duration of each of the 10 iterations, a final average, and the specified number of threads.

### Executing on the cluster:
- `qsub blas/blas.pbs`
- `qsub serial/serial.pbs`
- `qsub openmp/openmp.pbs`
- `qsub openmp2/openmp2.pbs`
- `qsub mpi/mpi.pbs`
- `qsub mpi2/mpi2.pbs`

The output is saved in the respective `.o` folder.
PBS files:
- sets cluster parameters
- loads modules
- automatically navigates to `/home/$(whoami)/ODMM`
- tests configurations with compilation end executions.

### Make your test
In order to perform specific tests on the cluster you can edit the respective PBS file, or clone it, and edit under line 24 (`# Run code`).