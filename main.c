#include "./permutedMultiplications/permutedMultiplications.h"

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>


// int N = 500; /* number of rows and column of matrices */
int nthreads = 4;
int chunk = 10;
double start, end, time_serial, time_parallel;

void Generate_matrix(char *prompt, double** mat, int N);
int Equal_matrixes(double** mat1, double** mat2, int N);
void Print_matrix(char *prompt, double** mat, int execNumber, int N);
double** allocate_matrix(int N);
void executeOneStepMultiplications(int execNumber,double** a,double** b,double** c,double** c2, int N);

void free_matrices(double** a, double** b, double** c, double** c2, int N) {
    for (int i = 0; i < N; i++) {
        free(a[i]);
        free(b[i]);
        free(c[i]);
        free(c2[i]);
    }
    free(a);
    free(b);
    free(c);
    free(c2);
}

int bestSerialAlgorithm = 1; // 1,2,3,4,5,6 for ijk,ikj,kij,jik,jki
double bestSerialTime = 10000000;
double bestParallelTime = 10000000;
double bestParallelBlockedTime = 10000000;
int bestParallelAlgorithm = 1; // 1,2,3,4,5,6 for ijk,ikj,kij,jik,jki

const char* algorithmNames[] = {
    "Unknown",  
    "ijk",      // Index 1
    "ikj",      // Index 2
    "kij",      // Index 3
    "jik",      // Index 4
    "jki"       // Index 5
};

int getNextMultipleDiv(int startPoint, int number)
{
    if (startPoint <= 0 || number <= 0)
    {
        return 100;
    }

    for (int i = startPoint * 2; i <= number; i += startPoint)
    {
        if (number % i == 0)
        {
            return i; 
        }
    }

    return 100;// return 100 to end the iteration by failing the iteration verification inside for
}

int getNextDiv(int startPoint, int number)
{
    if (startPoint <= 0 || number <= 0)
        return 100; 

    for (int i = startPoint; i <= number; i++)
    {
        if (number % i == 0)
            return i; // Found the next divisor
    }
    return 100; // return 100 to end the iteration by failing the iteration verification inside for
}


int main(int argc, char *argv[])
{
    //by running root@TMLS488W:/mnt/c/Users/uig60438# getconf LEVEL1_DCACHE_LINESIZE
    //I got 64 bytes cache block size

    //
    //int L1_size = 64;//because we do not have a function that will work on windows and linux in the same way, for now it is hardcoded, 
    //there is a need of an if based implemenation that will check the operational system and after that will return the L1 size
    //double size is 8 bytes, so in a cache line we can fit 8 values, we may choose the Q value by getting the 8 numbers multiples:
    // the best will be to find the multiples of 8 that will be divisors of N to use in the best way the cache line  
    // also consider that chunk size may also impact the speedup
    // to best fit the chunk should be also a multiple of 8
    for(int N = 500; N <= 3000; N += 500){

        double** a = allocate_matrix(N);
        double** b = allocate_matrix(N);
        double** c = allocate_matrix(N);
        double** c2 = allocate_matrix(N);

        printf("--->Matrix size N=%d \n", N);

        Generate_matrix("Generating matrix a ...", a, N);

        Generate_matrix("Generating matrix b ...", b, N);


    #ifdef DEBUG
        Print_matrix("Print matrix A: ...", a, -1, N);
        Print_matrix("Print matrix B: ...", b, -1, N);
    #endif

        for(int i = 1;i < 7;i++){
            executeOneStepMultiplications(i,a,b,c,c2, N);
        }

        printf("BEST SERIAL TIME %lf AND ALGORITHM: %s\n", bestSerialTime, algorithmNames[bestSerialAlgorithm]);    
        printf("BEST PARALLEL TIME %lf AND ALGORITHM: %s\n", bestParallelTime, algorithmNames[bestParallelAlgorithm]);

        printf("Starting computations for blocked multiplication with N size: %d\n",N);

        start = omp_get_wtime();
        c = serialBlockedMultiplication(a,b,N);
    #ifdef DEBUG
        Print_matrix("Serial blocked multiplication result: ...", c, -1, N);
    #endif
        end = omp_get_wtime();
        time_serial = end - start;
        printf("Serial time of serial blocked multiplication %lf seconds \n", time_serial);

        double bestSpeedup = 0;
        int bestChunkSize = 8;
        int bestBlockSize = 8;
        int bestThreadNumber = 1;
        // printf("HSH%d\n",getNextDiv(8*2,N));
        // printf("HGDSHSA%d\n",getNextDiv(1 + 1,(N / 8) * (N / 8)));
        //block size must be a multiple of 8 and a divisor of N

        for (int i = 8; i <= 32; i = getNextMultipleDiv(i, N)) // used to change block size
        {
            for (int j = 1; j <= 32; j = getNextDiv(j + 1,(N / i) * (N / i)))// chunk size
            // must be a div of iterations number inside multiplication
            // iter number: (N / block size) ^ 2
            {
                for (int k = 2; k <= 6; k += 2) //thread number
                {
                    // printf("The program started to compute with %d block size\n",i);
                    // printf("%d chunk size\n",j);
                    // printf("%d thread number\n",k);
                    // int iterationOuterFor = N / i;
                    // printf("Start working parallel with %d threads ... \n", nthreads);
                    start = omp_get_wtime();

                    // printf("CHUNK:%d\n",j);
                    c2 = parallelBlockedMultiplication(a,b,N,k, j,i);

                    end = omp_get_wtime();
                    time_parallel = (end - start);

                #ifdef DEBUG
                    Print_matrix("Parallel result V", c2, -1, N);
                #endif

                    // printf("Parallel time %lf seconds \n", time_parallel);
                    double localSpeedup = time_serial / time_parallel;
                    // printf("Speedup Block Size(%d)= %2.2lf\n",i, localSpeedup);
                    if (!Equal_matrixes(c, c2, N))
                        printf("Attention! Serial and Parallel Result not the same ! \n");
                    if(localSpeedup > bestSpeedup){
                        bestSpeedup = localSpeedup;
                        bestChunkSize = j;
                        bestBlockSize = i;
                        bestThreadNumber = k;
                        bestParallelBlockedTime = time_parallel;
                    }
                }
            }
        }

        printf("BLOCKED PARALLEL DATA:\n");

        printf("BEST THREAD NUMBER:%d || BEST CHUNK:%d || BLOCK SIZE: %d || PARALLEL TIME: %lf || SPEEDUP: %lf\n",bestThreadNumber,bestChunkSize,bestBlockSize,bestParallelBlockedTime,bestSpeedup);

        free_matrices(a, b, c, c2, N);
    }
    
    return 0;
}

void executeOneStepMultiplications(int execNumber,double** a,double** b,double** c,double** c2, int N){
printf("Start working serial V%d ... \n",execNumber);
    start = omp_get_wtime();

    switch (execNumber) {
        case 1:
            c = serial_multiply_v1(a, b, N);
            break;
        case 2:
            c = serial_multiply_v2(a, b, N);
            break;
        case 3:
            c = serial_multiply_v3(a, b, N);
            break;
        case 4:
            c = serial_multiply_v4(a, b, N);
            break;
        case 5:
            c = serial_multiply_v5(a, b, N);
            break;
        case 6:
            c = serial_multiply_v6(a, b, N);
            break;
        default:
            printf("Invalid execution version: %d\n", execNumber);
            return;
    }

    end = omp_get_wtime();
    time_serial = end - start;

    if(time_serial < bestSerialTime){
        bestSerialTime = time_serial;
        bestSerialAlgorithm = execNumber;
    }

#ifdef DEBUG
    Print_matrix("Serial result ...", c, execNumber, N);
#endif
    printf("Serial time V%d %lf seconds \n", execNumber, time_serial);

    printf("Start working parallel V%d with %d threads ... \n", execNumber, nthreads);
    start = omp_get_wtime();

    switch (execNumber) {
        case 1:
            c2 = parallel_multiply_v1(a,b,N,nthreads, chunk);
            break;
        case 2:
            c2 = parallel_multiply_v2(a,b,N,nthreads, chunk);
            break;
        case 3:
            c2 = parallel_multiply_v3(a,b,N,nthreads, chunk);
            break;
        case 4:
            c2 = parallel_multiply_v4(a,b,N,nthreads, chunk);
            break;
        case 5:
            c2 = parallel_multiply_v5(a,b,N,nthreads, chunk);
            break;
        case 6:
            c2 = parallel_multiply_v6(a,b,N,nthreads, chunk);
            break;
        default:
            printf("Invalid execution version: %d\n", execNumber);
            return;
    }

    end = omp_get_wtime();
    time_parallel = (end - start);

    if(time_parallel < bestParallelTime){
        bestParallelTime = time_parallel;
        bestParallelAlgorithm = execNumber;
    }

#ifdef DEBUG
    
    Print_matrix("Parallel result V", c2, execNumber, N);
#endif

    printf("Parallel time V%d %lf seconds \n", execNumber, time_parallel);
    printf("Speedup = %2.2lf\n", time_serial / time_parallel);
    if (!Equal_matrixes(c, c2, N))
        printf("Attention! Serial and Parallel V%d Result not the same ! \n",execNumber);


    printf("\n*********** V%d\n",execNumber);
}

void Generate_matrix(char *prompt, double** mat, int N)
{
    int i, j;

    srand(time(NULL));
    printf("%s\n", prompt);
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            mat[i][j] = rand() % 3000;
}

void Print_matrix(char *title, double** mat, int execNumber, int N)
{
    int i, j;

    char buffer[100];  
    if (execNumber != -1) {
        snprintf(buffer, sizeof(buffer), "%s%d... ", title, execNumber);
        printf("%s\n", buffer); 
    } else {
        printf("%s\n", title);  
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
            printf("%4.1f ", mat[i][j]);
        printf("\n");
    }
}

int Equal_matrixes(double** mat1, double** mat2, int N)
{
    int i, j;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
            if (fabs(mat1[i][j] - mat2[i][j]) > EPSILON)
            {
                printf("i:%d||j:%d\n",i,j);
                return 0;
            }
    }
    return 1;
}

double** allocate_matrix(int N) {
    double** matrix = malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++) {
        matrix[i] = malloc(N * sizeof(double));
    }
    return matrix;
}