#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "permutedMultiplications.h"




double** initialize_matrix(int N, int nthreads) {
    double** c2 = (double**)malloc(N * sizeof(double*));
    #pragma omp parallel num_threads(nthreads) default(none) shared(c2, N)
    {
        #pragma omp for
        for (int o = 0; o < N; o++) {
            c2[o] = (double*)calloc(N, sizeof(double)); // Allocate and initialize to zero
        }
    }
    return c2;
}

/*
matrix multiplication parallel version
input: global vars a, b
output: global var c2
*/
double** parallel_multiply_v1(double** a,double** b,int N,int nthreads, int chunk)
{
    int i, j, k;
    // double temp;
    double** c2 = initialize_matrix(N, nthreads); 

    #pragma omp parallel num_threads(nthreads), default(none), private(i, j, k), shared(a, N, b, c2, chunk)
    {
        #pragma omp for schedule(static, chunk)
            for (i = 0; i < N; i++){
                for (j = 0; j < N; j++)
                {
                    c2[i][j] = 0;
                    for (k = 0; k < N; k++)
                        c2[i][j] += a[i][k] * b[k][j];
                    // c2[i][j] = temp;
                }
            }
    }
    return c2;
}

/*
matrix multiplication parallel version (ikj)
input: global vars a, b
output: global var c2
*/
double** parallel_multiply_v2(double** a,double** b,int N,int nthreads, int chunk)
{
    int i, k;
    double aik;
    double** c2 = initialize_matrix(N, nthreads); 

    #pragma omp parallel num_threads(nthreads) default(none) private(i, k,aik) shared(a, b, c2,N, chunk)
    {
        double* temp = (double*)calloc(N, sizeof(double));
        #pragma omp for schedule(static, chunk)
            for (i = 0; i < N; i++){
                for (k = 0; k < N; k++)
                {
                    aik = a[i][k];
                    for (int j = 0; j < N; j++)
                    {
                        temp[j] += aik * b[k][j];
                    }
                }
                #pragma omp critical
                for (int o = 0; o < N; o++) {
                    c2[i][o] += temp[o];
                }
                memset(temp, 0, N * sizeof(double));
            }
        free(temp);
    }
    return c2;
}

//kij
// double** parallel_multiply_v3(double** a,double** b,int N,int nthreads, int chunk) {

//     int i, j, k;
//     double kij;
//     double** c2 = initialize_matrix(N, nthreads); 

//     #pragma omp parallel num_threads(nthreads) default(none) private(k, i, j, kij) shared(a, b, c2, N, chunk)
//     {
//         #pragma omp for schedule(static, chunk)
//         for (k = 0; k < N; k++) {
//             for (i = 0; i < N; i++) {
//                 kij = a[i][k]; 
//                 // double temp = 0;
//                 for (j = 0; j < N; j++) {
//                     #pragma omp atomic 
//                     c2[i][j] += kij * b[k][j];
//                 }

//             }
//         }
//     }
//     return c2;
// }
double** parallel_multiply_v3(double** a,double** b,int N,int nthreads, int chunk) {

    int i, j, k;
    double kij;
    double** c2 = initialize_matrix(N, nthreads); 

    #pragma omp parallel num_threads(nthreads) default(none) private(k, i, j, kij) shared(a, b, c2, N, chunk)
    {
        double* temp = (double*)calloc(N, sizeof(double));
        #pragma omp for schedule(static, chunk)
        for (k = 0; k < N; k++) {
            for (i = 0; i < N; i++) {
                kij = a[i][k]; 
                // double temp = 0;
                for (j = 0; j < N; j++) {
                    temp[j] += kij * b[k][j];
                }
            #pragma omp critical
            for (int o = 0; o < N; o++) {
                c2[i][o] += temp[o];
            }
            memset(temp, 0, N * sizeof(double));

            }
        }
        free(temp);
    }
    return c2;
}

//kji
// double** parallel_multiply_v4(double** a,double** b,int N,int nthreads, int chunk) {

//     int k, j;
//     double kji;
//     double** c2 = initialize_matrix(N, nthreads); 

//     #pragma omp parallel num_threads(nthreads) default(none) private(k, j, kji) shared(a, b, c2, N, chunk)
//     {
//         #pragma omp for schedule(static, chunk)
//             for (k = 0; k < N; k++) {
//                 for (j = 0; j < N; j++) {
//                     kji = b[k][j];
//                     for (int i = 0; i < N; i++) {
//                         #pragma omp atomic 
//                             c2[i][j] += kji * a[i][k];
//                     }
//                 }
//             }
//     }
//     return c2;
// }

double** parallel_multiply_v4(double** a,double** b,int N,int nthreads, int chunk) {

    int k, j, i;
    double kji;
    double** c2 = initialize_matrix(N, nthreads); 

    #pragma omp parallel num_threads(nthreads) default(none) private(k, j, i, kji) shared(a, b, c2, N, chunk)
    {
        double* temp = (double*)calloc(N, sizeof(double));
        #pragma omp for schedule(static, chunk)
            for ( k = 0; k < N; k++) {
                for ( j = 0; j < N; j++) {
                    kji = b[k][j];
                    for ( i = 0; i < N; i++) {
                            temp[i] += kji * a[i][k];
                    }

                    #pragma omp critical
                    for (int o = 0; o < N; o++) {
                        c2[o][j] += temp[o];
                    }
                    memset(temp, 0, N * sizeof(double));
                }
            }
        free(temp);
    }
    return c2;
}

//kji
double** parallel_multiply_v5(double** a,double** b,int N,int nthreads, int chunk) {

    int i, j;
    double** c2 = initialize_matrix(N, nthreads); 

    #pragma omp parallel num_threads(nthreads) default(none) private(i, j) shared(a, b, c2, N, chunk)
    {
        #pragma omp for schedule(static, chunk)
            for ( j = 0; j < N; j++) {
                for ( i = 0; i < N; i++) {
                    double sum = 0.0;
                    for (int k = 0; k < N; k++) {
                        sum += a[i][k] * b[k][j];
                    }
                    c2[i][j] = sum;
                }
            }
    }
    return c2;
}

//jki
double** parallel_multiply_v6(double** a,double** b,int N,int nthreads, int chunk) {

    int i, j;
    double jki;
    double** c2 = initialize_matrix(N, nthreads); 

    #pragma omp parallel num_threads(nthreads) default(none) private(i, j, jki) shared(a, b, c2, N, chunk)
    {
        double* temp = (double*)calloc(N, sizeof(double));
        #pragma omp for schedule(static, chunk)
            for ( j = 0; j < N; j++) {
                for (int k = 0; k < N; k++) {
                    jki = b[k][j];
                    for ( i = 0; i < N; i++) {
                        temp[i] += jki * a[i][k];
                    }
                }

                #pragma omp critical
                for (i = 0; i < N; i++) {
                    c2[i][j] += temp[i];
                }
                memset(temp, 0, N * sizeof(double));
            }
        free(temp);
    }
    return c2;
}

double** parallelBlockedMultiplication(double** a, double** b, int N,int nthreads, int chunk, int Q) {
    int i, j, k, ii, jj, kk;
    double** c = initialize_matrix(N, nthreads);


    //choosen because we are litterally breaking the first matrix in N / Q groups of rows, so each thread will take by one group and one group on columns
    #pragma omp parallel for collapse(2) private(i, j, k, ii, jj, kk) shared(a, b, c) schedule(static, chunk)
    //used to iterate trough first matrix rows, here we are gona consider rows taken by blocks of Q
    for (ii = 0; ii < N; ii += Q) {
    //used to iterate trough second matrix by columns, again considering an increment of Q
        for (jj = 0; jj < N; jj += Q) {
            //iterating trough results matrix
            for (kk = 0; kk < N; kk += Q) {
                //used to take first each row that were iterated before
                for (i = ii; i < ii + Q && i < N; i++) {
                    // used to iterate through columns by Q increment
                    for (j = jj; j < jj + Q && j < N; j++) {  
                        double sum = 0.0;
                        //iterating sub block of result array
                        for (k = kk; k < kk + Q && k < N; k++) { 
                            sum += a[i][k] * b[k][j];
                        }

                        c[i][j] += sum;
                    }
                }
            }
        }
    }

    return c;
}



