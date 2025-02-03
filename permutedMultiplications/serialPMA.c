#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "permutedMultiplications.h"

/*(i,j,k)
matrix multiplication - classical serial version
input: global vars a, b
output: global var c
*/
double** serial_multiply_v1(double** a,double** b,int N)
{
    double** c = (double**)malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++) {
        c[i] = (double*)malloc(N * sizeof(double));
    }

    int i, j, k;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
        {
            c[i][j] = 0;
            for (k = 0; k < N; k++)
                c[i][j] += a[i][k] * b[k][j];
        }

    return c;
}

/*
matrix multiplication serial version (ikj)
input: global vars a, b
output: global var c
*/
double** serial_multiply_v2(double** a,double** b,int N)
{
    int i, j, k;
    double aik;
    double** c = (double**)malloc(N * sizeof(double*));
    
    for (int i = 0; i < N; i++) {
        c[i] = (double*)malloc(N * sizeof(double));
    }

    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            c[i][j] = 0;

    for (i = 0; i < N; i++)
        for (k = 0; k < N; k++)
        {
            aik = a[i][k];
            for (int j = 0; j < N; j++)
            {
                c[i][j] += aik * b[k][j];
            }
        }
    
    return c;
}

//kij
double** serial_multiply_v3(double** a,double** b,int N) {

    double** c = (double**)malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++) {
        c[i] = (double*)malloc(N * sizeof(double));
        for (int j = 0; j < N; j++) {
            c[i][j] = 0.0; // Set each element to 0
        }
    }

    for (int k = 0; k < N; k++) {
        for (int i = 0; i < N; i++) {
            double kij = a[i][k];
            for (int j = 0; j < N; j++) {
                c[i][j] += kij * b[k][j];
            }
        }
    }
    return c;
}

//kji
double** serial_multiply_v4(double** a,double** b,int N) {

    double** c = (double**)malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++) {
        c[i] = (double*)malloc(N * sizeof(double));
        for (int j = 0; j < N; j++) {
            c[i][j] = 0.0; // Set each element to 0
        }
    }

    for (int k = 0; k < N; k++) {
        for (int j = 0; j < N; j++) {
            double kji = b[k][j];
            for (int i = 0; i < N; i++) {
                c[i][j] += a[i][k] * kji;
            }
        }
    }

    return c;
}

//kji
double** serial_multiply_v5(double** a,double** b,int N) {

    double** c = (double**)malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++) {
        c[i] = (double*)malloc(N * sizeof(double));
        for (int j = 0; j < N; j++) {
            c[i][j] = 0.0; // Set each element to 0
        }
    }

    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            double sum = 0.0;
            for (int k = 0; k < N; k++) {
                sum += a[i][k] * b[k][j];
            }
            c[i][j] = sum;
        }
    }

    return c;
}

//jki
double** serial_multiply_v6(double** a,double** b,int N) {

    double** c = (double**)malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++) {
        c[i] = (double*)malloc(N * sizeof(double));
        for (int j = 0; j < N; j++) {
            c[i][j] = 0.0; // Set each element to 0
        }
    }

    for (int j = 0; j < N; j++) {
        for (int k = 0; k < N; k++) {
            double jki = b[k][j];
            for (int i = 0; i < N; i++) {
                c[i][j] += a[i][k] * jki;
            }
        }
    }

    return c;
}

double** serialBlockedMultiplication(double** a, double** b, int N) {

    int i, j, k, ii, jj, kk;
    double** c = (double**)malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++) {
        c[i] = (double*)calloc(N,sizeof(double));
    }

    // for (int i = 0; i < N; i++) {
    //     for (int j = 0; j < N; j++) {
    //         c[i][j] = 0.0;  
    //     }
    // }

    int Q = 8;  

    for (ii = 0; ii < N; ii += Q) {
        for (jj = 0; jj < N; jj += Q) {

            for (kk = 0; kk < N; kk += Q) {
                for (i = ii; i < ii + Q && i < N; i++) {       
                    for (j = jj; j < jj + Q && j < N; j++) {  
                        for (k = kk; k < kk + Q && k < N; k++) { 
                            c[i][j] += a[i][k] * b[k][j];
                        }
                    }
                }
            }
        }
    }

    return c;
}