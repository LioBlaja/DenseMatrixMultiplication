#ifndef HEADER_H
#define HEADER_H

// #define EPSILON 0.000001
#define EPSILON 0.000001  



double** serial_multiply_v1(double**,double**,int);
double** serial_multiply_v2(double**,double**,int);
double** serial_multiply_v3(double**,double**,int);
double** serial_multiply_v4(double**,double**,int);
double** serial_multiply_v5(double**,double**,int);
double** serial_multiply_v6(double**,double**,int);

double** parallel_multiply_v1(double**,double**,int,int nthreads, int chunk);
double** parallel_multiply_v2(double**,double**,int,int nthreads, int chunk);
double** parallel_multiply_v3(double**,double**,int,int nthreads, int chunk);
double** parallel_multiply_v4(double**,double**,int,int nthreads, int chunk);
double** parallel_multiply_v5(double**,double**,int,int nthreads, int chunk);
double** parallel_multiply_v6(double**,double**,int,int nthreads, int chunk);

double** serialBlockedMultiplication(double**, double**,int);
double** parallelBlockedMultiplication(double**, double**, int, int, int,int);

#endif // HEADER_H
