#include "matrix_io.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Доступ к A как к плотной матрице n×n */
static double Aat(const double *A, int n, int i, int j)
{
    return A[(long)i * (long)n + (long)j];
}

static void Aset(double *A, int n, int i, int j, double v)
{
    A[(long)i * (long)n + (long)j] = v;
}

int read_matrix_from_file(const char *filename, int n, double *A)
{
    FILE *fp;
    int i, j;
    double x;
    
    if (filename == NULL || A == NULL || n <= 0) return -1;
    
    fp = fopen(filename, "r");
    if (fp == NULL) return -1;
    
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (fscanf(fp, "%lf", &x) != 1) {
                fclose(fp);
                return -2;
            }
            Aset(A, n, i, j, x);
        }
    }
    
    fclose(fp);
    return 0;
}

/*
f(k,n,i,j), i,j = 1..n
k=1: n - max{i,j} + 1
k=2: 2 if i=j; -1 if |i-j|=1; 0 otherwise
k=3: 1 if i=j<n; i if j=n; j if i=n; 0 otherwise
k=4: 1/(i+j-1)
*/
double f_formula(int k, int n, int i, int j)
{
    int mx, diff;
    
    mx = (i > j) ? i : j;
    
    if (k == 1) {
        return (double)(n - mx + 1);
    }
    if (k == 2) {
        diff = i - j;
        if (diff < 0) diff = -diff;
        if (i == j) return 2.0;
        if (diff == 1) return -1.0;
        return 0.0;
    }
    if (k == 3) {
        if (i == j && i < n) return 1.0;
        if (j == n) return (double)i;
        if (i == n) return (double)j;
        return 0.0;
    }
    if (k == 4) {
        return 1.0 / (double)(i + j - 1);
    }
    
    return 0.0;
}

void fill_matrix_by_formula(int n, int k, double *A)
{
    int i0, j0, i, j;
    
    if (A == NULL || n <= 0) return;
    
    for (i0 = 0; i0 < n; ++i0) {
        for (j0 = 0; j0 < n; ++j0) {
            i = i0 + 1;  /* 1-based для формул */
            j = j0 + 1;
            Aset(A, n, i0, j0, f_formula(k, n, i, j));
        }
    }
}

void print_matrix(int l, int n, int m, const double *A)
{
    int rr, cc, i, j;
    
    if (A == NULL || n <= 0) return;
    
    rr = (l < m) ? l : m;
    cc = (n < m) ? n : m;
    
    for (i = 0; i < rr; ++i) {
        for (j = 0; j < cc; ++j) {
            printf(" %10.3e", Aat(A, n, i, j));
        }
        printf("\n");
    }
}

double trace_matrix(int n, const double *A)
{
    int i;
    double tr = 0.0;
    
    if (A == NULL || n <= 0) return 0.0;
    
    for (i = 0; i < n; ++i) {
        tr += Aat(A, n, i, i);
    }
    return tr;
}

double frobenius_norm(int n, const double *A)
{
    int i, j;
    double s = 0.0, v;
    
    if (A == NULL || n <= 0) return 0.0;
    
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            v = Aat(A, n, i, j);
            s += v * v;
        }
    }
    return sqrt(s);
}