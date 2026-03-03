#include "matrix_io.h"
#include "eigen_algorithms.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

/* Получение времени в секундах */
static double wall_time_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* Парсинг целого числа */
static int parse_int(const char *s, int *out)
{
    char *end = NULL;
    long v = strtol(s, &end, 10);
    if (!end || *end != '\0') return -1;
    if (v < -2147483647L || v > 2147483647L) return -1;
    *out = (int)v;
    return 0;
}

/* Парсинг вещественного числа */
static int parse_double(const char *s, double *out)
{
    char *end = NULL;
    double v = strtod(s, &end);
    if (!end || *end != '\0') return -1;
    *out = v;
    return 0;
}

int main(int argc, char **argv)
{
    int n, m, k_formula;
    double eps;
    double *A;
    double *diag;
    double *off;
    double *eig;
    int rc;
    long double sum;
    long double sum2;
    double res1;
    double res2;
    double trA;
    double normA;
    double t_start;
    double t1;  /* Время приведения к трёхдиагональному виду */
    double t2;  /* Время нахождения собственных значений */
    int its;    /* Общее число итераций */
    char *filename;
    int i;
    
    /* Проверка аргументов командной строки */
    if (argc != 5 && argc != 6) {
        fprintf(stderr, "Usage: %s n m eps k [filename if k==0]\n", argv[0]);
        return 1;
    }
    
    if (parse_int(argv[1], &n) != 0 || n <= 0) {
        fprintf(stderr, "Error: invalid n\n");
        return 1;
    }
    if (parse_int(argv[2], &m) != 0 || m <= 0) {
        fprintf(stderr, "Error: invalid m\n");
        return 1;
    }
    if (parse_double(argv[3], &eps) != 0 || eps <= 0.0) {
        fprintf(stderr, "Error: invalid eps\n");
        return 1;
    }
    if (parse_int(argv[4], &k_formula) != 0 || k_formula < 0 || k_formula > 4) {
        fprintf(stderr, "Error: invalid k (0..4)\n");
        return 1;
    }
    
    if (k_formula == 0 && argc != 6) {
        fprintf(stderr, "Error: filename is required when k==0\n");
        return 1;
    }
    if (k_formula != 0 && argc != 5) {
        fprintf(stderr, "Error: filename must be absent when k!=0\n");
        return 1;
    }
    
    /* Память: A = n², плюс O(n) */
    A = (double*)malloc((size_t)n * (size_t)n * sizeof(double));
    diag = (double*)malloc((size_t)n * sizeof(double));
    off = (double*)malloc((size_t)(n > 1 ? n - 1 : 1) * sizeof(double));
    eig = (double*)malloc((size_t)n * sizeof(double));
    
    if (!A || !diag || !off || !eig) {
        fprintf(stderr, "Error: not enough memory\n");
        free(A); free(diag); free(off); free(eig);
        return 1;
    }
    
    /* Инициализация матрицы */
    rc = 0;
    if (k_formula == 0) {
        rc = read_matrix_from_file(argv[5], n, A);
        if (rc != 0) {
            fprintf(stderr, "Error: cannot read matrix from file (code %d)\n", rc);
            free(A); free(diag); free(off); free(eig);
            return 2;
        }
    } else {
        fill_matrix_by_formula(n, k_formula, A);
    }
    
    /* Печать исходной матрицы */
    print_matrix(n, n, m, A);
    
    /* Инварианты исходной матрицы */
    trA = trace_matrix(n, A);
    normA = frobenius_norm(n, A);
    if (fabs(normA) < eps) normA = 1.0;
    
    /* ===== T1: Приведение к трехдиагональному виду ===== */
    t_start = wall_time_sec();
    rc = tridiagonalize_symmetric(n, A, diag, off);
    t1 = wall_time_sec() - t_start;
    
    if (rc != 0) {
        fprintf(stderr, "Error: tridiagonalization failed\n");
        free(A); free(diag); free(off); free(eig);
        return 3;
    }
    
    /* ===== T2: Нахождение всех собственных значений бисекцией + Sturm ===== */
    t_start = wall_time_sec();
    its = bisection_all_eigenvalues(n, diag, off, eps, eig, &its);
    t2 = wall_time_sec() - t_start;
    
    /* Печать собственных значений как матрица 1×n */
    print_matrix(1, n, m, eig);
    
    /* ===== Вычисление невязок ===== */
    sum = 0.0L;
    sum2 = 0.0L;
    for (i = 0; i < n; ++i) {
        long double v = (long double)eig[i];
        sum += v;
        sum2 += v * v;
    }
    
    /* Residual1 = |trace(A) - Σλi| / ||A||_F */
    res1 = fabs((double)(trA - (double)sum)) / normA;
    
    /* Residual2 = | ||A||_F - sqrt(Σλi²) | / ||A||_F */
    res2 = fabs(normA - sqrt((double)sum2)) / normA;
    
    /* ===== ВЫВОД ОТЧЁТА (ТОЧНЫЙ ФОРМАТ) ===== */
    printf("%s: Residual1=%e Residual2=%e Iterations=%d Iterations1=%d Elapsed1=%.2f Elapsed2=%.2f\n",
           argv[0], res1, res2, its, (n > 0) ? (its / n) : 0, t1, t2);
    
    free(A);
    free(diag);
    free(off);
    free(eig);
    
    return 0;
}