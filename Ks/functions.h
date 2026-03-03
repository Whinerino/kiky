#ifndef FUNCTIONS
#define FUNCTIONS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define SWAP_DOUBLE(a, b) do { double _tmp = (a); (a) = (b); (b) = _tmp; } while (0)
#define SWAP_INT(a, b) do { int _tmp = (a); (a) = (b); (b) = _tmp; } while (0)

double f(int s, int n, int i, int j);

/* возвращает индекс элемента в A по индексам матрицы */
int index_block_matrix(int l, int n, int m, int i, int j);

int get_matrix(int argc, int n, int m, int s, double *A, char *file_name);

void get_vector(int n, int m, double *A, double *b);

void get_block(int n, int m, int i, int j, double *A, double *C);

void put_block(int n, int m, int i, int j, double *A, double *D);

void print_matrix(int l, int n, int m, int r, double *A);

double calculation_r1(int n, int m, double *A, double *x, double *b);

double calculation_r2(int n, double *x);

double norm_matrix(int m, double *B);

int get_inverse(int m, double norm_A, int *moved_columns, double *B, double *C);

void multiplication_blocks(int p, int q, int r, double *B, double *C, double *R);

int main_block(int n, int m, int i, double norm_A, double *A, double *B, double *C, int *s);

void swap_columns_blocks(int n, int m, int j1, int j2, double *A);

void get_vector_block(int n, int m, int i, double *b, double *C);

void put_vector_block(int n, int m, int i, double *b, double *D);

void subtraction_block(int m, int l, double *C, double *D);

int solution_system(int n, int m, double *A, double *b, double *x,
                    int *s, int *moved_columns_blocks, double *C, double *D, double *R);

#endif /* FUNCTIONS */