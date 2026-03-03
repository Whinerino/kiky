#ifndef MATRIX_IO_H
#define MATRIX_IO_H

/* Чтение матрицы из файла */
int read_matrix_from_file(const char *filename, int n, double *A);

/* Формула для элемента матрицы */
double f_formula(int k, int n, int i, int j);

/* Заполнение матрицы по формуле */
void fill_matrix_by_formula(int n, int k, double *A);

/* Вывод матрицы на экран */
void print_matrix(int l, int n, int m, const double *A);

/* След матрицы */
double trace_matrix(int n, const double *A);

/* Норма Фробениуса матрицы */
double frobenius_norm(int n, const double *A);

#endif /* MATRIX_IO_H */