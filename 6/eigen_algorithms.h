#ifndef EIGEN_ALGORITHMS_H
#define EIGEN_ALGORITHMS_H

/* Приведение симметричной матрицы к трёхдиагональному виду методом вращений */
int tridiagonalize_symmetric(int n, double *A, double *diag, double *off);

/* Число собственных значений < x (Sturm count через рекуррентные формулы) */
int sturm_count(int n, const double *diag, const double *off, double x);

/* Границы спектра трёхдиагональной матрицы */
void tridiag_spectrum_bounds(int n, const double *diag, const double *off,
                             double *lo, double *hi);

/* Нахождение k-го собственного значения методом бисекции */
int bisection_kth_eigenvalue(int n, const double *diag, const double *off,
                             double eps, int k, double *eigenvalue, int *iterations);

/* Нахождение всех собственных значений (для невязок) */
int bisection_all_eigenvalues(int n, const double *diag, const double *off,
                              double eps, double *eig, int *total_iterations);

#endif /* EIGEN_ALGORITHMS_H */