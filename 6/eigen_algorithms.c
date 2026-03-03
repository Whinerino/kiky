#include "eigen_algorithms.h"
#include <math.h>
#include <stdlib.h>


/* Надёжный sqrt(a^2+b^2) без hypot (C90) */
static double safe_hypot(double a, double b)
{
    double aa = fabs(a);
    double bb = fabs(b);
    if (aa < bb) {
        double tmp = aa; aa = bb; bb = tmp;
    }
    if (aa < 1e-300) return 0.0;
    double r = bb / aa;
    return aa * sqrt(1.0 + r * r);
}

/* Приведение симметричной матрицы к трехдиагональному виду вращениями Гивенса */
int tridiagonalize_symmetric(int n, double *A, double *diag, double *off)
{
    int i, j, k;
    double x, y, r, c, s;
    double tmp1, tmp2;
    
    if (n <= 0 || A == NULL || diag == NULL) return -1;
    if (n > 1 && off == NULL) return -1;
    
    for (j = 0; j < n - 2; j++) {
        for (i = n - 1; i > j + 1; i--) {
            x = A[(i-1) * n + j];
            y = A[i * n + j];
            
            if (fabs(y) < 1e-15) {
                A[i * n + j] = 0.0;
                A[j * n + i] = 0.0;
                continue;
            }
            
            r = safe_hypot(x, y);
            if (r < 1e-15) continue;
            
            c = x / r;
            s = -y / r;
            
            /* Левое умножение: строки i-1 и i */
            for (k = 0; k < n; k++) {
                tmp1 = A[(i-1) * n + k];
                tmp2 = A[i * n + k];
                A[(i-1) * n + k] = c * tmp1 - s * tmp2;
                A[i * n + k] = s * tmp1 + c * tmp2;
            }
            
            /* Правое умножение: столбцы i-1 и i */
            for (k = 0; k < n; k++) {
                tmp1 = A[k * n + (i-1)];
                tmp2 = A[k * n + i];
                A[k * n + (i-1)] = c * tmp1 - s * tmp2;
                A[k * n + i] = s * tmp1 + c * tmp2;
            }
            
            /* Зануление */
            A[i * n + j] = 0.0;
            A[j * n + i] = 0.0;
        }
    }
    
    /* Извлекаем диагональ и поддиагональ */
    for (i = 0; i < n; i++) {
        diag[i] = A[i * n + i];
        if (i < n - 1) {
            off[i] = fabs(A[(i+1) * n + i]);
        }
    }
    
    return 0;
}

/* Число собственных < x (Sturm count через рекуррентные формулы) */
int sturm_count(int n, const double *diag, const double *off, double x)
{
    int i, count;
    double p0, p1, p2;
    double e;
    const double TINY = 1e-18;
    const double BIG = 1e100;
    
    if (n <= 0 || diag == NULL) return 0;
    
    count = 0;
    
    /* Первый минор */
    p1 = diag[0] - x;
    if (fabs(p1) < TINY) {
        p1 = (p1 < 0.0) ? -TINY : TINY;
    }
    if (p1 < 0.0) count++;
    
    /* Рекуррентная формула для остальных миноров */
    p0 = 1.0;
    for (i = 1; i < n; i++) {
        e = (off != NULL && i > 0) ? off[i-1] : 0.0;
        
        /* p2 = (d_i - x) * p1 - e^2 * p0 */
        p2 = (diag[i] - x) * p1 - e * e * p0;
        
        /* Защита от переполнения */
        if (fabs(p2) > BIG) {
            p0 /= BIG;
            p1 /= BIG;
            p2 /= BIG;
        }
        
        /* Защита от деления на 0 */
        if (fabs(p2) < TINY) {
            p2 = (p2 < 0.0) ? -TINY : TINY;
        }
        
        /* Перемена знака */
        if ((p1 < 0.0 && p2 >= 0.0) || (p1 >= 0.0 && p2 < 0.0)) {
            count++;
        }
        
        p0 = p1;
        p1 = p2;
    }
    
    return count;
}

/* Границы спектра (теорема Гершгорина) */
void tridiag_spectrum_bounds(int n, const double *diag, const double *off,
                             double *lo, double *hi)
{
    int i;
    double mn, mx, r, center, a, b;
    
    if (n <= 0 || diag == NULL || lo == NULL || hi == NULL) {
        if (lo) *lo = 0.0;
        if (hi) *hi = 0.0;
        return;
    }
    
    mn = diag[0];
    mx = diag[0];
    
    for (i = 0; i < n; i++) {
        r = 0.0;
        if (i > 0 && off != NULL) r += fabs(off[i-1]);
        if (i < n - 1 && off != NULL) r += fabs(off[i]);
        
        center = diag[i];
        a = center - r;
        b = center + r;
        
        if (a < mn) mn = a;
        if (b > mx) mx = b;
    }
    
    double pad = 1e-12 * (fabs(mn) + fabs(mx) + 1.0);
    *lo = mn - pad;
    *hi = mx + pad;
}

/* Нахождение k-го собственного значения методом бисекции */
int bisection_kth_eigenvalue(int n, const double *diag, const double *off,
                             double eps, int k, double *eigenvalue, int *iterations)
{
    double left, right, mid;
    int count, target, it;
    const int ITMAX = 200;
    double width, scale;
    
    if (n <= 0 || diag == NULL || eigenvalue == NULL || iterations == NULL) {
        if (iterations) *iterations = 0;
        return -1;
    }
    
    tridiag_spectrum_bounds(n, diag, off, &left, &right);
    
    /* k-е по величине = (n-k+1)-е по возрастанию */
    target = n - k + 1;
    if (target < 1) target = 1;
    if (target > n) target = n;
    
    *iterations = 0;
    
    for (it = 0; it < ITMAX; it++) {
        mid = 0.5 * (left + right);
        count = sturm_count(n, diag, off, mid);
        (*iterations)++;
        
        if (count < target) {
            left = mid;
        } else {
            right = mid;
        }
        
        width = fabs(right - left);
        scale = fabs(mid) + 1.0;
        if (width <= eps * scale) break;
    }
    
    *eigenvalue = 0.5 * (left + right);
    return 0;
}

/* Нахождение всех собственных значений (для невязок) */
int bisection_all_eigenvalues(int n, const double *diag, const double *off,
                              double eps, double *eig, int *total_iterations)
{
    double left, right, mid, a, b;
    int k, it, count, its;
    const int ITMAX = 200;
    double width, scale;
    
    if (n <= 0 || diag == NULL || eig == NULL) return -1;
    
    tridiag_spectrum_bounds(n, diag, off, &left, &right);
    its = 0;
    
    for (k = 0; k < n; k++) {
        a = left;
        b = right;
        
        for (it = 0; it < ITMAX; it++) {
            mid = 0.5 * (a + b);
            count = sturm_count(n, diag, off, mid);
            its++;
            
            if (count <= k) {
                a = mid;
            } else {
                b = mid;
            }
            
            width = fabs(b - a);
            scale = fabs(mid) + 1.0;
            if (width <= eps * scale) break;
        }
        
        eig[k] = 0.5 * (a + b);
    }
    
    if (total_iterations) *total_iterations = its;
    return 0;
}