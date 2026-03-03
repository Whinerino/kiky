#include "functions.h"

double f(int s, int n, int i, int j)
{
    if (s == 1)
    {
        return n - MAX(i, j) + 1;
    }
    else if (s == 2)
    {
        return MAX(i, j);
    }
    else if (s == 3)
    {
        return abs(i - j);
    }
    return 1.0 / (i + j - 1); /* s == 4 */
}

/* возвращает индекс элемента в A по индексам матрицы */
int index_block_matrix(int l, int n, int m, int i, int j)
{
    int kl = l / m;
    int kn = n / m;
    int pl = l % m;
    int pn = n % m;
    int h_block = (i < kl * m) ? m : pl;
    int v_block = (j < kn * m) ? m : pn;
    int i_block = i / m;
    int j_block = j / m;
    int i_in_block = i % m;
    int j_in_block = j % m;
    
    int ij = i_block * n * m + j_block * m * h_block +
             i_in_block * v_block + j_in_block;
    return ij;
}

int get_matrix(int argc, int n, int m, int s, double *A, char *file_name)
{
    int i, j, ij, c;
    
    if (argc == 6 && s == 0)
    {
        FILE *file = fopen(file_name, "r");
        if (!file)
        {
            return 1;
        }
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                ij = index_block_matrix(n, n, m, i, j);
                if (fscanf(file, "%lf", &A[ij]) != 1)
                {
                    fclose(file);
                    return 2;
                }
            }
            while ((c = fgetc(file)) != '\n' && c != EOF)
            {
                if (c != ' ' && c != '\t' && c != '\r')
                {
                    fclose(file);
                    return 3;
                }
            }
        }
        while ((c = fgetc(file)) != EOF)
        {
            if (c != '\n' && c != ' ' && c != '\t' && c != '\r')
            {
                fclose(file);
                return 3;
            }
        }
        fclose(file);
        return 0;
    }
    else if (argc == 5 && 1 <= s && s <= 4)
    {
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                ij = index_block_matrix(n, n, m, i, j);
                A[ij] = f(s, n, i + 1, j + 1);
            }
        }
        return 0;
    }
    return 4;
}

void get_vector(int n, int m, double *A, double *b)
{
    int i, ij, t;
    for (i = 0; i < n; i++)
    {
        b[i] = 0;
        for (t = 0; t <= (n - 1) / 2; t++)
        {
            ij = index_block_matrix(n, n, m, i, 2 * t);
            b[i] += A[ij];
        }
    }
}

void get_block(int n, int m, int i, int j, double *A, double *C)
{
    int s, t;
    int k = n / m;
    int l = n % m;
    int h = (i < k ? m : l);
    int v = (j < k ? m : l);
    double *Aij = A + i * n * m + j * m * h;
    for (s = 0; s < h; s++)
    {
        for (t = 0; t < v; t++)
        {
            C[s * v + t] = Aij[s * v + t];
        }
    }
}

void put_block(int n, int m, int i, int j, double *A, double *D)
{
    int s, t;
    int k = n / m;
    int l = n % m;
    int h = (i < k ? m : l);
    int v = (j < k ? m : l);
    double *Aij = A + i * n * m + j * m * h;
    for (s = 0; s < h; s++)
    {
        for (t = 0; t < v; t++)
        {
            Aij[s * v + t] = D[s * v + t];
        }
    }
}

void print_matrix(int l, int n, int m, int r, double *A)
{
    int i, j, ij;
    for (i = 0; i < MIN(r, l); i++)
    {
        for (j = 0; j < MIN(r, n); j++)
        {
            ij = index_block_matrix(l, n, m, i, j);
            printf(" %10.3e", A[ij]);
        }
        printf("\n");
    }
}

double calculation_r1(int n, int m, double *A, double *x, double *b)
{
    int i, j, ij;
    double norm_dif = 0, norm_b = 0, mult;
    for (i = 0; i < n; i++)
    {
        mult = 0;
        for (j = 0; j < n; j++)
        {
            ij = index_block_matrix(n, n, m, i, j);
            mult += A[ij] * x[j];
        }
        norm_dif += fabs(mult - b[i]);
    }
    for (i = 0; i < n; i++)
    {
        norm_b += fabs(b[i]);
    }
    return norm_dif / norm_b;
}

double calculation_r2(int n, double *x)
{
    int i;
    double norm_dif = 0, norm_x1 = (n + 1) / 2;
    for (i = 0; i < n; i++)
    {
        norm_dif += fabs(x[i] - ((i + 1) % 2));
    }
    return norm_dif / norm_x1;
}

double norm_matrix(int m, double *B)
{
    int i, j;
    double max_sum = 0, sum;
    for (i = 0; i < m; i++)
    {
        sum = 0;
        for (j = 0; j < m; ++j)
        {
            sum += fabs(B[i * m + j]);
        }
        if (sum > max_sum)
        {
            max_sum = sum;
        }
    }
    return max_sum;
}

int get_inverse(int m, double norm_A, int *moved_columns, double *B, double *C)
{
    int i, j, t, max_j;
    double max_elem, elem_i;
    double eps = norm_A * 1e-15;
    
    for (i = 0; i < m * m; i++)
    {
        if ((i / m) == (i % m))
        {
            C[i] = 1;
        }
        else
        {
            C[i] = 0;
        }
    }
    for (i = 0; i < m; i++)
    {
        max_elem = B[i * m + i];
        max_j = i;
        for (j = i; j < m; j++)
        {
            if (fabs(B[i * m + j]) > fabs(max_elem))
            {
                max_elem = B[i * m + j];
                max_j = j;
            }
        }
        if (fabs(max_elem) < eps)
        {
            return 1;
        }
        moved_columns[i] = max_j;
        if (i != max_j)
        {
            for (j = 0; j < m; j++)
            {
                SWAP_DOUBLE(B[j * m + i], B[j * m + max_j]);
            }
        }
        for (j = 0; j < m; j++)
        {
            C[i * m + j] /= max_elem;
            if (j > i)
            {
                B[i * m + j] /= max_elem;
            }
            else if (j == i)
            {
                B[i * m + j] = 1;
            }
        }
        for (t = 0; t < m; t++)
        {
            if (t != i)
            {
                elem_i = B[t * m + i];
                for (j = 0; j < m; j++)
                {
                    C[t * m + j] -= (elem_i * C[i * m + j]);
                    if (j > i)
                    {
                        B[t * m + j] -= (elem_i * B[i * m + j]);
                    }
                }
                B[t * m + i] = 0;
            }
        }
    }
    for (i = m - 1; i >= 0; i--)
    {
        t = moved_columns[i];
        if (t != i)
        {
            for (j = 0; j < m; j++)
            {
                SWAP_DOUBLE(C[t * m + j], C[i * m + j]);
            }
        }
    }
    return 0;
}

void multiplication_blocks(int p, int q, int r, double *B, double *C, double *R)
{
    int i, j, t;
    int p3 = p % 3;
    int r3 = r % 3;
    double s00, s01, s02, s10, s11, s12, s20, s21, s22;
    for (i = 0; i < p3; i++)
    {
        for (j = 0; j < r3; j++)
        {
            s00 = 0;
            for (t = 0; t < q; t++)
            {
                s00 += B[i * q + t] * C[t * r + j];
            }
            R[i * r + j] = s00;
        }
        for (; j < r; j += 3)
        {
            s00 = 0; s01 = 0; s02 = 0;
            for (t = 0; t < q; t++)
            {
                s00 += B[i * q + t] * C[t * r + j];
                s01 += B[i * q + t] * C[t * r + j + 1];
                s02 += B[i * q + t] * C[t * r + j + 2];
            }
            R[i * r + j] = s00;
            R[i * r + j + 1] = s01;
            R[i * r + j + 2] = s02;
        }
    }
    for (; i < p; i += 3)
    {
        for (j = 0; j < r3; j++)
        {
            s00 = 0; s10 = 0; s20 = 0;
            for (t = 0; t < q; t++)
            {
                s00 += B[i * q + t] * C[t * r + j];
                s10 += B[(i + 1) * q + t] * C[t * r + j];
                s20 += B[(i + 2) * q + t] * C[t * r + j];
            }
            R[i * r + j] = s00;
            R[(i + 1) * r + j] = s10;
            R[(i + 2) * r + j] = s20;
        }
        for (; j < r; j += 3)
        {
            s00 = 0; s01 = 0; s02 = 0;
            s10 = 0; s11 = 0; s12 = 0;
            s20 = 0; s21 = 0; s22 = 0;
            for (t = 0; t < q; t++)
            {
                s00 += B[i * q + t] * C[t * r + j];
                s01 += B[i * q + t] * C[t * r + j + 1];
                s02 += B[i * q + t] * C[t * r + j + 2];
                s10 += B[(i + 1) * q + t] * C[t * r + j];
                s11 += B[(i + 1) * q + t] * C[t * r + j + 1];
                s12 += B[(i + 1) * q + t] * C[t * r + j + 2];
                s20 += B[(i + 2) * q + t] * C[t * r + j];
                s21 += B[(i + 2) * q + t] * C[t * r + j + 1];
                s22 += B[(i + 2) * q + t] * C[t * r + j + 2];
            }
            R[i * r + j] = s00;
            R[i * r + j + 1] = s01;
            R[i * r + j + 2] = s02;
            R[(i + 1) * r + j] = s10;
            R[(i + 1) * r + j + 1] = s11;
            R[(i + 1) * r + j + 2] = s12;
            R[(i + 2) * r + j] = s20;
            R[(i + 2) * r + j + 1] = s21;
            R[(i + 2) * r + j + 2] = s22;
        }
    }
}

int main_block(int n, int m, int i, double norm_A, double *A, double *B, double *C, int *s)
{
    int j;
    int k = n / m;
    int min_j = -1;
    double norm_block;
    double min_norm = 0;
    for (j = i; j < k; j++)
    {
        get_block(n, m, i, j, A, B);
        if (get_inverse(m, norm_A, s, B, C) == 0)
        {
            norm_block = norm_matrix(m, C);
            if (min_j == -1 || norm_block < min_norm)
            {
                min_norm = norm_block;
                min_j = j;
            }
        }
    }
    return min_j;
}

void swap_columns_blocks(int n, int m, int j1, int j2, double *A)
{
    int i, t, ij1, ij2;
    for (i = 0; i < n; i++)
    {
        for (t = 0; t < m; t++)
        {
            ij1 = index_block_matrix(n, n, m, i, j1 * m + t);
            ij2 = index_block_matrix(n, n, m, i, j2 * m + t);
            SWAP_DOUBLE(A[ij1], A[ij2]);
        }
    }
}

void get_vector_block(int n, int m, int i, double *b, double *C)
{
    int s;
    int k = n / m;
    int l = n % m;
    int h = (i < k ? m : l);
    for (s = 0; s < h; s++)
    {
        C[s] = b[i * m + s];
    }
}

void put_vector_block(int n, int m, int i, double *b, double *D)
{
    int s;
    int k = n / m;
    int l = n % m;
    int h = (i < k ? m : l);
    for (s = 0; s < h; s++)
    {
        b[i * m + s] = D[s];
    }
}

void subtraction_block(int m, int l, double *C, double *D)
{
    int i;
    for (i = 0; i < m * l; i++)
    {
        C[i] -= D[i];
    }
}

int solution_system(int n, int m, double *A, double *b, double *x,
                    int *s, int *moved_columns_blocks, double *C, double *D, double *R)
{
    int i, j, main_j, v, h, t, i1;
    int k = n / m;
    int l = n % m;
    double norm_A = norm_matrix(n, A);
    if (norm_A < 1e-64)
    {
        return 1;
    }
    for (i = 0; i < k; i++)
    {
        main_j = main_block(n, m, i, norm_A, A, C, D, s);
        if (main_j == -1)
        {
            return 1;
        }
        get_block(n, m, i, main_j, A, C);
        get_inverse(m, norm_A, s, C, D);
        for (j = i; j < k || (l != 0 && j == k); j++)
        {
            if (j != main_j)
            {
                v = (j < k ? m : l);
                get_block(n, m, i, j, A, C);
                multiplication_blocks(m, m, v, D, C, R);
                put_block(n, m, i, j, A, R);
            }
        }
        get_vector_block(n, m, i, b, C);
        multiplication_blocks(m, m, 1, D, C, R);
        put_vector_block(n, m, i, b, R);
        swap_columns_blocks(n, m, i, main_j, A);
        moved_columns_blocks[i] = main_j;
        for (t = 0; t < k || (l != 0 && t == k); t++)
        {
            if (t != i)
            {
                h = (t < k ? m : l);
                for (j = i + 1; j < k || (l != 0 && j == k); j++)
                {
                    v = (j < k ? m : l);
                    get_block(n, m, t, i, A, C);
                    get_block(n, m, i, j, A, D);
                    multiplication_blocks(h, m, v, C, D, R);
                    get_block(n, m, t, j, A, C);
                    subtraction_block(h, v, C, R);
                    put_block(n, m, t, j, A, C);
                }
            }
        }
        for (t = 0; t < k || (l != 0 && t == k); t++)
        {
            if (t != i)
            {
                h = (t < k ? m : l);
                get_block(n, m, t, i, A, C);
                get_vector_block(n, m, i, b, D);
                multiplication_blocks(h, m, 1, C, D, R);
                get_vector_block(n, m, t, b, C);
                subtraction_block(h, 1, C, R);
                put_vector_block(n, m, t, b, C);
            }
        }
    }
    if (l != 0)
    {
        get_block(n, m, k, k, A, C);
        if (get_inverse(l, norm_A, s, C, D) != 0)
        {
            return 1;
        }
        get_vector_block(n, m, k, b, C);
        multiplication_blocks(l, l, 1, D, C, R);
        put_vector_block(n, m, k, b, R);
        get_vector_block(n, m, k, b, C);
        for (t = 0; t < k; t++)
        {
            get_block(n, m, t, k, A, D);
            multiplication_blocks(m, l, 1, D, C, R);
            get_vector_block(n, m, t, b, D);
            subtraction_block(m, 1, D, R);
            put_vector_block(n, m, t, b, D);
        }
    }
    for (i = k - 1; i >= 0; i--)
    {
        if ((i1 = moved_columns_blocks[i]) != i)
        {
            for (j = 0; j < m; j++)
            {
                SWAP_DOUBLE(b[i * m + j], b[i1 * m + j]);
            }
        }
    }
    for (i = 0; i < n; i++)
    {
        x[i] = b[i];
    }
    return 0;
}