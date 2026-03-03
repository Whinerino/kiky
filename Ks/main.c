#include "functions.h"

int main(int argc, char *argv[])
{
    int n, m, r, s, k;
    int result = 0;
    int solution_result = 0;
    double r1 = -1;
    double r2 = -1;
    double solution_time = 0;
    double discrepancy_time = 0;
    char *file_name = NULL;
    double *A, *b, *x;
    int *h, *g;
    double *C, *D, *R;

    if (!(
        (argc == 5 || argc == 6) &&
        (sscanf(argv[1], "%d", &n) == 1) && (n > 0) &&
        (sscanf(argv[2], "%d", &m) == 1) && (m > 0) &&
        (sscanf(argv[3], "%d", &r) == 1) && (r > 0) &&
        (sscanf(argv[4], "%d", &s) == 1) &&
        (s == 0 || s == 1 || s == 2 || s == 3 || s == 4)
    ))
    {
        printf("Usage = %s n m r s\n", argv[0]);
        return 1;
    }

    if (s == 0)
    {
        if (argc == 5)
        {
            printf("Usage: %s n m r s [file_name]\n", argv[0]);
            return 1;
        }
        file_name = argv[5];
    }
    else if (argc == 6)
    {
        printf("Usage: %s n m r s\n", argv[0]);
        return 1;
    }

    A = (double *)malloc(n * n * sizeof(double));
    result = get_matrix(argc, n, m, s, A, file_name);
    if (result == 1)
    {
        printf("The file did not open.\n");
        free(A);
        return 1;
    }
    else if (result == 2)
    {
        printf("The element did not read.\n");
        free(A);
        return 1;
    }
    else if (result == 3)
    {
        printf("There are extra elements in the file.\n");
        free(A);
        return 1;
    }
    else if (result != 0)
    {
        printf("Something is wrong with getting the matrix.\n");
        free(A);
        return 1;
    }

    print_matrix(n, n, m, r, A);
    printf("\n");

    b = (double *)malloc(n * sizeof(double));
    get_vector(n, m, A, b);
    x = (double *)malloc(n * sizeof(double));

    k = n / m;
    h = (int *)malloc(m * sizeof(int));
    g = (int *)malloc(k * sizeof(int));
    C = (double *)malloc(m * m * sizeof(double));
    D = (double *)malloc(m * m * sizeof(double));
    R = (double *)malloc(m * m * sizeof(double));

    solution_time = clock();
    solution_result = solution_system(n, m, A, b, x, h, g, C, D, R);
    solution_time = (clock() - solution_time) / (double)CLOCKS_PER_SEC;

    result = get_matrix(argc, n, m, s, A, file_name);
    if (result != 0)
    {
        printf("Something is wrong with getting the matrix.\n");
        free(A); free(b); free(x); free(h); free(g); free(C); free(D); free(R);
        return 1;
    }

    get_vector(n, m, A, b);

    if (solution_result == 0)
    {
        discrepancy_time = clock();
        r1 = calculation_r1(n, m, A, x, b);
        r2 = calculation_r2(n, x);
        discrepancy_time = (clock() - discrepancy_time) / (double)CLOCKS_PER_SEC;
    }

    if (solution_result == 0)
    {
        printf("Solution:\n");
        print_matrix(1, n, m, r, x);
        printf("\n");
    }
    else
    {
        printf("The program could not find a solution.\n");
    }

    /* ИСПРАВЛЕННЫЙ ФОРМАТ: без пробелов вокруг : и = */
    printf("%s: Task=%d Res1=%e Res2=%e T1=%.2f T2=%.2f S=%d N=%d M=%d\n",
           argv[0], 16, r1, r2, solution_time, discrepancy_time, s, n, m);

    free(A);
    free(b);
    free(x);
    free(h);
    free(g);
    free(C);
    free(D);
    free(R);

    return 0;
}