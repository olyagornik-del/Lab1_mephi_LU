#include <stdio.h>
#include <stdlib.h>
#include "matrices.h"
#include "metrics.h"
#include "solvers.h"
#include "timing.h"
#include "experiments.h"

/* 4.1 - время решения одной системы разными методами */
void experiment_single_system(const int *sizes, int count) {
    printf("\n=== 4.1 Сравнение времени решения одной системы ===\n");
    printf("%-5s %-11s %-11s %-11s %-11s %-11s %-13s %-13s\n",
           "n", "Gauss", "Gauss piv", "LU decomp", "LU solve",
           "LU total", "Res (pivot)", "Res (LU)");

    for (int s = 0; s < count; s++) {
        int n = sizes[s];
        double **a = random_matrix(n);
        double *b = random_vector(n);
        double *x = alloc_vector(n);
        double t0;

        t0 = now_seconds();
        int ok_gauss = gaussian_no_pivot(a, b, n, x);
        double t_gauss = now_seconds() - t0;

        t0 = now_seconds();
        gaussian_partial_pivot(a, b, n, x);
        double t_pivot = now_seconds() - t0;
        double res_pivot = residual_norm(a, x, b, n);   /* x сейчас от метода с выбором */

        double **l = alloc_matrix(n);
        double **u = alloc_matrix(n);

        t0 = now_seconds();
        lu_decomposition(a, n, l, u);
        double t_decomp = now_seconds() - t0;

        t0 = now_seconds();
        solve_with_lu(l, u, b, n, x);
        double t_solve = now_seconds() - t0;
        double res_lu = residual_norm(a, x, b, n);

        printf("%-5d ", n);
        if (ok_gauss == 0)
            printf("%-11.6f ", t_gauss);
        else
            printf("%-11s ", "error");
        printf("%-11.6f %-11.6f %-11.6f %-11.6f %-13.3e %-13.3e\n",
               t_pivot, t_decomp, t_solve, t_decomp + t_solve, res_pivot, res_lu);

        free_matrix(a, n);
        free_matrix(l, n);
        free_matrix(u, n);
        free(b);
        free(x);
    }
}

/* 4.2 - выгода LU при нескольких правых частях */
void experiment_multiple_rhs(int n, const int *k_values, int count) {
    printf("\n=== 4.2 Экономия времени при множественных правых частях ===\n");
    printf("%-6s %-18s %-11s %-15s %-11s\n",
           "k", "Gauss pivot total", "LU decomp", "LU all solves", "LU total");

    double **a = random_matrix(n);
    double **l = alloc_matrix(n);
    double **u = alloc_matrix(n);
    double *x = alloc_vector(n);

    double t0 = now_seconds();
    lu_decomposition(a, n, l, u);
    double t_decomp = now_seconds() - t0;

    for (int s = 0; s < count; s++) {
        int k = k_values[s];

        /* k случайных правых частей */
        double **rhs = malloc(k * sizeof(double *));
        for (int i = 0; i < k; i++)
            rhs[i] = random_vector(n);

        /* Гаусс: каждую систему решаем заново целиком */
        t0 = now_seconds();
        for (int i = 0; i < k; i++)
            gaussian_partial_pivot(a, rhs[i], n, x);
        double t_gauss = now_seconds() - t0;

        /* LU: разложение уже готово, для каждой части - только подстановки */
        t0 = now_seconds();
        for (int i = 0; i < k; i++)
            solve_with_lu(l, u, rhs[i], n, x);
        double t_lu = now_seconds() - t0;

        printf("%-6d %-18.6f %-11.6f %-15.6f %-11.6f\n",
               k, t_gauss, t_decomp, t_lu, t_decomp + t_lu);

        for (int i = 0; i < k; i++)
            free(rhs[i]);
        free(rhs);
    }

    free_matrix(a, n);
    free_matrix(l, n);
    free_matrix(u, n);
    free(x);
}

/* печать одной строки результатов для эксперимента 4.3 */
static void print_hilbert_row(int n, const char *method, int ok,
                              double rel_err, double res) {
    if (ok == 0)
        printf("%-4d %-13s %-13.3e %-13.3e\n", n, method, rel_err, res);
    else
        printf("%-4d %-13s %-13s %-13s\n", n, method, "error", "-");
}

/* 4.3 - точность на матрицах Гильберта */
void experiment_hilbert(const int *sizes, int count) {
    printf("\n=== 4.3 Проверка точности на матрицах Гильберта ===\n");
    printf("%-4s %-13s %-13s %-13s\n", "n", "Method", "Rel error", "Residual");

    for (int s = 0; s < count; s++) {
        int n = sizes[s];
        double **h = hilbert_matrix(n);
        double *x_exact = alloc_vector(n);
        for (int i = 0; i < n; i++)
            x_exact[i] = 1.0;
        double *b = mat_vec(h, x_exact, n);    /* b = H * x_exact */
        double *x = alloc_vector(n);

        int ok = gaussian_no_pivot(h, b, n, x);
        print_hilbert_row(n, "Gauss", ok,
                          relative_error(x, x_exact, n), residual_norm(h, x, b, n));

        ok = gaussian_partial_pivot(h, b, n, x);
        print_hilbert_row(n, "Gauss pivot", ok,
                          relative_error(x, x_exact, n), residual_norm(h, x, b, n));

        double **l = alloc_matrix(n);
        double **u = alloc_matrix(n);
        ok = lu_decomposition(h, n, l, u);
        if (ok == 0)
            ok = solve_with_lu(l, u, b, n, x);
        print_hilbert_row(n, "LU", ok,
                          relative_error(x, x_exact, n), residual_norm(h, x, b, n));

        free_matrix(l, n);
        free_matrix(u, n);
        free_matrix(h, n);
        free(x_exact);
        free(b);
        free(x);
    }
}
