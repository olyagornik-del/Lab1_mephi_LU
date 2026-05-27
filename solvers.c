#include <stdlib.h>
#include <math.h>
#include "matrices.h"
#include "solvers.h"

#define EPS 1e-12

// решаем Ux = y, идём снизу вверх 
int backward_substitution(double **u, double *y, int n, double *x) {
    for (int i = n - 1; i >= 0; i--) {
        double s = 0.0;
        for (int j = i + 1; j < n; j++)
            s += u[i][j] * x[j];
        if (fabs(u[i][i]) < EPS)
            return -1;
        x[i] = (y[i] - s) / u[i][i];
    }
    return 0;
}

// решаем Ly = b, идём сверху вниз 
int forward_substitution(double **l, double *b, int n, double *y) {
    for (int i = 0; i < n; i++) {
        double s = 0.0;
        for (int j = 0; j < i; j++)
            s += l[i][j] * y[j];
        if (fabs(l[i][i]) < EPS)
            return -1;
        y[i] = (b[i] - s) / l[i][i];
    }
    return 0;
}

// классический метод Гаусса (без выбора ведущего элемента) 
int gaussian_no_pivot(double **a_in, double *b_in, int n, double *x) {
    double **a = copy_matrix(a_in, n);
    double *b = copy_vector(b_in, n);
    int status = 0;

    for (int k = 0; k < n - 1; k++) {
        if (fabs(a[k][k]) < EPS) { status = -1; break; }
        for (int i = k + 1; i < n; i++) {
            double factor = a[i][k] / a[k][k];
            a[i][k] = 0.0;
            for (int j = k + 1; j < n; j++)  // пересчитываем остальные элементы строки 
                a[i][j] -= factor * a[k][j];
            b[i] -= factor * b[k];
        }
    }

    if (status == 0)
        status = backward_substitution(a, b, n, x);

    free_matrix(a, n);
    free(b);
    return status;
}

// метод Гаусса с частичным выбором ведущего элемента по столбцу 
int gaussian_partial_pivot(double **a_in, double *b_in, int n, double *x) {
    double **a = copy_matrix(a_in, n);
    double *b = copy_vector(b_in, n);
    int status = 0;

    for (int k = 0; k < n - 1; k++) {
        int pivot_row = k;                       // ищем макс эл в столбце 
        for (int i = k + 1; i < n; i++)
            if (fabs(a[i][k]) > fabs(a[pivot_row][k]))
                pivot_row = i;

        if (fabs(a[pivot_row][k]) < EPS) { status = -1; break; }

        if (pivot_row != k) {                    // меняем строки местами 
            double *trow = a[k]; a[k] = a[pivot_row]; a[pivot_row] = trow;
            double tb = b[k]; b[k] = b[pivot_row]; b[pivot_row] = tb;
        }

        for (int i = k + 1; i < n; i++) {
            double factor = a[i][k] / a[k][k];
            a[i][k] = 0.0;
            for (int j = k + 1; j < n; j++)
                a[i][j] -= factor * a[k][j];
            b[i] -= factor * b[k];
        }
    }

    if (status == 0)
        status = backward_substitution(a, b, n, x);

    free_matrix(a, n);
    free(b);
    return status;
}

// LU-разложение без перестановок: A = L * U 
int lu_decomposition(double **a, int n, double **l, double **u) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            l[i][j] = (i == j) ? 1.0 : 0.0;
            u[i][j] = 0.0;
        }

    for (int k = 0; k < n; k++) {
        for (int j = k; j < n; j++) {
            double s = 0.0;
            for (int t = 0; t < k; t++)
                s += l[k][t] * u[t][j];
            u[k][j] = a[k][j] - s;
            // тк a_kj = l_k0*u_0j + ... + l_kk-1*u_k-1j + l_kk*u_kj,
            // а l_kk = 1, то u_kj = a_kj минус накопленная сумма
        }

        if (fabs(u[k][k]) < EPS)
            return -1;

        for (int i = k + 1; i < n; i++) {
            double s = 0.0;
            for (int t = 0; t < k; t++)
                s += l[i][t] * u[t][k];
            l[i][k] = (a[i][k] - s) / u[k][k];
        }
    }
    return 0;
}

int solve_with_lu(double **l, double **u, double *b, int n, double *x) {
    double *y = alloc_vector(n);
    int status = forward_substitution(l, b, n, y);
    if (status == 0)
        status = backward_substitution(u, y, n, x);
    free(y);
    return status;
}
