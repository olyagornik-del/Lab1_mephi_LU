#include <stdlib.h>
#include "matrices.h"
#include "metrics.h"

// невязка - проверяет, хорошо ли выполняется уравнение Ax = b
double residual_norm(double **a, double *x, double *b, int n) {
    double *ax = mat_vec(a, x, n);
    double *diff = vector_sub(ax, b, n);
    double r = norm2(diff, n);
    free(ax);
    free(diff);
    return r;
}

// относительная погрешность - насколько решение близко к точному
double relative_error(double *x_approx, double *x_exact, int n) {
    double *diff = vector_sub(x_approx, x_exact, n);
    double num = norm2(diff, n);
    double den = norm2(x_exact, n);
    free(diff);
    if (den == 0.0)
        return num;
    return num / den;
}
