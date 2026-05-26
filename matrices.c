#include <stdlib.h>
#include <math.h>
#include "matrices.h"

double **alloc_matrix(int n) {
    double **a = malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
        a[i] = malloc(n * sizeof(double));
    return a;
}

void free_matrix(double **a, int n) {
    for (int i = 0; i < n; i++)
        free(a[i]);
    free(a);
}

double *alloc_vector(int n) {
    return malloc(n * sizeof(double));
}

double **copy_matrix(double **a, int n) {
    double **r = alloc_matrix(n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            r[i][j] = a[i][j];
    return r;
}

double *copy_vector(double *b, int n) {
    double *r = alloc_vector(n);
    for (int i = 0; i < n; i++)
        r[i] = b[i];
    return r;
}

/* случайное число из [-1, 1] */
static double rand_uniform(void) {
    return 2.0 * rand() / RAND_MAX - 1.0;
}

double **random_matrix(int n) {
    double **a = alloc_matrix(n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            a[i][j] = rand_uniform();
    return a;
}

double *random_vector(int n) {
    double *v = alloc_vector(n);
    for (int i = 0; i < n; i++)
        v[i] = rand_uniform();
    return v;
}

/* Hij = 1 / (i + j - 1) при нумерации с 1; у нас индексы с 0 */
double **hilbert_matrix(int n) {
    double **h = alloc_matrix(n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            h[i][j] = 1.0 / (i + j + 1);
    return h;
}

/* умножение матрицы на вектор */
double *mat_vec(double **a, double *x, int n) {
    double *r = alloc_vector(n);
    for (int i = 0; i < n; i++) {
        double s = 0.0;
        for (int j = 0; j < n; j++)
            s += a[i][j] * x[j];
        r[i] = s;
    }
    return r;
}

double *vector_sub(double *a, double *b, int n) {
    double *r = alloc_vector(n);
    for (int i = 0; i < n; i++)
        r[i] = a[i] - b[i];
    return r;
}

/* длинна вектора */
double norm2(double *v, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++)
        s += v[i] * v[i];
    return sqrt(s);
}
