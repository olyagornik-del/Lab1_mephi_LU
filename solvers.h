#ifndef SOLVERS_H
#define SOLVERS_H

/* Все функции возвращают 0 при успехе и -1, если матрица вырождена
   (слишком маленький ведущий элемент). Результат пишется в заранее
   выделенные массивы x / y / L / U. */

int backward_substitution(double **u, double *y, int n, double *x);
int forward_substitution(double **l, double *b, int n, double *y);

int gaussian_no_pivot(double **a, double *b, int n, double *x);
int gaussian_partial_pivot(double **a, double *b, int n, double *x);

int lu_decomposition(double **a, int n, double **l, double **u);
int solve_with_lu(double **l, double **u, double *b, int n, double *x);

#endif
