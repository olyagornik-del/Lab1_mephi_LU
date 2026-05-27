#ifndef MATRICES_H
#define MATRICES_H

// выделение и освобождение памяти 
double **alloc_matrix(int n);
void free_matrix(double **a, int n);
double *alloc_vector(int n);

double **copy_matrix(double **a, int n);
double *copy_vector(double *b, int n);

// генерация тестовых данных 
double **random_matrix(int n);
double *random_vector(int n);
double **hilbert_matrix(int n);

// базовые операции (результат - новый массив, его нужно освободить) 
double *mat_vec(double **a, double *x, int n);
double *vector_sub(double *a, double *b, int n);
double norm2(double *v, int n);

#endif
