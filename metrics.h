#ifndef METRICS_H
#define METRICS_H

// невязка ‖Ax - b‖
double residual_norm(double **a, double *x, double *b, int n);

// относительная погрешность ‖x_approx - x_exact‖ / ‖x_exact‖
double relative_error(double *x_approx, double *x_exact, int n);

#endif
