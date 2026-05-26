#ifndef EXPERIMENTS_H
#define EXPERIMENTS_H

void experiment_single_system(const int *sizes, int count);
void experiment_multiple_rhs(int n, const int *k_values, int count);
void experiment_hilbert(const int *sizes, int count);

#endif
