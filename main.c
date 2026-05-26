#include <stdio.h>
#include <stdlib.h>
#include "experiments.h"

int main(void) {
    srand(42);  /* фиксированный seed для воспроизводимости */

    printf("Сравнение методов решения СЛАУ\n");
    printf("Все эксперименты воспроизводимы: используется srand(42)\n");

    /* если компьютер не тянет большие размеры, можно уменьшить */
    int single_sizes[]  = {100, 200, 500};
    int rhs_counts[]    = {1, 10, 100};
    int hilbert_sizes[] = {5, 10, 15};

    experiment_single_system(single_sizes, 3);
    experiment_multiple_rhs(500, rhs_counts, 3);
    experiment_hilbert(hilbert_sizes, 3);

    return 0;
}
