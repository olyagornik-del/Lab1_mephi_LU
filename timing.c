#define _POSIX_C_SOURCE 199309L
#include <time.h>
#include "timing.h"

double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec / 1e9;
}
