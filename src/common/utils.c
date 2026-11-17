#define _POSIX_C_SOURCE 200809L
#include "utils.h"
#include <time.h>
#include <math.h>



/*
    * Função get_time_sec: Obtém o tempo atual do sistema em segundos com precisão de nanossegundos.

    * Retorna:
        - Tempo monotônico decorrido em segundos (double).
*/
double get_time_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/*
    * Função compute_max_error: Calcula o maior erro absoluto entre a solução obtida e x* = 1.

    * Parâmetros:
        - x: Vetor solução calculado pelo solver.
        - N: Dimensão do vetor.

    * Retorna:
        - Erro absoluto máximo max_i |x[i] - 1.0|.
*/
double compute_max_error(const double *x, int N) {
    double max_err = 0.0;
    for (int i = 0; i < N; i++) {
        double err = fabs(x[i] - 1.0);
        if (err > max_err) {
            max_err = err;
        }
    }
    return max_err;
}