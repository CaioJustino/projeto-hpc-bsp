#ifndef UTILS_H
#define UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

/*
    * Função get_time_sec: Obtém o tempo atual do sistema em segundos com precisão de nanossegundos.

    * Retorna:
        - Tempo monotônico decorrido em segundos (double).
*/
double get_time_sec(void);

/*
    * Função compute_max_error: Calcula o maior erro absoluto entre a solução obtida e x* = 1.

    * Parâmetros:
        - x: Vetor solução calculado pelo solver.
        - N: Dimensão do vetor.

    * Retorna:
        - Erro absoluto máximo max_i |x[i] - 1.0|.
*/
double compute_max_error(const double *x, int N);

#ifdef __cplusplus
}
#endif

#endif /* UTILS_H */