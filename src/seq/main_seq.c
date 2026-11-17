#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "csr_matrix.h"
#include "utils.h"

#define DEFAULT_MAX_ITER 5000
#define DEFAULT_NUM_RUNS 10
#define DEFAULT_TOL 0.0



/*
    * Função gauss_jacobi_seq: Executa o método iterativo de Gauss-Jacobi de forma sequencial.

    * Parâmetros:
        - A: Ponteiro para a matriz esparsa CSR.
        - b: Vetor de termos independentes.
        - x: Vetor solução inicial (atualizado com o resultado final ao término).
        - next: Vetor auxiliar para armazenar os valores da próxima iteração.
        - max_iter: Número máximo de iterações permitidas.
        - tol: Tolerância mínima para o critério de parada (se <= 0.0, executa todas as iterações).
        - final_diff: Ponteiro de saída que recebe o resíduo final |x^(k+1) - x^(k)|.

    * Retorna:
        - Quantidade total de iterações executadas.
*/
int gauss_jacobi_seq(const CSRMatrix *A, const double *b, double *x, double *next,
                     int max_iter, double tol, double *final_diff) {
    int N = A->num_rows;
    double max_diff = 0.0;
    int iter = 0;

    while (iter < max_iter) {
        max_diff = 0.0;

        for (int i = 0; i < N; i++) {
            double sum = 0.0;

            for (int k = A->row_ptr[i]; k < A->row_ptr[i + 1]; k++) {
                int j = A->col_idx[k];
                if (i != j) {
                    sum += A->values[k] * x[j];
                }
            }

            next[i] = (b[i] - sum) / A->diag[i];

            double diff = fabs(next[i] - x[i]);
            if (diff > max_diff) {
                max_diff = diff;
            }
        }

        for (int i = 0; i < N; i++) {
            x[i] = next[i];
        }
        iter++;

        if (tol > 0.0 && max_diff < tol) {
            break;
        }
    }

    *final_diff = max_diff;
    return iter;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Uso: %s <arquivo.mtx> [num_execucoes] [max_iter] [tol]\n", argv[0]);
        return 1;
    }

    int num_runs = (argc >= 3 && atoi(argv[2]) > 0) ? atoi(argv[2]) : DEFAULT_NUM_RUNS;
    int max_iter = (argc >= 4 && atoi(argv[3]) > 0) ? atoi(argv[3]) : DEFAULT_MAX_ITER;
    double tol   = (argc >= 5) ? atof(argv[4]) : DEFAULT_TOL;

    printf("=================================================================\n");
    printf("INICIALIZAÇÃO DO MÉTODO DE GAUSS-JACOBI (SEQUENCIAL)\n");
    printf("=================================================================\n");
    printf("Matriz Alvo: %s | Repetições: %d | Teto de Iterações: %d\n", argv[1], num_runs, max_iter);
    if (tol > 0.0) {
        printf("Critério de Parada: Tolerância residual <= %.1e\n", tol);
    } else {
        printf("Critério de Parada: Carga Fixa (%d iterações completas)\n", max_iter);
    }
    fflush(stdout);

    /*
        1. Carregamento e Validação da Matriz Esparsa (Medição de I/O + Construção CSR).
    */
    CSRMatrix A;
    double t_load_start = get_time_sec();
    if (load_matrix_market_csr(argv[1], &A) != 0) {
        return 1;
    }
    double t_load_end = get_time_sec();
    double load_time  = t_load_end - t_load_start;

    int N = A.num_rows;
    printf("Dimensão: %d x %d | Não-zeros (NNZ): %d\n", N, A.num_cols, A.nnz);
    printf("Tempo de Carregamento (I/O + CSR): %.6f s\n\n", load_time);
    fflush(stdout);

    if (N <= 0) {
        fprintf(stderr, "Erro: dimensão inválida da matriz (N = %d)!\n", N);
        free_csr(&A);
        return 1;
    }

    for (int i = 0; i < N; i++) {
        if (fabs(A.diag[i]) < 1e-15) {
            fprintf(stderr, "Erro: Elemento zero na diagonal na linha %d!\n", i);
            free_csr(&A);
            return 1;
        }
    }

    /*
        2. Alocação e Inicialização dos Vetores.
    */
    size_t n_size     = (size_t)N;
    double *b         = (double *)malloc(n_size * sizeof(double));
    double *x         = (double *)malloc(n_size * sizeof(double));
    double *next      = (double *)malloc(n_size * sizeof(double));
    double *run_times = (double *)calloc((size_t)num_runs, sizeof(double));

    if (!b || !x || !next || !run_times) {
        fprintf(stderr, "Erro de alocação de memória para os vetores!\n");
        free(b);
        free(x);
        free(next);
        free(run_times);
        free_csr(&A);
        return 1;
    }

    initialize_rhs(&A, b);

    /*
        3. Execução do Benchmark (Loop de Repetições do Solver).
    */
    int last_iter         = 0;
    double last_max_diff  = 0.0;
    double last_sol_error = 0.0;

    for (int run = 0; run < num_runs; run++) {
        memset(x, 0, n_size * sizeof(double));

        double t_start = get_time_sec();
        last_iter = gauss_jacobi_seq(&A, b, x, next, max_iter, tol, &last_max_diff);
        double t_end = get_time_sec();

        run_times[run] = t_end - t_start;
        last_sol_error = compute_max_error(x, N);

        printf("  [Execução %2d/%d] Tempo: %.6f s | Iterações: %d | Resíduo: %.4e\n",
               run + 1, num_runs, run_times[run], last_iter, last_max_diff);
        fflush(stdout);
    }

    /*
        4. Cálculo das Estatísticas Consolidadas.
    */
    double sum_time = 0.0;
    double min_time = run_times[0];
    double max_time = run_times[0];

    for (int r = 0; r < num_runs; r++) {
        sum_time += run_times[r];
        if (run_times[r] < min_time) min_time = run_times[r];
        if (run_times[r] > max_time) max_time = run_times[r];
    }

    double mean_time = sum_time / num_runs;
    double variance  = 0.0;

    if (num_runs > 1) {
        for (int r = 0; r < num_runs; r++) {
            double diff = run_times[r] - mean_time;
            variance += diff * diff;
        }
        variance /= (num_runs - 1);
    }

    double std_time     = sqrt(variance);
    double total_flops  = (double)last_iter * (2.0 * A.nnz + 2.0 * N);
    double mean_mflops  = (mean_time > 0.0) ? (total_flops / (mean_time * 1e6)) : 0.0;
    double mean_iter_ms = (mean_time > 0.0 && last_iter > 0) ? ((mean_time / last_iter) * 1000.0) : 0.0;

    /*
        5. Exibição dos Resultados Consolidados (Separação Clara entre I/O e Solver).
    */
    char status_str[64];
    if (tol > 0.0) {
        snprintf(status_str, sizeof(status_str), "%s", (last_max_diff < tol) ? "CONVERGIU" : "LIMITE DE ITERAÇÕES");
    } else {
        snprintf(status_str, sizeof(status_str), "CONCLUÍDO (%d ITERAÇÕES)", last_iter);
    }

    printf("\nRESULTADOS CONSOLIDADOS (SEQUENCIAL - %d EXECUÇÕES)\n", num_runs);
    printf("Matriz avaliada: %s\n", argv[1]);
    printf("Status da convergência: %s\n", status_str);
    printf("Dimensão da matriz (N): %d\n", N);
    printf("Elementos não-nulos: %d (Densidade: %.4f%%)\n", A.nnz, ((double)A.nnz / ((double)N * N)) * 100.0);
    if (tol > 0.0) {
        printf("Tolerância estipulada: %.1e\n", tol);
    } else {
        printf("Tolerância estipulada: Desativada (Carga Fixa)\n");
    }
    printf("Iterações executadas: %d / %d\n", last_iter, max_iter);
    printf("Resíduo final (|x-x_old|): %.10e (%.12f)\n", last_max_diff, last_max_diff);
    printf("Erro absoluto máx (|x-1|): %.10e (%.12f)\n", last_sol_error, last_sol_error);

    printf("\n--- TEMPOS DE EXECUÇÃO ---\n");
    printf("- Tempo de Carregamento (I/O + CSR): %.6f s\n", load_time);
    printf("- Tempo Médio do Solver (%dx):        %.6f s (± %.6f s)\n", num_runs, mean_time, std_time);
    printf("- Tempo Mínimo / Máximo do Solver:   %.6f s / %.6f s\n", min_time, max_time);
    printf("- Tempo Acumulado do Solver (%dx):   %.6f s\n", num_runs, sum_time);
    printf("- Tempo Total (I/O + Solver Médio):  %.6f s\n", load_time + mean_time);
    printf("- Tempo Total Acumulado (I/O + %dx): %.6f s\n", num_runs, load_time + sum_time);
    printf("- Tempo Médio por Iteração:          %.6f ms\n", mean_iter_ms);
    printf("- Throughput Médio Estimado:         %.2f MFlops/s\n\n", mean_mflops);
    fflush(stdout);

    /*
        6. Salvamento dos Resultados Consolidados em CSV.
    */
    FILE *f_csv = fopen("results/results_seq.csv", "r");
    int write_header = (f_csv == NULL);
    if (f_csv) {
        fclose(f_csv);
    }

    f_csv = fopen("results/results_seq.csv", "a");
    if (f_csv) {
        if (write_header) {
            fprintf(f_csv, "matriz,dimensao,nnz,densidade_pct,repeticoes,iteracoes,max_iter,convergencia,"
                           "tempo_io_csr_s,tempo_solver_medio_s,desvio_padrao_s,tempo_solver_min_s,tempo_solver_max_s,"
                           "tempo_solver_acumulado_s,tempo_total_s,tempo_medio_iter_ms,mflops,residuo_final,erro_max_solucao\n");
        }
        fprintf(f_csv, "%s,%d,%d,%.6f,%d,%d,%d,%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.2f,%.10e,%.10e\n",
                argv[1], N, A.nnz, ((double)A.nnz / ((double)N * N)) * 100.0,
                num_runs, last_iter, max_iter, status_str,
                load_time, mean_time, std_time, min_time, max_time,
                sum_time, load_time + mean_time, mean_iter_ms, mean_mflops, last_max_diff, last_sol_error);
        fclose(f_csv);
    }

    /*
        7. Liberação de Recursos.
    */
    free(run_times);
    free(b);
    free(x);
    free(next);
    free_csr(&A);

    return 0;
}
