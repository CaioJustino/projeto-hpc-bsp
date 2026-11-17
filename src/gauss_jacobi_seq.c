#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define DEFAULT_MAX_ITER 5000
#define DEFAULT_NUM_RUNS 10
#define TOL 1e-6



/*
    * Estrutura CSRMatrix: Representa uma matriz esparsa no formato Compressed Sparse Row (CSR).

    * Campos:
        - num_rows: Número de linhas da matriz.
        - num_cols: Número de colunas da matriz.
        - nnz: Quantidade total de elementos não-nulos armazenados.
        - row_ptr: Vetor de deslocamentos que indica onde cada linha começa em col_idx e values.
        - col_idx: Vetor com os índices de coluna de cada elemento não-nulo.
        - values: Vetor com os valores numéricos não-nulos.
        - diag: Vetor auxiliar com os elementos da diagonal principal para acesso direto.
*/
typedef struct {
    int num_rows;
    int num_cols;
    int nnz;
    int *row_ptr;
    int *col_idx;
    double *values;
    double *diag;
} CSRMatrix;

/*
    * Função get_time_sec: Obtém o tempo atual do sistema em segundos com precisão de nanossegundos.

    * Retorna:
        - Tempo monotônico decorrido em segundos (double).
*/
static inline double get_time_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/*
    * Função load_matrix_market_csr: Lê um arquivo Matrix Market (.mtx) e converte para o formato CSR.

    * Parâmetros:
        - filename: Caminho do arquivo .mtx de entrada.
        - mat: Ponteiro para a estrutura CSRMatrix que será preenchida.

    * Retorna:
        - 0 em caso de sucesso, ou -1 se ocorrer erro de leitura ou alocação.
*/
int load_matrix_market_csr(const char *filename, CSRMatrix *mat) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Erro ao abrir o arquivo: %s\n", filename);
        return -1;
    }

    /*
        1. Leitura do cabeçalho e identificação de simetria
    */
    char line[1024];
    if (!fgets(line, sizeof(line), f)) {
        fclose(f);
        return -1;
    }
    int is_symmetric = (strstr(line, "symmetric") != NULL);

    while (fgets(line, sizeof(line), f)) {
        if (line[0] != '%') break;
    }

    int rows, cols, entries;
    if (sscanf(line, "%d %d %d", &rows, &cols, &entries) != 3) {
        fprintf(stderr, "Erro ao ler as dimensões da matriz.\n");
        fclose(f);
        return -1;
    }

    /*
        2. Leitura temporária das coordenadas (formato COO)
    */
    int max_entries = is_symmetric ? entries * 2 : entries;
    int *coo_r    = (int *)malloc((size_t)max_entries * sizeof(int));
    int *coo_c    = (int *)malloc((size_t)max_entries * sizeof(int));
    double *coo_v = (double *)malloc((size_t)max_entries * sizeof(double));

    if (!coo_r || !coo_c || !coo_v) {
        fprintf(stderr, "Erro de alocação temporária para leitura COO.\n");
        free(coo_r);
        free(coo_c);
        free(coo_v);
        fclose(f);
        return -1;
    }

    int actual_nnz = 0;
    for (int k = 0; k < entries; k++) {
        int r, c;
        double v;
        if (fscanf(f, "%d %d %lf", &r, &c, &v) != 3) break;
        r--;
        c--;

        coo_r[actual_nnz] = r;
        coo_c[actual_nnz] = c;
        coo_v[actual_nnz] = v;
        actual_nnz++;

        /* Espelha elementos fora da diagonal quando a matriz é simétrica */
        if (is_symmetric && r != c) {
            coo_r[actual_nnz] = c;
            coo_c[actual_nnz] = r;
            coo_v[actual_nnz] = v;
            actual_nnz++;
        }
    }
    fclose(f);

    /*
        3. Conversão de COO para CSR
    */
    mat->num_rows = rows;
    mat->num_cols = cols;
    mat->nnz      = actual_nnz;
    mat->row_ptr  = (int *)calloc((size_t)(rows + 1), sizeof(int));
    mat->col_idx  = (int *)malloc((size_t)actual_nnz * sizeof(int));
    mat->values   = (double *)malloc((size_t)actual_nnz * sizeof(double));
    mat->diag     = (double *)calloc((size_t)rows, sizeof(double));

    for (int k = 0; k < actual_nnz; k++) {
        mat->row_ptr[coo_r[k] + 1]++;
    }

    for (int i = 0; i < rows; i++) {
        mat->row_ptr[i + 1] += mat->row_ptr[i];
    }

    int *current_pos = (int *)malloc((size_t)rows * sizeof(int));
    memcpy(current_pos, mat->row_ptr, (size_t)rows * sizeof(int));

    for (int k = 0; k < actual_nnz; k++) {
        int r = coo_r[k];
        int dest = current_pos[r]++;

        mat->col_idx[dest] = coo_c[k];
        mat->values[dest]  = coo_v[k];

        if (r == coo_c[k]) {
            mat->diag[r] = coo_v[k];
        }
    }

    free(coo_r);
    free(coo_c);
    free(coo_v);
    free(current_pos);
    return 0;
}

/*
    * Função free_csr: Libera a memória alocada dinamicamente para a matriz CSR.

    * Parâmetros:
        - mat: Ponteiro para a estrutura CSRMatrix a ser liberada.
*/
void free_csr(CSRMatrix *mat) {
    free(mat->row_ptr);
    free(mat->col_idx);
    free(mat->values);
    free(mat->diag);
}

/*
    * Função initialize_rhs: Constrói o vetor independente b somando os elementos de cada linha,
      garantindo que a solução analítica esperada seja o vetor unitário x* = [1, 1, ..., 1]^T.

    * Parâmetros:
        - A: Ponteiro para a matriz esparsa CSR.
        - b: Ponteiro para o vetor de termos independentes a ser preenchido.
*/
void initialize_rhs(const CSRMatrix *A, double *b) {
    for (int i = 0; i < A->num_rows; i++) {
        double row_sum = 0.0;
        for (int k = A->row_ptr[i]; k < A->row_ptr[i + 1]; k++) {
            row_sum += A->values[k];
        }
        b[i] = row_sum;
    }
}

/*
    * Função gauss_jacobi_seq: Executa o método iterativo de Gauss-Jacobi de forma sequencial.

    * Parâmetros:
        - A: Ponteiro para a matriz esparsa CSR.
        - b: Vetor de termos independentes.
        - x: Vetor solução inicial (atualizado com o resultado final ao término).
        - next: Vetor auxiliar para armazenar os valores da próxima iteração.
        - max_iter: Número máximo de iterações permitidas.
        - tol: Tolerância mínima para o critério de parada pelo resíduo máximo.
        - final_diff: Ponteiro de saída que recebe o resíduo final |x^(k+1) - x^(k)|.

    * Retorna:
        - Quantidade total de iterações executadas.
*/
int gauss_jacobi_seq(const CSRMatrix *A, const double *b, double *x, double *next,
                     int max_iter, double tol, double *final_diff) {
    int N = A->num_rows;
    double *cur_x    = x;
    double *cur_next = next;
    double max_diff  = 0.0;
    int iter         = 0;

    while (iter < max_iter) {
        max_diff = 0.0;

        for (int i = 0; i < N; i++) {
            double sum = 0.0;

            for (int k = A->row_ptr[i]; k < A->row_ptr[i + 1]; k++) {
                int j = A->col_idx[k];
                if (i != j) {
                    sum += A->values[k] * cur_x[j];
                }
            }

            cur_next[i] = (b[i] - sum) / A->diag[i];

            double diff = fabs(cur_next[i] - cur_x[i]);
            if (diff > max_diff) {
                max_diff = diff;
            }
        }

        /* Troca os ponteiros dos vetores sem copiar memória */
        double *temp = cur_x;
        cur_x = cur_next;
        cur_next = temp;
        iter++;

        if (max_diff < tol) {
            break;
        }
    }

    /* Garante que o resultado final fique armazenado no vetor original x */
    if (cur_x != x) {
        memcpy(x, cur_x, (size_t)N * sizeof(double));
    }

    *final_diff = max_diff;
    return iter;
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

int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Uso: %s <arquivo.mtx> [num_execucoes] [max_iter]\n", argv[0]);
        return 1;
    }

    int num_runs = (argc >= 3 && atoi(argv[2]) > 0) ? atoi(argv[2]) : DEFAULT_NUM_RUNS;
    int max_iter = (argc >= 4 && atoi(argv[3]) > 0) ? atoi(argv[3]) : DEFAULT_MAX_ITER;

    printf("=================================================================\n");
    printf("INICIALIZAÇÃO DO MÉTODO DE GAUSS-JACOBI (SEQUENCIAL - BENCHMARK)\n");
    printf("=================================================================\n");
    printf("Matriz alvo: %s | Repetições: %d | Teto de iterações: %d\n", argv[1], num_runs, max_iter);
    fflush(stdout);

    /*
        1. Carregamento e Validação da Matriz Esparsa
    */
    CSRMatrix A;
    if (load_matrix_market_csr(argv[1], &A) != 0) {
        return 1;
    }

    int N = A.num_rows;
    printf("Dimensão: %d x %d | Não-zeros (NNZ): %d\n\n", N, A.num_cols, A.nnz);
    fflush(stdout);

    if (N <= 0) {
        fprintf(stderr, "Erro: dimensão inválida da matriz (N = %d).\n", N);
        free_csr(&A);
        return 1;
    }

    for (int i = 0; i < N; i++) {
        if (fabs(A.diag[i]) < 1e-15) {
            fprintf(stderr, "Erro: Elemento zero na diagonal na linha %d.\n", i);
            free_csr(&A);
            return 1;
        }
    }

    /*
        2. Alocação e Inicialização dos Vetores
    */
    size_t n_size     = (size_t)N;
    double *b         = (double *)malloc(n_size * sizeof(double));
    double *x         = (double *)malloc(n_size * sizeof(double));
    double *next      = (double *)malloc(n_size * sizeof(double));
    double *run_times = (double *)malloc((size_t)num_runs * sizeof(double));

    if (!b || !x || !next || !run_times) {
        fprintf(stderr, "Erro de alocação de memória para os vetores.\n");
        free(b);
        free(x);
        free(next);
        free(run_times);
        free_csr(&A);
        return 1;
    }

    initialize_rhs(&A, b);

    /*
        3. Execução do Benchmark (Loop de Repetições)
    */
    int last_iter         = 0;
    double last_max_diff  = 0.0;
    double last_sol_error = 0.0;

    for (int run = 0; run < num_runs; run++) {
        memset(x, 0, n_size * sizeof(double));

        double t_start = get_time_sec();
        last_iter = gauss_jacobi_seq(&A, b, x, next, max_iter, TOL, &last_max_diff);
        double t_end = get_time_sec();

        run_times[run] = t_end - t_start;
        last_sol_error = compute_max_error(x, N);

        printf("  [Execução %2d/%d] Tempo: %.6f s | Iterações: %d | Resíduo: %.4e\n",
               run + 1, num_runs, run_times[run], last_iter, last_max_diff);
        fflush(stdout);
    }

    /*
        4. Cálculo das Estatísticas Consolidadas
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
        5. Exibição dos Resultados Consolidados
    */
    printf("\nRESULTADOS CONSOLIDADOS (SEQUENCIAL - %d EXECUÇÕES)\n", num_runs);
    printf("\nMatriz avaliada: %s\n", argv[1]);
    printf("Status da convergência: %s\n", (last_max_diff < TOL) ? "CONVERGIU" : "LIMITE DE ITERAÇÕES");
    printf("Dimensão da matriz (N): %d\n", N);
    printf("Elementos não-nulos: %d (Densidade: %.4f%%)\n", A.nnz, ((double)A.nnz / ((double)N * N)) * 100.0);
    printf("Tolerância estipulada: %.1e\n", TOL);
    printf("Iterações executadas: %d / %d\n", last_iter, max_iter);
    printf("Resíduo final (|x-x_old|): %.10e (%.12f)\n", last_max_diff, last_max_diff);
    printf("Erro absoluto máx (|x-1|): %.10e (%.12f)\n", last_sol_error, last_sol_error);
    printf("Tempo médio de execução: %.6f s (± %.6f s)\n", mean_time, std_time);
    printf("Tempo mínimo / máximo: %.6f s / %.6f s\n", min_time, max_time);
    printf("Tempo total acumulado (%dx): %.6f s\n", num_runs, sum_time);
    printf("Tempo médio por iteração: %.6f ms\n", mean_iter_ms);
    printf("Throughput médio estimado: %.2f MFlops/s\n\n", mean_mflops);
    fflush(stdout);

    /*
        6. Liberação de Recursos
    */
    free(run_times);
    free(b);
    free(x);
    free(next);
    free_csr(&A);

    return 0;
}