#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_ITER 50000
#define TOL 1e-6



/*
- Estrutura CSRMatrix: Representa uma matriz esparsa no formato Compressed Sparse Row (CSR).
- Armazena ponteiros de linha, índices de coluna, valores não-nulos e acesso rápido à diagonal.
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
- Função get_time_sec: Obtém o tempo atual com precisão de nanossegundos.
- Utiliza CLOCK_MONOTONIC para evitar distorções por ajustes do relógio do sistema.

- Returns:
- Tempo decorrido em segundos como ponto flutuante de precisão dupla.
*/
static inline double get_time_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/*
- Função load_matrix_market_csr: Lê um arquivo no formato Matrix Market (.mtx) e constrói a estrutura CSR.
- Trata matrizes simétricas expandindo os pares fora da diagonal e extrai os coeficientes diagonais.

- Params:
- filename: caminho para o arquivo .mtx.
- mat: ponteiro para a estrutura CSRMatrix de destino.

- Returns:
- 0 em caso de sucesso, ou -1 em caso de erro de I/O ou formatação.
*/
int load_matrix_market_csr(const char *filename, CSRMatrix *mat) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Erro ao abrir o arquivo: %s\n", filename);
        return -1;
    }

    char line[1024];
    int is_symmetric = 0;

    if (fgets(line, sizeof(line), f) == NULL) {
        fclose(f);
        return -1;
    }
    if (strstr(line, "symmetric") != NULL) {
        is_symmetric = 1;
    }

    while (fgets(line, sizeof(line), f)) {
        if (line[0] != '%') break;
    }

    int rows, cols, entries;
    if (sscanf(line, "%d %d %d", &rows, &cols, &entries) != 3) {
        fprintf(stderr, "Erro ao processar as dimensoes da matriz.\n");
        fclose(f);
        return -1;
    }
    
    int max_entries = is_symmetric ? entries * 2 : entries;
    int *coo_r = (int *)malloc(max_entries * sizeof(int));
    int *coo_c = (int *)malloc(max_entries * sizeof(int));
    double *coo_v = (double *)malloc(max_entries * sizeof(double));

    int actual_nnz = 0;
    for (int k = 0; k < entries; k++) {
        int r, c;
        double v;
        if (fscanf(f, "%d %d %lf", &r, &c, &v) != 3) break;
        r--; c--;

        coo_r[actual_nnz] = r;
        coo_c[actual_nnz] = c;
        coo_v[actual_nnz] = v;
        actual_nnz++;

        if (is_symmetric && r != c) {
            coo_r[actual_nnz] = c;
            coo_c[actual_nnz] = r;
            coo_v[actual_nnz] = v;
            actual_nnz++;
        }
    }
    fclose(f);

    mat->num_rows = rows;
    mat->num_cols = cols;
    mat->nnz = actual_nnz;
    mat->row_ptr = (int *)calloc(rows + 1, sizeof(int));
    mat->col_idx = (int *)malloc(actual_nnz * sizeof(int));
    mat->values  = (double *)malloc(actual_nnz * sizeof(double));
    mat->diag    = (double *)malloc(rows * sizeof(double));

    /* Passo 1: Contagem do número de elementos não-nulos por linha */
    for (int k = 0; k < actual_nnz; k++) {
        mat->row_ptr[coo_r[k] + 1]++;
    }

    /* Passo 2: Soma cumulativa de prefixos para definir os deslocamentos das linhas */
    for (int i = 0; i < rows; i++) {
        mat->row_ptr[i + 1] += mat->row_ptr[i];
    }

    /* Passo 3: Distribuição dos elementos nas posições definitivas dos vetores CSR */
    int *current_pos = (int *)malloc(rows * sizeof(int));
    memcpy(current_pos, mat->row_ptr, rows * sizeof(int));

    for (int k = 0; k < actual_nnz; k++) {
        int r = coo_r[k];
        int dest = current_pos[r]++;
        mat->col_idx[dest] = coo_c[k];
        mat->values[dest] = coo_v[k];
        if (coo_r[k] == coo_c[k]) {
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
- Função free_csr: Desaloca toda a memória dinâmica associada a uma estrutura CSRMatrix.

- Params:
- mat: ponteiro para a estrutura CSRMatrix a ser liberada.
*/
void free_csr(CSRMatrix *mat) {
    if (mat->row_ptr) free(mat->row_ptr);
    if (mat->col_idx) free(mat->col_idx);
    if (mat->values)  free(mat->values);
    if (mat->diag)    free(mat->diag);
}



int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Uso: %s <arquivo.mtx>\n", argv[0]);
        return 1;
    }

    printf("INICIALIZAÇÃO DO MÉTODO DE GAUSS-JACOBI (SEQUENCIAL)\n");

    printf("\nCarregando matriz: %s\n", argv[1]);
    CSRMatrix A;
    if (load_matrix_market_csr(argv[1], &A) != 0) {
        return 1;
    }

    int N = A.num_rows;
    printf("Dimensão: %d x %d | Não-zeros (NNZ): %d\n", N, A.num_cols, A.nnz);

    if (N <= 0) {
        fprintf(stderr, "Erro: dimensão inválida da matriz (N = %d).\n", N);
        free_csr(&A);
        return 1;
    }

    for (int i = 0; i < N; i++) {
        if (fabs(A.diag[i]) < 1e-15) {
            fprintf(stderr, "Erro: Elemento zero encontrado na diagonal na linha %d.\n", i);
            free_csr(&A);
            return 1;
        }
    }

    /*
        1. Alocação dos Vetores de Estado e Termos Independentes
    */
    size_t n_size = (size_t)N;
    double *b    = (double *)malloc(n_size * sizeof(double));
    double *x    = (double *)calloc(n_size, sizeof(double));
    double *next = (double *)malloc(n_size * sizeof(double));
    
    if (!b || !x || !next) {
        fprintf(stderr, "Erro de alocação de memória para vetores.\n");
        free_csr(&A);
        return 1;
    }

    for (int i = 0; i < N; i++) {
        double row_sum = 0.0;
        for (int k = A.row_ptr[i]; k < A.row_ptr[i + 1]; k++) {
            row_sum += A.values[k];
        }
        b[i] = row_sum;
    }

    /*
        2. Configuração de Cronometragem
    */
    double t_start = get_time_sec();

    int iter;
    double max_diff = 0.0;

    /*
        3. Loop Iterativo (Modelo BSP Sequencial)
    */
    for (iter = 0; iter < MAX_ITER; iter++) {
        max_diff = 0.0;

        for (int i = 0; i < N; i++) {
            double sum = 0.0;
            int start = A.row_ptr[i];
            int end   = A.row_ptr[i + 1];

            for (int k = start; k < end; k++) {
                int j = A.col_idx[k];
                if (i != j) {
                    sum += A.values[k] * x[j];
                }
            }

            next[i] = (b[i] - sum) / A.diag[i];

            double diff = fabs(next[i] - x[i]);
            if (diff > max_diff) {
                max_diff = diff;
            }
        }

        memcpy(x, next, N * sizeof(double));

        if (max_diff < TOL) {
            iter++;
            break;
        }
    }

    double t_end = get_time_sec();
    double elapsed_time = t_end - t_start;

    double max_sol_error = 0.0;
    for (int i = 0; i < N; i++) {
        double err = fabs(x[i] - 1.0);
        if (err > max_sol_error) {
            max_sol_error = err;
        }
    }

    /*
        4. Exibição dos Resultados e Métricas Detalhadas
    */
    printf("\nRESULTADOS (SEQUENCIAL)\n");
    printf("\nStatus da convergência: %s\n", (max_diff < TOL) ? "CONVERGIU" : "LIMITE DE ITERAÇÕES");
    printf("Dimensão da matriz (N): %d\n", N);
    printf("Elementos não-nulos: %d (Densidade: %.4f%%)\n", A.nnz, ((double)A.nnz / ((double)N * N)) * 100.0);
    printf("Tolerância estipulada: %.1e\n", TOL);
    printf("Iterações executadas: %d / %d\n", iter, MAX_ITER);
    printf("Resíduo final (|x-x_old|): %.10e (%.12f)\n", max_diff, max_diff);
    printf("Erro absoluto máx (|x-1|): %.10e (%.12f)\n", max_sol_error, max_sol_error);
    printf("Tempo total de execução: %.6f s\n", elapsed_time);
    if (elapsed_time > 0.0) {
        printf("Tempo médio por iteração: %.6f ms\n", (elapsed_time / iter) * 1000.0);
        printf("Throughput estimado: %.2f MFlops/s\n", ((double)iter * (2.0 * A.nnz + 2.0 * N)) / (elapsed_time * 1e6));
    }

    free(b);
    free(x);
    free(next);
    free_csr(&A);

    return 0;
}