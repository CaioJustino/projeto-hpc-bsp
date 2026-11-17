#include "csr_matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



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
        1. Leitura do cabeçalho e identificação de simetria.
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
        2. Leitura temporária das coordenadas (formato COO).
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

        if (is_symmetric && r != c) {
            coo_r[actual_nnz] = c;
            coo_r[actual_nnz] = r;
            coo_v[actual_nnz] = v;
            actual_nnz++;
        }
    }
    fclose(f);

    /*
        3. Conversão de COO para CSR.
    */
    mat->num_rows = rows;
    mat->num_cols = cols;
    mat->nnz      = actual_nnz;
    mat->row_ptr  = (int *)calloc((size_t)(rows + 1), sizeof(int));
    mat->col_idx  = (int *)malloc((size_t)actual_nnz * sizeof(int));
    mat->values   = (double *)malloc((size_t)actual_nnz * sizeof(double));
    mat->diag     = (double *)calloc((size_t)rows, sizeof(double));

    if (!mat->row_ptr || !mat->col_idx || !mat->values || !mat->diag) {
        fprintf(stderr, "Erro de alocação para estrutura CSR.\n");
        free(coo_r);
        free(coo_c);
        free(coo_v);
        return -1;
    }

    for (int k = 0; k < actual_nnz; k++) {
        mat->row_ptr[coo_r[k] + 1]++;
    }

    for (int i = 0; i < rows; i++) {
        mat->row_ptr[i + 1] += mat->row_ptr[i];
    }

    int *current_pos = (int *)malloc((size_t)rows * sizeof(int));
    if (!current_pos) {
        fprintf(stderr, "Erro de alocação temporária em current_pos.\n");
        free(coo_r);
        free(coo_c);
        free(coo_v);
        return -1;
    }
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
    if (!mat) return;
    free(mat->row_ptr);
    free(mat->col_idx);
    free(mat->values);
    free(mat->diag);
    mat->row_ptr = NULL;
    mat->col_idx = NULL;
    mat->values = NULL;
    mat->diag = NULL;
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