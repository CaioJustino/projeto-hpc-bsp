#ifndef CSR_MATRIX_H
#define CSR_MATRIX_H

#ifdef __cplusplus
extern "C" {
#endif

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
    * Função load_matrix_market_csr: Lê um arquivo Matrix Market (.mtx) e converte para o formato CSR.

    * Parâmetros:
        - filename: Caminho do arquivo .mtx de entrada.
        - mat: Ponteiro para a estrutura CSRMatrix que será preenchida.

    * Retorna:
        - 0 em caso de sucesso, ou -1 se ocorrer erro de leitura ou alocação.
*/
int load_matrix_market_csr(const char *filename, CSRMatrix *mat);

/*
    * Função free_csr: Libera a memória alocada dinamicamente para a matriz CSR.

    * Parâmetros:
        - mat: Ponteiro para a estrutura CSRMatrix a ser liberada.
*/
void free_csr(CSRMatrix *mat);

/*
    * Função initialize_rhs: Constrói o vetor independente b somando os elementos de cada linha,
      garantindo que a solução analítica esperada seja o vetor unitário x* = [1, 1, ..., 1]^T.

    * Parâmetros:
        - A: Ponteiro para a matriz esparsa CSR.
        - b: Ponteiro para o vetor de termos independentes a ser preenchido.
*/
void initialize_rhs(const CSRMatrix *A, double *b);

#ifdef __cplusplus
}
#endif

#endif /* CSR_MATRIX_H */