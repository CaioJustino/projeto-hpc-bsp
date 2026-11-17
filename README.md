# Paralelização do Método de Gauss-Jacobi (BSP - HPC)

O presente projeto designa a implementação modular e análise de desempenho do método iterativo de **Gauss-Jacobi** para sistemas lineares esparsos de grande porte ($Ax = b$), utilizando o modelo de computação paralela **BSP (*Bulk Synchronous Parallel*)** no formato **CSR (*Compressed Sparse Row*)**.

Projeto desenvolvido para a disciplina de **Computação de Alto Desempenho (CAD)** com execução e avaliação experimental no supercomputador do **NPAD/UFRN**.

* **Download dos Datasets (Matrizes):** [Google Drive - Pasta de Matrizes (.mtx)](https://drive.google.com/drive/folders/1scfvr7DDkiVlgQ6_RSmYkfFa606cfFsS?usp=drive_link)

---

## Sumário

- [Visão Geral do Problema](#visão-geral-do-problema)
- [Estrutura do Projeto](#estrutura-do-projeto)
- [Conjunto de Dados (Suíte de Matrizes)](#conjunto-de-dados-suíte-de-matrizes)
- [Compilação](#compilação)
- [Execução no NPAD (SLURM)](#execução-no-npad-slurm)

---

## Visão Geral do Problema

O método de **Gauss-Jacobi** calcula a aproximação de cada incógnita $x_i^{(k+1)}$ exclusivamente a partir dos valores da iteração imediatamente anterior $x^{(k)}$:

$$x_i^{(k+1)} = \frac{1}{A_{ii}} \left( b_i - \sum_{j \neq i} A_{ij} x_j^{(k)} \right)$$

Essa formulação confere **elevado paralelismo intrínseco**, pois o cálculo de cada linha é independente durante o superpasso de computação, sem condições de corrida.

### Características da Implementação:

* **Formato CSR (*Compressed Sparse Row*):** Armazena apenas os valores não-nulos (`values`, `col_idx`, `row_ptr`) e um vetor direto da diagonal principal (`diag`), reduzindo drasticamente o uso de memória em relação a matrizes densas.

* **Validação Analítica Exata ($x^* = [1, 1, \dots, 1]^T$):** O vetor independente $b$ é gerado pela soma das linhas de $A$ ($b_i = \sum_j A_{ij}$), permitindo calcular em tempo real o erro analítico real $\|x - 1\|_\infty$ e verificar a corretude numérica de cada versão.

* **Critério de Parada:** Parada por tolerância de resíduo máximo $\|x^{(k+1)} - x^{(k)}\|_\infty < 10^{-6}$ (`TOL`) ou limite estipulado de iterações (`MAX_ITER`).

---

## Estrutura do Projeto

```text
projeto-hpc-bsp/
├── Makefile
├── README.md
├── .gitignore
├── data/
│   ├── bcsstk21.mtx
│   ├── memplus.mtx
│   ├── s3dkq4m2.mtx
│   └── thermal2.mtx
├── include/
│   ├── csr_matrix.h
│   └── utils.h
├── src/
│   ├── common/
│   │   ├── csr_matrix.c
│   │   └── utils.c
│   ├── seq/
│   │   └── main_seq.c
│   ├── mpi_omp/
│   ├── omp_gpu/
│   └── cuda/
├── bin/
├── scripts/
│   ├── job_seq.sh
│   ├── job_mpi_omp.sh
│   ├── job_omp_gpu.sh
│   └── job_cuda.sh
├── output/
├── errors/
└── results/
```

---

## Conjunto de Dados (Suíte de Matrizes)

As matrizes utilizadas no projeto provêm dos repositórios *NIST Matrix Market* e *SuiteSparse*:

| Matriz | Dimensão ($N$) | Elementos Não-Nulos ($NNZ$) | Densidade | Aplicação |
| :--- | :---: | :---: | :---: | :--- |
| `bcsstk21.mtx` | 3.600 | 26.600 | 0,205% | Rigidez estrutural (validação rápida) |
| `memplus.mtx` | 17.758 | 126.150 | 0,040% | Simulação de circuito de memória |
| `s3dkq4m2.mtx` | 90.449 | 4.427.725 | 0,054% | Análise de elementos finitos (cilindro) |
| `thermal2.mtx` | 1.228.045 | 8.580.313 | 0,0005% | Problema térmico não-estruturado em grande escala |

### Download das Matrizes:
1. Baixe os arquivos `.mtx` diretamente da pasta do Google Drive:
   👉 **[Google Drive - Matrizes do Projeto](https://drive.google.com/drive/folders/1scfvr7DDkiVlgQ6_RSmYkfFa606cfFsS?usp=drive_link)**
2. Mova ou descompacte os arquivos dentro do diretório `data/` na raiz do projeto:
   ```bash
   mkdir -p data
   ```

---

## Compilação

O projeto conta com compilação centralizada via `Makefile`.

```bash
# Compilar versão sequencial (baseline)
make seq

# Compilar versão híbrida (MPI + OpenMP)
make mpi_omp

# Compilar versão OpenMP Offloading (GPU)
make omp_gpu

# Compilar versão CUDA Nativo (GPU)
make cuda

# Compilar todas as versões disponíveis
make all

# Limpar executáveis gerados
make clean
```

---

## Execução no NPAD (SLURM)

A execução dos experimentos no supercomputador é gerenciada pelo **SLURM**. Cada versão possui seu respectivo script de submissão dentro do diretório `scripts/`, configurado com a partição e os recursos adequados para seu hardware-alvo (CPU ou nós com aceleradores GPU).

### 1. Submissão dos Jobs

* **Versão Sequencial (Baseline):**
  ```bash
  sbatch scripts/job_seq.sh
  ```

* **Versão Híbrida (MPI + OpenMP):**
  ```bash
  sbatch scripts/job_mpi_omp.sh
  ```

* **Versão OpenMP Offloading (GPU):**
  ```bash
  sbatch scripts/job_omp_gpu.sh
  ```

* **Versão CUDA Nativo (GPU):**
  ```bash
  sbatch scripts/job_cuda.sh
  ```

---

### 2. Gerenciamento e Monitoramento da Fila

```bash
# Verificar status dos jobs em execução na fila
squeue -u $USER

# Cancelar um job em execução
scancel <JOB_ID>

# Acompanhar a saída de um job em tempo real
tail -f output/<JOB_ID>.out
```

---

**Discentes:** Caio Justino, João Lucas e Liriel Felix.