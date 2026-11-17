# Compiladores
CC      = gcc
MPICC   = mpicc
NVCC    = nvcc

# Flags de compilação (otimização padrão, sem vetorização agressiva ou unroll)
CFLAGS  = -O2 -Wall -Wextra -Iinclude
LDFLAGS = -lm

# Diretórios do projeto
BIN_DIR = bin
SRC_DIR = src

# Fontes comuns
COMMON_SRC = $(SRC_DIR)/common/csr_matrix.c $(SRC_DIR)/common/utils.c

# Versão Sequencial
SEQ_SRC = $(SRC_DIR)/seq/main_seq.c
SEQ_BIN = $(BIN_DIR)/gauss_jacobi_seq

.PHONY: all seq dirs clean

all: seq

dirs:
	@mkdir -p $(BIN_DIR) output errors results

seq: dirs
	$(CC) $(CFLAGS) $(SEQ_SRC) $(COMMON_SRC) -o $(SEQ_BIN) $(LDFLAGS)

clean:
	rm -rf $(BIN_DIR)