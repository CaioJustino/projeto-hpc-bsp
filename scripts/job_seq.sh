#!/bin/bash
#SBATCH --job-name=gauss_jacobi_baseline
#SBATCH --partition=amd-512
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=output/baseline_%j.out
#SBATCH --error=errors/baseline_%j.err
#SBATCH --chdir=/home/cvajustino/projeto-hpc-bsp

mkdir -p output errors

# Compilação
gcc -O3 -march=native -funroll-loops src/gauss_jacobi_seq.c -o src/gauss_jacobi_seq -lm

# Execução
MATRIZ="data/thermal2.mtx"

./src/gauss_jacobi_seq "$MATRIZ"