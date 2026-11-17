#!/bin/bash
#SBATCH --job-name=gauss_jacobi_seq
#SBATCH --partition=amd-512
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --output=output/%j.out
#SBATCH --error=errors/%j.err
#SBATCH --chdir=/home/cvajustino/projeto-hpc-bsp

mkdir -p output errors results bin

make seq

# Parâmetros do Experimento
NUM_RUNS=20
MAX_ITER=5000

# Suíte de Matrizes
MATRIZES=(
    "data/bcsstk21.mtx"
    "data/memplus.mtx"
    "data/s3dkq4m2.mtx"
    "data/thermal2.mtx"
)

echo "EXECUÇÃO (SEQUENCIAL)"
echo ""
echo "- Data/Hora (Início): $(date)"

for MATRIZ in "${MATRIZES[@]}"; do
    if [ -f "$MATRIZ" ]; then
        ./bin/gauss_jacobi_seq "$MATRIZ" "$NUM_RUNS" "$MAX_ITER"
    else
        echo "AVISO: Arquivo $MATRIZ não encontrado!" >&2
    fi
done

echo "- Data/Hora (Final): $(date)"
