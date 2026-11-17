#!/bin/bash
#SBATCH --job-name=gauss_jacobi_seq
#SBATCH --partition=amd-512
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --output=output/%j.out
#SBATCH --error=errors/%j.err
#SBATCH --chdir=/home/cvajustino/projeto-hpc-bsp

mkdir -p output errors results bin

make seq

# Parâmetros Globais do Experimento
NUM_RUNS=10
TOL=0.0  # Carga fixa para benchmark sustentado (executa todas as iterações)

# Suíte de Matrizes e Carga de Iterações (Calibrado para durar ~1 hora no total com 10 repetições)
# - bcsstk21: 50.000 iterações (~25s total)
# - memplus:  100.000 iterações (~2,5min total)
# - s3dkq4m2: 20.000 iterações (~20min total)
# - thermal2: 10.000 iterações (~40min total)
# Tempo total estimado: ~63 minutos
EXPERIMENTOS=(
    "data/bcsstk21.mtx 50000"
    "data/memplus.mtx 100000"
    "data/s3dkq4m2.mtx 20000"
    "data/thermal2.mtx 10000"
)

echo "EXECUÇÃO (SEQUENCIAL - CARGA SUSTENTADA ~1H)"
echo ""
echo "- Data/Hora (Início): $(date)"

for EXP in "${EXPERIMENTOS[@]}"; do
    read -r MATRIZ MAX_ITER <<< "$EXP"
    if [ -f "$MATRIZ" ]; then
        ./bin/gauss_jacobi_seq "$MATRIZ" "$NUM_RUNS" "$MAX_ITER" "$TOL"
    else
        echo "AVISO: Arquivo $MATRIZ não encontrado!" >&2
    fi
done

echo "- Data/Hora (Final): $(date)"
