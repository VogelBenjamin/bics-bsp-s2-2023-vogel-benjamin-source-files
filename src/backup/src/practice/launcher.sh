#!/bin/bash -l
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 4
#SBATCH --time=0-00:05:00
#SBATCH -p batch
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
./a.out > output.txt
