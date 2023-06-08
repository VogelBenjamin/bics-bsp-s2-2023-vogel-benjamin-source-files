#!/bin/bash -l
#SBATCH --ntasks-per-node=1
#SBATCH -c 2
#SBATCH --time=0-00:05:00
#SBATCH -p batch

export OMP_NUM_THREADS=1

echo "Take 1" >> timeOutPut.txt

{ time ./pH proton.txt 10 1 0 1 0 ;} 2>> timeOutPut.txt

export OMP_NUM_THREADS=2

echo "Take 2" >> timeOutPut.txt

{ time ./pH proton.txt 10 1 0 1 0 ;} 2>> timeOutPut.txt


