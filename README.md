# bics-bsp-s2-2023-vogel-benjamin-source-files
This repository contains all the content concerning the BSP-S2 of student Vogel Benjamin at uni.lu

## Compilation
gcc -o "pHA" particles.c electromagnetismImp.c linearAlgebraImp.c
testingImp.c -fopenmp -lm -g

## Execution
Recommended execution on hpc cluster using as many threads as there are
particles to simulate.

Basic structure of the command:
./pHA 'particle.txt' '\# of particles' 'simulation time' '3 init speed components' '3 electric field components' '3 magnetic field components' ('-info' 'name of file in which info is stored')

Example:
./pHA proton.txt 128 10 10 10 10 50 50 50 0 0 1 -info demonstartionFile

Batch script:
- 1 Node
- 128 cores
- 30 min execution time max

to be executed on aion supercomputer -- sbatch execP.sh