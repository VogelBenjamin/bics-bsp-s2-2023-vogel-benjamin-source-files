#!/bin/bash

# this script executes multiple different
# version of the particle simulation
# and stores the results in a file

mkdir outPut

gcc -fopenmp -o "pH" particles.c linearAlgebraImp.c electromagnetismImp.c testingImp.c -lm -g

./pH proton.txt 1 1 0 1 0 > outPut/sEOut.txt

./pH proton.txt 1 1 1 0 1 > outPut/sMOut.txt

./pH proton.txt 2 1 0 1 0 > outPut/mEOut.txt

./pH proton.txt 10 10 9 3 7 > outPut/mCOut.txt