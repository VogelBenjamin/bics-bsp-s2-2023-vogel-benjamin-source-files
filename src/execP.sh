#!/bin/bash
#SBATCH -N 1
#SBATCH --time=0-00:30:00
#SBATCH --exclusive


# this script executes multiple different
# version of the particle simulation
# and stores the results in a file

gcc -fopnmp -o "pH" particles.c linearAlgebraImp.c electromagnetismImp.c testingImp.c -lm -g

./pH proton.txt 5 5 0 0 0 1 1 1 0 0 0 -info protonElField

./pH naPlus.txt 5 5 0 0 0 1 1 1 0 0 0 -info sodiumElField

./pH cuPP.txt 5 5 0 0 0 1 1 1 0 0 0 -info copperElField

./pH proton.txt 5 5 1 0 0 0 0 0 1 1 1 -info protonMagField

./pH naPlus.txt 5 5 1 0 0 0 0 0 1 1 1 -info sodiumMagField

./pH cuPP.txt 5 5 0 1 0 0 0 0 1 1 1 -info copperMagField

./pH proton.txt 20 5 0 0 0 0 0 0 0 0 0 -info protonBulk

./pH naPlus.txt 20 5 0 0 0 0 0 0 0 0 0 -info sodiumBulk

./pH cuPP.txt 20 5 0 0 0 0 0 0 0 0 0 -info copperBulk

for i in 1 2 4 8 16 28
do

echo "Next Step" >> outPut/timeCompare.txt

echo "Small execution time:" >> outPut/timeCompare.txt

{ time ./pH proton.txt 10 10 9 3 7 ;} 2>> outPut/timeCompare.txt

echo "Larger execution time:" >> outPut/timeCompare.txt

{ time ./pH proton.txt 50 10 9 3 7 ;} 2>> outPut/timeCompare.txt
done