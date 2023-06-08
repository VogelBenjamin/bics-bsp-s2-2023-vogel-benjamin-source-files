#!/bin/bash
#SBATCH -N 1
#SBATCH -c 128
#SBATCH --time=0-01:00:00



# this script executes multiple different
# version of the particle simulation
# and stores the results in a file

module purge

module load toolchain/foss/2020b

module load compiler/AOCC/3.1.0-GCCcore-10.2.0

gcc -fopenmp -o "pH" particles.c linearAlgebraImp.c electromagnetismImp.c testingImp.c -lm -g
export OMP_NUM_THREADS=128
./pHA cuPP.txt 128 10 10 10 10 50 50 50 0 0 10 -info 128visual
#./pH proton.txt 5 5 0 0 0 1 1 1 0 0 0 -info protonElField

#./pH naPlus.txt 5 5 0 0 0 1 1 1 0 0 0 -info sodiumElField

#./pH cuPP.txt 5 5 0 0 0 1 1 1 0 0 0 -info copperElField

#./pH proton.txt 10 3 1 0 0 0 0 0 1 1 1 -info protonMagField

#./pH naPlus.txt 10 3 1 0 0 0 0 0 1 1 1 -info sodiumMagFieldAccurate



#./pH cuPP.txt 5 5 0 1 0 0 0 0 1 1 1 -info copperMagField

#./pH proton.txt 20 5 0 0 0 0 0 0 0 0 0 -info protonBulk

#./pH naPlus.txt 20 5 0 0 0 0 0 0 0 0 0 -info sodiumBulk

#./pH cuPP.txt 20 5 0 0 0 0 0 0 0 0 0 -info copperBulk

#echo "Single Thread" >> outPut/timeCompare.txt

#echo "Small execution time:" >> outPut/timeCompare.txt

#{ time ./pH proton.txt 28 10 0 0 0 3 2 4 7 4 2 ;} 2>> outPut/timeCompare.txt

#echo "Larger execution time:" >> outPut/timeCompare.txt

#{ time ./pH proton.txt 16 10 0 0 0 3 2 4 7 4 2 ;} 2>> outPut/timeCompare.txt

#icc -qopenmp -o "pH" particles.c linearAlgebraImp.c electromagnetismImp.c testingImp.c -lm -g

#export OMP_NUMTHREADS=20

#vtune -collect hotspots ./pH proton.txt 20 1 0 0 0 0 0 0 0 0 0

#for I in 1,2,4,8,12,16,24
#do
#export OMP_NUM_THREADS=$I
#{ time ./pH proton.txt 28 10 0 0 0 3 2 4 7 4 2 ;} 2>> outPut/timeCompare.txt
#done
#echo "Multiple Threads" >> outPut/timeCompare.txt

#echo "Small execution time:" >> outPut/timeCompare.txt

#{ time ./pH proton.txt 20 10 0 0 0 3 2 4 7 4 2 ;} 2>> outPut/timeCompare.txt

#export OMP_NUM_THREADS=4

#{ time ./pH proton.txt 20 10 0 0 0 3 2 4 7 4 2 ;} 2>> outPut/timeCompare.txt

#export OMP_NUM_THREADS=28

#echo "Larger execution time:" >> outPut/timeCompare.txt

#{ time ./pH proton.txt  10 0 0 0 3 2 4 7 4 2 ;} 2>> outPut/timeCompare.txt
