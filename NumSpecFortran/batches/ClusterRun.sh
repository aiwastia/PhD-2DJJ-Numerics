#!/bin/bash -l
#SBATCH -D /home/aiwastia/Desktop/git/PhD-2DJJ-Numerics/NumSpecFortran/batches
#SBATCH --partition=long
#SBATCH --job-name=1massSOC
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=300:00:00

#SBATCH --output make.out

./JJ_gap infile_JJ_bands_1 infile_JJ_2
