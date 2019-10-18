#!/bin/bash -l
#SBATCH -D /home/aiwastia/Desktop/git/PhD-2DJJ-Numerics/NumSpecFortran/newMesh
#SBATCH --partition=long
#SBATCH --job-name=diffmass
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=90:00:00

#SBATCH --output make.out

./JJ_bands infile_JJ_bands_1 infile_JJ_2
