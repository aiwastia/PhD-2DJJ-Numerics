#!/bin/bash -l
#SBATCH -D /home/aiwastia/Desktop/NumSpecFortran/2mass
#SBATCH --partition=long
#SBATCH --job-name=1massSOC
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=90:00:00

#SBATCH --output make.out

make
