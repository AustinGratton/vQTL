#!/bin/bash
#
#
#SBATCH -J SOVinit
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user=azg5169@uncw.edu
#SBATCH --mail-type=all 
#SBATCH -p skx-normal
#SBATCH -t 14:00:00
#SBATCH -A TG-MCB180066
#SBATCH -o job_%j_%N.out
#------------------------------------------------------
mkdir -p output_init
Rscript --verbose ./SOVinit.R > ./output_init.Rout
