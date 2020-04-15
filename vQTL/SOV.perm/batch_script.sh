#!/bin/bash
#
#
#SBATCH -J SOVperm
#SBATCH -N 16
#SBATCH -n 16
#SBATCH --mail-user=azg5169@uncw.edu
#SBATCH --mail-type=all 
#SBATCH -p skx-normal
#SBATCH -t 14:00:00
#SBATCH -A TG-MCB180066
#SBATCH -o job_%j_%N.out
#------------------------------------------------------
mkdir -p output_perm
Rscript --verbose ./SOVperm.R > ./output_perm.Rout
