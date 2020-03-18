#!/bin/bash
#
#
#SBATCH -J minday
#SBATCH -N 4
#SBATCH -n 4
#SBATCH --mail-user=azg5169@uncw.edu
#SBATCH --mail-type=all 
#SBATCH -p skx-normal
#SBATCH -t 4:00:00
#SBATCH -A TG-MCB180066
#SBATCH -o job_%j_%N.out
#------------------------------------------------------
mkdir -p output_perm
Rscript --verbose ./SOVperm.R > ./output_perm.Rout
