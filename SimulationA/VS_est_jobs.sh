#!/bin/bash
#SBATCH -N 1
#SBATCH -t 6:00:00

module load openmpi/gnu
module load R/3.3.1
module load eb
module load intel/2016b
module load fortran/intel
module load mkl

# Remove blocking directory, if present
# rm -r /home/haslbeck/R/x86_64-pc-linux-gnu-library/3.2/00LOCK-mgm2
#export R_LIBS=$HOME/rpackages:$R_LIBS

cp -r "$HOME"/tvvarSim2019 "$TMPDIR"
cd "$TMPDIR"/tvvarSim2019

Rscript --vanilla VarSim_est.R iter

cp -r ./*.RDS "$HOME"/tvvarSim2019/output



