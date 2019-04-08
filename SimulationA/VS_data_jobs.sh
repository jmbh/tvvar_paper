#!/bin/bash
#SBATCH -N 1
#SBATCH -t 00:30:00

module load openmpi/gnu
module load R/3.3.1
module load eb
module load intel/2016b
module load fortran/intel
module load mkl

#export R_LIBS=$HOME/rpackages:$R_LIBS

cp -r "$HOME"/tvvarSim2019 "$TMPDIR"
cd "$TMPDIR"/tvvarSim2019

Rscript --vanilla VarSim_data.R iter

cp -r ./*.RDS "$HOME"/tvvarSim2019/output



