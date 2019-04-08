#!/bin/bash
#SBATCH -N 1
#SBATCH -t 07:00:00

module load openmpi/gnu
module load R/3.3.1
module load eb
module load intel/2016b
module load fortran/intel
module load mkl


#export R_LIBS=$HOME/rpackages:$R_LIBS

cp -r "$HOME"/VARsim_UT "$TMPDIR"
cd "$TMPDIR"/VARsim_UT

Rscript --vanilla VS_UT_est.R iter

cp -r ./*.RDS "$HOME"/VARsim_UT/output



