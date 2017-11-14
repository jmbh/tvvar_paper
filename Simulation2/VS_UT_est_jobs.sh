#!/bin/bash

#PBS -N VARsim_UT_Nov17
#PBS -lnodes=1
#PBS -lwalltime=5:00:00

module load openmpi/gnu
module load R/3.2.3
module load c/intel
module load fortran/intel


#export R_LIBS=$HOME/rpackages:$R_LIBS

cp -r "$HOME"/VARsim_UT "$TMPDIR"
cd "$TMPDIR"/VARsim_UT

Rscript --vanilla VS_UT_est.R iter

cp -r ./*.RDS "$HOME"/VARsim_UT/output



