#!/bin/bash

#PBS -N VARsim_data_Oct17_1
#PBS -lnodes=1
#PBS -lwalltime=00:40:00

module load openmpi/gnu
module load R/3.2.3
module load c/intel
module load fortran/intel

#export R_LIBS=$HOME/rpackages:$R_LIBS

cp -r "$HOME"/varSim "$TMPDIR"
cd "$TMPDIR"/varSim

Rscript --vanilla VarSim_data.R iter

cp -r ./*.RDS "$HOME"/varSim/output



