#!/bin/bash

#PBS -N VARsim_est_Oct17_1
#PBS -lnodes=1
#PBS -lwalltime=12:00:00

module load openmpi/gnu
module load R/3.2.3
module load c/intel
module load fortran/intel

# Remove blocking directory, if present
# rm -r /home/haslbeck/R/x86_64-pc-linux-gnu-library/3.2/00LOCK-mgm2


#export R_LIBS=$HOME/rpackages:$R_LIBS

cp -r "$HOME"/varSim "$TMPDIR"
cd "$TMPDIR"/varSim

Rscript --vanilla VarSim_est.R iter

cp -r ./*.RDS "$HOME"/varSim/output



