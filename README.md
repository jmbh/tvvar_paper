# tvvar_paper

Code Archive to reproduce all results in the paper https://arxiv.org/abs/1711.05204.

## Simulation 1: Random Graph

The simulation is split in two parts. In a first step the data is generated and in a second step the original models are recovered with the 5 methods. This results in 100 .RDS-files (one for each iteration) and 100 .RDS-files containing the estimates of the 5 methods.

These R-files allow to run a single iteration:

- VarSim_data.R generates data for one iteration (outputs a .RDS-data file)
- VarSim_est.R takes the RDS.-data file from VarSim_data.R as input and estimates the model with all 5 methods on the data; outputs a .RDS-file containing a list with all estimates
- VarSim_aux.R Functions used in VarSim_data.R to generate time-varying VAR models and sample observations from them
- VarSim_aux_Laura.R Functions used in VarSim_est.R to estimate VAR models with the GLM method
- VarSim_aux_tvvarGAM.R Functions to used in VARsim_est.R to estimate time-varying VAR models with the spline-smoothing method

The batch-files were used to run the 100 iterations in parallel on the UvA LISA supercomputer. All these scripts do is to call both VarSim_data.R and VarSim_est.R with a different iteration. The iteration number is used as a random seed and hence any setup (even locally) that does that replicates our simulation results exactly.

- VS_data_jobs.sh batch script sent to each node to generate data
- VS_data_sunmit.sh batch script to submit VS_data_jobs.sh to 100 nodes with iteration (seeds) 1:100
- VS_est_jobs.sh batch script sent to each node to estimate models
- VS_est_sunmit.sh batch script to submit VS_est_jobs.sh to 100 nodes with iteration (seeds) 1:100

The size of the 2*100 output files is 21GB and is therefore not included in this archive.

The R-file VarSim_eval.R contains the code to reproduce all figures reported in the paper from the 100 data files and 100 estimation files obtained from the above simulation files.

## Simulation 2: UT-Model

Similarly to Simulation 1, the simulation is split in data generation and estimation. In Simulation 2 much less data is generated and we therefore we generated data with the file VS_UT_data.R locally. The output of this R-file are 100 .RDS-data files similarly to Simulation 1.

The estimation is done in parallel in the same way as in Simulation 1:

- VS_UT_est.R takes the 100 .RDS-data files as input, estimates the models and outputs 100 .RDS output files
- VS_UT_est_jobs.sh is the batch script sent to each of the 100 nodes
- VS_UT_est_submit.sh is the batch script used to send VS_UT_est_jobs.sh to 100 nodes with iteration (seeds) 1:100
- VS_UT_aux.R contains functions called by VS_UT_data.R
- VarSim_aux_Laura.R is the same file as described in Simulation 1
- VarSim_aux_tvvarGAM.R is the same file as described in Simulation 1

The size of the 2*100 output files is 5.3GB and is therefore not included in this archive.

The code in VS_UT_eval.R preprocesses the 100 output files and reproduces Figure 10 in the paper.

## Tutorial

- tutorial_mgm.R contains the code (shown in the paper) to estimate a time-varying VAR model with the KS(L1) method, plus the code to reproduce Figure 11
- tutorial_tvvarGAM.R contains the code (shown in the paper) to estimate a time-varying VAR model with the GAM(st) method, plus the code to reproduce Figure 12


## Illustration Figures

This folder contains R-code to reproduce all illustration figures except Figure 4. Figure 4 requires input from data files and is therefore contained in the R-file VarSim_eval.R (Simulation 1 folder).

- IF_Figure1_Intro.R code for Figure 1
- IF_Figure2_GAM.R code for Figure 2
- IF_Figure3_KS.R code for Figure 3
- IF_Figure9_Simulation2.R code for Figure 9
