# jonashaslbeck@gmail.com; March 2019

# ----------------------------------------------------------------------------------
# ----------------------- 1) Take iteration from function call ---------------------
# ----------------------------------------------------------------------------------

#!/usr/bin/env Rscript
iter <- commandArgs(trailingOnly=TRUE)
iter <- as.numeric(iter)

# if (length(iter)==0){
#   iter <- 0
# }


# ----------------------------------------------------------------------------------
# ----------------------- 2) Load Packages, Source Aux Functions -------------------
# ----------------------------------------------------------------------------------

# ----- Source and Load -----

library(foreach)
library(doParallel)

# --- Laura's method ---
library(mgcv)
library(quantmod)
library(MASS)

source('VarSim_aux_Laura.R') # Estimation function GLM
source('VarSim_aux_tvvarGAM.R') # Estimation function GAM (instead of package, which caused trouble)

# --- Jonas ---
# library(devtools)
# install_github('jmbh/mgm')
library(mgm)


# ----------------------------------------------------------------------------------
# -------------------------------- 3) Load Data ------------------------------------
# ----------------------------------------------------------------------------------

l_data <- readRDS(file = paste0('output/VarSimData_Iter_', iter, '.RDS'))


# ----------------------------------------------------------------------------------
# --------------------------- 4) Estimate three Models -----------------------------
# ----------------------------------------------------------------------------------

nClust <- 12 # We run the n-variations in parallel

# Fire up foreach
cl <- makeCluster(nClust)
registerDoParallel(cl)

# start clustered loops
outlist <- foreach(n=1:12,
                   .packages=c("mgm", "quantmod", "mgcv", "MASS"),
                   .verbose=FALSE) %dopar% {
                     
                     out_sp <- vector('list', length = 3) # Storage
                     out_p <- vector('list', length = 3)
                     out_sp <- lapply(out_sp, function(x) out_p)
                     
                     
                     timer_tot <- proc.time()[3]
                     
                     for(sp in 3) { # loop over sparsity variations
                       
                       for(p in 2) { # loop over p-variations
                         
                         data_n <- l_data[[n]][[sp]][[p]]$Data
                         
                         nCases <- nrow(data_n)
                         nNodes <- ncol(data_n)
                         
                         l_n <- list('mvar' = NULL, # model objects
                                     'glmvar' = NULL,
                                     'tvmvar' = NULL,
                                     'tvmvar_unreg' = NULL,
                                     'gam' = NULL,
                                     'bwErrors' = NULL,
                                     'bwErrors_unreg' = NULL,
                                     'mvar.time' = NULL, # time object
                                     'glmvar.time' = NULL,
                                     'tvmvar.time' = NULL,
                                     'tvmvar.time_unreg' = NULL,
                                     'gam.time' = NULL,
                                     'bw.time' = NULL,
                                     'bw.time_unreg' = NULL)
                         
                         
                         # 1) ----- Estimate mVAR (lasso) -----
                         
                         timer <- proc.time()[3]
                         
                         l_n$mvar <- mvar(data = data_n,
                                          type = rep('g', nNodes),
                                          level = rep(1, nNodes), 
                                          lags = 1,
                                          lambdaSel = 'CV',
                                          pbar = FALSE,
                                          saveModels = FALSE)
                         
                         l_n$mvar.time <- proc.time()[3] - timer
                         
                         
                         # 2) ----- Estimate VAR (standard glm) -----
                         
                         timer <- proc.time()[3]
                         
                         l_n$glmvar <- Estimate_var_adj(data_n)
                         
                         l_n$glmvar.time <- proc.time()[3] - timer
                         
                         
                         
                         # 3) ----- Estimate tv-mVAR -----
                         
                         # Select Bandwidth using bwSelect
                         timer <- proc.time()[3]
                         bwSeq_full <- seq(0.01, .5, length = 15) # candidate bandwidth sequence
                         if(n <= 4) bwSeq <- bwSeq_full[6:15] # restrict candidate sequence depending on n
                         if(n > 4) bwSeq <- bwSeq_full[1:10]
                         
                         bwObj <- bwSelect(data = data_n,
                                           type = rep('g', nNodes),
                                           level = rep(1, nNodes),
                                           bwSeq = bwSeq,
                                           bwFolds = 1,
                                           bwFoldsize = round((nCases * .2)^(2/3)),
                                           modeltype = 'mvar', 
                                           lambdaSel = 'CV',
                                           lags = 1,
                                           saveModels = TRUE, 
                                           pbar = FALSE) # otherwise prediction is not possible
                         
                         # Delete model objects to save disc space
                         bwObj$bwModels <- NULL
                         
                         # Save error as function of bw
                         l_n$bwErrors <- bwObj$meanError
                         
                         # Time
                         l_n$bw.time <- proc.time()[3] - timer
                         
                         # Select bw that minimizes error
                         bw.min <- bwSeq[which.min(unlist(bwObj$meanError))]
                         
                         
                         # Use that bandwidth to fit model
                         timer <- proc.time()[3]
                         
                         l_n$tvmvar <- tvmvar(data = data_n,
                                              type = rep('g', nNodes),
                                              level = rep(1, nNodes), 
                                              lags = 1,
                                              lambdaSel = 'CV', 
                                              bandwidth = bw.min, 
                                              estpoints = seq(0, 1, length = 20), # estimation points
                                              pbar = FALSE)
                         
                         # Delete model objects to save disc space
                         l_n$tvmvar$tvmodels <- NULL
                         
                         l_n$tvmvar.time <- proc.time()[3] - timer
                         
                         
                         
                         # 4) ----- Estimate tv-mVAR [unregularized] -----
                         
                         # Select Bandwidth using bwSelect
                         timer <- proc.time()[3]
                         bwSeq_full <- seq(0.01, .5, length = 15) # candidate bandwidth sequence
                         if(n <= 4) bwSeq <- bwSeq_full[6:15] # restrict candidate sequence depending on n
                         if(n > 4) bwSeq <- bwSeq_full[1:10]
                         
                         bwObj_unreg <- bwSelect(data = data_n,
                                                 type = rep('g', nNodes),
                                                 level = rep(1, nNodes),
                                                 bwSeq = bwSeq,
                                                 bwFolds = 1,
                                                 bwFoldsize = round((nCases * .2)^(2/3)),
                                                 modeltype = 'mvar', 
                                                 lambdaSel = 'EBIC',
                                                 lambdaSeq = 0, # unregularized also for bwSelection
                                                 lags = 1,
                                                 saveModels = TRUE, 
                                                 pbar = FALSE, 
                                                 thresholds = "none") # otherwise prediction is not possible
                         
                         # Delete model objects to save disc space
                         bwObj_unreg$bwModels <- NULL
                         
                         # Save error as function of bw
                         l_n$bwErrors_unreg <- bwObj_unreg$meanError
                         
                         # Time
                         l_n$bw.time_unreg <- proc.time()[3] - timer
                         
                         # Select bw that minimizes error
                         bw.min <- bwSeq[which.min(unlist(bwObj_unreg$meanError))]
                         
                         
                         # Use that bandwidth to fit model
                         timer <- proc.time()[3]
                         
                         l_n$tvmvar_unreg <- tvmvar(data = data_n,
                                                    type = rep('g', nNodes),
                                                    level = rep(1, nNodes), 
                                                    lags = 1,
                                                    lambdaSel = 'EBIC',
                                                    lambdaSeq = 0,
                                                    bandwidth = bw.min, 
                                                    estpoints = seq(0, 1, length = 20), # estimation points
                                                    pbar = FALSE, 
                                                    thresholds = "none")
                         
                         # Delete model objects to save disc space
                         l_n$tvmvar_unreg$tvmodels <- NULL
                         
                         l_n$tvmvar.time_unreg <- proc.time()[3] - timer
                         
                         
                         
                         # 5) ----- Estimate GAM -----
                         
                         timer <- proc.time()[3]
                         
                         ## Compute number of Knots used
                         # As large as possible, but min=3, max=10
                         # we substract a small constant to avoid an equality in situations where no rounding is necessary, e.g. p=5, n=30 (new: July 2017)
                         nKnots <- max(3, min(floor(nCases/(nNodes+1) - 0.00001), 10))
                         
                         # Check whether GAM method can be used (only cases where not possible: n=25 for p=5 and n=25, 50 for p=10)
                         knot_check <- (nKnots * (nNodes+1)) < nCases
                         
                         if(knot_check) {
                           
                           
                           # # Old implementation from source file
                           # l_n$gam <- GAM_TVVAR(data_n,
                           #                      np = nNodes,
                           #                      N = nCases,
                           #                      K = nKnots)
                           
                           
                           # New implementation from tvvarGAM package
                           tvvarGAM_out <- tvvarGAM(data = data_n, 
                                                    nb = nKnots,
                                                    plot = FALSE, 
                                                    estimates = TRUE,
                                                    thresholding = FALSE,
                                                    pbar = FALSE)
                           
                           
                           # change export format so the old processing code still works
                           l_n$gam <- array(0, dim=c(nNodes+1, nNodes, nCases, 3))
                           l_n$gam[, , , 1] <- tvvarGAM_out$Results_GAM$CI_low
                           l_n$gam[, , , 2] <- tvvarGAM_out$Results_GAM$Estimate
                           l_n$gam[, , , 3] <- tvvarGAM_out$Results_GAM$CI_high
                           
                           
                         } else {
                           
                           l_n$gam <- array(0, dim=c(nNodes+1, nNodes, nCases, 3)) # same data structure to make evaluation easier
                           
                         }
                         
                         l_n$gam.time <- proc.time()[3] - timer
                         
                         
                         
                         # --- Save ---
                         
                         out_sp[[sp]][[p]] <- l_n
                         
                       } # end: for (p)
                       
                     } # end: for (sp)
                     
                     return(out_sp)
                     
                     total_time <- proc.time()[3] - timer_tot
                     
                     
                   } # end: foreach (n)


# Close down foreach
stopCluster(cl)


# ----------------------------------------------------------------------------------
# ---------------------------------- 5) Export -------------------------------------
# ----------------------------------------------------------------------------------

saveRDS(outlist, file = paste0('VarSimEst_Iter_', iter, '.RDS'))





