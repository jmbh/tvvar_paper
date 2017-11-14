# jonashaslbeck@gmail.com, November 2017


# ----------------------------------------------------------------------------------
# ----------------------- 1) Take iteration from command line call -----------------
# ----------------------------------------------------------------------------------

#!/usr/bin/env Rscript
iter <- commandArgs(trailingOnly=TRUE)
# if (length(iter)==0){
#   iter <- 0
# }

# iter <- 1


# ----------------------------------------------------------------------------------
# ----------------------- 2) Load Packages, Source Aux Functions -------------------
# ----------------------------------------------------------------------------------

# ----- Install Packages -----

# if(!require(glmnet)) install.packages('glmnet', repos="http://cran.rstudio.com/", dependencies = TRUE)
# if(!require(quantmod)) install.packages('quantmod', repos="http://cran.rstudio.com/", dependencies = TRUE)
# if(!require(stringr)) install.packages('stringr', repos="http://cran.rstudio.com/", dependencies = TRUE)
# if(!require(matrixcalc)) install.packages('matrixcalc', repos="http://cran.rstudio.com/", dependencies = TRUE)
# if(!require(foreach)) install.packages('foreach', repos="http://cran.rstudio.com/", dependencies = TRUE)
# if(!require(doParallel)) install.packages('doParallel', repos="http://cran.rstudio.com/", dependencies = TRUE)
# if(!require(mgcv)) install.packages('mgcv', repos="http://cran.rstudio.com/", dependencies = TRUE)
# if(!require(MASS)) install.packages('mvtnorm', repos="http://cran.rstudio.com/", dependencies = TRUE)

# ----- Source and Load -----

library(foreach)
library(doParallel)

# Laura
library(mgcv)
library(quantmod)
library(MASS)

source('VarSim_aux_Laura.R') # Estimation function GAM


# Jonas
library(devtools)
install_github('jmbh/mgm')
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
                     
                     for(sp in 1:3) {
                       
                       for(p in 1:3) {
                         
                         data_n <- l_data[[n]][[sp]][[p]]$Data
                         
                         nCases <- nrow(data_n)
                         nNodes <- ncol(data_n)
                         
                         l_n <- list('mvar' = NULL,
                                     'glmvar' = NULL,
                                     'tvmvar' = NULL,
                                     'gam' = NULL,
                                     'bwObj' = NULL,
                                     'mvar.time' = NULL,
                                     'glmvar.time' = NULL,
                                     'tvmvar.time' = NULL,
                                     'gam.time' = NULL,
                                     'bwObj.time' = NULL)
                         
                         
                         # a) ----- Estimate mVAR (lasso) -----
                         
                         timer <- proc.time()[3]
                         
                         l_n$mvar <- mvar(data = data_n,
                                          type = rep('g', nNodes),
                                          level = rep(1, nNodes), 
                                          lags = 1,
                                          lambdaSel = 'CV',
                                          pbar = FALSE,
                                          saveModels = FALSE)
                         
                         l_n$mvar.time <- proc.time()[3] - timer
                         
                         
                         # c) ----- Estimate VAR (standard glm) -----
                         
                         timer <- proc.time()[3]
                         
                         l_n$glmvar <- Estimate_var_adj(data_n)
                         
                         l_n$glmvar.time <- proc.time()[3] - timer
                         
                         
                         
                         # c) ----- Estimate tv-mVAR -----
                         
                         # Select Bandwidth using bwSelect
                         timer <- proc.time()[3]
                         bwSeq_full <- seq(0.01, .5, length = 15)
                         if(n <= 4) bwSeq <- bwSeq_full[6:15]
                         if(n > 4) bwSeq <- bwSeq_full[1:10]
                         
                         l_n$bwObj <- bwSelect(data = data_n,
                                               type = rep('g', nNodes),
                                               level = rep(1, nNodes),
                                               bwSeq = bwSeq,
                                               bwFolds = 1,
                                               bwFoldsize = round((nCases * .2)^(2/3)),
                                               modeltype = 'mvar', 
                                               lambdaSel = 'CV',
                                               lags = 1,
                                               pbar = T,
                                               saveModels = TRUE) # otherwise prediction is not possible
                         
                         # Delete model objects to save disc space
                         l_n$bwObj$bwModels <- NULL

                         l_n$bwObj.time <- proc.time()[3] - timer
                         
                         bw.min <- bwSeq[which.min(unlist(l_n$bwObj$meanError))]
                         
                         
                         # Use that bandwidth to fit model
                         timer <- proc.time()[3]
                         
                         l_n$tvmvar <- tvmvar(data = data_n,
                                              type = rep('g', nNodes),
                                              level = rep(1, nNodes), 
                                              lags = 1,
                                              lambdaSel = 'CV', 
                                              bandwidth = bw.min, 
                                              estpoints = seq(0, nCases-1, length = 20), # estimation points
                                              pbar = F)
                         
                         # Delete model objects to save disc space
                         l_n$tvmvar$tvmodels <- NULL
                         
                         l_n$tvmvar.time <- proc.time()[3] - timer
                         
                         
                         
                         # d) ----- Estimate GAM -----
                         
                         timer <- proc.time()[3]
                         
                         ## Compute number of Knots used
                         # As large as possible, but min=3, max=10
                         # we substract a small constant to avoid an equality in situations where no rounding is necessary, e.g. p=5, n=30 (new: July 2017)
                         nKnots <- max(3, min(floor(nCases/(nNodes+1) - 0.00001), 10))
                         
                         # Check whether GAM method can be used (only cases where not possible: n=25 for p=5 and n=25, 50 for p=10)
                         knot_check <- (nKnots * (nNodes+1)) < nCases

                         if(knot_check) {
                           
                           l_n$gam <- GAM_TVVAR(data_n, 
                                                np = nNodes, 
                                                N = nCases,  # Number of estimated time points; for now has to be equal to number of cases
                                                K = nKnots)
                           
                         } else {
                           
                           l_n$gam <- array(0, dim=c(nNodes+1, nNodes, nCases, 3)) # same data structure to make evaluation easier
                           
                         }
                         
                         l_n$gam.time <- proc.time()[3] - timer
                         
                         
                         
                         # --- Save ---
                         
                         out_sp[[sp]][[p]] <- l_n
                         
                       } # end: for (p)
                       
                     } # end: for (sp)
                     
                     return(out_sp)
                     
                   } # end: foreach (n)


# Close down foreach
stopCluster(cl)



# ----------------------------------------------------------------------------------
# ---------------------------------- 5) Export -------------------------------------
# ----------------------------------------------------------------------------------

saveRDS(outlist, file = paste0('VarSimEst_Iter_', iter, '.RDS'))





