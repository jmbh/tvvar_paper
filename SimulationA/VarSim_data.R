# jonashaslbeck@gmail.com; March 2019

# ----------------------------------------------------------------------------------
# ----------------------- 1) Take iteration from function call ---------------------
# ----------------------------------------------------------------------------------

#!/usr/bin/env Rscript
iter <- commandArgs(trailingOnly=TRUE)
# if (length(iter)==0){
#   iter <- 0
# }

iter <- as.numeric(iter)

# ----------------------------------------------------------------------------------
# ----------------------- 2) Load Packages, Source Aux Functions -------------------
# ----------------------------------------------------------------------------------

# foreach
library(foreach)
library(doParallel)

# Source Aux Functions
source('VarSim_aux.R')


# ----------------------------------------------------------------------------------
# --------------------------- 3) Sample Graphs and Data ----------------------------
# ----------------------------------------------------------------------------------

# a) Set Variations

v_p <- c(5, 10, 20)
v_sparse <- c(.05, .1, .2)

n_seq_log <- seq(3, 7.5, length = 12)
n_seq <- round(exp(n_seq_log))

theta <- .35 # selected so we get model matrices within constrain with reasonable amount of iterations
sigma <- .1

nClust <- 12 # running n variation in parallel uses cores most efficiently


# Fire up foreach
cl <- makeCluster(nClust)
registerDoParallel(cl)

# start clustered loops
outlist <- foreach(n=1:12,
                   .verbose=FALSE) %dopar% {
                     
                     # time <- proc.time()[3]
                     
                     # Storage
                     out_p <- out_sp <- vector('list', length = 3) # Storage
                     out_sp <- lapply(out_sp, function(x) out_p)
                     
                     for(sp in 1:3) { 
                       
                       for(p in 1:3) {
                         
                         print(p)
                         
                         l_p <- list('Graph' = NULL,
                                     'Data' = NULL)
                         
                         l_p$Graph <- GenerateGraph(p = v_p[p],
                                                    sp = v_sparse[sp],
                                                    theta = theta,
                                                    N = n_seq[n],
                                                    seed = iter)
                         
                         th <- matrix(0, nrow = n_seq[n], ncol = v_p[p]) # Define 'time-varying' Means as constant zero
                         
                         l_p$Data <- VARData(graph = l_p$Graph$G,
                                             th = th,
                                             sd = sigma,
                                             seed = iter)
                         
                         out_sp[[sp]][[p]] <- l_p
                         
                       }
                     }
                     
                     # proc.time()[3] - time
                     
                     return(out_sp)
                     
                   }


# Close down foreach
stopCluster(cl)


# ----------------------------------------------------------------------------------
# ---------------------------------- 4) Export -------------------------------------
# ----------------------------------------------------------------------------------

saveRDS(outlist, file = paste0('VarSimData_Iter_', iter, '.RDS'))




