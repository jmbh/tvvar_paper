
# ----------------------------------------------------------------------------------
# ----------------------- 1) Take iteration from function call ---------------------
# ----------------------------------------------------------------------------------

#!/usr/bin/env Rscript
iter <- commandArgs(trailingOnly=TRUE)
# if (length(iter)==0){
#   iter <- 0
# }

# iter <- 1

for(iter in 1:100) {


# ----------------------------------------------------------------------------------
# ----------------------- 2) Load Packages, Source Aux Functions -------------------
# ----------------------------------------------------------------------------------

# foreach
library(foreach)
library(doParallel)

# Source Aux Functions
source('VS_UT_aux.R')


# ----------------------------------------------------------------------------------
# --------------------------- 3) Sample Graphs and Data ----------------------------
# ----------------------------------------------------------------------------------

# n-variations
n_seq_log <- seq(3, 7.5, length = 12)
n_seq <- round(exp(n_seq_log))

# Fixed
theta <- .35 # selected so we get model matrices within constrain with reasonable amount of iterations
sigma <- .1
p <- 20

nClust <- 12 # running n variation in parallel uses cores most efficiently


# Fire up foreach
cl <- makeCluster(nClust)
registerDoParallel(cl)

# start clustered loops
outlist <- foreach(n=1:12,
                   .verbose=FALSE) %dopar% {
                     
                         l_p <- list('Graph' = NULL,
                                     'Data' = NULL)
                         
                         l_p$Graph  <- GenerateGraph_fix(theta = .35,
                                                  N = n_seq[n],
                                                  p = 20,
                                                  seed = iter,
                                                  max_check=100)

                         th <- matrix(0, 
                                      nrow = n_seq[n], 
                                      ncol = p) # Define 'time-varying' Means as constant zero
                         
                         l_p$Data <- VARData(graph = l_p$Graph$G,
                                             th = th,
                                             sd = sigma,
                                             seed = iter)

                     
                     # proc.time()[3] - time
                     
                     return(l_p)
                     
                   }


# Close down foreach
stopCluster(cl)


# ----------------------------------------------------------------------------------
# ---------------------------------- 4) Export -------------------------------------
# ----------------------------------------------------------------------------------

saveRDS(outlist, file = paste0('output/VarSimData_Iter_', iter, '.RDS'))

print(iter)

} # end for: iter

