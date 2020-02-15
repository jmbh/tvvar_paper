# jonashaslbeck@gmail.com; March 2019

# Define Dirs
simDir <- '...' #  SPECIFY MAIN DIRECTORY
figDir <- paste0(simDir, 'figures/')
figDir_paper <- paste0(simDir, 'figures_paper/')
dir_files <- paste0(simDir, "output/")


# DEVEL:
figDir_paper <- paste0(simDir, "figures/")

# ----------------------------------------------------------------------------------
# --------------------------------- 1) Load Data -----------------------------------
# ----------------------------------------------------------------------------------

v_files <- list.files(dir_files)
v_files <- v_files[order(v_files)] # make sure that order in vector corresponds to iter (seed)

# Load Data
v_files_data <- v_files[grepl('Data', v_files)]
n_files <- length(v_files_data)
l_data <- list()
for(i in 1:n_files) {
  l_data[[i]] <- readRDS(paste0(dir_files, v_files_data[i]))
  print(i)  
}

# Load Estimates
v_files_est <- v_files[grepl('Est', v_files)]
n_files <- length(v_files_est)
l_est <- list()
for(i in 1:n_files) {
  l_est[[i]] <- readRDS(paste0(dir_files, v_files_est[i]))
  print(i)  
}


# ----------------------------------------------------------------------------------
# --------------------------------- 2) Preprocess ----------------------------------
# ----------------------------------------------------------------------------------

## Data Structure Strategy:
# For each condition (n x p x sp) and estimation method we have
# a matrix with columns: iteration, row, col, type, est1, est2, ..., est20
# - row/col indicate the edge in the parameter matrix
# - type is the type of the edge (0, 1, 2, ..., 8)

n <- 12
sp <- 3 # P(e) = 0.2
p <- 2 # p = 10
i <- 1

# example true type graph:
X <- l_data[[i]][[n]][[sp]][[p]]$Graph$Gind
sum(X!=0) - p


# ---------- 2.1) Get Estimates ----------

#### Aux function to get estimates ###

varSim_PP <- function(l_data, 
                      l_est, 
                      degree,
                      pbar = TRUE) {
  
  
  
  n_seq_log <- seq(3, 7.5, length = 12)
  n_seq <- round(exp(n_seq_log))
  
  # Create Storage Object
  nest6 <- vector('list', length = 7) # true, var, mvar, tvvar, tvvar_unreg, gam and thresholded gam
  nest3 <- vector('list', length = 3) # for sp/p
  nest8 <- vector('list', length = 12) # for n
  nest_n <- vector('list', length = n_files)
  
  l_sp <- nest3
  l_sp <- lapply(l_sp, function(x) nest6)
  l_p <- nest3
  l_p <- lapply(l_sp, function(x) l_sp)
  l_n <- nest8
  l_n <- lapply(l_n, function(x) l_p)
  
  for(n in 1:12) {
    for(sp in 3) {
      for(p in 2) {
        
        # save matrix for each estimation type
        l_Plist_true <- list()
        l_Plist_mvar <- list()
        l_Plist_varGLM <- list()
        l_Plist_tvmvar <- list()
        l_Plist_tvmvar_unreg <- list()
        l_Plist_gam <- list()
        l_Plist_gam_th <- list()
        
        for(i in 1:n_files) {
          
          # 1) ----- True Graph -----
          
          # Indicator Matrix
          Gind <- l_data[[i]][[n]][[sp]][[p]]$Graph$Gind
          nNodes <- nrow(Gind)
          
          # true graph, intrapolated to 20 time points
          thin <- round(seq(1, n_seq[n], length = 20))
          tg <- l_data[[i]][[n]][[sp]][[p]]$Graph$G[, , thin]
          
          # Take rowsums, to be able to subset for degree (needed for Simulation B)
          RS <- rowSums(Gind!=0)
          nodes_within_degree <- sum(RS %in% degree)
          
          # Storage for local matrix
          i_mat <- matrix(NA, ncol=24, nrow= nodes_within_degree * nNodes)
          
          browser()
          
          # Loop over all edge parameters
          counter <- 1
          for(rr in 1:nNodes) { # loop rows
            for(cc in 1:nNodes) { # loop columns
              
              if(RS[rr] %in% degree) {
                i_mat[counter, 1] <- i # edge id
                i_mat[counter, 2] <- rr # row
                i_mat[counter, 3] <- cc # column
                i_mat[counter, 4] <- Gind[rr, cc] # type
                i_mat[counter, 5:24] <- tg[rr, cc, 1:20] # true edge, thinned to 20 points
                counter <- counter + 1
              } # end if: degree
              
            } # end for: cc
          } # end for: rr
          
          l_Plist_true[[i]] <- i_mat
          
          
          # 2) ----- Stationary mVAR Estimate -----
          
          # Get out of Object
          wadj <- l_est[[i]][[n]][[sp]][[p]]$mvar$wadj[, , 1]
          signs <- l_est[[i]][[n]][[sp]][[p]]$mvar$signs[, , 1]
          
          # Add sign information
          wadj_s <- wadj
          wadj_s[!is.na(signs)] <- wadj_s[!is.na(signs)] * signs[!is.na(signs)]
          
          # Strech into array with 20 time points
          a_wadj <- array(dim=c(nrow(wadj), nrow(wadj), 20))
          a_wadj[ , , ] <- wadj_s
          
          # Storage for local matrix
          i_mat <- matrix(NA, ncol=24, nrow= nodes_within_degree * nNodes)
          
          # Loop over all edge parameters
          counter <- 1
          for(rr in 1:nNodes) { # loop rows
            for(cc in 1:nNodes) { # loop columns
              
              if(RS[rr] %in% degree) {
                i_mat[counter, 1] <- i
                i_mat[counter, 2] <- rr
                i_mat[counter, 3] <- cc
                i_mat[counter, 4] <- Gind[rr, cc]
                i_mat[counter, 5:24] <- a_wadj[rr, cc, 1:20]
                counter <- counter + 1
              } # end if: degree
              
            } # end for: cc
          } # end for: rr
          
          l_Plist_mvar[[i]] <- i_mat
          
          
          # 3) ----- Stationary GLM VAR Estimate -----
          
          # Get out of Object
          wadj <- l_est[[i]][[n]][[sp]][[p]]$glmvar[,-1]
          
          # Strech into array with 20 time points
          a_wadj <- array(dim=c(nrow(wadj), nrow(wadj), 20))
          a_wadj[,,] <- wadj
          
          # Storage for local matrix
          i_mat <- matrix(NA, ncol=24, nrow= nodes_within_degree * nNodes)
          
          # Loop over all edge parameters
          counter <- 1
          for(rr in 1:nNodes) { # loop rows
            for(cc in 1:nNodes) { # loop columns
              
              if(RS[rr] %in% degree) {
                i_mat[counter, 1] <- i
                i_mat[counter, 2] <- rr
                i_mat[counter, 3] <- cc
                i_mat[counter, 4] <- Gind[rr, cc]
                i_mat[counter, 5:24] <- a_wadj[rr, cc, 1:20]
                counter <- counter + 1
              } # end if: degree
              
            } # end for: cc
          } # end for: rr
          
          l_Plist_varGLM[[i]] <- i_mat
          
          
          # 4) ----- Time-varying mVAR Estimate -----
          
          # Get out of Object
          wadj <- l_est[[i]][[n]][[sp]][[p]]$tvmvar$wadj[, , 1, ]
          signs <- l_est[[i]][[n]][[sp]][[p]]$tvmvar$signs[, , 1, ]
          
          # Add sign information
          wadj_s <- wadj
          wadj_s[!is.na(signs)] <- wadj_s[!is.na(signs)] * signs[!is.na(signs)]
          
          # Storage for local matrix
          i_mat <- matrix(NA, ncol=24, nrow= nodes_within_degree * nNodes)
          
          # Loop over all edge parameters
          counter <- 1
          for(rr in 1:nNodes) { # loop rows
            for(cc in 1:nNodes) { # loop columns
              
              if(RS[rr] %in% degree) {
                i_mat[counter, 1] <- i
                i_mat[counter, 2] <- rr
                i_mat[counter, 3] <- cc
                i_mat[counter, 4] <- Gind[rr, cc]
                i_mat[counter, 5:24] <- wadj_s[rr, cc, 1:20]
                counter <- counter + 1
              } # end if: degree
              
            } # end for: cc
          } # end for: rr
          
          l_Plist_tvmvar[[i]] <- i_mat
          
          
          # 5) ----- Time-varying mVAR Estimate [unregularized] -----
          
          # Get out of Object
          wadj <- l_est[[i]][[n]][[sp]][[p]]$tvmvar_unreg$wadj[, , 1, ]
          signs <- l_est[[i]][[n]][[sp]][[p]]$tvmvar_unreg$signs[, , 1, ]
          
          # Add sign information
          wadj_s <- wadj
          wadj_s[!is.na(signs)] <- wadj_s[!is.na(signs)] * signs[!is.na(signs)]
          
          # Storage for local matrix
          i_mat <- matrix(NA, ncol=24, nrow= nodes_within_degree * nNodes)
          
          # Loop over all edge parameters
          counter <- 1
          for(rr in 1:nNodes) { # loop rows
            for(cc in 1:nNodes) { # loop columns
              
              if(RS[rr] %in% degree) {
                i_mat[counter, 1] <- i
                i_mat[counter, 2] <- rr
                i_mat[counter, 3] <- cc
                i_mat[counter, 4] <- Gind[rr, cc]
                i_mat[counter, 5:24] <- wadj_s[rr, cc, 1:20]
                counter <- counter + 1
              } # end if: degree
              
            } # end for: cc
          } # end for: rr
          
          l_Plist_tvmvar_unreg[[i]] <- i_mat
          
          
          # 6) ----- GAM VAR Estimate -----
          
          # Get out of Object
          wadj <- l_est[[i]][[n]][[sp]][[p]]$gam
          
          # Thin down to 20 dime points
          n_est <- dim(wadj)[3]
          wadj_20 <- wadj[,,round(seq(1, n_est, length=20)), 2] # new object: estimate at 4th dim 2
          
          # Delete intercepts and transpose
          wadj_20_t <- wadj_20[-1,,]
          for(tt in 1:20) wadj_20_t[, , tt] <- t(wadj_20[,,tt][-1,])
          
          # Storage for local matrix
          i_mat <- matrix(NA, ncol=24, nrow= nodes_within_degree * nNodes)
          
          # Loop over all edge parameters
          counter <- 1
          for(rr in 1:nNodes) { # loop rows
            for(cc in 1:nNodes) { # loop columns
              
              if(RS[rr] %in% degree) {
                i_mat[counter, 1] <- i
                i_mat[counter, 2] <- rr
                i_mat[counter, 3] <- cc
                i_mat[counter, 4] <- Gind[rr, cc]
                i_mat[counter, 5:24] <- wadj_20_t[rr, cc, 1:20]
                counter <- counter + 1
              } # end if: degree
              
            } # end for: cc
          } # end for: rr
          
          l_Plist_gam[[i]] <- i_mat
          
          
          # 7) ----- GAM VAR Estimate - thresholded 95% CI -----
          
          # Get out of Object
          wadj <- l_est[[i]][[n]][[sp]][[p]]$gam
          
          # Thin down to 20 dime points
          n_est <- dim(wadj)[3]
          wadj_20 <- wadj[,,round(seq(1, n_est, length=20)), 2] # new object: estimate at 4th dim 2
          wadj_lower_CI <- wadj[,,round(seq(1, n_est, length=20)), 1] # this was incorrect
          
          # Threshold
          wadj_20[wadj_lower_CI<0] <- 0
          
          # Delete intercepts and transpose
          wadj_20_t <- wadj_20[-1,,]
          for(tt in 1:20) wadj_20_t[, , tt] <- t(wadj_20[,,tt][-1,])
          
          # Storage for local matrix
          i_mat <- matrix(NA, ncol=24, nrow= nodes_within_degree * nNodes)
          
          # Loop over all edge parameters
          counter <- 1
          for(rr in 1:nNodes) { # loop rows
            for(cc in 1:nNodes) { # loop columns
              
              if(RS[rr] %in% degree) {
                i_mat[counter, 1] <- i
                i_mat[counter, 2] <- rr
                i_mat[counter, 3] <- cc
                i_mat[counter, 4] <- Gind[rr, cc]
                i_mat[counter, 5:24] <- wadj_20_t[rr, cc, 1:20]
                counter <- counter + 1
              } # end if: degree
              
            } # end for: cc
          } # end for: rr
          
          l_Plist_gam_th[[i]] <- i_mat
          
          
          print(paste0("n = ", n, " sp = ", sp, " p = ", p, " i = ", i))
          
        } # end for: iterations

        # Collapse Across Iterations & fill into l_n
        l_n[[n]][[sp]][[p]][[1]] <- do.call(rbind, l_Plist_true)
        l_n[[n]][[sp]][[p]][[2]] <- do.call(rbind, l_Plist_mvar)
        l_n[[n]][[sp]][[p]][[3]] <- do.call(rbind, l_Plist_varGLM)
        l_n[[n]][[sp]][[p]][[4]] <- do.call(rbind, l_Plist_tvmvar)
        l_n[[n]][[sp]][[p]][[5]] <- do.call(rbind, l_Plist_tvmvar_unreg)
        l_n[[n]][[sp]][[p]][[6]] <- do.call(rbind, l_Plist_gam)
        l_n[[n]][[sp]][[p]][[7]] <- do.call(rbind, l_Plist_gam_th)
        
      } # end for: p
    } # end for: sp
    if(pbar) print(n)
  } # end for: n
  
  return(l_n)
  
} # end of preprocessing function




#### Call function for different degrees ###

# a) For all degrees (main analysis)
l_n <- varSim_PP(l_data = l_data, 
                 l_est = l_est,
                 degree = 0:20,
                 pbar = TRUE)

saveRDS(l_n, 'files/SimA_preprocessed.RDS')
# l_n <- readRDS(paste0(simDir, '/files/SimA_preprocessed.RDS'))

# Checking Example
j <- 10
n <- 12
l_n[[n]][[3]][[2]][[1]][j, 5:24] # true
l_n[[n]][[3]][[2]][[2]][j, 5:24] # var lasso
l_n[[n]][[3]][[2]][[3]][j, 5:24] # var glm
l_n[[n]][[3]][[2]][[4]][j, 5:24] # tvmvar
l_n[[n]][[3]][[2]][[5]][j, 5:24] # tvmvar [unregularized]
l_n[[n]][[3]][[2]][[6]][j, 5:24] # tv gam
l_n[[n]][[3]][[2]][[7]][j, 5:24] # tv gam thresholded



# ---------- 2.2) Get Computation Time ----------

# set up storage
nest3 <- vector('list', length = 3)
nest8 <- vector('list', length = 12)
nest3 <- lapply(nest3, function(x) matrix(NA, nrow = n_files, ncol=11))
nest3a <- lapply(nest3, function(x) nest3)
nest8 <- lapply(nest8, function(x) nest3a)

# Loop
for(sp in 3) {
  for(p in 2) {
    for(n in 1:12) {
      for(i in 1:n_files) {
        
        c_row <- rep(NA, 11)
        c_row[1:4] <- c(i, n, sp, p)
        
        for(est in 1:4) {
          c_row[5] <- l_est[[i]][[n]][[sp]][[p]]$mvar.time
          c_row[6] <- l_est[[i]][[n]][[sp]][[p]]$glmvar.time
          c_row[7] <- l_est[[i]][[n]][[sp]][[p]]$tvmvar.time
          c_row[8] <- l_est[[i]][[n]][[sp]][[p]]$tvmvar.time_unreg
          c_row[9] <- l_est[[i]][[n]][[sp]][[p]]$gam.time
          c_row[10] <- l_est[[i]][[n]][[sp]][[p]]$bw.time
          c_row[11] <- l_est[[i]][[n]][[sp]][[p]]$bw.time_unreg
        }
        nest8[[n]][[sp]][[p]][i,] <- c_row
        
      }
    }
  }
}

# Collapse all together

time_collapse <- lapply(nest8, function(x) {
  
  do.call(rbind,lapply(x, function(y) {
    
    do.call(rbind, y)
    
  }))
  
})

time_collapse2 <- do.call(rbind, time_collapse)
time_collapse2 <- as.data.frame(time_collapse2)
colnames(time_collapse2) <- c('iter', 'n', 'sp', 'p',
                              'mvar.time', 'glmvar.time', 
                              'tvmvar.time', 'tvmvar.time_unreg',
                              'gam.time', 
                              'bw.time', 
                              'bw.time_unreg')

time_collapse2_noNA <- na.omit(time_collapse2)

saveRDS(time_collapse2_noNA, paste0(simDir,"files/varSim_CompTime.RDS"))
# time_collapse2_noNA <- readRDS(paste0(simDir,"files/varSim_CompTime.RDS"))


# ---------- 2.3) Compute Sparsity and Precision ----------

# Copy structure
a_SP <- array(dim = c(12, 3, 3, 7, 2))

for(n in 1:12) {
  for(sp in 3) {
    for(p in 2) {
      for(est in 2:7) {
        
        structure_true <- l_n[[n]][[sp]][[p]][[1]][, 5:24] != 0
        structure_est <- l_n[[n]][[sp]][[p]][[est]][, 5:24] != 0
        
        sen <- mean(structure_true[structure_true==1] == structure_est[structure_true==1])
        pre <- mean(structure_true[structure_true==0] == structure_est[structure_true==0])
        
        # Save
        a_SP[n, sp, p, est, 1] <- sen
        a_SP[n, sp, p, est, 2] <- pre
        
      }
    }
  }
}

saveRDS(a_SP, paste0(simDir,"files/varSim_StructureRecovery.RDS"))
# a_SP <- readRDS(paste0(simDir,"files/varSim_StructureRecovery.RDS"))


# --------------------------------------------------------------------------------------------
# --------------------------------- 3) Figures: Pretty for Paper -----------------------------
# --------------------------------------------------------------------------------------------


# ---------- Simulation Setup: The 8 Edge types ----------


# Take time-varying parameter functions from VarSim_aux.R :


edgetypes <- list()
N <- 1000
x <- seq(0, theta, length = N)
k <- 15
theta <- .35


edgetypes[[1]] <- rep(theta, N) # Edge 1: Constant
edgetypes[[2]] <- seq(0, theta, length = N) # Edge 2: Linear Increase
edgetypes[[3]] <- seq(theta, 0, length = N) # Edge 3: Linear Decrease
edgetypes[[4]] <- theta / (1 + exp (-k * (x - theta/2))) # Edge 4: Sigmoid Increase
edgetypes[[5]] <- theta / (1 + exp (k * (x - theta/2))) # Edge 5: Sigmoid Decrease
edgetypes[[6]] <- c(rep(theta, ceiling(N/2)), rep(0, floor(N/2))) # Edge 6: Step fuction increase
edgetypes[[7]] <- c(rep(0, floor(N/2)), rep(theta, ceiling(N/2))) # Edge 7: Step fuction decrease
edgetypes[[8]] <- rep(0, N) # Edge 8: Zero Constant


v_texts <- c("(a) Constant Nonzero",
             "(b) Linear Increase", 
             "(c) Linear Decrease", 
             "(d) Sigmoid Increase", 
             "(e) Sigmoid Decrease", 
             "(f) Step Function Up", 
             "(g) Step Function Down", 
             "(h) Constant Zero")

# Some plotting options
lwd <- 2

pdf(paste0(figDir_paper, 'Fig_8Types.pdf'), width = 12, height = 6)

par(mfrow=c(2,4))

# Loop over 8 types
for(type in 1:8) {
  
  # Setup plotting area
  plot.new()
  par(mar=c(3, 5, 3, 1))
  plot.window(xlim = c(0, N), ylim = c(-.1, .4))
  
  axis(2, seq(-.1, .4, length=6), las=2)
  axis(1, labels=FALSE)
  mtext(text = v_texts[type], side = 3, line=1, cex=1.2)

  # Continuous functions
  if(type %in% c(1:5, 8)) {
    lines(1:1000, edgetypes[[type]], lwd=lwd)
  }
  
  # Step functions
  if(type %in% 6:7) {
    lines(1:500, edgetypes[[type]][1:500], lwd=lwd)
    lines(501:1000, edgetypes[[type]][501:1000], lwd=lwd)
    segments(500.5, 0, 500.5, theta, lty=2)
  }
  
  # Legends and titles
  if(type %in% c(1, 5))  title(ylab = 'Parameter Value', cex.lab=1.4)
  if(type %in% 5:8)  title(xlab = 'Time', line = 1.5, cex.lab=1.4)
  
  text(18, .55, paste0('(', letters[counter], ')'), cex = 2)

  
}

dev.off()


# ---------- Figure A: Average Absolute Error ----------


plotFigureA <- function(CI=FALSE) {
  

# Define some vars
p_seq <- c(5, 10, 20)
sp_seq <- c(.05, .1, .2)
n_seq_log <- seq(3, 7.5, length = 12)
n_seq <- round(exp(n_seq_log))
cols <- RColorBrewer::brewer.pal(6, 'Paired')

# ----- Some Settings -----

# Fix p / np
p <- 2
sp <- 3

# define jittering, to aviod exactly overlapping lines
j <- .075
jitter <- c(3*j, -2*j, 2*j, -j, j)

# ----- Set up layout -----
mat <- matrix(1:6, ncol=2, byrow = T)
lo <- layout(mat)
# layout.show(lo)

# ----- Plot Results -----

v_titles <- c("(a) Constant Nonzero", 
              "(b) Linear", 
              "(c) Sigmoid", 
              "(d) Step Function", 
              "(e) Constant Zero")


for(ti in 1:5) {
  
  type <- c(1, 2, 4, 7, 0)[ti]
  
  # Setup plotting area
  plot.new()
  par(mar = c(4,4,1,2))
  ymax <- .4
  plot.window(xlim = c(1, 12), ylim = c(0, ymax))
  axis(1, at = 1:12, labels = n_seq, cex.axis = .8, las=2)
  axis(2, round(seq(0, ymax, length = 5), 2), las = 2, cex.axis = .8)
  
  # mtext(text = v_titles[ti], side = 3, line=1)
  text(x = 6.5, y=.4, v_titles[ti], cex=1.3)
  
  # Aggregate and plot
  l_n_errors <- list()
  l_n_errors_qt25 <- list()
  l_n_errors_qt75 <- list()
  
  for(est in 2:7) {
    
    for(n in 1:12) {
      
      data_sp_p <- l_n[[n]][[sp]][[p]]
      

      # Get True and Estimated Matrices
    
      
      # Estimated Matrix
      # Using if statements to combine symmetrical results
      if(type==1) {
        ss_est <- data_sp_p[[est]][data_sp_p[[est]][, 4] == 1, ][,5:24]
        ss_true <- data_sp_p[[1]][data_sp_p[[1]][, 4] == 1, ][,5:24]

      }
        
      if(type==0) {
        ss_est <- data_sp_p[[est]][data_sp_p[[est]][, 4] == 0, ][,5:24]
        ss_true <- data_sp_p[[1]][data_sp_p[[1]][, 4] == 0, ][,5:24]
      }
        
      if(type==2) {
        ss_est1 <- data_sp_p[[est]][data_sp_p[[est]][, 4] == 2, ][,5:24]
        ss_est2 <- data_sp_p[[est]][data_sp_p[[est]][, 4] == 3, ][,24:5]
        ss_est <- rbind(ss_est1, ss_est2)
        
        ss_true1 <- data_sp_p[[1]][data_sp_p[[1]][, 4] == 2, ][,5:24]
        ss_true2 <- data_sp_p[[1]][data_sp_p[[1]][, 4] == 3, ][,24:5]
        ss_true <- rbind(ss_true1, ss_true2)
      }
      if(type==4) {
        ss_est1 <- data_sp_p[[est]][data_sp_p[[est]][, 4] == 4, ][,5:24]
        ss_est2 <- data_sp_p[[est]][data_sp_p[[est]][, 4] == 5, ][,24:5]
        ss_est <- rbind(ss_est1, ss_est2)
        
        ss_true1 <- data_sp_p[[1]][data_sp_p[[1]][, 4] == 4, ][,5:24]
        ss_true2 <- data_sp_p[[1]][data_sp_p[[1]][, 4] == 5, ][,24:5]
        ss_true <- rbind(ss_true1, ss_true2)
        
      }
      
      if(type==7) {
        ss_est1 <- data_sp_p[[est]][data_sp_p[[est]][, 4] == 6, ][,5:24]
        ss_est2 <- data_sp_p[[est]][data_sp_p[[est]][, 4] == 7, ][,24:5]
        ss_est <- rbind(ss_est1, ss_est2)
        
        ss_true1 <- data_sp_p[[1]][data_sp_p[[1]][, 4] == 6, ][,5:24]
        ss_true2 <- data_sp_p[[1]][data_sp_p[[1]][, 4] == 7, ][,24:5]
        ss_true <- rbind(ss_true1, ss_true2)
        
      }
      
        
      # Compute Compute Absolute Error
      ss_error <- ss_true - ss_est
      ss_error_pos <- abs(ss_error)
      ss_error_pos_cm <- colMeans(ss_error_pos)
      ss_error_pt <- ss_error_pos_cm
      
      l_n_errors_qt25[[n]] <- apply(ss_error_pos, 2, function(x) quantile(x, p=0.25))
      l_n_errors_qt75[[n]] <- apply(ss_error_pos, 2, function(x) quantile(x, p=0.75))
      
      l_n_errors[[n]] <- ss_error_pt
      
      # Set to zero if tvGAM wasn't estimated
      if(est > 5) {
        if(n < 4 & p == 3) l_n_errors[[n]] <- NA
        if(n < 3 & p == 2) l_n_errors[[n]] <- NA
      }
      
      if(est > 5) {
        if(n < 3 & p == 2) {
          l_n_errors_qt25[[n]] <- NA
          l_n_errors_qt75[[n]] <- NA
        } 
      }
      
    } # end for: n
    
    # display means
    means <- unlist(lapply(l_n_errors, mean))
    # points(1:12, means, col = cols[c(1,2,3,6,4,5)][est-1], pch=20, cex = 1)
    lines(1:12, means, col = cols[est-1], pch=20, lwd=2, lty = est-1)
    
    # Display CIs
    if(CI) {
      

      # Lower CI
      mean_25 <- unlist(lapply(l_n_errors_qt25, mean)) # average across 20 estimation points
      lines(1:12, mean_25, lty = est-1, col = cols[est-1])
      # points(1:12, mean_25, lty=2, col = cols[c(1,2,3,6,5,4)][est-1], pch=21)
      
      # Upper CI
      mean_75 <- unlist(lapply(l_n_errors_qt75, mean)) # average across 20 estimation points
      lines(1:12, mean_75, lty = est-1, col = cols[est-1])
      # points(1:12, mean_75, lty=2, col = cols[c(1,2,3,6,5,4)][est-1], pch=21)
      
    }
    
  } # end for: est
  
  if(type == 0) title(xlab = 'Number of Time Points')
  if(type %in% c(1,4,0)) title(ylab = 'Mean Absolute Error')
  
} # end for: types


# ----- Plot Legend -----

plot.new()
plot.window(xlim=c(0, 1), ylim=c(0, 1))

legend("center",
       # c('GLM(L1)', 'GLM', 'KS(L1)', 'KS', 'GAM', 'GAM(st)'),
       c('GLM(L1)', 'GLM', 'KS(L1)', 'KS', 'GAM', 'GAM(st)'),
       col = cols, 
       bg = 'white',
       bty = "n",
       cex = 1.5, 
       lty=1:6, 
       lwd = rep(2, 6))




} # end of plotting Function

sc <- .9
pdf(paste0(figDir_paper, 'Fig_Sim_A_abserror.pdf'), width = 6.5*sc, height = 8*sc)
plotFigureA(CI = FALSE)
dev.off()


height = 9
pdf(paste0(figDir_paper, 'Fig_Sim_A_abserror_withCIs.pdf'), width = 6.5*sc, height = 8*sc)
plotFigureA(CI = TRUE)
dev.off()




# ---------- Figure B: For 2 types (linear up, constant up) and n = {103, 530, 1803} ----------


plotFigureB <- function(CI=FALSE) {
  
  cols <- RColorBrewer::brewer.pal(6, 'Paired')
    
  # Setup Layout
  mat <- cbind(c(7,8,9),
               c(1,2,3),
               c(4,5,6))
  
  mat <- rbind(c(10,11,12), mat)
  
  lo <- layout(mat,
               widths = c(1.2, 5.5, 5.5),
               heights = c(.8, 6, 6, 6))
  
  # Fix parameters
  p <- 2
  sp <- 3
  
  letter_c <- 1
  v_letters <- paste0("(",letters[1:9], ")")
  
  
  for(type in 1:2) {
    
    for(n in c(5, 9, 12)) {
      
      
      # Subset data
      data_sp_p <- l_n[[n]][[sp]][[p]]
      
      plot.new()
      par(mar=c(4,4,.5,.5))
      plot.window(xlim=c(1,20), ylim=c(-.2, .8))
      
      # true parameter function
      lines(data_sp_p[[1]][data_sp_p[[1]][,4] == type, ][1, 5:24], lty=1)
      
      
      # define jittering, to aviod exactly overlapping lines
      j <- .15
      jitter <- c(3*j, -2*j, 2*j, -j, j, 2*j)
      
      # Loop over estimation methods (1=true)
      for(est in 2:7) {
        
        # take subset for given estimation method
        ss_est <- data_sp_p[[est]][data_sp_p[[1]][,4] == type, ]
        ss_est <- ss_est[, -c(1:4)]
        
        # Set to zero if tvGAM wasn't estimated
        if(est == 5) {
          if(n == 1 & p > 1) ss_est <- matrix(rep(NA, 20), nrow=1)
          if(n == 2 & p > 2) ss_est <- matrix(rep(NA, 20), nrow=1)
        }
        
        
        # display means (on top)
        means <- colMeans(ss_est)
        # points(1:20 + jitter[est-1], means, col = cols[c(1,2,3,6,4,5)][est-1], pch=20, cex = 1)
        lines(1:20, means, col = cols[est-1], pch=20, lwd=1.5, lty=est-1)
        
        # display SDs
        
        
        if(CI) {
          
          qtls <- apply(ss_est, 2, function(x) quantile(x, probs = c(.10, .90)))

          lines(1:20, qtls[1, ], col = cols[est-1], pch=20, lwd=.9, lty=est-1)
          lines(1:20, qtls[2, ], col = cols[est-1], pch=20, lwd=.9, lty=est-1)

          # segments(x0 = 1:20 + jitter[est-1],
          #          y0 = qtls[1, ],
          #          x1 = 1:20 + jitter[est-1],
          #          y1 = qtls[2, ],
          #          col = cols[est-1])
          
        } 
        
      
      } # end for: est
      
      
      # Axis
      axis(1, c(1, 5, 10, 15, 20))
      axis(2, c(-.2, 0, .2, .4, .6, .8), las=2)
      
      
      # Legend
      if(type == 2 & n == 9) {
        
        # plot.new()
        # plot.window(xlim=c(0, 1), ylim=c(0, 1))
        
        legend("topleft",
               c('GLM(L1)', 'GLM', 'KS(L1)', 'KS', 'GAM', 'GAM(st)'),
               col = cols, 
               lty = 1:6,
               bg = 'white',
               bty = "n",
               cex = 1)

      }
      
      if(type == 1)     title(ylab = 'Parameter Value')
      if(n == 12) title(xlab = 'Estimation Points')
      
      
      text(19, .75, v_letters[letter_c], cex = 1.3)
      
      letter_c <- letter_c + 1
      
    }
    
  }
  
  
  # add N variation
  
  cex_metaT <- 1.7
  
  for(n_ylab in c(103, 530, 1803)) {
    
    par(mar=c(0,0,0,0))
    plot.new()
    plot.window(xlim=c(-.01, .01), ylim=c(-.01, .01))
    # rect(-.01, -.01, .01, .01, col="red")
    text(0, 0, paste0('         n = ', n_ylab), srt=90, cex=cex_metaT)
    
  }
  
  
  # Add title for edge types
  plot.new() # empty
  
  plot.new()
  par(mar=c(0,0,0,0))
  plot.window(xlim=c(-.01, .01), ylim=c(-.01, .01))
  text(0, 0, '          Constant Nonzero', cex=cex_metaT)
  
  plot.new()
  par(mar=c(0,0,0,0))
  plot.window(xlim=c(-.01, .01), ylim=c(-.01, .01))
  text(0, 0, '     Linear Increase', cex=cex_metaT)
  
} # eoF


sc <- .9
pdf(paste0(figDir_paper, 'Fig_Sim_B_2Examples.pdf'), width = 6.5*sc, height = 8*sc)
plotFigureB(CI=FALSE)
dev.off()


sc <- .9
pdf(paste0(figDir_paper, 'Fig_Sim_B_2Examples_withCIs.pdf'), width = 6.5*sc, height = 8*sc)
plotFigureB(CI=TRUE)
dev.off()





# ---------- Figure C: Precision/Sensitivity and total absolute error ----------


n_seq_log <- seq(3, 7.5, length = 12)
n_seq <- round(exp(n_seq_log))

height = 4
pdf(paste0(figDir_paper, 'Fig_Sim_C_Recov_PoolAv.pdf'), width = height * 2, height = height)

cols <- RColorBrewer::brewer.pal(6, 'Paired')

p <- 2
sp <- 3

# Set up Plotting area
par(mfrow = c(1,2),
    mar = c(4,3,3,2))


# Set un-identified variations to zero
# a_SP <- readRDS("varSim_StructureRecovery.RDS")
a_SP_new <- a_SP
a_SP_new[1:2, 3, 2, 6:7, ] <- NA

# Sensitivity
plot.new()
plot.window(xlim=c(1, 12), ylim=c(0, 1))

est_all <- 7

for(est in 2:est_all) {
  lines(a_SP_new[, sp, p, est, 1], col = cols[est-1], lwd=2, lty=est-1)
  # points(a_SP_new[, sp, p, est, 1], col = cols[est-1], pch=20)
}
axis(1, at = 1:12, label = n_seq, cex.axis = .9, las=2)
title(xlab = 'Number of Time Points')
axis(2, c(0, .5, 1), las=2)
# title(ylab = 'Sensitivity', cex = 1.5)
mtext("Sensitivity", side = 3, line=1, cex=1.2)

# Precision
plot.new()
plot.window(xlim=c(1, 12), ylim=c(0, 1))
for(est in 2:est_all) {
  lines(a_SP_new[, sp, p, est, 2], col = cols[est-1], lwd=2, lty=est-1)
  # points(a_SP_new[, sp, p, est, 2], col = cols[est-1], pch=20)
}
axis(1, at = 1:12, label = n_seq, cex.axis = .9, las=2)
title(xlab = 'Number of Time points')
axis(2, c(0, .5, 1), las=2)
# title(ylab = 'Precision', cex = 1.5)
mtext("Precision", side = 3, line=1, cex=1.2)

# Legend
legend("left",
       c('GLM(L1)', 'GLM', 'KS(L1)', 'KS', 'GAM', 'GAM(st)'),
       col = cols, 
       lty = 1:6,
       bg = 'white',
       bty = "n",
       cex = 1)


dev.off()





# ---------- Figure D: Computational Cost ----------

TC <- as.data.frame(time_collapse2_noNA)

# Plotting
pdf(paste0(figDir_paper, 'Fig_Sim_D_CompCost.pdf'), width = 6, height = 4)

par(mfrow=c(1,1), mar = c(4, 4, 1, 1))

maxtime <- 15
p <- 2 # fix p

# Setup
plot.new()
plot.window(xlim=c(1,12), ylim=c(0,maxtime))
# box()
axis(1, n_seq, at = 1:12, cex = .50)
axis(2, c(1, 5, 10,15), las = 2)


# Compute averages
n_vec_mvar <- rep(NA, 12)
n_vec_bw <- rep(NA, 12)
n_vec_GAMst <- rep(NA, 12)

for(n in 1:12) {
  n_vec_mvar[n] <- mean(TC[TC$n == n & TC$p == p, ]$tvmvar.time)
  n_vec_bw[n] <- mean(TC[TC$n == n & TC$p == p, ]$bw.time)
  n_vec_GAMst[n] <- mean(TC[TC$n == n & TC$p == p, ]$gam.time)
}
n_vec_mvar <- n_vec_mvar / 60
n_vec_bw <- n_vec_bw / 60
n_vec_GAMst <- n_vec_GAMst / 60

# Plot lines
lines(1:12, n_vec_mvar, lty = 1, lwd=2)
lines(1:12, n_vec_bw, lty = 2, lwd=2)
lines(1:12, n_vec_GAMst, lty = 3, lwd=2)

# Titles
title(xlab = "Number of time points n")
title(ylab = "Computational cost (minutes)")

# Legend
legend("topleft", 
       c("KS(L1)", 
         "Bandwidth Selection for KS(L1)", 
         "GAM(st)"), 
       lty = 1:3, lwd = rep(2, 3), 
       bty = "n")

dev.off()


