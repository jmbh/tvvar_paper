# jonashaslbeck@gmail.com; March 2019

# Define Dirs
simDir <- '...' #  SPECIFY MAIN DIRECTORY
figDir_paper <- paste0(simDir, 'figures_paper/')
dataDir_data <- paste0(simDir, "output_data/") # simulated data friles
dataDir_est <- paste0(simDir, "output_est/") # simulated estimation files

# ----------------------------------------------------------------------------------
# --------------------------------- 1) Load Data -----------------------------------
# ----------------------------------------------------------------------------------


# Data
v_files_data <- list.files(dataDir_data)
n_files <- length(v_files_data)
l_data <- list()
for(i in 1:n_files) {
  l_data[[i]] <- readRDS(paste0(dataDir_data, v_files_data[i]))
  print(i)  
}

# Estimates
v_files_est <- list.files(dataDir_est)
n_files <- length(v_files_est)
l_est <- list()
for(i in 1:n_files) {
  l_est[[i]] <- readRDS(paste0(dataDir_est, v_files_est[i]))
  print(i)  
}


# ----------------------------------------------------------------------------------
# --------------------------------- 2) Preprocessing Aux Function ------------------
# ----------------------------------------------------------------------------------

VS_UT_PP <- function(l_data, 
                     l_est, 
                     pbar = TRUE)
{
  
  n_seq_log <- seq(3, 7.5, length = 12)
  n_seq <- round(exp(n_seq_log))
  
  # ---------- Create Output Object ----------
  
  out_est <- vector("list", length = 6)
  out_n <- vector("list", length = 12)
  out_n <- lapply(out_n, function(x) out_est)
  
  
  # ---------- Loop over i & n; retrieve results of each estimator ----------
  
  for(n in 1:12) {
    
    # Storage for lists over i
    l_parList_true <- list()
    l_parList_mvar <- list()
    l_parList_varGLM <- list()
    l_parList_tvmvar <- list()
    l_parList_tvmvar_unreg <- list()
    
    l_parList_gam <- list()
    l_parList_gam_th <- list()
    
    for(i in 1:n_files) {
      
      # 1) ----- True Graph -------------------------------------
      
      # Indicator Matrix
      Gind <- l_data[[i]][[n]]$Graph$G_type
      nNodes <- nrow(Gind)
      
      # true graph, intrapolated to 20 time points
      thin <- round(seq(1, n_seq[n], length = 20))
      tg <- l_data[[i]][[n]]$Graph$G[, , thin]
      
      # Storage for local matrix
      i_mat <- matrix(NA, ncol = 24, nrow = nNodes^2)
      
      # Loop over all edge parameters
      counter <- 1
      for(rr in 1:nNodes) { # loop rows
        for(cc in 1:nNodes) { # loop columns
          i_mat[counter, 1] <- i
          i_mat[counter, 2] <- rr
          i_mat[counter, 3] <- cc
          i_mat[counter, 4] <- Gind[rr, cc]
          i_mat[counter, 5:24] <- tg[rr, cc, 1:20]
          counter <- counter + 1
        } # end for: cc
      } # end for: rr
      
      l_parList_true[[i]] <- i_mat
      
      
      # 2) ----- mVar Estimate -------------------------------------
      
      # Get out of Object
      wadj <- l_est[[i]][[n]]$mvar$wadj[, , 1]
      signs <- l_est[[i]][[n]]$mvar$signs[, , 1]
      
      # Add sign information
      wadj_s <- wadj
      wadj_s[!is.na(signs)] <- wadj_s[!is.na(signs)] * signs[!is.na(signs)]
      
      # Strech into array with 20 time points
      a_wadj <- array(dim=c(nrow(wadj), nrow(wadj), 20))
      a_wadj[ , , ] <- wadj_s
      
      # Storage for local matrix
      i_mat <- matrix(NA, ncol = 24, nrow = nNodes^2)
      
      # Loop over all edge parameters
      counter <- 1
      for(rr in 1:nNodes) { # loop rows
        for(cc in 1:nNodes) { # loop columns
          i_mat[counter, 1] <- i
          i_mat[counter, 2] <- rr
          i_mat[counter, 3] <- cc
          i_mat[counter, 4] <- Gind[rr, cc]
          i_mat[counter, 5:24] <- a_wadj[rr, cc, 1:20]
          counter <- counter + 1
        } # end for: cc
      } # end for: rr
      
      l_parList_mvar[[i]] <- i_mat
      
      
      # 3) ----- GLM Estimate -------------------------------------
      
      # Get out of Object
      wadj <- l_est[[i]][[n]]$glmvar[,-1]
      
      # Strech into array with 20 time points
      a_wadj <- array(dim=c(nrow(wadj), nrow(wadj), 20))
      a_wadj[,,] <- wadj
      
      # Storage for local matrix
      i_mat <- matrix(NA, ncol = 24, nrow = nNodes^2)
      
      if(n > 1) { # for n = 1 (N=20), this model is not identified
        
        # Loop over all edge parameters
        counter <- 1
        for(rr in 1:nNodes) { # loop rows
          for(cc in 1:nNodes) { # loop columns
            i_mat[counter, 1] <- i
            i_mat[counter, 2] <- rr
            i_mat[counter, 3] <- cc
            i_mat[counter, 4] <- Gind[rr, cc]
            i_mat[counter, 5:24] <- a_wadj[rr, cc, 1:20]
            counter <- counter + 1
          } # end for: cc
        } # end for: rr
        
      }
      
      l_parList_varGLM[[i]] <- i_mat
      
      
      # 4) ----- tv mVar Estimate -------------------------------------     
      
      # Get out of Object
      wadj <- l_est[[i]][[n]]$tvmvar$wadj[, , 1, ]
      signs <- l_est[[i]][[n]]$tvmvar$signs[, , 1, ]
      
      # Add sign information
      wadj_s <- wadj
      wadj_s[!is.na(signs)] <- wadj_s[!is.na(signs)] * signs[!is.na(signs)]
      
      # Storage for local matrix
      i_mat <- matrix(NA, ncol =24, nrow = nNodes^2)
      
      # Loop over all edge parameters
      counter <- 1
      for(rr in 1:nNodes) { # loop rows
        for(cc in 1:nNodes) { # loop columns
          i_mat[counter, 1] <- i
          i_mat[counter, 2] <- rr
          i_mat[counter, 3] <- cc
          i_mat[counter, 4] <- Gind[rr, cc]
          i_mat[counter, 5:24] <- wadj_s[rr, cc, 1:20]
          counter <- counter + 1
        } # end for: cc
      } # end for: rr
      
      l_parList_tvmvar[[i]] <- i_mat
      
      
      # 5) ----- tv mVar Estimate -------------------------------------     
      
      # Get out of Object
      wadj <- l_est[[i]][[n]]$tvmvar$wadj[, , 1, ]
      signs <- l_est[[i]][[n]]$tvmvar$signs[, , 1, ]
      
      # Add sign information
      wadj_s <- wadj
      wadj_s[!is.na(signs)] <- wadj_s[!is.na(signs)] * signs[!is.na(signs)]
      
      # Storage for local matrix
      i_mat <- matrix(NA, ncol =24, nrow = nNodes^2)
      
      # Loop over all edge parameters
      counter <- 1
      for(rr in 1:nNodes) { # loop rows
        for(cc in 1:nNodes) { # loop columns
          i_mat[counter, 1] <- i
          i_mat[counter, 2] <- rr
          i_mat[counter, 3] <- cc
          i_mat[counter, 4] <- Gind[rr, cc]
          i_mat[counter, 5:24] <- wadj_s[rr, cc, 1:20]
          counter <- counter + 1
        } # end for: cc
      } # end for: rr
      
      l_parList_tvmvar_unreg[[i]] <- i_mat
      
      
      # 6) ----- tv GAM Estimate -------------------------------------     
      
      # Get out of Object
      wadj <- l_est[[i]][[n]]$gam
      
      # Thin down to 20 dime points
      n_est <- dim(wadj)[3] # = always equal to N
      wadj_20 <- wadj[,,round(seq(1, n_est, length=20)), 2] # new object: estimate at 4th dim 2
      
      # Delete intercepts and transpose
      wadj_20_t <- wadj_20[-1 , ,]
      for(tt in 1:20) wadj_20_t[, , tt] <- t(wadj_20[,,tt][-1,])
      
      # Storage for local matrix
      i_mat <- matrix(NA, ncol = 24, nrow = nNodes^2)
      
      # Loop over all edge parameters
      counter <- 1
      for(rr in 1:nNodes) { # loop rows
        for(cc in 1:nNodes) { # loop columns
          i_mat[counter, 1] <- i
          i_mat[counter, 2] <- rr
          i_mat[counter, 3] <- cc
          i_mat[counter, 4] <- Gind[rr, cc]
          i_mat[counter, 5:24] <- wadj_20_t[rr, cc, 1:20]
          counter <- counter + 1
        } # end for: cc
      } # end for: rr
      
      l_parList_gam[[i]] <- i_mat
      
      
      # 7) ----- tv GAM (st) Estimate -------------------------------------     
      
      # Get out of Object
      wadj <- l_est[[i]][[n]]$gam
      
      # Thin down to 20 dime points
      n_est <- dim(wadj)[3]
      wadj_20 <- wadj[,,round(seq(1, n_est, length=20)), 2] # new object: estimate at 4th dim 2
      wadj_lower_CI <- wadj[,,round(seq(1, n_est, length=20)), 1] # lower CI sufficient, since all true parameters are positive
      
      # Threshold
      wadj_20[wadj_lower_CI < 0] <- 0
      
      # Delete intercepts and transpose
      wadj_20_t <- wadj_20[-1,,]
      for(tt in 1:20) wadj_20_t[, , tt] <- t(wadj_20[,,tt][-1,])
      
      # Storage for local matrix
      i_mat <- matrix(NA, ncol = 24, nrow = nNodes^2)
      
      # Loop over all edge parameters
      counter <- 1
      for(rr in 1:nNodes) { # loop rows
        for(cc in 1:nNodes) { # loop columns
          i_mat[counter, 1] <- i
          i_mat[counter, 2] <- rr
          i_mat[counter, 3] <- cc
          i_mat[counter, 4] <- Gind[rr, cc]
          i_mat[counter, 5:24] <- wadj_20_t[rr, cc, 1:20]
          counter <- counter + 1
        } # end for: cc
      } # end for: rr
      
      l_parList_gam_th[[i]] <- i_mat
      
    } # end for: i
    
    out_n[[n]][[1]] <- do.call(rbind, l_parList_true)
    out_n[[n]][[2]] <- do.call(rbind, l_parList_mvar)
    out_n[[n]][[3]] <- do.call(rbind, l_parList_varGLM)
    out_n[[n]][[4]] <- do.call(rbind, l_parList_tvmvar)
    out_n[[n]][[5]] <- do.call(rbind, l_parList_tvmvar_unreg)
    out_n[[n]][[6]] <- do.call(rbind, l_parList_gam)
    out_n[[n]][[7]] <- do.call(rbind, l_parList_gam_th)
    
  } # end for: n
  
  if(pbar) print(n)
  
  return(out_n)
  
} # EoF


# ----------------------------------------------------------------------------------
# --------------------------------- 3) Preprocess ----------------------------------
# ----------------------------------------------------------------------------------


VS_UT_prep <- VS_UT_PP(l_data = l_data,
                       l_est = l_est,
                       pbar = TRUE)

saveRDS(VS_UT_prep, file="files/VS_UT_prep.RDS")
VS_UT_prep <- readRDS(file="SimulationB/files/VS_UT_prep.RDS")

# ----------------------------------------------------------------------------------
# --------------------------------- 4) Figure --------------------------------------
# ----------------------------------------------------------------------------------

# ----- Aux functions -----

f_text <- function(text, 
                   srt = 1,
                   cex = 1, 
                   mar = c(0,0,0,0)) 
{
  
  par(mar = mar)
  plot.new()
  plot.window(xlim=c(-1, 1), ylim=c(-1,1))
  
  text(0,0, text, cex = cex, srt = srt)  
  
}


# ----- Compute aux variables -----

# Define some vars
n_seq_log <- seq(3, 7.5, length = 12)
n_seq <- round(exp(n_seq_log))
cols <- RColorBrewer::brewer.pal(6, 'Paired')

# define jittering, to aviod exactly overlapping lines
j <- .075
jitter <- c(3*j, -2*j, 2*j, -j, j, j*2)

# ----- Plotting -----

pdf(paste0(figDir_paper, 'Fig_Sim_VS_TU_condDegrees.pdf'), width = 8, height = 10)

# --- Set up layout ---

par(mar = c(0,0,0,0))

lo <- matrix(c(1, 2, 3, 4,
               5, 11, 17, 23,
               6, 12, 18, 24, 
               7, 13, 19, 25,
               8, 14, 20, 26, 
               9, 15, 21, 27,
               10,16, 22, 28), nrow = 7, byrow = T)

layout(mat = lo, 
       widths = c(1,1,1,1), 
       heights = c(.3,1,1,1,1,1,1))


# A) --- Top Row: degrees ---

# Select degrees
d_seq <- c(1, 10, 20)

plot.new()
for(d in 1:3) f_text(paste0("        Indegree = ", d_seq[d]), cex = 1.5)


# margins for all remaining panels ..
mar <- c(3,4,1,1)

# B) --- Left column: types of parameters ---


# -- Define true parameter functions --

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

titles <- c("Constant nonzero",
            "Linear increase", 
            "Sigmoid increase",
            "Step function",
            "Constant zero")


# -- Plotting --

# All together
plot.new()
par(mar = mar)
plot.window(xlim = c(0, N), ylim = c(-.1, .4))
lwd <- 1
axis(2, seq(-.1, .4, length=6), las=2)
axis(1, labels=FALSE)

for(type in 1:8) {
  # Continuous functions
  if(type %in% c(1:5, 8)) {
    lines(1:1000, edgetypes[[type]], lwd=lwd, col="grey")
  }
  
  # Step functions
  if(type %in% 6:7) {
    lines(1:500, edgetypes[[type]][1:500], lwd=lwd, col="grey")
    lines(501:1000, edgetypes[[type]][501:1000], lwd=lwd, col="grey")
    segments(500.5, 0, 500.5, theta, lty=2, col="grey")
  }
}

# Legends and titles
title(ylab = 'Parameter Value', cex.lab = 1)
title(xlab = 'Time', line = .8, cex.lab = 1)
mtext("Averaged across Types", side = 3, line = .3, cex = .8)


# Seperately

count <- 1
for(type in c(1, 2, 4, 6, 8)) {
  
  plot.new()
  par(mar = mar)
  plot.window(xlim = c(0, N), ylim = c(-.1, .4))
  lwd <- 1
  axis(2, seq(-.1, .4, length=6), las=2)
  axis(1, labels=FALSE)
  
  
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
  title(ylab = 'Parameter Value', cex.lab = 1)
  title(xlab = 'Time', line = .8, cex.lab = 1)
  mtext(titles[count], side = 3, line = .3, cex = .8)
  count <- count + 1
}




# C) --- Results ---

par(mar = c(0,0,0,0))

letter_c <- 1
l_types <- list()
l_types[[1]] <- c(1, 2, 4, 7, 0)
l_types[[2]] <- 1
l_types[[3]] <- 2:3 # collapse across symmetric functions
l_types[[4]] <- 4:5
l_types[[5]] <- 6:7
l_types[[6]] <- 0

lo <- matrix(c(1, 2, 3, 4,
               5, 11, 17, 23,
               6, 12, 18, 24, 
               7, 13, 19, 25,
               8, 14, 20, 26, 
               9, 15, 21, 27,
               10,16, 22, 28), nrow = 7, byrow = T)


# Loop over times (as in Figure A)
v_letters <- paste0("(",letters[1:18], ")")
v_letters_ord <- c(1, 4, 7, 10, 13, 16, 
                   2, 5, 8, 11, 14, 17,
                   3, 6, 9, 12, 15, 18)

# Loop over degree
for(d in d_seq) {
  
  for(type in 1:6) {
    
    plot.new()
    par(mar = mar)
    ymax <- .5
    plot.window(xlim = c(1, 12), ylim = c(0, ymax))
    
    l_n_errors <- list()
    
    for(est in 2:7) {
      
      for(n in 1:12) {
        
        # Set to zero if tvGAM/GLM wasn't estimated
        if(est > 4 & n < 4) {
          
          l_n_errors[[n]] <- NA
          
        } else if(est == 3 & n == 1) {
          
          l_n_errors[[n]] <- NA
          
        } else {
          
          data <- VS_UT_prep[[n]]
          
          d_inv <- abs(d - 21) # degree = abs(node number - 21)
          
          ss_true <- data[[1]][data[[1]][,4] %in% l_types[[type]] & data[[1]][,2] == d_inv, ][,5:24] # True matrix (est==1)
          ss_est <- data[[est]][data[[est]][,4] %in% l_types[[type]] & data[[est]][,2] == d_inv, ][,5:24] # Estimated Matrix
          
          
          
          
          # Compute Absolute Error
          ss_error <- ss_true - ss_est
          ss_error_pos <- abs(ss_error)
          ss_error_pos_cm <- colMeans(ss_error_pos)
          ss_error_pt <- ss_error_pos_cm
          
          l_n_errors[[n]] <- ss_error_pt # still errors for each time point
          
        }
        
      } # end for: n
      
      # display means
      means <- unlist(lapply(l_n_errors, mean))
      # means <- unlist(lapply(l_n_errors, sd))
      
      if(d==d_seq[3] & type == 6) {
        legend("left",
               c('GLM(L1)', 'GLM', 'KS(L1)', 'KS', 'GAM', 'GAM(st)'),
               col = cols, 
               lty = 1:6,
               bg = 'white',
               bty = "n",
               cex = 1.3, 
               lwd = rep(1.5, 6))
      } else {
        # points(1:12, means, col = cols[est-1], pch=20, cex = 1)
        lines(1:12, means, col = cols[est-1], pch=20, lty=est-1, lwd=1.5)  
      }
      
    } # end for: est
    
    if(d==d_seq[3] & type == 6) {
      
    } else {
      if(type == 0) title(xlab = 'Number of Time points')
      if(d == 1) title(ylab = 'Mean Absolute Error')
      axis(1, at = 1:12, labels = n_seq, cex.axis = .65, las =2)
      axis(2, round(seq(0, ymax, length = 5), 2), las = 2, cex.axis = .8)
      text(11.5, .48, v_letters[v_letters_ord[letter_c]], cex = 1.3)
    }
    
    
    letter_c <- letter_c + 1
    
  } # end for: types
  
} # end for: degree


dev.off()


