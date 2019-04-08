
# ----------------------------------------------------------------------------------
# ---------------------------- Graph Generating Function ---------------------------
# ----------------------------------------------------------------------------------

GenerateGraph_fix <- function(theta, 
                              N,
                              p, 
                              seed, 
                              verbose=FALSE,
                              max_check=1000) 
  
  
{
  
  # ---------- Define fixed Graph ----------
  
  G_init <- matrix(0, p, p)
  diag(G_init) <- 1
  G_init[upper.tri(G_init)] <- 1
  
  
  # ---------- Define time-varying parameter functions ----------
  
  # Define Graph Types
  edgetypes <- list()
  x <- seq(0, theta, length = N)
  k <- 15
  
  edgetypes[[1]] <- rep(theta, N) # Edge 1: Constant
  edgetypes[[2]] <- seq(0, theta, length = N) # Edge 2: Linear Increase
  edgetypes[[3]] <- seq(theta, 0, length = N) # Edge 3: Linear Decrease
  edgetypes[[4]] <- theta / (1 + exp (-k * (x - theta/2))) # Edge 4: Sigmoid Increase
  edgetypes[[5]] <- theta / (1 + exp (k * (x - theta/2))) # Edge 5: Sigmoid Decrease
  edgetypes[[6]] <- c(rep(theta, ceiling(N/2)), rep(0, floor(N/2))) # Edge 6: Step fuction increase
  edgetypes[[7]] <- c(rep(0, floor(N/2)), rep(theta, ceiling(N/2))) # Edge 7: Step fuction decrease
  
  
  # ---------- Assign parameters to Graph ----------
  
  # get indicator for which entries are 1
  ind_present <- which(G_init == 1, arr.ind = TRUE)
  
  # ---- 2) Assign the Type of Edges -----
  
  G_type <- matrix(0, p, p) # Storage: summary information
  G <- array(0, dim = c(p, p, N)) # the 3D time-varying VAR array
  
  for(i in 1:nrow(ind_present)) {
    type_i <- sample(1:7, size = 1)
    G[ind_present[i, 1], ind_present[i, 2], ] <- edgetypes[[type_i]]
    G_type[ind_present[i, 1], ind_present[i, 2]] <- type_i
  }
  
  
  # ---------- Check whether VAR matrix is stable ----------
  
  # Check condition for each time step
  G_check <- apply(G, 3, function(x) {
    eig_x <- eigen(x)
    sum(abs(eig_x$values) >= 1) > 0 # TRUE if violated
  })
  
  counter <- 1
  if(sum(G_check) == 0) {
    check <- 1 # break the loop
  } else {
    counter <- counter + 1
  }
  
  # To avoid infinite loop: break after 1000 loops:
  
  if(counter > max_check) stop('No graph within the constrained region found in 1000 iterations.')
  if(verbose) print(counter)
  
  
  # ---------- Return tv-VAR array ----------
  
  outlist <- list('G_init' = G_init,
                  'G_type' = G_type,
                  "G" = G,
                  'counter' = counter)
  
  return(outlist)
  
} # end of Function


# GGf <- GenerateGraph_fix(theta = .35,
#                          N = 1000,
#                          p = 20,
#                          seed = 1,
#                          verbose=TRUE,
#                          max_check=100)


# ----------------------------------------------------------------------------------
# --------------------- Generate Data from VAR-(1) model ---------------------------
# ----------------------------------------------------------------------------------

VARData <- function(graph,
                    th,
                    sd,
                    seed) {
  
  set.seed(seed)
  
  N <- dim(graph)[3]
  p <- dim(graph)[1]
  
  data <- matrix(NA, N, p)
  data[1,] <- rnorm(p, 0, sd)
  
  
  for(t in 2:N) {
    for(node in 1:p) {
      data[t, node] <- th[t, node] + sum(data[t-1,] * graph[node,,t]) + rnorm(1, mean = 0, sd = sd)
    }
  }
  
  return(data)
  
}
