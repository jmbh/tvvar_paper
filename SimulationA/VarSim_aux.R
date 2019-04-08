

# ----------------------------------------------------------------------------------
# ---------------------------- Graph Generating Function ---------------------------
# ----------------------------------------------------------------------------------


# ---------- Graph Generator: 5 Edge Types ----------

GenerateGraph <- function(p,
                          sp,
                          theta,
                          N,
                          seed, 
                          verbose=FALSE, 
                          max_check=10000) {

  set.seed(seed)

  # ---------- Define Graph Types ----------

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


  # ---------- Generate and Assign Edges ----------

  check <- 0
  counter <- 1

  while(check == 0) {
    
    # ---- 1) Generate Graph Structure -----
    
    # What are we doing here?
    # a) Set all autocorrelations in the diagonal to be present
    # b) sample size off-diagonal elements, reflecting the specified P(e)

    # Generate Empty Graph
    n_edges <- p^2 - p # minus diagonal
    n_upper_diag <- p*(p-1)/2 # same for lower diag, of course
    size <- round(n_edges*sp)
    v_e <- sample(2:n_edges, size = size, replace = FALSE)

    # Fill in 2d Dummy matrix
    G_aux <- matrix(NA, p, p)
    diag(G_aux) <- 1
    G_aux[upper.tri(G_aux)] <- 2:(n_upper_diag+1) # not using 1, so I can subset using unequal 3 lines down
    G_aux[lower.tri(G_aux)] <- (n_upper_diag+2):(n_upper_diag*2+1)
    G_aux[which(G_aux %in% v_e, arr.ind = TRUE)] <- 1
    G_aux[G_aux != 1] <- 0

    # get indicator for which entries are 1
    ind_present <- which(G_aux == 1, arr.ind = TRUE)

    # ---- 2) Assign the Type of Edges -----
    
    Gind <- matrix(0, p, p) # Storage: summary information
    G <- array(0, dim = c(p, p, N)) # the 3D time-varying VAR array
    
    for(i in 1:nrow(ind_present)) {
      type_i <- sample(1:7, size = 1)
      G[ind_present[i, 1], ind_present[i, 2], ] <- edgetypes[[type_i]]
      Gind[ind_present[i, 1], ind_present[i, 2]] <- type_i
    }

    # ---------- Check Eigenvalue Condition ----------

    # Check condition for each time step
    G_check <- apply(G, 3, function(x) {
      eig_x <- eigen(x)
      sum(abs(eig_x$values) >= 1) > 0 # TRUE if violated
    })

    if(sum(G_check) == 0) {
      check <- 1 # break the loop
    } else {
      counter <- counter + 1
    }

    # To avoid infinite loop: break after 1000 loops:

    if(counter > max_check) stop('No graph within the constrained region found in 1000 iterations.')
    if(verbose) print(counter)

  }




  # ---------- Export ----------

  outlist <- list('G' = G,
                  'Gind' = Gind,
                  'counter' = counter)

  return(outlist)

} # eoF: GenerateGraph


# 
# seed <- 11
# 
# G <- GenerateGraph(p = 10,
#                    sp = .05,
#                    theta = .35,
#                    N = 10,
#                    seed = seed,
#                    verbose = TRUE)
# 
# G$Gind
# g <- G$G[, , 1]
# 
# G2 <- GenerateGraph(p = 10,
#                    sp = .3,
#                    theta = .35,
#                    N = 100,
#                    seed = seed,
#                    verbose = TRUE)
# 
# G2$Gind
# 
# g2 <- G2$G[, , 1]
# 
# par(mfrow=c(1,2))
# qgraph(g, asize = 5)
# qgraph(g2, asize = 5)



# 
# EV <- eigen(G$G[,,1])$values
# abs(EV)
# 
# any(class(EV) == "complex")
# 
# class(EV[1])
# abs(eigen(G$G[,,1])$values)
# 
# 
# G_check <- apply(G, 3, function(x) {
#   eig_x <- eigen(x)
#   sum(abs(eig_x$values) >= 1) > 0 # TRUE if violated
# })
# 

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


# 
# data <- VARData(graph = G$G,
#                 th = matrix(0, 1000, 20), 
#                 sd = 1, 
#                 seed = 1)
# 
# 
# # check convergence
# plot.new()
# plot.window(xlim = c(1,1000), ylim=c(-4.5, 4.5))
# for(i in 1:20) lines(data[,i], type='l', col = i)
# axis(2)
# 
# 
# 
# 


