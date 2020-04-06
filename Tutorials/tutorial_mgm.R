# jonashaslbeck@gmail.com, February 2020

library(devtools)
# install_github("jmbh/mgm")

figDir <- "" # Specify directory in which figures should be saved
library(qgraph)

# ----------------------------------------------------------------------------------
# ----------------------- 1) Data Preparation --------------------------------------
# ----------------------------------------------------------------------------------

library(mgm) # 1.2-8

# Subset Mood variables
mood_data <- as.matrix(symptom_data$data[, 1:12])
mood_labels <- symptom_data$colnames[1:12]
colnames(mood_data) <- mood_labels
time_data <- symptom_data$data_time

dim(mood_data)
head(mood_data[,1:7])
head(time_data)

# ----------------------------------------------------------------------------------
# ----------------------- 2) Estimation --------------------------------------------
# ----------------------------------------------------------------------------------

# Bandwidth selection (in Appendix A)
bwSeq <- seq(0.01, 1, length = 10)

set.seed(1)
bw_object <- bwSelect(data = mood_data,
                      type = rep("g", 12),
                      level = rep(1, 12),
                      bwSeq = bwSeq,
                      bwFolds = 1,
                      bwFoldsize = 10,
                      modeltype = "mvar",
                      lags = 1,
                      scale = TRUE,
                      timepoints = time_data$time_norm,
                      beepvar = time_data$beepno,
                      dayvar = time_data$dayno,
                      pbar = TRUE)

bandwidth <- bwSeq[which.min(bw_object$meanError)]
bandwidth

bandwidth <- .34

# Estimate Model on Full Dataset
set.seed(1)
tvvar_obj <- tvmvar(data = mood_data,
                    type = rep("g", 12),
                    level = rep(1, 12), 
                    lambdaSel = "CV",
                    timepoints = time_data$time_norm, 
                    estpoints = seq(0, 1, length = 20), 
                    bandwidth = bandwidth,
                    lags = 1,
                    beepvar = time_data$beepno,
                    dayvar = time_data$dayno,
                    scale = TRUE,
                    pbar = TRUE)

# Check on how much data was used
tvvar_obj


# ----------------------------------------------------------------------------------
# ----------------------- 3) Reliability -------------------------------------------
# ----------------------------------------------------------------------------------

t1 <- proc.time()[3]

res_obj <- resample(object = tvvar_obj, 
                    data = mood_data, 
                    nB = 50, 
                    blocks = 10,
                    seeds = 2:51, 
                    quantiles = c(.05, .95))

# saveRDS(res_obj, file="res_obj_nB50_bw34.RDS")
# res_obj <- readRDS(file="res_obj_nB50_bw34.RDS")

proc.time()[3] - t1 # Note that this takes a while

res_obj$bootQuantiles[2, 2, 1, , 1]
res_obj$bootQuantiles[2, 2, 1, , 2]


# ----------------------------------------------------------------------------------
# ----------------------- 4) Predictability -------------------------------------------
# ----------------------------------------------------------------------------------


pred_obj <- predict(object = tvvar_obj, 
                    data = mood_data, 
                    errorCon = c("R2", "RMSE"),
                    tvMethod = "weighted", 
                    consec = time_data$beepno)

pred_obj$errors
pred_obj$tverrors


# ----------------------------------------------------------------------------------
# ----------------------- 5) Visualization -----------------------------------------
# ----------------------------------------------------------------------------------


# ----- Aux Functions ------

# Function for timeline
f_timeline_new <- function(length = .15, 
                           gap = .005, 
                           mar = c(0,0,0,0), 
                           ylim = c(-.1,.1), 
                           ytext = -.1,
                           cex = 1) {
  
  
  par(mar=mar)
  plot.new()
  plot.window(xlim=c(0,1), ylim=ylim)
  # box()
  
  # arrows
  p_weeks <- c(4,14,4,12)
  bor_end <- c(0,cumsum(p_weeks)/sum(p_weeks))
  for(i in 1:4) {
    arrows(bor_end[i]+gap, 0, bor_end[i+1]-gap, code=3, length=length, lwd=1.5)
  }
  
  # text
  t_lengths <- p_weeks / sum(p_weeks)
  midpoints <- bor_end[-1] - t_lengths/2
  
  text(midpoints, rep(ytext, 4), c("Baseline (4w)",
                                   "Double-blind period (14w)",
                                   "Postassessment (4w)",
                                   "Additional Postassessment (12w)"),
       cex = cex)
  
  # change of medication
  points(c(42,98) / (sum(p_weeks)*7), rep(0,2), pch=20, cex=1.5)
  
}


# ----- Preprocessing  ------

# Compute mean movel over time to create decent layout
mean_wadj <- apply(tvvar_obj$wadj[, , 1, ], 1:2, mean)

par_ests <- tvvar_obj$wadj
ind_negative <- which(tvvar_obj$signs == -1, arr.ind = T)
par_ests[ind_negative] <- par_ests[ind_negative] * -1

# Find parameters with highest SD
wadj_ws <- tvvar_obj$wadj
wadj_ws[tvvar_obj$edgecolor=="red"] <- wadj_ws[tvvar_obj$edgecolor=="red"] * -1
parm_sds <- apply(wadj_ws, 1:2, sd)
parm_sds_mat <- matrix(NA, 12^2, 3)
counter <- 1
for(i in 1:12) {
  for(j in 1:12) {
    parm_sds_mat[counter, ] <- c(i, j, parm_sds[i, j]) 
    counter <- counter + 1
  }
}

parm_sds_mat_ord <- parm_sds_mat[order(parm_sds_mat[, 3], decreasing = TRUE), ]
head(parm_sds_mat_ord) # six most time-varying parameters


# ----- Plotting ------

library(qgraph)

pdf(paste0(figDir, "Fig_Application_mgm.pdf"), width = 8, height = 7)

# 1) Define Layout

lmat <- matrix(c(1, 2, 3,
                 4, 4, 4,
                 5, 5, 5), ncol=3, byrow = T)
lo <- layout(lmat, 
             heights = c(.7,.1, .6), 
             widths = c(1, 1, 1))


# 2) Two Network Plots

# Get layout of mean graph
Q <- qgraph(t(mean_wadj), DoNotPlot=TRUE)
saveRDS(Q$layout, "Tutorials/files/layout_mgm.RDS")

# Plot graph at selected fixed time points
tpSelect <- c(2, 10, 18)

# Switch to colorblind scheme
tvvar_obj$edgecolor[, , , ][tvvar_obj$edgecolor[, , , ] == "darkgreen"] <- c("darkblue")
lty_array <- array(1, dim=c(12, 12, 1, 20))
lty_array[tvvar_obj$edgecolor[, , , ] != "darkblue"] <- 2

for(tp in tpSelect) {
  qgraph(t(tvvar_obj$wadj[, , 1, tp]), 
         layout = Q$layout,
         edge.color = t(tvvar_obj$edgecolor[, , 1, tp]), 
         labels = mood_labels, 
         vsize = 13, 
         esize = 10,
         asize = 10, 
         mar = rep(5, 4), 
         minimum = 0, 
         maximum = .5, 
         lty = t(lty_array[, , 1, tp]),
         pie = pred_obj$tverrors[[tp]][, 3])
}

# 4) Timeline
f_timeline_new(length = .1, 
               mar=c(0, 4, 0, 1), 
               ylim = c(-1.2, .2), 
               ytext = -.9, 
               cex = 1)

# 5) Line-plots + CIs
plot.new()
par(mar = c(4,4,0,1))
plot.window(xlim=c(1, 20), ylim=c(-.25, .55))
axis(1, c(1, 5, 10, 15, 20), labels=T)
axis(2, c(-.25, 0, .25, .5), las=2)
abline(h = 0, col = "grey", lty=2)
title(xlab = "Estimation points", cex.lab = 1.2)
title(ylab = "Parameter estimate", cex.lab = 1.2)


head(parm_sds_mat_ord) # pick three highest
m_par_display <- matrix(c(1, 1, 
                          4, 12, 
                          10, 4), ncol = 2, byrow = T)

# Select colors
library(RColorBrewer)
library(scales)
cols <- brewer.pal(5, "Set1")[c(2,4,5)] # avoid red/green because used for edges in upper panel

for(i in 1:nrow(m_par_display)) {
  
  par_row <- m_par_display[i, ]
  
  ## Plot point estimates
  P1_pointest <- par_ests[par_row[1], par_row[2], 1, ]
  lines(1:20, P1_pointest, col = cols[i], lwd = 2, lty=i) 
  
  
  ## Plot uncertainty estimates [new shading]
  # Compute CIs
  CIs <- apply(res_obj$bootParameters[par_row[1], par_row[2], 1, , ], 1, function(x) {
    quantile(x, probs = c(.05, .95))
  } )
  
  # Plot shading
  polygon(x = c(1:20, 20:1), y = c(CIs[1,], rev(CIs[2,])), col=alpha(colour = cols[i], alpha = .3), border=FALSE)
  

  
} # end for: i

# Legend
legend_labels <- c(expression("Relaxed"["t-1"]  %->%  "Relaxed"["t"]),
                   expression("Satisfied"["t-1"]  %->%  "Strong"["t"]),
                   expression("Guilty"["t-1"]  %->%  "Satisfied"["t"]))

legend(1, .49, 
       legend_labels,
       col = cols, 
       lwd = 2, bty = "n", cex = 1.5, horiz=T, lty=1:3)


dev.off()


# ----------------------------------------------------------------------------------
# ----------------------- 4) Time-varying or not? A Bootstrap test -----------------
# ----------------------------------------------------------------------------------


# ----- A) Fit once time-varying on actual data  -----------------------------------

# Get RMSE from time-varying model
pred_emp <- predict(object = tvvar_obj, 
                    data = mood_data, 
                    tvMethod = "closestModel")

error_emp <- mean(pred_emp$errors$RMSE)



# ----- B) Fit stationary model to be able to simulate data ------------------------

var_obj <- mvar(data = mood_data,
                type = rep("g", 12),
                level = rep(1, 12), 
                lambdaSel = "CV",
                timepoints = time_data$time_norm, 
                lags = 1,
                beepvar = time_data$beepno,
                dayvar = time_data$dayno,
                scale = TRUE,
                pbar = TRUE)

# ----- C) Fit time-varying model to stationary data -------------------------------

library(mlVAR)

# Generate data
n <- 876 #numer of efficient data points, see above
nIter <- 100

l_data <- list()
l_model <- list()
l_meanerror <- rep(NA, nIter)

# get parameters out of model object
pars <- var_obj$wadj[, , 1]
pars[var_obj$edgecolor=="red"] <- pars[var_obj$edgecolor=="red"] * -1
intercepts <- unlist(var_obj$intercepts)

for(i in 1:nIter) {
  
  l_data[[i]] <- simulateVAR(pars = pars, 
                             means = intercepts, 
                             Nt = n, 
                             residuals = 1)  
  
  l_model[[i]] <- tvmvar(data = l_data[[i]],
                         type = rep("g", 12),
                         level = rep(1, 12), 
                         lambdaSel = "CV",
                         timepoints = 1:n, 
                         estpoints = seq(0, 1, length = 20), 
                         bandwidth = bandwidth,
                         lags = 1,
                         scale = TRUE,
                         pbar = FALSE, 
                         signInfo = FALSE)
  
  pred_res <- predict(object = l_model[[i]], 
                      data = l_data[[i]], 
                      tvMethod = "closestModel")
  l_meanerror[i] <- mean(pred_res$errors$RMSE)
  
  print(i)
  
} # end for: i


# saveRDS(l_model, file="l_model.RDS")
# saveRDS(l_meanerror, file="l_meanerror.RDS")


# ----- D) Evaluate ----------------------------------------------------------------

hist(l_meanerror, xlim=c(0.9, 1)) # sampling distribution under H0
abline(v = error_emp, col="red") # empirical error

# Conclusion: significant by essentially any alpha-level





