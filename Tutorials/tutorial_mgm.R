# jonashaslbeck@gmail.com, August 2017

remove.packages("mgm") # in case old version is installed
.rs.restartR()
library(devtools)
install_github("jmbh/mgm")

figDir <- "" # Specify directory in which figures should be saved
library(qgraph)

# ----------------------------------------------------------------------------------
# ----------------------- 1) Data Preparation --------------------------------------
# ----------------------------------------------------------------------------------

library(mgm) # 1.2-2

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
                    estpoints = seq(0, nrow(mood_data) - 1, length = 20), 
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
                    seeds = 1:50, 
                    quantiles = c(.05, .95))

saveRDS(res_obj, file="res_obj_nB50_newbw34.RDS")
res_obj <- readRDS(file="res_obj_nB50_newbw34.RDS")

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

mean_wadj <- apply(tvvar_obj$wadj[, , 1, ], 1:2, mean) # also used for mean layout below

par_ests <- tvvar_obj$wadj
ind_negative <- which(tvvar_obj$signs == -1, arr.ind = T)
par_ests[ind_negative] <- par_ests[ind_negative] * -1


# get largest effects
larg <- sort(as.numeric(mean_wadj), decreasing = T)
m_largPar <- matrix(NA, nrow = 20, ncol = 2)
for(i in 1:20) m_largPar[i, ] <- which(mean_wadj == larg[i], arr.ind = TRUE)


# ----- Plotting ------

pdf(paste0(figDir, "Fig_Application_mgm.pdf"), width = 8, height = 7)

# 1) Define Layout

lmat <- matrix(c(1, 2, 3,
                 4, 4, 4,
                 5, 5, 5), ncol=3, byrow = T)
lo <- layout(lmat, 
             heights = c(.7,.1, .85), 
             widths = c(1, 1, 1))
# layout.show(lo)


# 2) Two Network Plots

# Get layout of mean graph
Q <- qgraph(t(mean_wadj), DoNotPlot=TRUE)
saveRDS(Q$layout, "layout_mgm.RDS")

# Plot graph at selected fixed time points
tpSelect <- c(8, 15, 18)

for(tp in tpSelect) {
  qgraph(t(tvvar_obj$wadj[, , 1, tp]), 
         layout = Q$layout,
         edge.color = t(tvvar_obj$edgecolor[, , 1, tp]), 
         labels = mood_labels, 
         vsize = 13, 
         esize = 10,
         asize = 10, 
         mar = c(6, 6, 6, 6), 
         minimum = 0, 
         maximum = .45, 
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
par(mar = c(4,4,1,1))
plot.window(xlim=c(1, 20), ylim=c(-.5, .75))
axis(1, c(1, 5, 10, 15, 20), labels=T)
axis(2, c(-.5, -.25, 0, .25, .5, .75), las=2)
abline(h = 0, col = "grey", lty=2)
title(xlab = "Estimation points", cex.lab = 1.2)
title(ylab = "Parameter estimate", cex.lab = 1.2)

# Displayed parameters:
# Down (2) -> Down (2)
# Satisfied (4) -> Relaxed (1) 
# Satisfied (4) -> Down (2) 
# mood_labels
# m_largPar

# mood_labels[c(2,4)]

m_par_display <- matrix(c(2, 2, 
                          1, 4, 
                          2, 4), ncol = 2, byrow = T)

# Select colors
library(RColorBrewer)
cols <- brewer.pal(5, "Set1")[c(2,4,5)] # avoid red/green because used for edges in upper panel
v_jitter <- c(-.1, 0, .1)

for(i in 1:nrow(m_par_display)) {
  
  # Plot point estimates
  par_row <- m_par_display[i, ]
  P1_pointest <- par_ests[par_row[1], par_row[2], 1, ]
  points((1:20)+v_jitter[i], P1_pointest, col = cols[i], pch = 20, cex = 2) 
  lines((1:20)+v_jitter[i], P1_pointest, col = cols[i], lwd = 2) 
  
  # Plot uncertainty estimates
  CIs <- apply(res_obj$bootParameters[par_row[1], par_row[2], 1, , ], 1, function(x) {
    quantile(x, probs = c(.05, .95))
  } )
  CIs <- CIs #- P1_pointest # center them!
  
  segments((1:20)+v_jitter[i], CIs[1,]  , (1:20)+v_jitter[i], CIs[2,],
           col = cols[i],
           lwd = 2)
  
}


# Legend
legend_labels <- c(expression("Down"["t-1"]  %->%  "Down"["t"]),
                   expression("Satisfied"["t-1"]  %->%  "Relaxed"["t"]),
                   expression("Satisfied"["t-1"]  %->%  "Down"["t"]))

legend(1, .72, 
       legend_labels,
       col = cols, 
       lwd = 2, bty = "n", cex = 1.5, horiz=T)


dev.off()








