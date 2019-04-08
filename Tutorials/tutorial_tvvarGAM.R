# jonashaslbeck@gmail.com, March 2019

# remove.packages("tvvarGAM") # in case old version is installed
# .rs.restartR()
# library(devtools)
# install_github("LauraBringmann/tvvarGAM")

library(tvvarGAM)
library(qgraph) # plotting

# ----------------------------------------------------------------------------------
# ----------------------- 1) Data Preparation --------------------------------------
# ----------------------------------------------------------------------------------

library(mgm) # For symptom dataset

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

tvvargam_obj <- tvvarGAM(data = mood_data,
                         beepvar = time_data$beepno,
                         dayvar = time_data$dayno,
                         estimates = TRUE,
                         plot = FALSE,
                         nb = 20,
                         scale = TRUE)

# Save model object
saveRDS(tvvargam_obj, file="tvvargam20_obj.RDS")
tvvargam_obj <- readRDS(file="tvvargam20_obj.RDS")


# ----------------------------------------------------------------------------------
# ----------------------- 3) Visualization -----------------------------------------
# ----------------------------------------------------------------------------------

library(qgraph)

figDir <- "/Volumes/Macintosh HD 2/Dropbox/Time-Varying Networks/paper2_ComparingMethods/3_code/gamtvvar_sim/figures_paper/"

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
  #for(i in 1:3) abline(v=bor_end[i+1], lty=2, col="grey")
  
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

# ----- Preprocessing tvvarGAM() output ------

# Thin to 20 estimation points and delete intercept
timeL <- dim(tvvargam_obj$Results_GAM$Estimate)[3]
wadj_point <- tvvargam_obj$Results_GAM$Estimate[-1, , ]

ind_thin <- round(seq(1, timeL, length = 20))

wadj_point <- tvvargam_obj$Results_GAM$Estimate[-1, , ind_thin]
wadj_CIlow <- tvvargam_obj$Results_GAM$CI_low[-1, , ind_thin]
wadj_CIup <- tvvargam_obj$Results_GAM$CI_high[-1, , ind_thin]
wadj_point_unthresh <- wadj_point

## Code Laura: Start
# This is to only show the significant arrows, see Bringmann et al. 2017 Psychological Methods for more information.
unlist_this<-function(x){matrix(unlist(x$s.table[-1,4]),12,1)}
sign<-matrix(unlist(lapply(lapply(tvvargam_obj$model,summary),unlist_this)),12,12)

for(ii in 1:12){
  for(j in 1:12){
    if(sign[j,ii]<0.05)
    {wadj_point[j,ii,]<-wadj_point[j,ii,]}
    else{wadj_point[j,ii,]<-0}}}
wadj_point_thresh<-wadj_point
## Code Laura: End

ind_overlap <- sign(round(wadj_CIlow,2)) == sign(round(wadj_CIup,2)) # if signs are the same, there is no overlap with zero
wadj_point_thresh[!ind_overlap] <- 0


# ----- Plotting ------

pdf("Fig_Application_tvvarGAM_thresh.pdf", width = 8, height = 7)

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
Q <- list()
Q$layout <- readRDS("layout_mgm.RDS")


# Plot graph at selected fixed time points
tpSelect <- c(8, 15, 18)

# Edgecolor & lines
a_edgecolor <- array("darkblue", dim=c(12, 12, 20))
a_edgecolor[wadj_point < 0] <- "red"
a_lty <- array(1, dim=c(12, 12, 20))
a_lty[wadj_point < 0] <- 2

for(tp in tpSelect) {
  qgraph(wadj_point_thresh[, , tp], 
         layout = Q$layout,
         #         edge.color = t(tvvar_obj$edgecolor[, , 1, tp]), 
         labels = mood_labels, 
         vsize = 13, 
         esize = 10,
         asize = 10, 
         mar = c(6, 6, 6, 6), 
         minimum = 0, 
         maximum = .8, 
         edge.color = a_edgecolor[, , tp], 
         lty = a_lty[, , tp])
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
# doubt (11) -> doubt (11)
# satisfied (4) -> relaxed (1) 
# satisfied (4) -> guilty (9) 

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
  P1_pointest <- wadj_point_unthresh[par_row[1], par_row[2], ]
  points((1:20)+v_jitter[i], P1_pointest, col = cols[i], pch = 20, cex = 2) 
  lines((1:20)+v_jitter[i], P1_pointest, col = cols[i], lwd = 2) 
  
  # Center CI
  CI_low_par <- wadj_CIlow[par_row[1], par_row[2], ] 
  CI_up_par <- wadj_CIup[par_row[1], par_row[2], ] 
  
  
  segments((1:20)+v_jitter[i], CI_low_par  , (1:20)+v_jitter[i], CI_up_par, 
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


