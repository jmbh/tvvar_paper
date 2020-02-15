# jonashaslbeck@gmail.com; November 2017

library(qgraph)
library(mgm)

codeDir <- '' # Specify directory of code files
figDir <- "" # Specify directory in which Figures should be saved

source(paste0(codeDir, 'VarSim_aux.R'))


# ---------- Figure: Illustrating kernel-smoothing Idea ----------

# ----- Part A: weight plots -----
sd1 <- .2
sd2 <- .05

# Define weight functions
x <- seq(.1,1, length = 10000)
thin <- seq(1, 10000, length = 10)
e1 <- dnorm(x, .2, sd1)
e1n <- e1/max(e1)
e1s <- dnorm(x, .2, sd2)
e1sn <- e1s/max(e1s)
e2 <- dnorm(x, .3, sd1)
e2n <- e2/max(e2)
e2s <- dnorm(x, .3, sd2)
e2sn <- e2s/max(e2s)


# Plotting

pdf(paste0(figDir_paper, 'KS_illu_partA.pdf'), width = 5, height = 4)

plot.new()
par(mar=c(4,3.5,1,1))
plot.window(xlim=c(.1,1),
            ylim=c(0,1))


# box()
axis(1, seq(.1, 1, length = 10), 1:10)
axis(2, c(0, .5, 1), las=2)

title(xlab = 'Time points', line = 2.5)
title(ylab = bquote(paste("w"["t"["e"]],""["=3"] )), line=2.5)

# plot weight densities
lwd <- 2
lines(x, e2n, col='tomato', lwd = lwd)
lines(x, e2sn, col='lightblue', lty=1, lwd = lwd)
points(x[thin], e2n[thin], col='tomato', lwd = lwd, pch=20, cex=1.5)
points(x[thin], e2sn[thin], col='lightblue', lty=1, lwd = lwd, pch=20, cex=1.5)

text(.55, .85, "b = 0.2", col = "tomato", cex=1.5)
text(.55, .1, "b = 0.05", col = "lightblue", cex=1.5)

dev.off()

# For the latex part:
round(e2n[thin], 2)
round(e2sn[thin], 2)



# ----- Part B: The densities alone (to paste into Latex thingy) -----

pdf(paste0(figDir, 'KS_illu_partB.pdf'), width = 5, height = 2)
plot.new()
par(mar=c(4,4,1,1))
plot.window(xlim=c(.1,1),
            ylim=c(0,1))

lines(x, e2n, col='tomato', lwd = lwd)
lines(x, e2sn, col='lightblue', lty=1, lwd = lwd)
dev.off()
































