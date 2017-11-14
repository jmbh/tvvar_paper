# jonashaslbeck@gmail.com; November 2017

library(qgraph)
library(mgm)

codeDir <- '' # Specify directory of code files
figDir <- "" # Specify directory in which Figures should be saved

source(paste0(codeDir, 'VarSim_aux.R'))

# ---------- 1) Figure: Conceptual Example ----------

# Define Graph
p <- 3
Gw1 <- Gw2 <- Gw3 <- matrix(0, p, p)

r <- -.1
s <- .5

Gw1[1,1] <- .3 + r
Gw1[2,2] <- .4 + r
Gw1[3,3] <- .5 + r

Gw2[1,1] <- .3 + r
Gw2[2,2] <- .4 + r
Gw2[3,3] <- .5 + r
Gw2[3,1] <- .55 + s

Gw3[1,1] <- .3 + r
Gw3[2,2] <- .4 + r
Gw3[3,3] <- .5 + r
Gw3[2,1] <- .3 + s

l_G <- list(Gw1, Gw2, Gw3)

n_week <- 25

th <- matrix(0, nrow = n_week*3, ncol = p)
sd <- matrix(1, nrow = n_week*3, ncol = p)

tvGraph <- array(dim=c(p, p, n_week*3))
tvGraph[,, 1:n_week] <- Gw1
tvGraph[,, (n_week+1):(n_week*2)] <- Gw2
tvGraph[,, (n_week*2+1):(n_week*3)] <- Gw3


# Sample Data
data <- VARData(graph = tvGraph,
                   th = th,
                   sd = sd,
                   seed = 1)



# ----- Plot all Together -----

pdf(paste0(figDir, 'Fig1_illu_new.pdf'), width = 12, height = 8)

# Set up layout
mat <- rbind(c(1,1,1, 7),
             c(2,2,2, 7),
             c(3,3,3, 7),
             c(4,5,6, 8))

lo <- layout(mat, heights = c(1,1,1,2))
# layout.show(lo)

# Plot: Raw Time Series
for(i in 1:3) {
plot.new()
par(mar=rep(0,4))
plot.window(ylim = c(-4.2, 4.2), xlim=c(1,(n_week*3)))
# box()
lines(data[,i])
abline(v=c(n_week+.5, n_week*2+.5), lty = 2)
}


# Plot: Graphs
layout <- rbind(c(-.5, .5),
                c(0, .7),
                c(-.2, -.8))

vsize <- 25
esize <- 14
asize <- 15
m <- 18
labels <- c('Mood', 'Anxiety', 'Worrying')

for(i in 1:3) {
qgraph(l_G[[i]],
       layout=layout,
       mar = rep(m, 4),
       vsize = vsize,
       esize = esize,
       asize = asize,
       labels = labels,
       edge.color = 'black', 
       directed = TRUE, 
       minimum = 0, 
       maximum = 1.05)
title(paste0('Week ', i), line = -3, cex.main = 2)
}


# Plot Legend and y-axis
cex <- 2
plot.new()
plot.window(xlim=c(0,1), ylim=c(0,4))
x_left <- 0
text(x_left, 3.4, 'Depressed Mood', adj = 0, cex = cex)
text(x_left, 2, 'Anxiety', adj = 0, cex = cex)
text(x_left, .6,'Worrying', adj = 0, cex = cex)

# arrows(.05, 0, .05, 4, lwd=2, code=3)


# Plot Average Model

Gav <- (Gw1 + Gw2 + Gw3) / 3
qgraph(Gav,
       layout=layout,
       mar = rep(m, 4),
       vsize = vsize,
       esize = esize,
       asize = asize,
       labels = labels,
       edge.color = 'black', 
       minimum = 0, 
       maximum = 1.05)
title('Average', line = -3, cex.main = 2, col.main='black')



dev.off()

