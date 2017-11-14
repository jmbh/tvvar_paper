# jonashaslbeck@gmail.com; November 2017

p <- 6
m <- matrix(0, p, p)
m[upper.tri(m)] <- 1
diag(m) <- 1

library(qgraph)

pdf("VS_UT_SystemFigure.pdf", width = 4, height = 4)
qgraph(m, 
       layout="circle", 
       diag = TRUE)
dev.off()





