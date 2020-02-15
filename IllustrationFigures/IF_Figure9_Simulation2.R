# jonashaslbeck@gmail.com; September 2019

figDir_paper <- "/" # Define !!!

p <- 6
m <- matrix(0, p, p)
m[upper.tri(m)] <- 1
diag(m) <- 1

library(qgraph)

pdf(paste0(figDir_paper,"VS_UT_SystemFigure.pdf"), width = 4, height = 4)
qgraph(m, 
       layout="circle", 
       diag = TRUE, 
       vsize = 14,
       esize= 7, 
       asize=7,
       mar = rep(7, 4))
dev.off()





