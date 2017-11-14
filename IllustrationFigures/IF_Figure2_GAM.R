
figDir_paper <- "" # Specify directory in which figure should be saved

library(mgcv)

width <- 12
pdf(paste0(figDir_paper, 'Illustr_GAM.pdf'), width = width, height = width/2)


# --------- Set up plotting area ----------

par(mfrow=c(2,4),
    oma = c(0, 0, 0, 0),
    mar=c(4.5,6,1,1))


# --------- Code adapted from Laura ----------

set.seed("3021")
nT=30
P=20
t=seq(1:nT)
y=sin(2*pi*t/P)+rnorm(nT)*.3
#t=x
#y=wear
tt=1:nT

k=5

cex.lab <- 2


mod=gam(y~s(t,bs="tp",k=k))

plot(y,ylab=bquote(paste(y)),cex.lab=1.4,
     ylim = c(-1.5, 1.5),
     xlab = '(a)', 
     cex.lab = cex.lab, xaxt='n',pch=20,cex=1.5)


cc=matrix(mod$coefficients)# the estimated alpha coefficients
tt=seq(1,nT,.01)# values for prediction
newd=data.frame(t=tt)# transform to a data.frame
Xp=predict(mod,newd,type="lpmatrix") # matrix containing the values of the linear predictor

k=dim(Xp)[2]# shows that the prediction matrix Xp has indeed 6 basis functions
j=0

v_xlabs <- c('(b)', '(c)', '(d)', '(e)', '(f)')
for (i in c(1,5,2,4,3)){
  j=j+1
  plot(tt,Xp[,i],type="l", ylab=bquote(paste(R[.(j)](t))),cex.lab=1.4, 
       xlab = v_xlabs[j], 
       cex.lab = cex.lab, xaxt='n',ylim = c(-1.5, 1.5))
  
}

plot(tt,Xp[,1]*cc[1],cex.lab=1.4,type="l",ylab=bquote(paste("",paste(alpha[(k)]),R[(k)](t))),
     xlab = '(g)', 
     cex.lab = cex.lab, xaxt='n',ylim = c(-1.5, 1.5))
for (i in 1:k){
  lines(tt,Xp[,i]*cc[i],col="grey42")
}# Plot weigthed basis functions 

plot(mod,rug=FALSE,shift=coef(mod)[1],ylab=bquote(paste(hat(beta)[10][t])),cex.lab=1.4,
     xlab = '(h)', 
     cex.lab = cex.lab, xaxt='n',ylim = c(-1.5, 1.5))
points(t,y,pch=20,cex=1.5)



dev.off()


