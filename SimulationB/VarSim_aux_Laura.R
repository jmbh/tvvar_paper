####  Part 1: creating the data #############
#%##########################################%##
#######Bivariate process #################

nvar=2
cormatrix2=function(x,np){
  CC.int=matrix(x,np,np,byrow=TRUE)
  CC.int=CC.int+t(CC.int)
  diag(CC.int)<-1
  return(CC.int)   
}
#five different options
##1. The five generating functions (genfun) options##
#1A. Time invariant 
invariant<-function(N,MaxAbsValue)  #N is the sample size, and MaxAbsValue is the maximum absolute value of the function.  
{genfun=rep(NA,(N)) #Creating the function first with NAs (missing values).
genfun=rep(MaxAbsValue,(N)) #Here the actual invariant function is created.
return(genfun)
}

#1B. Linear
linear<-function(N,MaxAbsValue)  #N is the sample size, and MaxAbsValue is the maximum absolute value of the function.  
{genfun=rep(NA,(N)) #Creating the function first with NAs (missing values).
genfun=seq(0,MaxAbsValue,length.out=(N))  #Here the actual linear function is created.    

return(genfun)
}

#1C. Sine 
sine<-function(N,MaxAbsValue) #N is the sample size, and MaxAbsValue is the maximum absolute value of the function.  
{genfun=rep(NA,(N)) #Creating the function first with NAs (missing values).
tt=1:(N) #Defining a time parameter in order to create the sine function
genfun=MaxAbsValue*sin(2*pi*tt/(N)) #Here the actual sine function is created.              
return(genfun)
}


choose.coefaint<-function(FUN,N,MaxAbsValue,np){  
  FUNchoose=c(invariant,linear,sine)
  aint=matrix(NA,(N),np)
  for(i in 1:np){
    aint[,i]=FUNchoose[[FUN[i]]](N,MaxAbsValue[i])}
  return(aint)
}

#in order to be able to start stationairity check
stat.check=function(N,rho){
  WAAR=c()
  duizend=rep(NA,N)  
  for (t in 1:(N)){
    duizend[t]=ifelse (max(abs(eigen( matrix(rho[t,],np,np))$values))<1,1,0)
  }
  WAAR=all(duizend==1)  
  return(WAAR) 
}


choose.coefaint<-function(FUN,N,MaxAbsValue,np){  
  FUNchoose=c(invariant,linear,sine)
  aint=matrix(NA,N,np)
  for(i in 1:np){
    aint[,i]=FUNchoose[[FUN[i]]](N,MaxAbsValue[i])}
  return(aint)
}




choose.coefrho<-function(FUNr,N,MaxAbsValue,np){  #choose one of the above intercept simulations
  FUNchoose=c(invariant,linear,sine)
  rho=matrix(NA,(N),(np*np))
  for(i in 1:(np*np)){
    rho[,i]=FUNchoose[[FUNr[i]]](N,MaxAbsValue[i])
  }
  
  s.check=stat.check(N=N,rho=rho)
  conv=FALSE   
  while (!conv){
    if (s.check==TRUE) {conv=TRUE} 
    else{
      for(i in 1:(np*np)){
        rho[,i]=FUNchoose[[FUNr[i]]](N,MaxAbsValue[i])
      }
      s.check=stat.check(N=N,rho)
    }}
  return(rho)
}



creat.y<-function(aint,rho,N,np){
  aint=aint
  sigma2 <- cov.mat<-cormatrix2(0.1,np)
  y=matrix(0,N,np)
  vecsigma2=matrix(sigma2,np*np,1)
  pphi=kronecker(matrix(rho[1,],np,np,byrow=T),matrix(rho[1,],np,np,byrow=T))
  matrix(solve(diag(np*np)-pphi)%*%vecsigma2,np,np)
  
  y[1,]<-mvrnorm(1,mean=solve(diag(np)-matrix(rho[1,],np,np,byrow=T))%*%matrix(aint[1,],np,1),sigma=matrix(solve(diag(np*np)-pphi)%*%vecsigma2,np,np))
  
  for (t in 2:N){
    y[t,]=matrix(aint[t,],np,1)+matrix(rho[t-1,],np,np,byrow=T)%*%y[t-1,]+matrix(mvrnorm(1,sigma=sigma2),np,1) #t=>1  # u should be a martingal difference noise
    
  }
  
  return(list(y=y,aint=aint,rho=rho))
}


#estimating function
Estimate_gam<-function(y,N,np,K){ #y is the y variable,N is the number of time points, np=number of variables and k is the number of knots
  tt=1:N
  
  data1=matrix(0,N,(np*2),byrow=T)
  for (h in (np+1):(np*2)){
    data1[,(h-np)]=y[,(h-np)]# data wordt hier gelagged    
    
    data1[,h]=c(NA,y[1:(N-1),(h-np)])# data wordt hier gelagged    
  }
  
  data1=as.data.frame(data1)
  colnames(data1)=c(paste("y",1:np,sep=""),paste("y",1:np,"L",sep=""))
  
  coln=colnames(data1)[1:np]
  colnL=colnames(data1)[(np+1):(np*2)]
  allcol2=c()
  for(i in 1:np){allcol2[i]=paste("s(tt,by=",colnL[i],",k=K",")",sep="")}
  allcol3=paste(allcol2,collapse="+")
  
  
  model=list()
  for (j in 1:np){
    ff <- as.formula(paste(coln[j]," ~ ","s(tt,k=K)","+",allcol3))
    model[[j]]<-gam(ff,data=data1,seWithMean=TRUE)
    
  }
  
  return(model)
  
  }




Estimate_var<-function(y,
                       N,
                       np){ #x is the index, y is the y variable
  tt=1:N
  
  data1=matrix(0,N,(np*2),byrow=T)
  for (h in (np+1):(np*2)){
    data1[,(h-np)]=y[,(h-np)]# data wordt hier gelagged    
    
    data1[,h]=c(NA,y[1:(N-1),(h-np)])# data wordt hier gelagged    
  }
  
  data1=as.data.frame(data1)
  colnames(data1)=c(paste("y",1:np,sep=""),paste("y",1:np,"L",sep=""))
  
  coln=colnames(data1)[1:np]
  colnL=colnames(data1)[(np+1):(np*2)]
  
  allcol3=paste(colnL,collapse="+")
  
  
  model=list()
  for (j in 1:np){
    ff <- as.formula(paste(coln[j]," ~ ",allcol3))
    model[[j]]<-lm(ff,data=data1)
    
  }
  
  return(model)
  
  
}



GAM_TVVAR<-function(Data,np,N,K){
  Results_GAM<-array(NA,c(c(np+1,np),N,3))
  
  
  
  
  N=N
  #N=100 # 100 500 
  tt=1:N
  
  #%##########################################%###
  ####  Part 1: creating the data #############
  #%##########################################%##   
  #The functions in this part are created in the file (source code for simulation article 3.R)
  
  
  
  y<-Data
  
  
  #%##########################################%###
  ####  Part 2: ESTIMATING GAM##### #############
  #%##########################################%##   
  
  for (ii in 1:np){
    
    mod<-Estimate_gam(y,N,np,K)[[ii]] 
    
    mat_dat<-matrix(c(tt,rep(rep(1,N),np)),length(tt),np+1)
    coln_data<-paste("y",1:np,"L",sep="")
    coln_data_full<-c("tt",coln_data)
    colnames(mat_dat)<-coln_data_full
    newd<-as.data.frame(mat_dat)
    Xp=predict(mod,newd,type="lpmatrix",seWithMean = TRUE)
    kdim=dim(Xp)[2]/c(np+1) 
    newpre=predict(mod,new.data=newd,type="terms",se=TRUE) 
    Results_GAM[1,ii,1:N,2]<-Xp[,1:kdim]%*% coef(mod)[1:kdim]#basis functions intercept!
    
    Numbrep=1000
    modr<-mvrnorm(Numbrep,coef(mod),mod$Vp+diag((np+1)*kdim)*10^(-14))
    
    #The confidence intervals
    int.ci<-matrix(NA,N,Numbrep)
    for (m in 1:Numbrep){
      int.ci[,m]<-  Xp[,1:kdim]%*%modr[m,1:kdim]
    }
    
    Results_GAM[1,ii,1:N,1]<-apply(int.ci,1,quantile,c(.975)) 
    Results_GAM[1,ii,1:N,3]<-apply(int.ci,1,quantile,c(.025)) 
    
    
    for (j in 1:np){
      Results_GAM[j+1,ii,1:N,2]<-Xp[,(j*kdim+1):((j+1)*kdim)]%*% coef(mod)[(j*kdim+1):((j+1)*kdim)]
      
      #The confidence intervals
      phi.ci<-matrix(NA,N,Numbrep)
      for (m in 1:Numbrep){ 
        phi.ci[,m]<-Xp[,(j*kdim+1):((j+1)*kdim)]%*%modr[m,(j*kdim+1):((j+1)*kdim)]  
      }
      Results_GAM[j+1,ii,1:N,1]<-apply(phi.ci,1,quantile,c(.975))
      Results_GAM[j+1,ii,1:N,3]<-apply(phi.ci,1,quantile,c(.025))
    }
    
    
    
    
  }
  return(Results_GAM=Results_GAM)
}






# ---------- Addition Jonas; Feb 24th ----------


Estimate_var_adj <- function(y)
  
  
{ #x is the index, y is the y variable
  
  N <- nrow(y)
  np <- ncol(y)
  
  tt=1:N
  
  data1=matrix(0,N,(np*2),byrow=T)
  for (h in (np+1):(np*2)){
    data1[,(h-np)]=y[,(h-np)]# data wordt hier gelagged    
    
    data1[,h]=c(NA,y[1:(N-1),(h-np)])# data wordt hier gelagged    
  }
  
  data1=as.data.frame(data1)
  colnames(data1)=c(paste("y",1:np,sep=""),paste("y",1:np,"L",sep=""))
  
  coln=colnames(data1)[1:np]
  colnL=colnames(data1)[(np+1):(np*2)]
  
  allcol3=paste(colnL,collapse="+")
  
  
  model=list()
  for (j in 1:np){
    ff <- as.formula(paste(coln[j]," ~ ",allcol3))
    model[[j]]<-lm(ff,data=data1)
    
  }
  
  # Collapse coefficients into matrix
  mat_ret <- do.call(rbind, lapply(model, function(x) x$coefficients))
  
  return(mat_ret)
  
}











