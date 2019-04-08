



tvvarGAM <- function(data, # the n x p data matrix
                     nb = 10,
                     consec,
                     SIMdata=NULL,
                     simulated = FALSE,
                     plot = FALSE,
                     scale = FALSE,
                     beepvar,
                     dayvar,
                     estimates = FALSE,
                     tvvarOpt = "TVVAR",
                     thresholding = FALSE,
                     pbar){
  
  
  #---------- Input check ---------
  
  if(!(is.numeric(data))) stop('Object data has to be numeric.')
  ifelse(simulated==TRUE,ifelse(is.null(SIMdata)==FALSE,ifelse(plot==TRUE,TRUE,TRUE),stop("SIMdata had to be defined")),TRUE)
  if(tvvarOpt=="VAR" & estimates==TRUE & plot==FALSE) stop("You cannot get the estimates of a VAR model with this argument, please put estimates=FALSE")
  if(tvvarOpt=="VAR" & plot==TRUE) stop("With this function you cannot get a plot of the parameters of a VAR model, please put plot=FALSE")
  
  # --------- Fill in defaults ---------
  
  if(missing(pbar)) pbar <- TRUE
  if(missing(consec)) consec <- NULL
  if(missing(beepvar)) beepvar <- NULL
  if(missing(dayvar)) dayvar <- NULL
  
  # --------- Compute Aux Variables ---------
  
  # ----- Compute consec argument -----
  
  
  # Input checks (can only specify consec OR beepvar and dayvar)
  
  if(!is.null(consec) & !is.null(beepvar)) stop("Please specify the consecutiveness of measurements either via consec, or via dayvar and beepvar")
  if(!is.null(consec) & !is.null(dayvar)) stop("Please specify the consecutiveness of measurements either via consec, or via dayvar and beepvar")
  
  if(!is.null(dayvar)) if(is.null(beepvar)) stop("Argument beepvar not specified.")
  if(!is.null(beepvar)) if(is.null(dayvar)) stop("Argument dayvar not specified.")
  
  if(!is.null(beepvar) & !is.null(dayvar)) {
    
    consec <- mgm:::beepday2consec(beepvar = beepvar,
                                   dayvar = dayvar)
    
  }
  else{consec <- 1:nrow(data)}
  # --------- Compute Aux Variables ---------
  
  nt <- nrow(data)
  nv <- ncol(data)
  if(nv==1){if(is.numeric(data)==FALSE) stop('Object data has to be numeric')}
  else{if(all(apply(data[,],2,is.numeric))==FALSE) stop('Object data has to be numeric')}
  
  
  tt=1:nt
  
  # Define colnames, if not provided with data
  if(is.null(colnames(data))) colnames(data) <- paste0("X", 1:nv)
  
  coln<-colnames(data)#Defining the colnames
  if(is.null(coln)) stop("Colnames need to be defined")
  colnL=paste(coln,"L",sep="") #And the lagged colnames
  
  call <- list(data= data,nb=nb,consec=consec,
               SIMdata=SIMdata,
               simulated=simulated,
               plot=plot,
               scale=scale,
               estimates=estimates,
               tvvarOpt=tvvarOpt,
               thresholding=thresholding)
  
  #%##########################################%###
  ####  Part 1: creating the data #############
  #%##########################################%##
  #The functions in this part are created in the file (source code for simulation article 3.R)
  
  y <- data
  
  #%##########################################%###
  ####  Part 2: ESTIMATING GAM##### #############
  #%##########################################%##
  
  mod_all <- tvvarDATA(data = y,
                       nb = nb,
                       pbar = pbar,scale=scale)$model
  
  # --------- Case: estimates = TRUE ---------
  
  if(estimates==TRUE | plot==TRUE) {
    
    
    if(plot==TRUE & simulated==TRUE){
      par(mfrow=c(nv,(nv+1)))
      tt=1:nt
      k=0
      for (i in 1:nv){
        mod<-mod_all[[i]]
        
        
        for ( j in 1:(nv+1))
          
          if(j==1)
          {plot.gam(mod,seWithMean = TRUE,select=1,rug=F,shift=coef(mod)[1],xlab="Time",ylab=paste("Intercept of variable",coln[i],sep=""))
            lines(tt,SIMdata$aint[,i],col="red")
          }
        else {plot.gam(mod,seWithMean = TRUE,select=j,rug=F,ylim=c(-1,1),xlab="Time",ylab=paste(coln[i]," is regressed on ",colnL[j-1],sep=""))
          k=1+k
          
          lines(tt,SIMdata$rho[,k],col="red")
          
        }
      }
      
    }  else if (plot==TRUE) {  par(mfrow=c(nv,(nv+1)))
      for (i in 1:nv){
        mod<-mod_all[[i]]
        
        
        for ( j in 1:(nv+1))
          if(j==1)
          {plot.gam(mod,seWithMean = TRUE,select=1,rug=F,shift=coef(mod)[1],ylab=paste("Intercept of variable",coln[i],sep=""))
          }
        else {plot.gam(mod,seWithMean = TRUE,select=j,rug=F,ylim=c(-1,1),xlab="Time",ylab=paste(coln[i]," is regressed on ",colnL[j-1],sep=""))
          
          
        }
      }
      
    }
    
    
    
    if(estimates==TRUE){
      
      Results_GAM<-array(NA,c(c(nv+1,nv),nt,3))
      
      
      
      
      #-----------------thresholding---------------
      if(thresholding==TRUE){
        for (ii in 1:nv){
          
          mod <- mod_all[[ii]]
          
          mat_dat<-matrix(c(tt,rep(rep(1,nt),nv)),length(tt),nv+1)
          coln_Data<-colnL
          coln_Data_full<-c("tt",coln_Data)
          colnames(mat_dat)<-coln_Data_full
          newd<-as.data.frame(mat_dat)
          Xp=predict(mod,newd,type="lpmatrix",seWithMean = TRUE)
          kdim=dim(Xp)[2]/c(nv+1)
          newpre=predict(mod,new.data=newd,type="terms",se=TRUE)
          Results_GAM[1,ii,1:nt,2]<-Xp[,1:kdim]%*% coef(mod)[1:kdim]#basis functions intercept!
          
          
          Numbrep=5000
          modr<-mvrnorm(Numbrep,coef(mod),mod$Vp+diag((nv+1)*kdim)*10^(-30))
          
          #The confidence intervals
          
          
          
          
          
          int.ci<-matrix(NA,nt,Numbrep)
          for (m in 1:Numbrep){
            int.ci[,m]<-  Xp[,1:kdim]%*%modr[m,1:kdim]
          }
          
          
          
          
          Results_GAM[1,ii,1:nt,1]<-apply(int.ci,1,quantile,c(.975))
          Results_GAM[1,ii,1:nt,3]<-apply(int.ci,1,quantile,c(.025))
          
          for(x in 1:nt){
            ifelse((Results_GAM[1,ii,x,1]< 0 &&  Results_GAM[1,ii,x,3] > 0) || ( Results_GAM[1,ii,x,1] > 0 &&  Results_GAM[1,ii,x,3] < 0)==TRUE,Results_GAM[1,ii,x,2]<-0,Results_GAM[1,ii,x,2])
            
          }
          
          
          for (j in 1:nv){
            Results_GAM[j+1,ii,1:nt,2]<-Xp[,(j*kdim+1):((j+1)*kdim)]%*% coef(mod)[(j*kdim+1):((j+1)*kdim)]
            
            #The confidence intervals
            phi.ci<-matrix(NA,nt,Numbrep)
            for (m in 1:Numbrep){
              phi.ci[,m]<-Xp[,(j*kdim+1):((j+1)*kdim)]%*%modr[m,(j*kdim+1):((j+1)*kdim)]
            }
            Results_GAM[j+1,ii,1:nt,1]<-apply(phi.ci,1,quantile,c(.975))
            Results_GAM[j+1,ii,1:nt,3]<-apply(phi.ci,1,quantile,c(.025))
            
            for(x in 1:nt){
              ifelse((Results_GAM[j+1,ii,x,1]< 0 &&  Results_GAM[j+1,ii,x,3] > 0) || ( Results_GAM[j+1,ii,x,1] > 0 &&  Results_GAM[j+1,ii,x,3] < 0)==TRUE,Results_GAM[j+1,ii,x,2]<-0,Results_GAM[j+1,ii,x,2])
              
            }
            
            
          }
          
          
          
          Results<-list('Estimate'=Results_GAM[, , , 2],'CI_low'=Results_GAM[, , , 3],'CI_high'=Results_GAM[, , , 1])
          
          
          outlist<-list(call=call,Results_GAM=Results,model=mod_all)
          
        }
        
        
        return(outlist)
      }
      
      
      
      
      
      #---------non thresholding ------------------------
      else{
        for (ii in 1:nv){
          
          mod <- mod_all[[ii]]
          
          mat_dat<-matrix(c(tt,rep(rep(1,nt),nv)),length(tt),nv+1)
          coln_Data<-colnL
          coln_Data_full<-c("tt",coln_Data)
          colnames(mat_dat)<-coln_Data_full
          newd<-as.data.frame(mat_dat)
          Xp=predict(mod,newd,type="lpmatrix",seWithMean = TRUE)
          kdim=dim(Xp)[2]/c(nv+1)
          newpre=predict(mod,new.data=newd,type="terms",se=TRUE)
          Results_GAM[1,ii,1:nt,2]<-Xp[,1:kdim]%*% coef(mod)[1:kdim]#basis functions intercept!
          
          
          Numbrep=5000
          modr<-mvrnorm(Numbrep,coef(mod),mod$Vp+diag((nv+1)*kdim)*10^(-30))
          
          
          int.ci<-matrix(NA,nt,Numbrep)
          for (m in 1:Numbrep){
            int.ci[,m]<-  Xp[,1:kdim]%*%modr[m,1:kdim]
          }
          
          Results_GAM[1,ii,1:nt,1]<-apply(int.ci,1,quantile,c(.975))
          Results_GAM[1,ii,1:nt,3]<-apply(int.ci,1,quantile,c(.025))
          
          
          for (j in 1:nv){
            Results_GAM[j+1,ii,1:nt,2]<-Xp[,(j*kdim+1):((j+1)*kdim)]%*% coef(mod)[(j*kdim+1):((j+1)*kdim)]
            
            #The confidence intervals
            phi.ci<-matrix(NA,nt,Numbrep)
            for (m in 1:Numbrep){
              phi.ci[,m]<-Xp[,(j*kdim+1):((j+1)*kdim)]%*%modr[m,(j*kdim+1):((j+1)*kdim)]
            }
            Results_GAM[j+1,ii,1:nt,1]<-apply(phi.ci,1,quantile,c(.975))
            Results_GAM[j+1,ii,1:nt,3]<-apply(phi.ci,1,quantile,c(.025))
          }
          
          
          Results<-list('Estimate'=Results_GAM[, , , 2],'CI_low'=Results_GAM[, , , 3],'CI_high'=Results_GAM[, , , 1])
          
          
          outlist<-list(call=call,Results_GAM=Results,model=mod_all)
          
          
        }
        return(outlist)
        
        
      }
      
      
    }
    
    
    
    
    # Return Estimates
    
    
    
    
    
  }
  else{
    outlist<-list(call=call,model=mod_all)
    return(outlist)
  }
  
}


























tvvarDATA <- function(data, # the n x p data matrix
                      nb = 10,
                      consec,
                      beepvar,
                      dayvar,
                      scale,
                      tvvarOpt="TVVAR",
                      pbar){ #data is the data,nt is the number of time points, nv=number of variables and k is the number of knots
  
  # --------- Fill in defaults ---------
  
  if(missing(pbar)) pbar <- TRUE
  if(missing(pbar)) scale <- FALSE
  if(missing(consec)) consec <- NULL
  if(missing(beepvar)) beepvar <- NULL
  if(missing(dayvar)) dayvar <- NULL
  # --------- Compute Aux Variables ---------
  
  # ----- Compute consec argument -----
  
  
  # Input checks (can only specify consec OR beepvar and dayvar)
  
  if(!is.null(consec) & !is.null(beepvar)) stop("Please specify the consecutiveness of measurements either via consec, or via dayvar and beepvar")
  if(!is.null(consec) & !is.null(dayvar)) stop("Please specify the consecutiveness of measurements either via consec, or via dayvar and beepvar")
  
  if(!is.null(dayvar)) if(is.null(beepvar)) stop("Argument beepvar not specified.")
  if(!is.null(beepvar)) if(is.null(dayvar)) stop("Argument dayvar not specified.")
  
  if(!is.null(beepvar) & !is.null(dayvar)) {
    
    consec <- mgm:::beepday2consec(beepvar = beepvar,
                                   dayvar = dayvar)
    
  }
  else{consec <- 1:nrow(data)}
  # if: specification of consecutiveness via beepvar and dayvar
  
  
  
  nt <- nrow(data)
  nv <- ncol(data)
  tt=1:nt
  
  # Standardize or scale your data
  if(scale==TRUE){data<-scale(data) }
  # Use lagData() from the mgm package
  
  lagD_obj <- mgm:::lagData(data,
                            lags = 1,
                            consec = consec)
  
  
  # Back to the variable names:
  Data2 <- cbind(lagD_obj$data_response, lagD_obj$l_data_lags[[1]])
  Data2 <- rbind(rep(NA, ncol(Data2)), Data2) # to make the code below work
  Data2 <- as.data.frame(Data2)
  coln<-colnames(data)
  colnL=paste(coln,"L",sep="")
  colnames(Data2)=c(coln,colnL)
  Data1 <- Data2
  tt <- 1:(dim(Data1)[1])
  
  
  
  
  # --------- Case: tvvarOpt = TVVAR ---------
  if(tvvarOpt=="TVVAR") {
    allcol2=c()
    for(i in 1:nv){allcol2[i]=paste("s(tt,by=",colnL[i],",k=nb",")",sep="")}
    allcol3=paste(allcol2,collapse="+")
    
    
    # Progress bar
    if(pbar==TRUE) pb <- txtProgressBar(min = 0, max=nv, initial=0, char="-", style = 3)
    
    model=list()
    for (j in 1:nv){
      
      ff <- as.formula(paste(coln[j]," ~ ","s(tt,k=nb)","+",allcol3))
      model[[j]]<-gam(ff,data=Data1,seWithMean=TRUE)
      
      # Update Progress Bar
      if(pbar==TRUE) setTxtProgressBar(pb, j)
      
      
    }
    outlist<-list(model=model)
    
    
    
    
    return(outlist)
    
  }
  
  # --------- Case: tvvarOpt = "VAR" and thus a standard VAR is estimated ---------
  
  else{
    allcol2=c()
    for(i in 1:nv){allcol2[i]=paste(colnL[i],sep="")}
    allcol3=paste(allcol2,collapse="+")
    
    
    # Progress bar
    if(pbar==TRUE) pb <- txtProgressBar(min = 0, max=nv, initial=0, char="-", style = 3)
    
    model=list()
    for (j in 1:nv){
      
      ff <- as.formula(paste(coln[j]," ~ ","tt","+",allcol3))
      model[[j]]<-gam(ff,data=Data1,seWithMean=TRUE)
      
      
      # Update Progress Bar
      if(pbar==TRUE) setTxtProgressBar(pb, j)
      
      
    }
    
    
    
    outlist<-list(model=model)
    
    
    
    
    return(outlist)
    
  }
}















cormatrix2=function(x,nv){
  CC.int=matrix(x,nv,nv,byrow=TRUE)
  CC.int=CC.int+t(CC.int)
  diag(CC.int)<-1
  return(CC.int)
}
#five different options
##1. The three generating functions (genfun) options##
#1A. Time invariant
invariant<-function(nt,MaxAbsValue)  #nt is the sample size, and MaxAbsValue is the maximum absolute value of the function.
{genfun=rep(NA,(nt)) #Creating the function first with NAs (missing values).
genfun=rep(MaxAbsValue,(nt)) #Here the actual invariant function is created.
return(genfun)
}

#1B. Linear
linear<-function(nt,MaxAbsValue)  #nt is the sample size, and MaxAbsValue is the maximum absolute value of the function.
{genfun=rep(NA,(nt)) #Creating the function first with NAs (missing values).
genfun=seq(0,MaxAbsValue,length.out=(nt))  #Here the actual linear function is created.

return(genfun)
}

#1C. Sine
sine<-function(nt,MaxAbsValue) #nt is the sample size, and MaxAbsValue is the maximum absolute value of the function.
{genfun=rep(NA,(nt)) #Creating the function first with NAs (missing values).
tt=1:(nt) #Defining a time parameter in order to create the sine function
genfun=MaxAbsValue*sin(2*pi*tt/(nt)) #Here the actual sine function is created.
return(genfun)
}


choose.coefaint<-function(FUN,nt,MaxAbsValue,nv){
  FUNchoose=c(invariant,linear,sine)
  aint=matrix(NA,(nt),nv)
  for(i in 1:nv){
    aint[,i]=FUNchoose[[FUN[i]]](nt,MaxAbsValue[i])}
  return(aint)
}

#in order to be able to start stationairity check
stat.check=function(nt,rho){
  WAAR=c()
  duizend=rep(NA,nt)
  for (t in 1:(nt)){
    duizend[t]=ifelse (max(abs(eigen(matrix(rho[t,],nv,nv))$values))<1,1,0)
  }
  WAAR=all(duizend==1)
  return(WAAR)
}


choose.coefaint<-function(FUN,nt,MaxAbsValue,nv){
  FUNchoose=c(invariant,linear,sine)
  aint=matrix(NA,nt,nv)
  for(i in 1:nv){
    aint[,i]=FUNchoose[[FUN[i]]](nt,MaxAbsValue[i])}
  return(aint)
}




choose.coefrho<-function(FUNr,nt,MaxAbsValue,nv){
  FUNchoose=c(invariant,linear,sine)
  rho=matrix(NA,nt,nv*nv)
  for(i in 1:(nv*nv)){
    rho[,i]=FUNchoose[[FUNr[i]]](nt,MaxAbsValue[i])
  }
  
  s.check=stat.check(nt=nt,rho=rho)
  conv=FALSE
  while (!conv){
    if (s.check==TRUE) {conv=TRUE}
    else{
      for(i in 1:(nv*nv)){
        rho[,i]=FUNchoose[[FUNr[i]]](nt,MaxAbsValue[i])
      }
      s.check=stat.check(nt=nt,rho)
    }}
  return(rho)
}



tvvarSIM<-function(aint,rho,nt,nv){
  aint=aint
  sigma2 <- cov.mat<-cormatrix2(0.1,nv)
  y=matrix(0,nt,nv)
  vecsigma2=matrix(sigma2,nv*nv,1)
  pphi=kronecker(matrix(rho[1,],nv,nv,byrow=T),matrix(rho[1,],nv,nv,byrow=T))
  matrix(solve(diag(nv*nv)-pphi)%*%vecsigma2,nv,nv)
  
  y[1,]<-rmvnorm(1,mean=solve(diag(nv)-matrix(rho[1,],nv,nv,byrow=T))%*%matrix(aint[1,],nv,1),sigma=matrix(solve(diag(nv*nv)-pphi)%*%vecsigma2,nv,nv))
  
  for (t in 2:nt){
    y[t,]=matrix(aint[t,],nv,1)+matrix(rho[t-1,],nv,nv,byrow=T)%*%y[t-1,]+matrix(rmvnorm(1,sigma=sigma2),nv,1) #t=>1  # u should be a martingal difference noise
    
  }
  
  colnames(y)=paste("y",1:nv,sep="")
  return(list(y=y,aint=aint,rho=rho))
}
