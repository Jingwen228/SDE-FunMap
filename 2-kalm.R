
KalmanFilter <- function(parr,par.sigma,par.covar,y,TT,state0,delTtime=1){
  
  sigma <- par.sigma
  r <- parr[1]
  K <- parr[2]
 
  x0 <- t(state0)
  dimX <- length(x0)
  matA <- c(r-2*r*x0/K)
  SIG <- sigma
  
  dimY <- length(SIG)
  
  S <- diag(unlist(par.covar),dimY,dimY)
  tau <- TT[2]-TT[1]
  tmp <-  rbind( cbind(-matA  , SIG*SIG) ,
                 cbind( 0 , matA ))*tau
  Pint  <- matexp(tmp) 
  PHI0   <- t.default(Pint[(dimX+1),(2*dimX),drop=FALSE])
  P0    <- PHI0%*% Pint[1,(2*dimX),drop=FALSE] 
  
  dimT <- length(TT)
  Xp      <- array(NA,c(dimX,dimT))
  Xf      <- array(NA,c(dimX,dimT))
  Yp      <- array(NA,c(dimY,dimT))
  KfGain  <- array(NA,c(dimY,dimT))
  Pf      <- array(NA,c(dimX,dimT))
  Pp      <- array(NA,c(dimX,dimT))
  R       <- array(NA,c(dimY,dimT))
  
  
  Pp[,1] <- P0
  Xp[,1]  <- state0
  
  DIAGDIMY  <- diag(1,dimY)
  SIGSIGT <- SIG*SIG
  ZERODIMXDIMX <- matrix(0,nrow=dimX,ncol=dimX)
  
  
  matASIGSIGTzerosmatAT <-  rbind( cbind(-matA , SIGSIGT ) ,
                                   cbind( ZERODIMXDIMX , matA))
  fy <- 0
  for(k in 1:dimT){
    ObsIndex  <- which(!is.na(y[,k]))
    
    E         <- DIAGDIMY[ObsIndex,,drop=FALSE] 
    Yp[,k] <- E%*%Xp[,k,drop=FALSE]
    R[,k]  <- E%*%Pp[,k]%*%t.default(E) + E%*%S%*%t.default(E)
    
    
    if(length(ObsIndex)>0){ 
      KfGain[,k] <- Pp[,k]%*%t.default(E) %*% solve.default(R[, k]) 
      e       <- as.numeric(y[,k,drop=FALSE])-as.numeric(Yp[,k,drop=FALSE])
      KFg     <- KfGain[,k,drop=FALSE]
      Xf[,k]  <- Xp[,k,drop=FALSE] + KFg%*%e
      Pf[,k] <- Pp[,k] - KFg %*% R[,k] %*% t.default(KFg)
      tmpR        <-  R[,k,drop=FALSE]
      logdetR     <-  determinant.matrix(tmpR)$modulus 
      fy  <-  fy + ((-dimY/2)*log(2*pi)-logdetR/2 +t.default(e)%*%solve.default(tmpR)%*%e/2*(-1) )
      
    } else { 
      Xf[,k] <- Xp[,k]
      Pf[,k] <- Pp[,k]
    } 
    
    if(k==dimT) break
    
    tau   <- TT[k+1]-TT[k]   
    tmp   <- matexp(tau * matASIGSIGTzerosmatAT)
    PHI   <- t.default(tmp[(dimX+1),(2*dimX),drop=FALSE])
    
    IntExpASIG <- PHI %*% tmp[1,(dimX*2),drop=FALSE]
    Pp[,k+1] <- PHI %*% Pf[,k] %*% t.default(PHI) + IntExpASIG
    Xp[,k+1] <- PHI %*% Xf[,k]
    
  }
  
  return(list(Time=TT,Xp=Xp, Xf=Xf, Yp=Yp, KfGain=KfGain,Pf=Pf, Pp=Pp, R=R,fy=fy))
}



