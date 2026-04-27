test <- function(dat,interval=c(1,8305),Time){
 
  phen <- dat$pf1
  Time <- dat$Times
  geno <- dat$markerf
 
  
  nm <- dim(geno)[1]
  n1 <- interval[1]
  if(length(interval)==1){
    n2 <- interval[1]
  }else{ 
    n2 <- interval[2]}
  
  if(n2 >=nm){
    n2 <- nm
  }
  res <- matrix(NA,nrow=length(c(n1:n2)),ncol=18)
  
  for(i in n1:n2){
    
    SNP <- as.numeric(geno[i,])      
    NSNP <- as.character(unlist(matrix(SNP,1)))
    missing <- which(NSNP==9)
    if ( length(missing) > 0)
    {
      SNP1 <- NSNP[-(missing)]
      yy<- phen[-(missing), ]
    }else{
      SNP1 <- NSNP                                 
      yy <- phen
    }
    
    ndat <- dat
    ndat$pf1 <- yy
    
    h01 <- try(fun.H0(ndat),TRUE)
    
    if(class(h01) == "try-error") next;
    
    h02 <- try(fun.para(y11=yy,SNP1,init.par=as.numeric(h01$para),Time),TRUE)
    
    if (class(h02) == "try-error")
      h02 <- NA
    
    LR <- 2*(h01$value-h02$value)
    if(is.na(h01)||is.na(h02)){
      allpar <- c(LR,rep(NA,25))
    }else{
      allpar <- c(LR,h01$value,h02$value,h02$para)
      
    }
    
    cat("snp", i, "=", allpar, "\n");
    
   res[i-(n1-1),(1:length(allpar))] <- allpar
    
  }
  return(res)
}


fun.para <- function(y11,SNP1,init.par,Time){
  
  polar.Y <- y11
  index <- table(SNP1)
  Marker.type <- names(index)
  
  lower.x <- c()
  upper.x <- c()
  g.par <- c()
  Marker.index <- list()
  for(i in 1:length(Marker.type)){
    
    Marker.n <- which(SNP1==Marker.type[i])
    Marker.index[[i]] <- Marker.n
    g.par <- c(g.par,init.par[-c(1:3)])
    lower.x <- c(lower.x,dat$lower.par[-c(1:3)])
    upper.x<- c(upper.x,dat$upper.par[-c(1:3)])
  } 
  
  parinx <- c(init.par[1:3],g.par)
  lower <- c(dat$lower.par[1:3],lower.x)
  upper <- c(dat$upper.par[1:3],upper.x)
  rr <-  optim(parinx,cm.mle.fun ,polar.Y=polar.Y,Time=Time,Marker.index=Marker.index,ng=length(index), method="L-BFGS-B",
               lower = lower,upper = upper, control=list(maxit=5000))
  
  return(list(para=rr$par,value=rr$value))
  
}

cm.mle.fun <- function(newpar,polar.Y,Marker.index,Time,ng,lower,upper){
  
  len.cov <- 2
  
  par.covar1 <-  newpar[1:len.cov]
  par.sigma <- newpar[3]
  
  sigma <- SAD1_get_matrix( par.covar1, Time, options=list())
  s.par <-  newpar[-c(1:3)]
  
  len <- 0
  len.len <- 2
  A1 <- c()
  for(i in 1:ng){
    
    par.curve <- s.par[(len+1):(len+len.len)]
    
    y  <- polar.Y[(Marker.index[[i]]),]
    state0<- as.numeric(colMeans(y))[1]
    N  <- dim(y)[1];  M  <- dim(y)[2]
    my <- matrix(NA,N,M)
    
    for(j in 1:N){
      ny   <- matrix(y[j,],1,byrow = TRUE)
      my[j,] <- KalmanFilter(par.curve,par.sigma=par.sigma,par.covar=par.covar1[2],y=ny,TT=Time,state0,delTtime=1)$Yp
    }
    
    Y.delt <- y-my
    pv <- dmvnorm(Y.delt, rep(0, NCOL(sigma)), sigma,log = T)
    
    A1 <- c(A1,-sum(pv))
    len <- len +len.len 
    
  }
  A <- sum(A1)
  return(A);
  
}


