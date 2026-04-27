fun.H0 <- function(dat){
  
  phe.Y <- as.matrix(dat$pf1)
  ym <- as.numeric(colMeans(phe.Y))
  Time <- dat$Times
  
  parinx <- c(1.200133582,0.02478588,0.01 ,0.125995,10)
  rr <- optim(parinx,curve.mle,yy=phe.Y,time.std=Time,x0=ym[1],method="L-BFGS-B",
              lower = dat$lower.par,upper = dat$upper.par, control=list(maxit=1000))
  
  return(list(para=rr$par,value=rr$value))
}

curve.mle <- function(par,yy,time.std,x0){
  
  len.cov <- 2
  
  par.covar1 <-  par[1:len.cov]
  
  sigma <-SAD1_get_matrix( par.covar1, time.std, options=list())
  par.sigma <- par[3]
  par.curve <- par[-c(1:3)]
  
  N  <- dim(yy)[1];M <- dim(yy)[2]
  my <- matrix(NA,N,M)
  
  for(j in 1:N){
    ny <- matrix(yy[j,],1,byrow = TRUE)
    my[j,] <- KalmanFilter(par.curve,par.sigma=par.sigma,par.covar=par.covar1[2],y=ny,TT=time.std,state0=x0,delTtime=1)$Yp
    
  }
  
  Y.delt <- yy-my
  pv <- dmvnorm(Y.delt, rep(0, NCOL(sigma)), sigma,log = T)
  A <- -sum(pv);
  return(A)
}

SAD1_get_matrix <- function(par, times, options=list()){
  
  n <- ifelse(is.vector(times), length(times), NCOL(times))
  phi<- par[1]
  v2 <- par[2]
  tmp <- (1-phi^2)
  
  sigma <- array(1, dim=c(n,n))
  for(i in 1:n)
  {
    sigma[i,i:n] <- phi^( c(i:n) - i ) * (1-phi^(2*i ))/tmp
    sigma[i:n,i] <- sigma[i,i:n]
  }
  
  sigma <- sigma*abs(v2)
  return(sigma)
} 

get_con_param<-function(parm.id)
{
  for (e in commandArgs()) 
  {
    ta = strsplit(e,"=", fixed=TRUE);
    if(! is.na( ta[[1]][2])) 
    {
      temp = ta[[1]][2];
      if( ta[[1]][1] == parm.id)
        return (temp );
    }
  }
  return(NA);
}




