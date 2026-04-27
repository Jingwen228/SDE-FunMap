f.load <-function(geno="Map-Genotype.csv",
                  pheno1="height1.csv",
                  time="Taproot-Height-Time.csv"){
  ph1 <- read.csv(pheno1)[,-1]
  geno1 <- read.csv(geno)
  Time <- read.csv(time)[,-1]
  marker.info <- geno1[,1:3]
  marker <- geno1[,-c(1:3)]
  pname <- paste("X",ph1[,1],sep="")
  mname <- colnames(marker)
  rr <- intersect(pname,mname)
  
  index1 <- c();
  for(i in 1:length(rr)){
    index1 <-c(index1,which(pname==rr[i]))
  }
  
  pf1 <- ph1[index1,]
  
  TT <- Time[index1,]
  
  markerf <- as.matrix(marker[,index1])
  
  markerf[which(markerf=="--")] <- 9
  markerf[which(markerf=="hh")] <- 2
  markerf[which(markerf=="hk")] <- 1
  markerf[which(markerf=="kk")] <- 0
  markerf[which(markerf=="ll")] <- 2
  markerf[which(markerf=="lm")] <- 1
  markerf[which(markerf=="nn")] <- 2
  markerf[which(markerf=="np")] <- 1

  lower.par=c(0.1,0.001,-Inf,0.01,8.6);
   upper.par=c(1.5,4.218195,Inf,0.2,11)
  
  
  Times <- round(colMeans(TT),0)
  list(markerf=markerf,pf1=pf1[,-1],TT=TT,Times=Times,marker.info=marker.info,lower.par=lower.par,upper.par=upper.par)
}



