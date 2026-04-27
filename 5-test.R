rm(list = ls())

require("PSM")
library(mvtnorm)
library(PSM)
require("deSolve")
library(deSolve)

source("1-sde.load.R")
source("2-kalm.R")
source("3-sde.curve.R")
source("4-sde.est.R")



dat <- f.load(geno="Map-Genotype.csv",
              pheno1="height1.csv",
              time="Taproot-Height-Time.csv")

ret.H0 <-fun.H0(dat)
ret.H1 <- test(dat,interval=c(1,5),Time)
