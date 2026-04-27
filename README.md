SDE-FunMap is a model for mapping quantitative trait loci (QTLs) underlying complex traits based on stochastic differential equations (SDE). With hypothesis testing, this model can identify biologically significant QTLs associated with primary root growth that cannot be detected by deterministic models.

the main progress is as follows:


rm(list = ls()) require("PSM") library(mvtnorm) library(PSM) require("deSolve") library(deSolve)
source("1-sde.load.R") source("2-kalm.R") source("3-sde.curve.R") source("4-sde.est.R")


# load data and function

dat <- f.load(geno="Map-Genotype.csv",
              pheno1="height1.csv",
              time="Taproot-Height-Time.csv")

# Null hypothesis and Alternative hypothesis

ret.H0 <-fun.H0(dat)
ret.H1 <- test(dat,interval=c(1,5),Time)
