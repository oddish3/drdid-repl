#---------------------------------------------------
#      Standard and Stabalized-weighted ATT of DID
#      Repeated Cross Section case
#      All DID estimators
#---------------------------------------------------
#-----------------------------------------------------------------------------
# Startup - clear memory, load packages, and set parameters
# Clear memory
rm(list=ls())
#-----------------------------------------------------------------------------
# Basic parameters for the simulation - Doesn't change over setups
#nrep    <- 1000        # Monte Carlo replications
ncores  <- 14           # Number of cores to use in parallel
seed1   <- 1234         # Set initial seed (guaranteed reproducibility)
Xsi.ps <- .75           # pscore index (strength of common support)
nrep <- 10000           # Monte Carlo replications
#-----------------------------------------------------------------------------
# Set the Working Directory
address <- "set/here/your/working/directory/for/rc"
setwd(address)
#-----------------------------------------------------------------------------
# load the necessary libraries
# Install DRDID package
devtools::install_github("pedrohcgs/DRDID")
library(foreach)
library(doSNOW)
library(doRNG)
library(trust)
library(DRDID)
#-----------------------------------------------------------------------------
# Source my functions
#source(paste0(address, "/", "dgps.R", sep=""))
source(paste0(address, "/dgps_KS_rc.R", sep=""))
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Set seed
set.seed(seed1)
#-----------------------------------------------------------------------------
################################################################
# Run the simulations for all DGPs and all sample sizes (no bootstrap)
bboot <- F

for (nl in 1:1){
  for (nn in 2:3){
    for (dgp in 1:4){
      
      #set sample size
      if(nn==1) n <- 200
      if(nn==2) n <- 500
      if(nn==3) n <- 1000
      if(nn==4) n <- 5000
      #if(nn==5) n <- 10000
      
      # Select lambda
      if(nl==1) lambda <- 0.5
      if(nl==2) lambda <- 0.25
      if(nl==3) lambda <- 0.75
      
      #do the Monte Carlo
      source(paste0(address,"/","sim_rc.R",sep=""))
    }
  }
}

