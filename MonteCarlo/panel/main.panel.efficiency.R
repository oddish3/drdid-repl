#---------------------------------------------------
#      Standard and Stabalized-weighted ATT of DID
#      Panel Data case
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
nrep <- 1000           # Monte Carlo replications
#-----------------------------------------------------------------------------
# Set the Working Directory
address <- "set/here/your/working/directory/for/panel"
setwd(address)
#-----------------------------------------------------------------------------
# load the necessary libraries
# Install DRDID package
devtools::install_github("pedrohcgs/DRDID")
library(foreach)
library(doSNOW)
library(doRNG)
#library(trust)
#library(DRDID)
#-----------------------------------------------------------------------------
# Source my functions
#source(paste0(address, "/", "dgps.R", sep=""))
source(paste0(address, "/", "dgps.KS.efficiency.R", sep=""))
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Set seed
set.seed(seed1)
#-----------------------------------------------------------------------------
################################################################
# Run the simulations for all DGPs and all sample sizes (no bootstrap)
bboot <- F

for (dgp in 1:4){

#dgp=1
  #set sample size
  n <- 10000000
  #do the Monte Carlo
  source(paste0(address,"/","sim.panel.efficiency.R",sep=""))
}


