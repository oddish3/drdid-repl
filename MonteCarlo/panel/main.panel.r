#---------------------------------------------------
#      Standard and Stabalized-weighted ATT of DID
#      Panel Data case
#      All DID estimators
#
# Simulations based on R 3.6.2
# Updated: 01/20/2020
#---------------------------------------------------
#-----------------------------------------------------------------------------
# Startup - clear memory, load packages, and set parameters
# Clear memory
rm(list=ls())
#-----------------------------------------------------------------------------
# Basic parameters for the simulation - Doesn't change over setups
#nrep    <- 1000        # Monte Carlo replications
ncores  <- 12           # Number of cores to use in parallel
seed1   <- 1234         # Set initial seed (guaranteed reproducibility)
Xsi.ps <- .75           # pscore index (strength of common support)
nrep <- 1000           # Monte Carlo replications
#-----------------------------------------------------------------------------
# Set the Working Directory
address <- "~/Documents/R_folder/drdid-repl/MonteCarlo/panel"

setwd(address)
#-----------------------------------------------------------------------------
# load the necessary libraries
# Install DRDID package
#devtools::install_github("pedrohcgs/DRDID")
library(foreach)
library(doSNOW)
library(doRNG)
library(trust)
library(DRDID)
#-----------------------------------------------------------------------------
# Source my functions
#source(paste0(address, "/", "dgps.R", sep=""))
source(paste0(address, "/", "dgps.KS.R", sep=""))
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Set seed
set.seed(seed1)
#-----------------------------------------------------------------------------
################################################################
# Run the simulations for all DGPs and all sample sizes (no bootstrap)
bboot <- F

for (nn in 2:2){ # nn = 2
  for (dgp in 1:4){ # dgp = 1
    
    #set sample size
    if(nn==1) n <- 200
    if(nn==2) n <- 500
    if(nn==3) n <- 1000
    #if(nn==4) n <- 5000
    #if(nn==5) n <- 10000
    #do the Monte Carlo
    source(paste0(address,"/","sim.panel.R",sep=""))
  }
}

