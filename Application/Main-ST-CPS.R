#------------------------------------------------------------------------------------------------------
#      Application: DR ATT estimators based on DID
#      Smith and Todd, 2005 (JoE)
#      Selection effect
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# MODIFIED    :   August 24, 2019

# AUTHOR      :   Pedro H. C. Sant'Anna
# AFFILIATION :   Vanderbilt University
# EMAIL       :   pedro.h.santanna@vanderbilt.edu

# AUTHOR      :   Jun Zhao
# AFFILIATION :   Vanderbilt University
# EMAIL       :   jun.zhao@vanderbilt.edu
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# Startup - clear memory, load packages, set working directory, and import the dataset
# Clear memory
rm(list=ls())
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# Set the Working Directory
address <- "set/here/your/working/directory"
setwd(address)
#-----------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# Functions for the applications
# function to do Implement all DID estimators at once
source(paste0(address, "/", "all.did.estimators.R", sep=""))
# Function to implement all DID estimators for each subsample
source(paste0(address, "/", "all.did.subsample.R", sep=""))
# Function to generate Tables
source(paste0(address, "/", "out.table.lalonde.R", sep=""))
#------------------------------------------------------------------------------------------------------
# load the necessary libraries
devtools::install_github("pedrohcgs/DRDID")
library(DRDID)
library(sem)
library(foreign)
library(readstata13)
library(matlib)
library(trust)
#------------------------------------------------------------------------------------------------------
# import the data
ST_2005 <- read.dta13(paste(address,"/data/nsw.dta",sep = ""),
                      missing.type=T)
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# Create Subsamples
source(paste0(address, "/", "subsamples.R", sep=""))
#------------------------------------------------------------------------------------------------------
set.seed(123)
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#####################             #CPS Comparison group       ##############################
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# Evaluation bias
#------------------------------------------------------------------------------------------------------
# Lalonde sample
eval.cps.lalonde <- all.did.subsample(eval.lalonde.cps)
# DW sample
eval.cps.dw <- all.did.subsample(eval.dw.cps)
# Early RA sample
eval.cps.early <- all.did.subsample(eval.early.cps)
#------------------------------------------------------------------------------------------------------
# generate tables
eval.cps.lalonde.tab <- out.table.lalonde(eval.cps.lalonde, benchmark = 886)
eval.cps.dw.tab <- out.table.lalonde(eval.cps.dw, benchmark = 1794)
eval.cps.early.tab <- out.table.lalonde(eval.cps.early, benchmark = 2748)
# Put the tables together
cps.eval.bias <- cbind(eval.cps.lalonde.tab,
                       eval.cps.dw.tab,
                       eval.cps.early.tab)

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
save.image(paste0(address, "/", "results","/Smith-Todd-application.Rdata", sep=""))
#------------------------------------------------------------------------------------------------------
# Save CSV
out.eval <- paste0("ST-eval-CPS.csv")
out.eval <- paste0(address, "/", "results","/", out.eval, sep="")
write.csv(cps.eval.bias, file = out.eval)
#------------------------------------------------------------------------------------------------------
