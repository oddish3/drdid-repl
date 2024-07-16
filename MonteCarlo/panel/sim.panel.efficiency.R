#---------------------------------------------------
#      Perform the Monte Carlo for ATT of DID
#      Panel Data
#      Further-improved and Traditional DR-DID estimators
#      IPW, Regression and Two-Way Fixed Effects, too
#---------------------------------------------------
#-----------------------------------------------------------------------------
#Make cluster
cl <- makeCluster(ncores)
registerDoSNOW(cl)

#Set seed
base::set.seed(seed1)
iseed <- floor(seed1+dgp*711)
base::set.seed(iseed)

#-----------------------------------------------------------------------------
#Start the MONTE CARLO loop
a <- foreach (nn = 1:nrep) %dorng% {
  #---------------------------------------------------------------------------
  #                          Data generating
  #---------------------------------------------------------------------------
  #Generate the data
  eff.bound <- dgps(dgp,n)
  
  
  ##############################################################################
  ############                PUT ALL IN A TABLE        ########################
  ##############################################################################
  out <- eff.bound
  
  
  
  #Return the output
  return(out)
}

#-----------------------------------------------------------------------------
#Stop the cluster
stopCluster(cl)
#Put the Monte Carlo Results in an Array
mc <- array(unlist(a), dim = c(nrow(a[[1]]), ncol(a[[1]]), length(a)))
mc <- t(matrix(mc, 1, nrep))

#-----------------------------------------------------------------------------
# Mean in the Monte Carlo
mean.mc <- base::colMeans(mc, na.rm = T)

out1 <- paste0("eff.bound.dgp-",dgp,".csv")
out1 <- paste0(address, "/", "Results.panel","/", out1, sep="")
write.csv(mean.mc, file = out1)