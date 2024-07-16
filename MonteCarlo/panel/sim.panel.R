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
iseed <- floor(seed1+n*2598+dgp*711)
base::set.seed(iseed)

#-----------------------------------------------------------------------------
#Start the MONTE CARLO loop
a <- foreach (nn = 1:nrep) %dorng% {
  #---------------------------------------------------------------------------
  #                          Data generating
  #---------------------------------------------------------------------------
  #Generate the data
  data.did <- dgps(dgp,n)
  y0 <- data.did$y0
  y1 <- data.did$y1
  d <-  data.did$d
  x <-  data.did$x
  att.true <- data.did$att
  att.unf <- data.did$att.unf
  eff <- data.did$eff
  #---------------------------------------------------------------------------
  ############################################################################
  ##############          IMPROVED DR ESTIMATOR               ###############
  ############################################################################
  #---------------------------------------------------------------------------
  # Improved DR estimator
  dr.imp <- DRDID::drdid_imp_panel(y1, y0, d, x, boot = bboot)
  # Whether the CI covers the true ATT (coverage probability)
  cp_imp <- as.numeric((dr.imp$lci <= att.true) * (dr.imp$uci >= att.true))
  #Length of confidence interval
  len_imp <- dr.imp$uci - dr.imp$lci
  #---------------------------------------------------------------------------
  ############################################################################
  ##############               DR ESTIMATOR                  #################
  ##############            TRADITIONAL APPROACH             #################
  ############################################################################
  #---------------------------------------------------------------------------
  # Traditional DR
  #dr.tr <- dr.did.panel.tr.linear(y1, y0, d, x, x, boot = bboot, stabilized = T)
  dr.tr <- DRDID::drdid_panel(y1, y0, d, x, boot = bboot)
  # Whether the CI covers the true ATT (coverage probability)
  cp_tr <- as.numeric((dr.tr$lci <= att.true) * (dr.tr$uci >= att.true))
  #Length of confidence interval
  len_tr <- dr.tr$uci - dr.tr$lci
  #---------------------------------------------------------------------------
  #---------------------------------------------------------------------------
  ############################################################################
  ##############          Abadie's IPW ESTIMATOR             #################
  ############################################################################
  #---------------------------------------------------------------------------
  # IPW estimator: Abadie (2005) NOT standardized
  #abadie.ipw <- did.ipw.panel(y1, y0, d, x, boot = bboot)
  abadie.ipw <- DRDID::ipw_did_panel(y1, y0, d, x, boot = bboot)
  # Whether the CI covers the true ATT (coverage probability)
  cp_ipw <- as.numeric((abadie.ipw$lci <= att.true) * (abadie.ipw$uci >= att.true))
  #Length of confidence interval
  len_ipw <- abadie.ipw$uci - abadie.ipw$lci
  #---------------------------------------------------------------------------
  ############################################################################
  ##############    Abadie's IPW ESTIMATOR, Standardized     #################
  ############################################################################
  #---------------------------------------------------------------------------
  # IPW estimator: Abadie (2005) standardized
  #abadie.ipw <- did.ipw.panel.std(y1, y0, d, x, boot = bboot)
  abadie.std.ipw <- DRDID::std_ipw_did_panel(y1, y0, d, x, boot = bboot)
  # Whether the CI covers the true ATT (coverage probability)
  cp_std_ipw <- as.numeric((abadie.std.ipw$lci <= att.true) * (abadie.std.ipw$uci >= att.true))
  #Length of confidence interval
  len_std_ipw <- abadie.std.ipw$uci - abadie.std.ipw$lci
  #---------------------------------------------------------------------------
  ############################################################################
  ##############           Regression DID ESTIMATOR           #################
  ############################################################################
  #---------------------------------------------------------------------------
  # Regression-based DID estimator
  #reg <- did.reg.panel.linear(y1, y0, d, x, boot = bboot)
  reg <- DRDID::reg_did_panel(y1, y0, d, x, boot = bboot)
  # Whether the CI covers the true ATT (coverage probability)
  cp_reg <- as.numeric((reg$lci <= att.true) * (reg$uci >= att.true))
  #Length of confidence interval
  len_reg <- reg$uci - reg$lci
  #---------------------------------------------------------------------------
  ############################################################################
  ##############                TWFE ESTIMATOR              ##################
  ############################################################################
  #---------------------------------------------------------------------------
  # Two Way Fixed Effets
  #twfe <- did.twfe.panel(y1, y0, d, x, boot = bboot)
  twfe <-  DRDID::twfe_did_panel(y1, y0, d, x, boot = bboot)
  # Whether the CI covers the true ATT (coverage probability)
  cp_twfe <- as.numeric((twfe$lci <= att.true) * (twfe$uci >= att.true))
  #Length of confidence interval
  len_twfe <- twfe$uci - twfe$lci
  #---------------------------------------------------------------------------
  
  ##############################################################################
  ############                PUT ALL IN A TABLE        ########################
  ##############################################################################
  out <- cbind(
    # true ATT
    att.true,
    # unfeasible ATT
    att.unf,
    
    # Improved DR estimator
    dr.imp$ATT,  #3
    (dr.imp$se * sqrt(n))^2,   #4
    cp_imp,      #5
    len_imp,    #6
    
    #Traditional DR estimator
    dr.tr$ATT,  #7
    (dr.tr$se * sqrt(n))^2,   #8
    cp_tr,      #9
    len_tr,    #10
    
    #Abadie's IPW estimator (NOT normalized)
    abadie.ipw$ATT,  #11
    (abadie.ipw$se * sqrt(n))^2,   #12
    cp_ipw,      #13
    len_ipw,    #14
    
    #Regression estimator
    reg$ATT,  #15
    (reg$se * sqrt(n))^2,   #16
    cp_reg,      #17
    len_reg,    #18
    
    #TWFE estimator
    twfe$ATT,  #19
    (twfe$se * sqrt(n))^2,   #20
    cp_twfe,      #21
    len_twfe,    #22
    
    #Abadie's IPW estimator (normalized)
    abadie.std.ipw$ATT,  #23
    (abadie.std.ipw$se * sqrt(n))^2,   #24
    cp_std_ipw,      #25
    len_std_ipw,    #26
    
    # FLAG for Bias-reduced pscore
    dr.imp$ps.flag, #27
    
    # Semiparametric Efficiency bound
    eff
  )
  
  
  #Return the output
  return(out)
}

#-----------------------------------------------------------------------------
#Stop the cluster
stopCluster(cl)
#Put the Monte Carlo Results in an Array
mc <- array(unlist(a), dim = c(nrow(a[[1]]), ncol(a[[1]]), length(a)))
mc <- t(matrix(mc, 28, nrep))

#-----------------------------------------------------------------------------
# Mean in the Monte Carlo
mean.mc <- base::colMeans(mc, na.rm = T)
# Median in MC
median.mc <- base::apply(mc, 2, FUN = median, na.rm=T)
# Bias
bias.mc <- base::colMeans(mc - mc[,1], na.rm = T)

median.bias.mc <- base::apply(mc - mc[,1], 2, FUN = median, na.rm=T)


# Standard deviation
sd.mc <- (base::colMeans(mc^2, na.rm = T) - base::colMeans(mc, na.rm = T)^2)^0.5

# RMSE
rmse.mc <- base::colMeans((mc - mc[1])^2, na.rm = T)^0.5

#MAE
mae.mc <- base::colMeans(abs(mc - mc[1]), na.rm = T)

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Create output tables
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#Summary Table with all results
# Create vector of summary statistics for all estimators
unf.summary <- c(dgp,
                 n,
                 mean.mc[2], #average of the estimator
                 mean.mc[2] - mean.mc[1], # Bias
                 median.bias.mc[2],#Median Bias
                 rmse.mc[2], #RMSE
                 mae.mc[2], #MAE
                 sd.mc[2], #Monte Carlo Std deviation
                 
                 NA, #Average Asy. variance
                 NA, # Empirical Coverage
                 NA, #Length of 95% Conf. Int.
                 NA # Semiparametric efficiency Bound
)

dr.imp.summary <- c(dgp,
                    n,
                    mean.mc[3], #average of the estimator
                    mean.mc[3] - mean.mc[1], # Bias
                    median.bias.mc[3],#Median Bias
                    rmse.mc[3], #RMSE
                    mae.mc[3], #MAE
                    sd.mc[3], #Monte Carlo Std deviation
                    
                    #median.mc[4], #Average Asy. variance
                    mean.mc[4],
                    mean.mc[5], # Empirical Coverage
                    mean.mc[6], #Length of 95% Conf. Int.
                    mean.mc[28] # Semiparametric efficiency Bound
)

dr.tr.summary <- c(dgp,
                   n,
                   mean.mc[7], #average of the estimator
                   mean.mc[7] - mean.mc[1], # Bias
                   median.bias.mc[7],#Median Bias
                   rmse.mc[7], #RMSE
                   mae.mc[7], #MAE
                   sd.mc[7], #Monte Carlo Std deviation
                   
                   #median.mc[8], #Average Asy. variance
                   mean.mc[8],
                   mean.mc[9], # Empirical Coverage
                   mean.mc[10], #Length of 95% Conf. Int.
                   mean.mc[28] # Semiparametric efficiency Bound
)


ipw.summary <- c(dgp,
                 n,
                 mean.mc[11], #average of the estimator
                 mean.mc[11] - mean.mc[1], # Bias
                 median.bias.mc[11],#Median Bias
                 rmse.mc[11], #RMSE
                 mae.mc[11], #MAE
                 sd.mc[11], #Monte Carlo Std deviation
                 
                 #median.mc[12], #Average Asy. variance
                 mean.mc[12],
                 mean.mc[13], # Empirical Coverage
                 mean.mc[14], #Length of 95% Conf. Int.
                 mean.mc[28] # Semiparametric efficiency Bound
)



reg.summary <- c(dgp,
                 n,
                 mean.mc[15], #average of the estimator
                 mean.mc[15] - mean.mc[1], # Bias\
                 median.bias.mc[15],#Median Bias
                 rmse.mc[15], #RMSE
                 mae.mc[15], #MAE
                 sd.mc[15], #Monte Carlo Std deviation
                 
                 #median.mc[16], #Average Asy. variance
                 mean.mc[16],
                 mean.mc[17], # Empirical Coverage
                 mean.mc[18], #Length of 95% Conf. Int.
                 mean.mc[28] # Semiparametric efficiency Bound
)

twfe.summary <- c(dgp,
                  n,
                  mean.mc[19], #average of the estimator
                  mean.mc[19] - mean.mc[1], # Bias
                  median.bias.mc[19],#Median Bias
                  rmse.mc[19], #RMSE
                  mae.mc[19], #MAE
                  sd.mc[19], #Monte Carlo Std deviation
                  
                  #median.mc[20], #Average Asy. variance
                  mean.mc[20],
                  mean.mc[21], # Empirical Coverage
                  mean.mc[22], #Length of 95% Conf. Int.
                  mean.mc[28] # Semiparametric efficiency Bound
)


std.ipw.summary <- c(dgp,
                     n,
                     mean.mc[23], #average of the estimator
                     mean.mc[23] - mean.mc[1], # Bias
                     median.bias.mc[23],#Median Bias
                     rmse.mc[23], #RMSE
                     mae.mc[23], #MAE
                     sd.mc[23], #Monte Carlo Std deviation
                     
                     #median.mc[24], #Average Asy. variance
                     mean.mc[24],
                     mean.mc[25], # Empirical Coverage
                     mean.mc[26], #Length of 95% Conf. Int.
                     mean.mc[28] # Semiparametric efficiency Bound
)

#-----------------------------------------------------------------------------
# Create the table
mc.summary <- rbind(unf.summary,
                    dr.tr.summary,
                    dr.imp.summary,
                    ipw.summary, std.ipw.summary,
                    reg.summary, twfe.summary
                    
)

rownames(mc.summary) <- c("unf.att",
                          "DR.tr",
                          "DR.imp",
                          
                          "IPW",
                          "IPW-STD",
                          "REG",
                          'TWFE'
                          
)

colnames(mc.summary) <- c("DGP", "n", "Estimator", "Av.Bias", "Med.Bias","RMSE", "MAE", "MCSD", 
                          "Asy.Var", "Coverage", "Lenth-CI", "Sem. Eff. Bound")

out1 <- paste0("mc.summary",".dgp-",dgp,".n-",n,".csv")
out1 <- paste0(address, "/", "Results.panel","/", out1, sep="")
write.csv(mc.summary, file = out1)
#-----------------------------------------------------------------------------
# Save simulation results as RData
out1Data <- paste0(address,"/MC-results/", "mc.dgp-",dgp,".n-",n,".RData")
save.image(out1Data)
