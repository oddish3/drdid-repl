#---------------------------------------------------
#      Perform the Monte Carlo for ATT of DID
#      Repeated Cross Section
#      Improved and Traditional DR-DID estimators
#      IPW, Regression and Two-Way Fixed Effects, too
#---------------------------------------------------
#-----------------------------------------------------------------------------
#Make cluster
cl <- makeCluster(ncores)
registerDoSNOW(cl)

#Set seed
base::set.seed(seed1)
iseed <- floor(seed1+n*2598+dgp*711+lambda*7152)
base::set.seed(iseed)

#-----------------------------------------------------------------------------
#Start the MONTE CARLO loop
a <- foreach (nn = 1:nrep) %dorng% {
  #---------------------------------------------------------------------------
  #                          Data generating
  #---------------------------------------------------------------------------
  #Generate the data
  data.did <- dgps_rc(dgp, n, lambda)
  y <- data.did$y
  d <-  data.did$d
  x <-  data.did$x
  post <- data.did$post
  att.true <- data.did$att
  att.unf <- data.did$att.unf
  eff <- data.did$eff
  #---------------------------------------------------------------------------
  ############################################################################
  ##############            IMPROVED DR ESTIMATOR             ###############
  ############################################################################
  #---------------------------------------------------------------------------
  # Improved DR estimator - not locally efficient
  dr.imp1 <- DRDID::drdid_imp_rc1(y, post, d, x, boot = bboot)
  # Whether the CI covers the true ATT (coverage probability)
  cp_imp1 <- as.numeric((dr.imp1$lci <= att.true) * (dr.imp1$uci >= att.true))
  #Length of confidence interval
  len_imp1 <- dr.imp1$uci - dr.imp1$lci
  #---------------------------------------------------------------------------
  # Improved DR estimator - locally efficient
  dr.imp2 <- DRDID::drdid_imp_rc(y, post, d, x, boot = bboot)
  # Whether the CI covers the true ATT (coverage probability)
  cp_imp2 <- as.numeric((dr.imp2$lci <= att.true) * (dr.imp2$uci >= att.true))
  #Length of confidence interval
  len_imp2 <- dr.imp2$uci - dr.imp2$lci
  #---------------------------------------------------------------------------
  ############################################################################
  ##############               DR ESTIMATOR                  #################
  ##############            TRADITIONAL APPROACH             #################
  ############################################################################
  #---------------------------------------------------------------------------
  # Traditional DR - NOT LOCALLY EFFICIENT
  dr.tr1 <- DRDID::drdid_rc1(y, post, d, x, boot = bboot)
  # Whether the CI covers the true ATT (coverage probability)
  cp_tr1 <- as.numeric((dr.tr1$lci <= att.true) * (dr.tr1$uci >= att.true))
  #Length of confidence interval
  len_tr1 <- dr.tr1$uci - dr.tr1$lci
  #---------------------------------------------------------------------------
  # Traditional DR - LOCALLY EFFICIENT
  dr.tr2 <- DRDID::drdid_rc(y, post, d, x, boot = bboot)
  # Whether the CI covers the true ATT (coverage probability)
  cp_tr2 <- as.numeric((dr.tr2$lci <= att.true) * (dr.tr2$uci >= att.true))
  #Length of confidence interval
  len_tr2 <- dr.tr2$uci - dr.tr2$lci
  #---------------------------------------------------------------------------
  #---------------------------------------------------------------------------
  ############################################################################
  ##############          Abadie's IPW ESTIMATOR             #################
  ############################################################################
  #---------------------------------------------------------------------------
  # IPW estimator: Abadie (2005)
  abadie.ipw <- DRDID::ipw_did_rc(y, post, d, x, boot = bboot)
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
  abadie.std.ipw <- DRDID::std_ipw_did_rc(y, post, d, x, boot = bboot)
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
  reg <- DRDID::reg_did_rc(y, post, d, x, boot = bboot)
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
  twfe <-  DRDID::twfe_did_rc(y, post, d, x, boot = bboot)
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
    
    # Improved DR estimator (but not locally efficient)
    dr.imp1$ATT,  #3
    (dr.imp1$se * sqrt(n))^2,   #4
    cp_imp1,      #5
    len_imp1,    #6
    
    #Traditional DR estimator - Not locally efficient
    dr.tr1$ATT,  #7
    (dr.tr1$se * sqrt(n))^2,   #8
    cp_tr1,      #9
    len_tr1,    #10
    
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
    dr.imp1$ps.flag, #27
    
    # Improved and locally efficient DR estimator
    dr.imp2$ATT,  #28
    (dr.imp2$se * sqrt(n))^2,   #29
    cp_imp2,      #30
    len_imp2,    #31
    
    #Traditional DR estimator - locally efficient
    dr.tr2$ATT,  #32
    (dr.tr2$se * sqrt(n))^2,   #33
    cp_tr2,      #34
    len_tr2,   #35
    
    # Semiparametric Efficiency bound
    eff #36
    
   
  )
  
  
  #Return the output
  return(out)
}

#-----------------------------------------------------------------------------
#Stop the cluster
stopCluster(cl)
#Put the Monte Carlo Results in an Array
mc <- array(unlist(a), dim = c(nrow(a[[1]]), ncol(a[[1]]), length(a)))
mc <- t(matrix(mc, 36, nrep))

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
                 sd.mc[2], #Monte Carlo Std error
                 
                 NA, #Average Asy. variance
                 NA, # Empirical Coverage
                 NA, #Length of 95% Conf. Int.
                 NA # Semiparametric efficiency Bound

)

dr.imp1.summary <- c(dgp,
                   n,
                   mean.mc[3], #average of the estimator
                   mean.mc[3] - mean.mc[1], # Bias
                   median.bias.mc[3],#Median Bias
                   rmse.mc[3], #RMSE
                   mae.mc[3], #MAE
                   sd.mc[3], #Monte Carlo Std error
                   
                   #median.mc[4], #Average Asy. variance
                   mean.mc[4],
                   mean.mc[5], # Empirical Coverage
                   mean.mc[6], #Length of 95% Conf. Int.
                   mean.mc[36] # Semiparametric efficiency Bound
)

dr.tr1.summary <- c(dgp,
                   n,
                   mean.mc[7], #average of the estimator
                   mean.mc[7] - mean.mc[1], # Bias
                   median.bias.mc[7],#Median Bias
                   rmse.mc[7], #RMSE
                   mae.mc[7], #MAE
                   sd.mc[7], #Monte Carlo Std error
                   
                   #median.mc[8], #Average Asy. variance
                   mean.mc[8],
                   mean.mc[9], # Empirical Coverage
                   mean.mc[10], #Length of 95% Conf. Int.
                   mean.mc[36] # Semiparametric efficiency Bound
)


ipw.summary <- c(dgp,
                 n,
                 mean.mc[11], #average of the estimator
                 mean.mc[11] - mean.mc[1], # Bias
                 median.bias.mc[11],#Median Bias
                 rmse.mc[11], #RMSE
                 mae.mc[11], #MAE
                 sd.mc[11], #Monte Carlo Std error
                 
                 #median.mc[12], #Average Asy. variance
                 mean.mc[12],
                 mean.mc[13], # Empirical Coverage
                 mean.mc[14], #Length of 95% Conf. Int.
                 mean.mc[36] # Semiparametric efficiency Bound
)



reg.summary <- c(dgp,
                 n,
                 mean.mc[15], #average of the estimator
                 mean.mc[15] - mean.mc[1], # Bias
                 median.bias.mc[15],#Median Bias
                 rmse.mc[15], #RMSE
                 mae.mc[15], #MAE
                 sd.mc[15], #Monte Carlo Std error
                 
                 #median.mc[16], #Average Asy. variance
                 mean.mc[16],
                 mean.mc[17], # Empirical Coverage
                 mean.mc[18], #Length of 95% Conf. Int.
                 mean.mc[36] # Semiparametric efficiency Bound
)

twfe.summary <- c(dgp,
                  n,
                  mean.mc[19], #average of the estimator
                  mean.mc[19] - mean.mc[1], # Bias
                  median.bias.mc[19],#Median Bias
                  rmse.mc[19], #RMSE
                  mae.mc[19], #MAE
                  sd.mc[19], #Monte Carlo Std error
                  
                  #median.mc[20], #Average Asy. variance
                  mean.mc[20],
                  mean.mc[21], # Empirical Coverage
                  mean.mc[22], #Length of 95% Conf. Int.
                  mean.mc[36] # Semiparametric efficiency Bound
)


std.ipw.summary <- c(dgp,
                     n,
                     mean.mc[23], #average of the estimator
                     mean.mc[23] - mean.mc[1], # Bias
                     median.bias.mc[23],#Median Bias
                     rmse.mc[23], #RMSE
                     mae.mc[23], #MAE
                     sd.mc[23], #Monte Carlo Std error
                     
                     #median.mc[24], #Average Asy. variance
                     mean.mc[24],
                     mean.mc[25], # Empirical Coverage
                     mean.mc[26], #Length of 95% Conf. Int.
                     mean.mc[36] # Semiparametric efficiency Bound
)

dr.imp2.summary <- c(dgp,
                 n,
                 mean.mc[28], #average of the estimator
                 mean.mc[28] - mean.mc[1], # Bias
                 median.bias.mc[28],#Median Bias
                 rmse.mc[28], #RMSE
                 mae.mc[28], #MAE
                 sd.mc[28], #Monte Carlo Std error
                 
                 #median.mc[29], #Average Asy. variance
                 mean.mc[29],
                 mean.mc[30], # Empirical Coverage
                 mean.mc[31], #Length of 95% Conf. Int.
                 mean.mc[36] # Semiparametric efficiency Bound
)

dr.tr2.summary <- c(dgp,
                    n,
                    mean.mc[32], #average of the estimator
                    mean.mc[32] - mean.mc[1], # Bias
                    median.bias.mc[32],#Median Bias
                    rmse.mc[32], #RMSE
                    mae.mc[32], #MAE
                    sd.mc[32], #Monte Carlo Std error
                    
                    #median.mc[33], #Average Asy. variance
                    mean.mc[33],
                    mean.mc[34], # Empirical Coverage
                    mean.mc[35], #Length of 95% Conf. Int.
                    mean.mc[36] # Semiparametric efficiency Bound
)

#-----------------------------------------------------------------------------
# Create the table
mc.summary <- rbind(unf.summary,
                    twfe.summary,
                    reg.summary,
                    ipw.summary, std.ipw.summary,
                    dr.tr1.summary,
                    dr.tr2.summary,
                    dr.imp1.summary,
                    dr.imp2.summary
)

rownames(mc.summary) <- c("unf.att",
                          'TWFE',
                          "REG",
                          "IPW",
                          "IPW-STD",
                          "DR.tr1",
                          "DR.tr2",
                          "DR.imp1",
                          "DR.imp2"
)

colnames(mc.summary) <- c("DGP", "n", "Estimator", "Av.Bias", "Med.Bias","RMSE", "MAE", "MCSD", 
                          "Asy.Var", "Coverage", "Lenth-CI", "Sem. Eff. Bound")

out1 <- paste0("mc.summary-rc",".dgp-",dgp,".n-",n,".lambda-",lambda,".csv")
out1 <- paste0(address, "/", "Results.rc","/", out1, sep="")
write.csv(mc.summary, file = out1)
#-----------------------------------------------------------------------------
# Save simulation results as RData
out1Data <- paste0(address,"/MC-results/", "mc-rc.dgp-",dgp,".n-",n,".lambda-",lambda,".RData")
save.image(out1Data)
