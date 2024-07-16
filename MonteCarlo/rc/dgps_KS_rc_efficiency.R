#DGP's for "Doubly RObust Difference-in-Difference Estimators"
#08/25/2019 - Repeated Cross sectiion

# DGPs based on Kang and Schafer
#Researcher always observes Z
#-----------------------------------------------------------------------------
# Mean and Std deviation of Z's without truncation
mean.z1 <- exp(0.25/2)
sd.z1 <- sqrt((exp(0.25) - 1) * exp(0.25))

mean.z2 <- 10
sd.z2 <- 0.54164

mean.z3 <- 0.21887
sd.z3 <-   0.04453

mean.z4 <- 402
sd.z4 <-  56.63891
#-----------------------------------------------------------------------------

dgps_rc = function(dgp, n, lambda=0.5){
  
  # Gen covariates
  x1 <- stats::rnorm(n, mean = 0, sd = 1)
  x2 <- stats::rnorm(n, mean = 0, sd = 1)
  x3 <- stats::rnorm(n, mean = 0, sd = 1)
  x4 <- stats::rnorm(n, mean = 0, sd = 1)
  
  z1 <- exp(x1/2)
  z2 <- x2/(1 + exp(x1)) + 10
  z3 <- (x1 * x3/25 + 0.6)^3
  z4 <- (x1 + x4 + 20)^2
  
  z1 <- (z1 - mean.z1)/sd.z1
  z2 <- (z2 - mean.z2)/sd.z2
  z3 <- (z3 - mean.z3)/sd.z3
  z4 <- (z4 - mean.z4)/sd.z4
  
  x <- cbind(x1, x2, x3, x4)
  z <- cbind(z1, z2, z3, z4)
  
  post <- as.numeric(stats::runif(n) <= lambda)
  #-----------------------------------------------------------------------------
  # Always observe Z
  #-----------------------------------------------------------------------------
  # Both models are functions of Z
  if(dgp ==1) {
    # Gen groups
    # Propensity score
    pi <- stats::plogis(Xsi.ps * (- z1 + 0.5 * z2 - 0.25 * z3 - 0.1 * z4)) 
    #pi <- pmin(pi, 0.999)
    d  <- as.numeric(runif(n) <= pi)
    
    # Generate aux indexes for the potential outcomes
    index.lin <- 210 + 27.4*z1 + 13.7*(z2 + z3 + z4)
    index.unobs.het <- d * (index.lin)
    index.att <- 0#10 + 27.4*x1 + 13.7*(x2 + x3 + x4)
    
    #This is the key for consistency out outcome regression
    index.trend <- 210 + 27.4*z1 + 13.7*(z2 + z3 + z4)
    
    #v is the unobserved heterogeneity
    v <- stats::rnorm(n, mean = index.unobs.het, sd = 1)
    
    #Gen realized outcome at time 0
    y0 <- index.lin + v + stats::rnorm(n) 
    
    # gen outcomes at time 1
    # First let's generate potential outcomes: y_1_potential
    y10 <- index.lin + v + stats::rnorm(n, mean = 0, sd = 1) +#This is the baseline
      index.trend #this is for the trend based on X  
    
    y11 <- index.lin + v + stats::rnorm(n, mean = 0, sd = 1) +#This is the baseline
      index.trend + #this is for the trend based on X
      index.att # This is the treatment effects
    
    
    # Get unfeasible att
    att.unf <- (mean(d*y11) - mean(d*y10))/mean(d)
    
    # Gen realized outcome at time 1
    y1 <- d * y11 + (1 - d) * y10
    
    # observed outcome
    y <- post * y1 + (1 - post) * y0
    
    
    # DOUBLE CHECK!!! THIS IS WHERE I STOPPED!!
    term1 <- d * (((index.trend + index.att) - (index.trend))^2)
    term2 <- (d*post/(lambda^2)) * ((y - (index.lin+index.lin + index.trend + index.att))^2)
    
    term3 <- (d*(1 - post)/((1-lambda)^2)) * ((y - (index.lin + index.lin))^2)
   
    term4 <- (1-d)*post* ((pi/((1-pi)*lambda))^2) * (y - (index.lin + index.trend))^2
    term5 <- (1-d)*(1 - post)* ((pi/((1-pi)*lambda))^2) * (y - (index.lin ))^2
    
    efficiency.bound = mean(term1 + term2 + term3 + term4 + term5)/(mean(d)^2)
    
    
    
    
    return(efficiency.bound)
      
  }
  
  
  
  # Pscore depends on X but Regressions depend on Z
  if(dgp == 2) {
    # Gen groups
    # Propensity score
    pi <- stats::plogis( Xsi.ps * (- x1 + 0.5 * x2 - 0.25 * x3 - 0.1 * x4)) 
    #pi <- pmin(pi, 0.999)
    d  <- as.numeric(runif(n) <= pi)
    
    # Generate aux indexes for the potential outcomes
    index.lin <- 210 + 27.4*z1 + 13.7*(z2 + z3 + z4)
    index.unobs.het <- d * (index.lin)
    index.att <- 0 # 10 + 27.4*x1 + 13.7*(x2 + x3 + x4)
    
    #This is the key for consistency out outcome regression
    index.trend <- 210 + 27.4*z1 + 13.7*(z2 + z3 + z4)
    
    #v is the unobserved heterogeneity
    v <- stats::rnorm(n, mean = index.unobs.het, sd = 1)
    
    #Gen realized outcome at time 0
    y0 <- index.lin + v + stats::rnorm(n) 
    
    # gen outcomes at time 1
    # First let's generate potential outcomes: y_1_potential
    y10 <- index.lin + v + stats::rnorm(n, mean = 0, sd = 1) +#This is the baseline
      index.trend #this is for the trend based on X  
    
    y11 <- index.lin + v + stats::rnorm(n, mean = 0, sd = 1) +#This is the baseline
      index.trend + #this is for the trend based on X
      index.att # This is the treatment effects
    
    # Get unfeasible att
    att.unf <- (mean(d*y11) - mean(d*y10))/mean(d)
    
    # Gen realized outcome at time 1
    y1 <- d * y11 + (1 - d) * y10
    
    # observed outcome
    y <- post * y1 + (1 - post) * y0
    # DOUBLE CHECK!!! THIS IS WHERE I STOPPED!!
    term1 <- d * (((index.trend + index.att) - (index.trend))^2)
    term2 <- (d*post/(lambda^2)) * ((y - (index.lin+index.lin + index.trend + index.att))^2)
    
    term3 <- (d*(1 - post)/((1-lambda)^2)) * ((y - (index.lin + index.lin))^2)
    
    term4 <- (1-d)*post* ((pi/((1-pi)*lambda))^2) * (y - (index.lin + index.trend))^2
    term5 <- (1-d)*(1 - post)* ((pi/((1-pi)*lambda))^2) * (y - (index.lin ))^2
    
    efficiency.bound = mean(term1 + term2 + term3 + term4 + term5)/(mean(d)^2)
    
    
    
    
    return(efficiency.bound)
  }
  
  # Pscore depends on Z but Regressions depend on X
  if(dgp ==3) {
    
    # Gen groups
    # Propensity score
    pi <- stats::plogis(Xsi.ps * ( - z1 + 0.5 * z2 - 0.25 * z3 - 0.1 * z4)) 
    #pi <- pmin(pi, 0.999)
    d  <- as.numeric(runif(n) <= pi)
    
    # Generate aux indexes for the potential outcomes
    index.lin <- 210 + 27.4*x1 + 13.7*(x2 + x3 + x4)
    index.unobs.het <- d * (index.lin)
    index.att <- 0# 10 + 27.4*x1 + 13.7*(x2 + x3 + x4)
    
    #This is the key for consistency out outcome regression
    index.trend <- 210 + 27.4*x1 + 13.7*(x2 + x3 + x4)
    #v is the unobserved heterogeneity
    v <- stats::rnorm(n, mean = index.unobs.het, sd = 1)
    
    #Gen realized outcome at time 0
    y0 <- index.lin + v + stats::rnorm(n) 
    
    # gen outcomes at time 1
    # First let's generate potential outcomes: y_1_potential
    y10 <- index.lin + v + stats::rnorm(n, mean = 0, sd = 1) +#This is the baseline
      index.trend #this is for the trend based on X  
    
    y11 <- index.lin + v + stats::rnorm(n, mean = 0, sd = 1) +#This is the baseline
      index.trend + #this is for the trend based on X
      index.att # This is the treatment effects
    
    # Get unfeasible att
    att.unf <- (mean(d*y11) - mean(d*y10))/mean(d)
    
    # Gen realized outcome at time 1
    y1 <- d * y11 + (1 - d) * y10
    
    # observed outcome
    y <- post * y1 + (1 - post) * y0
    
    # DOUBLE CHECK!!! THIS IS WHERE I STOPPED!!
    term1 <- d * (((index.trend + index.att) - (index.trend))^2)
    term2 <- (d*post/(lambda^2)) * ((y - (index.lin+index.lin + index.trend + index.att))^2)
    
    term3 <- (d*(1 - post)/((1-lambda)^2)) * ((y - (index.lin + index.lin))^2)
    
    term4 <- (1-d)*post* ((pi/((1-pi)*lambda))^2) * (y - (index.lin + index.trend))^2
    term5 <- (1-d)*(1 - post)* ((pi/((1-pi)*lambda))^2) * (y - (index.lin ))^2
    
    efficiency.bound = mean(term1 + term2 + term3 + term4 + term5)/(mean(d)^2)
    
    
    
    
    return(efficiency.bound)
  }
  
  
  # Both regression and pscore depend on x
  # BOTH MODELS ARE MISSPECIFIED
  if(dgp == 4) {
    
    # Gen groups
    # Propensity score
    pi <- stats::plogis( Xsi.ps * (- x1 + 0.5 * x2 - 0.25 * x3 - 0.1 * x4)) 
    #pi <- pmin(pi, 0.999)
    d  <- as.numeric(runif(n) <= pi)
    
    
    # Generate aux indexes for the potential outcomes
    index.lin <- 210 + 27.4*x1 + 13.7*(x2 + x3 + x4)
    index.unobs.het <- d * (index.lin)
    index.att <- 0# 10 + 27.4*x1 + 13.7*(x2 + x3 + x4)
    
    #This is the key for consistency out outcome regression
    index.trend <- 210 + 27.4*x1 + 13.7*(x2 + x3 + x4)
    
    #v is the unobserved heterogeneity
    v <- stats::rnorm(n, mean = index.unobs.het, sd = 1)
    
    #Gen realized outcome at time 0
    y0 <- index.lin + v + stats::rnorm(n) 
    
    # gen outcomes at time 1
    # First let's generate potential outcomes: y_1_potential
    y10 <- index.lin + v + stats::rnorm(n, mean = 0, sd = 1) +#This is the baseline
      index.trend #this is for the trend based on X  
    
    y11 <- index.lin + v + stats::rnorm(n, mean = 0, sd = 1) +#This is the baseline
      index.trend + #this is for the trend based on X
      index.att # This is the treatment effects
    
    # Get unfeasible att
    att.unf <- (mean(d*y11) - mean(d*y10))/mean(d)
    
    # Gen realized outcome at time 1
    y1 <- d * y11 + (1 - d) * y10
    
    # observed outcome
    y <- post * y1 + (1 - post) * y0
    
    # DOUBLE CHECK!!! THIS IS WHERE I STOPPED!!
    term1 <- d * (((index.trend + index.att) - (index.trend))^2)
    term2 <- (d*post/(lambda^2)) * ((y - (index.lin+index.lin + index.trend + index.att))^2)
    
    term3 <- (d*(1 - post)/((1-lambda)^2)) * ((y - (index.lin + index.lin))^2)
    
    term4 <- (1-d)*post* ((pi/((1-pi)*lambda))^2) * (y - (index.lin + index.trend))^2
    term5 <- (1-d)*(1 - post)* ((pi/((1-pi)*lambda))^2) * (y - (index.lin ))^2
    
    efficiency.bound = mean(term1 + term2 + term3 + term4 + term5)/(mean(d)^2)
    
    
    
    
    return(efficiency.bound)
  }
  
  
}
