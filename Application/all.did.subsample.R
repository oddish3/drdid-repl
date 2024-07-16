all.did.subsample <- function(subsample){
  #------------------------------------------------------------------------------------------------------
  # Attach the subsample
  attach(subsample)
  #------------------------------------------------------------------------------------------------------
  # create all specifications' sets of controls, outcomes pre and post, treatment
  # The different specifications
  x.linear <- cbind(age, educ, black, married, nodegree, hisp, re74)
  x.dw <- cbind(age, educ, black, married, nodegree, hisp, re74, ze74, agesq, agecub, educsq, educre74)
  x.aug.dw <- cbind(age, educ, black, married, nodegree, hisp, re74, ze74, agesq, agecub, educsq, educre74, mare74, maze74)
  # outcomes
  y1 <- c(re78)
  y0 <- c(re75)
  # Group/Treatment dummy
  d <- c(treated2)
  #------------------------------------------------------------------------------------------------------
  # detach the subsample
  detach(subsample)
  #------------------------------------------------------------------------------------------------------
  # Compute all DID's with X linear
  did.linear <- all.did.estimators(y1, y0, d, x.linear)
  # Compute all DID's with X as in DW
  did.dw <- all.did.estimators(y1, y0, d, x.dw)
  # Compute all DID's with X as in Augmented DW
  did.aug.dw <- all.did.estimators(y1, y0, d, x.aug.dw)

  #------------------------------------------------------------------------------------------------------
  # put all outputs in a list
  all.did <- list( linear = did.linear,
                   dw = did.dw,
                   aug.dw = did.aug.dw)
  return(all.did)
} 