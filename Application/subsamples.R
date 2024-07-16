# Create All subsamples

#------------------------------------------------------------------------------------------------------
# create nonlinear covariate terms
attach(ST_2005)
ST_2005$ze74 <- as.numeric(re74==0)
ST_2005$ze75 <- as.numeric(re75==0)

ST_2005$agesq <- age^2
ST_2005$agecub <- age^3/1000

ST_2005$educsq <- educ^2
ST_2005$re74sq <- re74^2

ST_2005$educre74 <- educ*re74
ST_2005$ze74hisp <- as.numeric(re74==0)*hisp

ST_2005$mare74 <- married*re74
ST_2005$maze74 <- married*as.numeric(re74==0)

detach(ST_2005)
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# create all "Selection" subsamples to analyze Selection Bias
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# Create "selection" treatment dummy: 1 if in experimental sample, 0 if in non-experimental
# if treatment=NA, it is the non-experimental comparison group, treated whenever in exp group
ST_2005$treated2 <- ifelse(is.na(ST_2005$treated), 0 , 1)

eval.lalonde.cps <- subset(ST_2005, ST_2005$treated==0 | ST_2005$sample==2)
eval.dw.cps <- subset(ST_2005, (ST_2005$dwincl==1 & ST_2005$treated==0) | ST_2005$sample==2)
eval.early.cps <- subset(ST_2005, (ST_2005$early_ra==1 & ST_2005$treated==0) | ST_2005$sample==2)
#------------------------------------------------------------------------------------------------------
# Delete the variable treated2 from the data
ST_2005 <- ST_2005[,colnames(ST_2005) !="treated2"]
#------------------------------------------------------------------------------------------------------
detach()
