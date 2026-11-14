########################################
# Simulation for mass imputation
# Estimates: Gold standard estimate (from P sample if Y were observed)
# Author: Rebecca Andridge
# Last modified: 01/17/2026
########################################
library(tidyverse)
require(doParallel)
require(foreach)
library(survey)

# Set up parallel processing -----
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

# Loop through populations -----
start.time <- Sys.time()
for (POPNUM in 1:3)
{
  print(paste0("Pop:", POPNUM))
  
  # Load population data
  filename <- paste0("./simdata/pop",POPNUM,".RData")
  load(filename)
  
  # Stratum population sizes
  N1 <- sum(pop$STRATUM == 1)
  N2 <- sum(pop$STRATUM == 2)
  
  # Load matrices with lists of sampled units
  filename2 <- paste0("./simdata/Psamp_pop",POPNUM,".RData")
  load(filename2)
    
  NREPS <- ncol(Smat.P.SRS)
    
  # Run in parallel -----
  res <- foreach(i=1:NREPS, .combine=rbind, .inorder=FALSE, .packages = c("survey")) %dopar%
    {
      # SRS -----
      # Create sample using sampled indices
      samp.srs <- pop[Smat.P.SRS[,i],]
      samp.srs$WT <- wts$wt_srs
      samp.srs$N <- N1 + N2
      
      # (Gold Standard) Estimates from probability sample
      des <- survey::svydesign(~1, weights = ~WT, fpc = ~N, data = samp.srs)
      est <- survey::svymean(~Y, des)
      est_SRS <- as.numeric(coef(est))
      se_SRS <- as.numeric(sqrt(vcov(est)))
      ci <- confint(est)
      lb_SRS <- ci[1]
      ub_SRS <- ci[2]
      B_SRS <- as.numeric(coef(svyglm(Y ~ Z1 + Z2, des)))
      
      # Stratified sample -----
      # Create sample using sampled indices
      samp.str <- pop[Smat.P.STRAT[,i],]
      samp.str$WT <- wts$wt_strat
      samp.str$N <- ifelse(samp.str$STRATUM == 1, N1, N2)
      
      # (Gold Standard) Estimates from probability sample
      des <- survey::svydesign(~1, weights = ~WT, strata = ~STRATUM, fpc = ~N, data = samp.str)
      est <- survey::svymean(~Y, des)
      est_STR <- as.numeric(coef(est))
      se_STR <- as.numeric(sqrt(vcov(est)))
      ci <- confint(est)
      lb_STR <- ci[1]
      ub_STR <- ci[2]
      B_STR <- as.numeric(coef(svyglm(Y ~ Z1 + Z2, des)))
      
      # Return results -----
      estimates <- data.frame(rbind(c(est_SRS, se_SRS, lb_SRS, ub_SRS, B_SRS),
                                    c(est_STR, se_STR, lb_STR, ub_STR, B_STR)))
      names(estimates) <- c("est_P", "se_P", "lb_P", "ub_P", "B0_P", "B1_P", "B2_P")
      estimates$probtype <- c("SRS", "STR")
      estimates$rep <- i
      
      estimates
    }
  
  # Add pop numbers -----
  res$pop <- POPNUM
  
  # Save out results -----
  filename3 <- paste0("./simresults/gold_pop",POPNUM,".RData")
  save(res, file = filename3)
    
}
  
stop.time <- Sys.time()
# Deregister parallel processing
stopCluster(cl)
# check runtime
stop.time - start.time

# MAC: Time difference of 12.54511 secs