########################################
# Simulation for mass imputation
# Estimates: Gold standard and Naive NP estimate
# Author: Rebecca Andridge
# Last modified: 06/27/2025
########################################
library(tidyverse)
require(doParallel)
require(foreach)
# require(doRNG)

# Set up parallel processing -----
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

# Loop through populations -----
start.time <- Sys.time()
for (POPNUM in 1:3)
{
  print(paste0("Pop:", POPNUM))
  # Load population data -----
  filename <- paste0("./simdata/pop",POPNUM,".RData")
  load(filename)
  
  # Loop through selection mechanisms
  for (MECHNUM in 1:9)
  {
    print(paste0("Mech:", MECHNUM))
    # Load samples info -----
    filename2 <- paste0("./simdata/samp_pop",POPNUM,"mech",MECHNUM,".RData")
    load(filename2)
    NREPS <- ncol(Smat.NP)
    
    # Run in parallel -----
    res <- foreach(i=1:NREPS, .combine=rbind, .inorder=FALSE) %dopar%
      {
        # Create P and NP samples using sampled indices -----
        # prob sample
        samp.P <- pop[Smat.P[,i],]
        n_P <- nrow(samp.P)
        # nonprob sample
        samp.NP <- pop[Smat.NP[,i],]
        n_NP <- nrow(samp.NP)
        
        # (Gold Standard) Estimated mean from probability sample -----
        est_P <- mean(samp.P$Y)
        se_P <- sd(samp.P$Y)/sqrt(n_P)
        ci <- confint(lm(Y~1, samp.P))
        lb_P <- ci[1]
        ub_P <- ci[2]
        
        # (Naive Estimate) Estimated mean from nonprobability sample -----
        est_NP <- mean(samp.NP$Y)
        se_NP <- sd(samp.NP$Y)/sqrt(n_NP)
        ci <- confint(lm(Y~1, samp.NP))
        lb_NP <- ci[1]
        ub_NP <- ci[2]
        
        # Regression parameter estimates from Y|Z -----
        # probability sample
        B_P <- as.numeric(coef(lm(Y ~ Z1 + Z2, data = samp.P)))
        # nonprobability sample
        B_NP <- as.numeric(coef(lm(Y ~ Z1 + Z2, data = samp.NP)))
        
        # Return results -----
        data.frame(cbind(est_P, se_P, lb_P, ub_P,
                         est_NP, se_NP, lb_NP, ub_NP,
                         B0_P = B_P[1], B1_P = B_P[2], B2_P = B_P[3],
                         B0_NP = B_NP[1], B1_NP = B_NP[2], B2_NP = B_NP[3]))
      }
    
    # Add pop and mech numbers -----
    res$pop <- POPNUM
    res$mech <- MECHNUM
    
    # Save out results -----
    filename3 <- paste0("./simresults/gold_naive_pop",POPNUM,"mech",MECHNUM,".RData")
    save(res, file = filename3)
  }
}
  
stop.time <- Sys.time()
# Deregister parallel processing
stopCluster(cl)
# check runtime
stop.time - start.time
#  Time difference of 23.59679 secs

