########################################
# Simulation for mass imputation
# Estimates: Naive NP estimate
# Author: Rebecca Andridge
# Last modified: 01/17/2026
########################################
library(tidyverse)
require(doParallel)
require(foreach)

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
  
  # Loop through relative sample sizes -----
  for (SAMPSIZES in 1:2)
  {
    print(paste0("RelSize:", SAMPSIZES))
    
    # Loop through NP selection mechanisms
    for (MECHNUM in 1:9)
    {
      print(paste0("Selection:", MECHNUM))
      
      # Load matrices with lists of sampled units
      filename2 <- paste0("./simdata/NPsamp_pop",POPNUM,"mech",MECHNUM,"relsizes",SAMPSIZES,".RData")
      load(filename2)
      
      NREPS <- ncol(Smat.NP)
      
      # Run in parallel -----
      res <- foreach(i=1:NREPS, .combine=rbind, .inorder=FALSE) %dopar%
        {
          # NP sample -----
          # Create sample using sampled indices
          samp.NP <- pop[Smat.NP[,i],]
          n_NP <- nrow(samp.NP)
          
          # (Naive Estimate) Estimates from nonprobability sample -----
          est_NP <- mean(samp.NP$Y)
          se_NP <- sd(samp.NP$Y)/sqrt(n_NP)
          ci <- confint(lm(Y~1, samp.NP))
          lb_NP <- ci[1]
          ub_NP <- ci[2]
          B_NP <- as.numeric(coef(lm(Y ~ Z1 + Z2, data = samp.NP)))
          
          # Return results -----
          estimates <- data.frame(cbind(est_NP, se_NP, lb_NP, ub_NP, 
                                        B0_NP = B_NP[1], B1_NP = B_NP[2], B2_NP = B_NP[3]))
          estimates$rep <- i
          
          estimates
        }
      
      # Add pop and selection mechanism and sample sizes numbers -----
      res$pop <- POPNUM
      res$mech <- MECHNUM
      res$sampsize <- SAMPSIZES
      
      # Save out results -----
      filename3 <- paste0("./simresults/naive_pop",POPNUM,"mech",MECHNUM,"relsizes",SAMPSIZES,".RData")
      save(res, file = filename3)
    }
  }
}

stop.time <- Sys.time()
# Deregister parallel processing
stopCluster(cl)
# check runtime
stop.time - start.time
# Time difference of 38.58224 secs