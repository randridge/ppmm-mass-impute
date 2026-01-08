########################################
# Simulation for mass imputation
# Estimates: Chen DR/MI/IPW approach
# Author: Rebecca Andridge
# Last modified: 06/27/2025
########################################
library(tidyverse)
require(doParallel)
require(foreach)
require(doRNG)
library(nonprobsvy)

# POPNUM <- 3
# MECHNUM <- 9

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
    res <- foreach(i=1:NREPS, .combine=rbind, .inorder=FALSE, .packages = c("nonprobsvy")) %dopar%
      {
        # Create P and NP samples using sampled indices -----
        # P sample (not including Y)
        samp.P <- pop[Smat.P[,i],c("Z1","Z2","A")]
        n_P <- nrow(samp.P)
        # NP sample
        samp.NP <- pop[Smat.NP[,i],]
        n_NP <- nrow(samp.NP)
        
        # Nonprobsvy methods -----
        # survey design object for reference sample
        # add weights for probability sample
        samp.P$W <- 100000/n_P
        ref.des <- svydesign(ids = ~1, weights = ~W, data = samp.P)
        # Doubly robust
        chen_dr <- nonprobsvy::nonprob(selection = ~ Z1 + Z2 + A, 
                                       outcome = Y ~ Z1 + Z2 + A, 
                                       data = samp.NP, 
                                       svydesign = ref.des, 
                                       method_outcome = "glm", 
                                       family_outcome = "gaussian")
        # Mass imputation
        chen_mi <- nonprobsvy::nonprob(outcome = Y ~ Z1 + Z2 + A, 
                                       data = samp.NP, 
                                       svydesign = ref.des, 
                                       method_outcome = "glm", 
                                       family_outcome = "gaussian")
        # IPW
        chen_ipw <- nonprobsvy::nonprob(selection = ~ Z1 + Z2 + A, 
                                        target = ~ Y, 
                                        data = samp.NP, 
                                        svydesign = ref.des)
        # estimated mean of Y from each method
        est_chendr  <- as.numeric(chen_dr$output[1])
        est_chenmi  <- as.numeric(chen_mi$output[1])
        est_chenipw <- as.numeric(chen_ipw$output[1])
        # standard errors
        se_chendr  <- as.numeric(chen_dr$output[2])
        se_chenmi  <- as.numeric(chen_mi$output[2])
        se_chenipw <- as.numeric(chen_ipw$output[2])
        # 95% CI
        lb_chendr  <- as.numeric(chen_dr$confidence_interval[1])
        ub_chendr  <- as.numeric(chen_dr$confidence_interval[2])
        lb_chenmi  <- as.numeric(chen_mi$confidence_interval[1])
        ub_chenmi  <- as.numeric(chen_mi$confidence_interval[2])
        lb_chenipw <- as.numeric(chen_ipw$confidence_interval[1])
        ub_chenipw <- as.numeric(chen_ipw$confidence_interval[2])
        
        # Return results -----
        data.frame(cbind(est_chendr,  se_chendr,  lb_chendr,  ub_chendr,
                         est_chenmi,  se_chenmi,  lb_chenmi,  ub_chenmi,
                         est_chenipw, se_chenipw, lb_chenipw, ub_chenipw))
      }
    
    # Add pop and mech numbers -----
    res$pop <- POPNUM
    res$mech <- MECHNUM
    
    # Save out results -----
    filename3 <- paste0("./simresults/chen_pop",POPNUM,"mech",MECHNUM,".RData")
    save(res, file = filename3)
  }
}
stop.time <- Sys.time()
# Deregister parallel processing
stopCluster(cl)
# check runtime
stop.time - start.time
# Time difference of 2.124186 mins
